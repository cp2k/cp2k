/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include "../../offload/offload_runtime.h"
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_NNP)

#include <stdio.h>

#include "nnp_gpu_internal.h"
#include "nnp_gpu_state.h"

#if defined(_OMP_H)
#error "OpenMP should not be used in .cu files to accommodate HIP."
#endif

/*
 * On-device cell-list walk and the per-step host-to-device uploads that feed
 * it. The walk kernel reproduces the host routine
 * nnp_compute_neighbors_cell_list: for every atom it visits the bins of its
 * local stencil, follows each bin's linked list of periodic images, applies
 * the same self-skip, image-match and cutoff filters as the host, and appends
 * the surviving neighbours into the per-(atom, group) walk buffers. The result
 * is the exact same neighbour SET the host produces; only the append order
 * within a set differs (see the kernel banner). The upload launchers copy the
 * coordinates, image table, cell-list bins and the species-pair routing table
 * the kernel reads, all onto the single NNP work stream so downstream kernels
 * observe the data through stream ordering.
 */

/*******************************************************************************
 * \brief Device helper for the 1-based linear bin index, matching the host
 *        routine nnp_cell_list_linear_index. bx, by, bz are 1-based bin
 *        coordinates; nbx, nby are the bin counts along x and y.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
__device__ __forceinline__ int nnp_walk_lin_idx_1based(int bx, int by, int bz,
                                                       int nbx, int nby) {
  return bx + nbx * ((by - 1) + nby * (bz - 1));
}

/*******************************************************************************
 * \brief Device version of the host routine nnp_cell_list_exact_image_match.
 *        perd is 0/1 per dimension (periodic flag); exact_pbc_copies bounds
 *        the accepted shift along each periodic dimension. Returns 1 to accept
 *        the image, 0 to reject it. coord_scaled is laid out as
 *        [3 * num_atoms_global], xyz-interleaved, matching Fortran
 *        column-major (3, n) storage when uploaded as a flat array. i_global
 *        and j_global are 0-based atom indices.
 *
 *        The Fortran NINT rounds half away from zero, whereas the CUDA
 *        intrinsic __double2int_rn rounds half to even. The two disagree only
 *        when the coordinate difference lands exactly on a half-integer, which
 *        continuous coordinates essentially never do; a portable
 *        round-half-away-from-zero is used here so the device result is
 *        bit-exact with the host.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
__device__ __forceinline__ int
nnp_walk_image_match(const double *coord_scaled, int i_global, int j_global,
                     int sx, int sy, int sz, const int *perd,
                     const int *exact_pbc_copies) {
  int shifts[3] = {sx, sy, sz};
  for (int d = 0; d < 3; d++) {
    if (perd[d] == 1) {
      double diff =
          coord_scaled[3 * i_global + d] - coord_scaled[3 * j_global + d];
      double rounded = (diff >= 0.0) ? floor(diff + 0.5) : -floor(-diff + 0.5);
      int minimal_shift = (int)rounded;
      int delta = minimal_shift - shifts[d];
      if (delta < 0)
        delta = -delta;
      if (delta > exact_pbc_copies[d])
        return 0;
    } else if (shifts[d] != 0) {
      return 0;
    }
  }
  return 1;
}

/*******************************************************************************
 * \brief On-device cell-list walk kernel. One warp handles one
 *        (atom_in_e_slot, e_idx) pair: blockIdx.x selects the 0-based atom
 *        slot within element e, blockIdx.y selects the 0-based element index,
 *        and the 32 lanes of the block share the work for that central atom.
 *        The kernel reproduces the host routine nnp_compute_neighbors_cell_list
 *        and appends the surviving neighbours into the per-(atom, group) flat
 *        buffers walk_*_ind_dev / walk_*_image_idx_dev, tracking per-(slot, s,
 *        e) counts in walk_n_*_pa_dev.
 *
 *        Walk math (line-for-line with the host):
 *          - Bin index from an INT-and-clamp on the central atom's primary
 *            coordinate (nnp_cell_list_bin_from_position), with a small width
 *            floor to avoid dividing by a degenerate bin width.
 *          - Linear bin index bx + nbx * ((by - 1) + nby * (bz - 1)) via
 *            nnp_walk_lin_idx_1based.
 *          - image_match: NINT(coord_scaled_i - coord_scaled_j) minus the
 *            image shift, bounded per periodic dimension by exact_pbc_copies
 *            (nnp_walk_image_match).
 *          - Displacement reconstruction:
 *              image_pos = coord(j) + image_translation(img),
 *              dr        = coord(i) - image_pos,
 *              norm      = sqrt(dr . dr).
 *          - Filtering: reject beyond the pair maximum cutoff, then per radial
 *            / ang1 / ang2 group with its per-(s, ele) cutoff.
 *
 *        Each lane strides over the flattened bin stencil and follows the
 *        linked-list chains; matches are appended with atomicAdd on the
 *        per-(atom, group) counter, so the neighbour append ORDER differs from
 *        the serial host walk. The downstream symmetry functions sum over the
 *        neighbour SET, so this only reorders the floating-point reduction and
 *        the results match the host. On a per-(slot, s, e) overflow the kernel
 *        sets overflow_flag_dev via atomicExch; the launcher reads the flag
 *        back and asks the caller to grow the buffers and retry.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
__global__ void nnp_cell_list_walk_kernel(
    int num_atoms_global, int n_ele, int n_atoms_in_e_max, int max_n_neigh_rad,
    int max_n_neigh_ang1, int max_n_neigh_ang2, int max_n_radgrp,
    int max_n_anggrp, int max_pair_grps, int nbx, int nby, int nbz,
    double bx_lo, double by_lo, double bz_lo, double bw_x, double bw_y,
    double bw_z, int span_x, int span_y, int span_z, int perd_x, int perd_y,
    int perd_z, int epbc_x, int epbc_y, int epbc_z,
    const double *__restrict__ coord_dev,
    const double *__restrict__ image_table_dev,
    const double *__restrict__ coord_scaled_dev,
    const int *__restrict__ ele_ind_global_dev,
    const int *__restrict__ image_atom_dev,
    const int *__restrict__ image_shift_dev, const int *__restrict__ head_dev,
    const int *__restrict__ next_dev,
    const int *__restrict__ atoms_global_per_e_dev,
    const int *__restrict__ n_atoms_in_e_dev,
    const int *__restrict__ pair_n_rad_dev,
    const int *__restrict__ pair_n_ang1_dev,
    const int *__restrict__ pair_n_ang2_dev,
    const double *__restrict__ pair_max_cutoff_dev,
    const int *__restrict__ pair_rad_groups_dev,
    const int *__restrict__ pair_ang1_groups_dev,
    const int *__restrict__ pair_ang2_groups_dev,
    const double *__restrict__ rad_cutoff_dev,
    const double *__restrict__ ang_cutoff_dev,
    int *__restrict__ walk_rad_ind_dev,
    int *__restrict__ walk_rad_image_idx_dev,
    int *__restrict__ walk_ang1_ind_dev,
    int *__restrict__ walk_ang1_image_idx_dev,
    int *__restrict__ walk_ang2_ind_dev,
    int *__restrict__ walk_ang2_image_idx_dev,
    int *__restrict__ walk_n_rad_pa_dev, int *__restrict__ walk_n_ang1_pa_dev,
    int *__restrict__ walk_n_ang2_pa_dev, int *__restrict__ overflow_flag_dev) {

  /* Block layout: blockIdx.x = atom_in_e_slot (0-based),
   *               blockIdx.y = e_idx (0-based). One warp per atom: every
   *               lane recomputes the (cheap) central-atom setup, then the
   *               stencil bins are split across threadIdx.x and matches are
   *               appended with atomicAdd on the per-(atom, group) counter. */
  const int slot = blockIdx.x;
  const int e0 = blockIdx.y;
  if (e0 >= n_ele)
    return;
  if (slot >= n_atoms_in_e_dev[e0])
    return;

  /* Pull i_global (1-based in host, 0-based on device for indexing). */
  const int i1 = atoms_global_per_e_dev[slot + n_atoms_in_e_max * e0];
  if (i1 <= 0)
    return;
  const int i0 = i1 - 1;
  const int i_ele1 = ele_ind_global_dev[i0]; /* 1-based */
  const int i_ele0 = i_ele1 - 1;

  const double cx = coord_dev[3 * i0 + 0];
  const double cy = coord_dev[3 * i0 + 1];
  const double cz = coord_dev[3 * i0 + 2];

  /* Bin index, matching host nnp_cell_list_bin_from_position. */
  const double rebuild_eps = 1.0e-12;
  const double bw_xs = (bw_x > rebuild_eps) ? bw_x : rebuild_eps;
  const double bw_ys = (bw_y > rebuild_eps) ? bw_y : rebuild_eps;
  const double bw_zs = (bw_z > rebuild_eps) ? bw_z : rebuild_eps;
  int bin_ix = (int)((cx - bx_lo) / bw_xs) + 1;
  int bin_iy = (int)((cy - by_lo) / bw_ys) + 1;
  int bin_iz = (int)((cz - bz_lo) / bw_zs) + 1;
  if (bin_ix < 1)
    bin_ix = 1;
  else if (bin_ix > nbx)
    bin_ix = nbx;
  if (bin_iy < 1)
    bin_iy = 1;
  else if (bin_iy > nby)
    bin_iy = nby;
  if (bin_iz < 1)
    bin_iz = 1;
  else if (bin_iz > nbz)
    bin_iz = nbz;

  const int bx_lo_i = max(1, bin_ix - span_x);
  const int bx_hi_i = min(nbx, bin_ix + span_x);
  const int by_lo_i = max(1, bin_iy - span_y);
  const int by_hi_i = min(nby, bin_iy + span_y);
  const int bz_lo_i = max(1, bin_iz - span_z);
  const int bz_hi_i = min(nbz, bin_iz + span_z);

  /* Bin walk: flatten the stencil and stride the bins across the warp.
   * Each lane scans its bins' chains; matches append via atomicAdd, so the
   * neighbour append order differs from the serial host walk order. The SF
   * sums over the neighbour SET, so this only reorders the FP reduction. */
  const int nbx_w = bx_hi_i - bx_lo_i + 1;
  const int nby_w = by_hi_i - by_lo_i + 1;
  const int nbz_w = bz_hi_i - bz_lo_i + 1;
  const int n_bins_w = nbx_w * nby_w * nbz_w;
  for (int bsel = threadIdx.x; bsel < n_bins_w; bsel += blockDim.x) {
    const int bx = bx_lo_i + (bsel % nbx_w);
    const int by = by_lo_i + ((bsel / nbx_w) % nby_w);
    const int bz = bz_lo_i + (bsel / (nbx_w * nby_w));
    const int lin = nnp_walk_lin_idx_1based(bx, by, bz, nbx, nby);
    int img1 = head_dev[lin - 1]; /* head is 1-based-indexed Fortran array */
    while (img1 > 0) {
      const int j1 = image_atom_dev[img1 - 1];
      const int j0 = j1 - 1;
      const int sx = image_shift_dev[3 * (img1 - 1) + 0];
      const int sy_ = image_shift_dev[3 * (img1 - 1) + 1];
      const int sz_ = image_shift_dev[3 * (img1 - 1) + 2];

      /* Self-skip at zero shift */
      if (j1 == i1 && sx == 0 && sy_ == 0 && sz_ == 0) {
        img1 = next_dev[img1 - 1];
        continue;
      }
      /* Image-match check */
      const int perd_arr[3] = {perd_x, perd_y, perd_z};
      const int epbc_arr[3] = {epbc_x, epbc_y, epbc_z};
      if (!nnp_walk_image_match(coord_scaled_dev, i0, j0, sx, sy_, sz_,
                                perd_arr, epbc_arr)) {
        img1 = next_dev[img1 - 1];
        continue;
      }

      const double image_x =
          coord_dev[3 * j0 + 0] + image_table_dev[3 * (img1 - 1) + 0];
      const double image_y =
          coord_dev[3 * j0 + 1] + image_table_dev[3 * (img1 - 1) + 1];
      const double image_z =
          coord_dev[3 * j0 + 2] + image_table_dev[3 * (img1 - 1) + 2];
      const double drx = cx - image_x;
      const double dry = cy - image_y;
      const double drz = cz - image_z;
      const double norm = sqrt(drx * drx + dry * dry + drz * drz);

      const int j_ele1 = ele_ind_global_dev[j0];
      const int j_ele0 = j_ele1 - 1;
      const int pair_idx = i_ele0 + n_ele * j_ele0;
      const double pair_max = pair_max_cutoff_dev[pair_idx];

      if (norm < pair_max) {
        /* Radial */
        const int n_rad_pair = pair_n_rad_dev[pair_idx];
        for (int g = 0; g < n_rad_pair; g++) {
          const int s1 = pair_rad_groups_dev[g + max_pair_grps * pair_idx];
          const int s0 = s1 - 1;
          const double cutoff = rad_cutoff_dev[s0 + max_n_radgrp * i_ele0];
          if (norm < cutoff) {
            /* 64-bit offsets, matching the reader side in nnp_gpu_pack.h:
             * the neighbour tables exceed 2^31 entries on large ranks. */
            const long long pa_off =
                (long long)slot + (long long)n_atoms_in_e_max *
                                      ((long long)s0 + (long long)max_n_radgrp *
                                                           (long long)i_ele0);
            const int k = atomicAdd(&walk_n_rad_pa_dev[pa_off], 1);
            if (k >= max_n_neigh_rad) {
              atomicExch(overflow_flag_dev, 1);
            } else {
              const long long out_off =
                  (long long)k + (long long)max_n_neigh_rad * pa_off;
              walk_rad_ind_dev[out_off] = j0;             /* 0-based */
              walk_rad_image_idx_dev[out_off] = img1 - 1; /* 0-based */
            }
          }
        }
        /* Ang1 */
        const int n_ang1_pair = pair_n_ang1_dev[pair_idx];
        for (int g = 0; g < n_ang1_pair; g++) {
          const int s1 = pair_ang1_groups_dev[g + max_pair_grps * pair_idx];
          const int s0 = s1 - 1;
          const double cutoff = ang_cutoff_dev[s0 + max_n_anggrp * i_ele0];
          if (norm < cutoff) {
            const long long pa_off =
                (long long)slot + (long long)n_atoms_in_e_max *
                                      ((long long)s0 + (long long)max_n_anggrp *
                                                           (long long)i_ele0);
            const int k = atomicAdd(&walk_n_ang1_pa_dev[pa_off], 1);
            if (k >= max_n_neigh_ang1) {
              atomicExch(overflow_flag_dev, 1);
            } else {
              const long long out_off =
                  (long long)k + (long long)max_n_neigh_ang1 * pa_off;
              walk_ang1_ind_dev[out_off] = j0;             /* 0-based */
              walk_ang1_image_idx_dev[out_off] = img1 - 1; /* 0-based */
            }
          }
        }
        /* Ang2 */
        const int n_ang2_pair = pair_n_ang2_dev[pair_idx];
        for (int g = 0; g < n_ang2_pair; g++) {
          const int s1 = pair_ang2_groups_dev[g + max_pair_grps * pair_idx];
          const int s0 = s1 - 1;
          const double cutoff = ang_cutoff_dev[s0 + max_n_anggrp * i_ele0];
          if (norm < cutoff) {
            const long long pa_off =
                (long long)slot + (long long)n_atoms_in_e_max *
                                      ((long long)s0 + (long long)max_n_anggrp *
                                                           (long long)i_ele0);
            const int k = atomicAdd(&walk_n_ang2_pa_dev[pa_off], 1);
            if (k >= max_n_neigh_ang2) {
              atomicExch(overflow_flag_dev, 1);
            } else {
              const long long out_off =
                  (long long)k + (long long)max_n_neigh_ang2 * pa_off;
              walk_ang2_ind_dev[out_off] = j0;             /* 0-based */
              walk_ang2_image_idx_dev[out_off] = img1 - 1; /* 0-based */
            }
          }
        }
      }

      img1 = next_dev[img1 - 1];
    }
  }
}

/*******************************************************************************
 * \brief Per-MD-step host-to-device copy of cell_list_cache%coord_primary into
 *        state->coord_dev. h_coord is a contiguous 3*num_atoms host buffer
 *        (xyz interleaved). The device buffer is grow-only via
 *        nnp_gpu_pbuf_ensure. The copy is queued on state->stream and the
 *        launcher returns immediately; downstream kernels on the same stream
 *        observe the updated positions in order.
 *
 *        This is the first GPU call of every MD step, so the observability step
 *        counter state->md_step is incremented here.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
extern "C" int nnp_gpu_upload_positions_step_c(int num_atoms,
                                               const double *h_coord) {
  if (num_atoms <= 0) {
    fprintf(stderr, "ERROR: nnp_gpu_upload_positions_step: bad num_atoms=%d\n",
            num_atoms);
    return -1;
  }
  nnp_gpu_state_t *state = nnp_gpu_state_get();
  if (state == NULL)
    return -1;

  /* First GPU call of every MD step. The observability step counter is
   * advanced here so that per-step device work can be attributed. */
  state->md_step++;

  const size_t bytes = sizeof(double) * 3 * (size_t)num_atoms;
  if (nnp_gpu_pbuf_ensure((void **)&state->coord_dev, &state->coord_dev_cap,
                          bytes))
    return -1;
  offloadMemcpyAsyncHtoD(state->coord_dev, h_coord, bytes, state->stream);
  return 0;
}

/*******************************************************************************
 * \brief Per-MD-step host-to-device copy of cell_list_cache%image_translation
 *        into state->image_table_dev. h_image_table is a contiguous
 *        3*n_images host buffer (xyz per image, the columns of the Fortran
 *        image_translation(3, n_images) array). The device buffer is grow-only
 *        and the copy is queued on state->stream. The image table changes only
 *        when the cell-list bins rebuild, so the driver calls this only on
 *        rebuild steps.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
extern "C" int nnp_gpu_upload_image_table_c(int n_images,
                                            const double *h_image_table) {
  if (n_images <= 0) {
    fprintf(stderr, "ERROR: nnp_gpu_upload_image_table: bad n_images=%d\n",
            n_images);
    return -1;
  }
  nnp_gpu_state_t *state = nnp_gpu_state_get();
  if (state == NULL)
    return -1;

  const size_t bytes = sizeof(double) * 3 * (size_t)n_images;
  if (nnp_gpu_pbuf_ensure((void **)&state->image_table_dev,
                          &state->image_table_dev_cap, bytes))
    return -1;
  offloadMemcpyAsyncHtoD(state->image_table_dev, h_image_table, bytes,
                         state->stream);
  return 0;
}

/*******************************************************************************
 * \brief Host-to-device copy of the cell-list bins (image_atom, image_shift,
 *        head, next) into state, read by the walk kernel. Sizes:
 *          image_atom:  n_images ints
 *          image_shift: 3 * n_images ints
 *          head:        n_cells ints
 *          next:        n_images ints
 *        The device buffers are grow-only via nnp_gpu_pbuf_ensure. The four
 *        copies are queued on state->stream and the stream is synchronized
 *        before returning: the caller passes kind-converted temporaries whose
 *        lifetime ends at the end of the call statement. Called whenever the
 *        cell-list is rebinned.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
extern "C" int nnp_gpu_upload_cell_list_c(int n_images, int n_cells,
                                          const int *h_image_atom,
                                          const int *h_image_shift,
                                          const int *h_head,
                                          const int *h_next) {
  if (n_images <= 0 || n_cells <= 0) {
    fprintf(stderr,
            "ERROR: nnp_gpu_upload_cell_list: bad sizes "
            "n_images=%d n_cells=%d\n",
            n_images, n_cells);
    return -1;
  }
  nnp_gpu_state_t *state = nnp_gpu_state_get();
  if (state == NULL)
    return -1;

  const size_t bytes_imgs = sizeof(int) * (size_t)n_images;
  const size_t bytes_shifts = sizeof(int) * 3 * (size_t)n_images;
  const size_t bytes_head = sizeof(int) * (size_t)n_cells;

  if (nnp_gpu_pbuf_ensure((void **)&state->cell_image_atom_dev,
                          &state->cell_image_atom_dev_cap, bytes_imgs))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->cell_image_shift_dev,
                          &state->cell_image_shift_dev_cap, bytes_shifts))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->cell_head_dev,
                          &state->cell_head_dev_cap, bytes_head))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->cell_next_dev,
                          &state->cell_next_dev_cap, bytes_imgs))
    return -1;

  offloadMemcpyAsyncHtoD(state->cell_image_atom_dev, h_image_atom, bytes_imgs,
                         state->stream);
  offloadMemcpyAsyncHtoD(state->cell_image_shift_dev, h_image_shift,
                         bytes_shifts, state->stream);
  offloadMemcpyAsyncHtoD(state->cell_head_dev, h_head, bytes_head,
                         state->stream);
  offloadMemcpyAsyncHtoD(state->cell_next_dev, h_next, bytes_imgs,
                         state->stream);
  /* Sync so the host buffers can go out of scope safely after return. */
  offloadStreamSynchronize(state->stream);
  return 0;
}

/*******************************************************************************
 * \brief Once-at-NNP-init host-to-device copy of the species-pair routing
 *        table. The caller passes host arrays already flat in the layouts the
 *        walk kernel expects. The upload is idempotent: a second call with the
 *        same shape detects state->pair_map_uploaded and returns 0 without
 *        re-uploading, and a call with a different shape is rejected with an
 *        error. On the first call the device buffers are grown, the tables are
 *        copied on state->stream, the stream is synchronized so the host
 *        buffers may be freed on return, and the recorded shape fields are set.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
extern "C" int nnp_gpu_upload_pair_map_c(
    int n_ele, int max_pair_grps, int max_n_radgrp, int max_n_anggrp,
    const int *h_n_rad, const int *h_n_ang1, const int *h_n_ang2,
    const double *h_pair_max_cutoff, const int *h_rad_groups,
    const int *h_ang1_groups, const int *h_ang2_groups,
    const double *h_rad_cutoff, const double *h_ang_cutoff) {
  if (n_ele <= 0 || max_pair_grps <= 0 || max_n_radgrp <= 0 ||
      max_n_anggrp <= 0) {
    fprintf(stderr,
            "ERROR: nnp_gpu_upload_pair_map: bad sizes "
            "n_ele=%d max_pair_grps=%d max_n_radgrp=%d max_n_anggrp=%d\n",
            n_ele, max_pair_grps, max_n_radgrp, max_n_anggrp);
    return -1;
  }
  nnp_gpu_state_t *state = nnp_gpu_state_get();
  if (state == NULL)
    return -1;

  if (state->pair_map_uploaded) {
    if (state->pair_map_n_ele == n_ele &&
        state->pair_map_max_pair_grps == max_pair_grps &&
        state->pair_map_max_n_radgrp == max_n_radgrp &&
        state->pair_map_max_n_anggrp == max_n_anggrp) {
      return 0;
    }
    fprintf(stderr,
            "ERROR: nnp_gpu_upload_pair_map: shape mismatch on second "
            "upload (n_ele %d->%d, max_pair_grps %d->%d, "
            "max_n_radgrp %d->%d, max_n_anggrp %d->%d)\n",
            state->pair_map_n_ele, n_ele, state->pair_map_max_pair_grps,
            max_pair_grps, state->pair_map_max_n_radgrp, max_n_radgrp,
            state->pair_map_max_n_anggrp, max_n_anggrp);
    return -1;
  }

  const size_t n_pairs = (size_t)n_ele * (size_t)n_ele;
  const size_t bytes_counts = sizeof(int) * n_pairs;
  const size_t bytes_max_cut = sizeof(double) * n_pairs;
  const size_t bytes_groups = sizeof(int) * (size_t)max_pair_grps * n_pairs;
  const size_t bytes_rcut =
      sizeof(double) * (size_t)max_n_radgrp * (size_t)n_ele;
  const size_t bytes_acut =
      sizeof(double) * (size_t)max_n_anggrp * (size_t)n_ele;

  if (nnp_gpu_pbuf_ensure((void **)&state->pair_n_rad_dev,
                          &state->pair_n_rad_dev_cap, bytes_counts))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pair_n_ang1_dev,
                          &state->pair_n_ang1_dev_cap, bytes_counts))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pair_n_ang2_dev,
                          &state->pair_n_ang2_dev_cap, bytes_counts))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pair_max_cutoff_dev,
                          &state->pair_max_cutoff_dev_cap, bytes_max_cut))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pair_rad_groups_dev,
                          &state->pair_rad_groups_dev_cap, bytes_groups))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pair_ang1_groups_dev,
                          &state->pair_ang1_groups_dev_cap, bytes_groups))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pair_ang2_groups_dev,
                          &state->pair_ang2_groups_dev_cap, bytes_groups))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->rad_cutoff_dev,
                          &state->rad_cutoff_dev_cap, bytes_rcut))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->ang_cutoff_dev,
                          &state->ang_cutoff_dev_cap, bytes_acut))
    return -1;

  offloadMemcpyAsyncHtoD(state->pair_n_rad_dev, h_n_rad, bytes_counts,
                         state->stream);
  offloadMemcpyAsyncHtoD(state->pair_n_ang1_dev, h_n_ang1, bytes_counts,
                         state->stream);
  offloadMemcpyAsyncHtoD(state->pair_n_ang2_dev, h_n_ang2, bytes_counts,
                         state->stream);
  offloadMemcpyAsyncHtoD(state->pair_max_cutoff_dev, h_pair_max_cutoff,
                         bytes_max_cut, state->stream);
  offloadMemcpyAsyncHtoD(state->pair_rad_groups_dev, h_rad_groups, bytes_groups,
                         state->stream);
  offloadMemcpyAsyncHtoD(state->pair_ang1_groups_dev, h_ang1_groups,
                         bytes_groups, state->stream);
  offloadMemcpyAsyncHtoD(state->pair_ang2_groups_dev, h_ang2_groups,
                         bytes_groups, state->stream);
  offloadMemcpyAsyncHtoD(state->rad_cutoff_dev, h_rad_cutoff, bytes_rcut,
                         state->stream);
  offloadMemcpyAsyncHtoD(state->ang_cutoff_dev, h_ang_cutoff, bytes_acut,
                         state->stream);

  /* Sync so the host buffers can go out of scope safely after return. */
  offloadStreamSynchronize(state->stream);

  state->pair_map_uploaded = 1;
  state->pair_map_n_ele = n_ele;
  state->pair_map_max_pair_grps = max_pair_grps;
  state->pair_map_max_n_radgrp = max_n_radgrp;
  state->pair_map_max_n_anggrp = max_n_anggrp;
  return 0;
}

/*******************************************************************************
 * \brief Launcher for nnp_cell_list_walk_kernel. Exported with C linkage so
 *        Fortran calls it via ISO_C_BINDING. Uploads the per-step support
 *        buffers (scaled coordinates, global element indices, the per-element
 *        atom lists and counts), grows and zeroes the per-(slot, s, e) walk
 *        output and count buffers, clears the overflow flag, and launches the
 *        walk with grid (n_atoms_in_e_max, n_ele) and one warp per block. It
 *        then synchronizes the stream and reads the overflow flag back: on
 *        overflow it returns -2 so the caller can grow the per-atom neighbour
 *        capacities and retry, otherwise it returns 0. The overflow flag lives
 *        in state->walk_overflow_dev, a small grow-only device buffer.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
extern "C" int nnp_gpu_cell_list_walk_c(
    int num_atoms_global, int n_ele, int n_images, int n_cells,
    int n_atoms_in_e_max, int max_n_neigh_rad, int max_n_neigh_ang1,
    int max_n_neigh_ang2, int max_n_radgrp, int max_n_anggrp, const int *h_nbin,
    const double *h_bin_lower, const double *h_bin_width, const int *h_bin_span,
    const int *h_perd, const int *h_exact_pbc_copies,
    const double *h_coord_scaled, const int *h_ele_ind_global,
    const int *h_atoms_global_per_e, const int *h_n_atoms_in_e) {
  if (num_atoms_global <= 0 || n_ele <= 0 || n_images <= 0 || n_cells <= 0 ||
      n_atoms_in_e_max <= 0 || max_n_neigh_rad <= 0 || max_n_neigh_ang1 <= 0 ||
      max_n_neigh_ang2 <= 0 || max_n_radgrp <= 0 || max_n_anggrp <= 0) {
    fprintf(stderr,
            "ERROR: nnp_gpu_cell_list_walk: bad sizes "
            "num_atoms_global=%d n_ele=%d n_images=%d n_cells=%d\n",
            num_atoms_global, n_ele, n_images, n_cells);
    return -1;
  }
  nnp_gpu_state_t *state = nnp_gpu_state_get();
  if (state == NULL)
    return -1;
  if (!state->pair_map_uploaded) {
    fprintf(stderr, "ERROR: nnp_gpu_cell_list_walk: pair_map not uploaded "
                    "(call nnp_gpu_upload_pair_map_c first)\n");
    return -1;
  }
  if (state->pair_map_max_pair_grps <= 0) {
    fprintf(stderr,
            "ERROR: nnp_gpu_cell_list_walk: pair_map_max_pair_grps=%d\n",
            state->pair_map_max_pair_grps);
    return -1;
  }
  /* The walk chases the cell-list chain and the image table, both of which are
   * uploaded by earlier launchers. Check them here rather than let a NULL or a
   * short buffer surface as an illegal access at the next stream sync, which
   * would blame the sync instead of the missing upload. */
  if (state->coord_dev == NULL || state->image_table_dev == NULL ||
      state->cell_head_dev == NULL || state->cell_next_dev == NULL ||
      state->cell_image_atom_dev == NULL ||
      state->cell_image_shift_dev == NULL) {
    fprintf(stderr, "ERROR: nnp_gpu_cell_list_walk: positions, image table or "
                    "cell-list chain not uploaded for this step\n");
    return -1;
  }
  if (state->cell_head_dev_cap < sizeof(int) * (size_t)n_cells ||
      state->image_table_dev_cap < sizeof(double) * 3 * (size_t)n_images ||
      state->cell_image_atom_dev_cap < sizeof(int) * (size_t)n_images ||
      state->cell_image_shift_dev_cap < sizeof(int) * 3 * (size_t)n_images) {
    fprintf(stderr,
            "ERROR: nnp_gpu_cell_list_walk: device cell list is sized for a "
            "smaller geometry than this call's n_cells=%d n_images=%d\n",
            n_cells, n_images);
    return -1;
  }
  const int max_pair_grps = state->pair_map_max_pair_grps;

  const size_t bytes_coord_scaled =
      sizeof(double) * 3 * (size_t)num_atoms_global;
  const size_t bytes_ele_ind = sizeof(int) * (size_t)num_atoms_global;
  const size_t bytes_atoms_pe =
      sizeof(int) * (size_t)n_atoms_in_e_max * (size_t)n_ele;
  const size_t bytes_n_in_e = sizeof(int) * (size_t)n_ele;
  const size_t pa_count_rad =
      (size_t)n_atoms_in_e_max * (size_t)max_n_radgrp * (size_t)n_ele;
  const size_t pa_count_ang =
      (size_t)n_atoms_in_e_max * (size_t)max_n_anggrp * (size_t)n_ele;
  const size_t bytes_rad_ind =
      sizeof(int) * (size_t)max_n_neigh_rad * pa_count_rad;
  const size_t bytes_ang1_ind =
      sizeof(int) * (size_t)max_n_neigh_ang1 * pa_count_ang;
  const size_t bytes_ang2_ind =
      sizeof(int) * (size_t)max_n_neigh_ang2 * pa_count_ang;
  const size_t bytes_n_rad_pa = sizeof(int) * pa_count_rad;
  const size_t bytes_n_ang1_pa = sizeof(int) * pa_count_ang;
  const size_t bytes_n_ang2_pa = sizeof(int) * pa_count_ang;
  int h_overflow = 0;

  /* Per-step support buffer uploads. */
  if (nnp_gpu_pbuf_ensure((void **)&state->coord_scaled_dev,
                          &state->coord_scaled_dev_cap, bytes_coord_scaled))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->ele_ind_global_dev,
                          &state->ele_ind_global_dev_cap, bytes_ele_ind))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->walk_atoms_global_dev,
                          &state->walk_atoms_global_dev_cap, bytes_atoms_pe))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->walk_n_atoms_in_e_dev,
                          &state->walk_n_atoms_in_e_dev_cap, bytes_n_in_e))
    return -1;

  offloadMemcpyAsyncHtoD(state->coord_scaled_dev, h_coord_scaled,
                         bytes_coord_scaled, state->stream);
  offloadMemcpyAsyncHtoD(state->ele_ind_global_dev, h_ele_ind_global,
                         bytes_ele_ind, state->stream);
  offloadMemcpyAsyncHtoD(state->walk_atoms_global_dev, h_atoms_global_per_e,
                         bytes_atoms_pe, state->stream);
  offloadMemcpyAsyncHtoD(state->walk_n_atoms_in_e_dev, h_n_atoms_in_e,
                         bytes_n_in_e, state->stream);

  /* Output buffer ensure + zero. */
  if (nnp_gpu_pbuf_ensure((void **)&state->walk_rad_ind_dev,
                          &state->walk_rad_ind_dev_cap, bytes_rad_ind))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->walk_rad_image_idx_dev,
                          &state->walk_rad_image_idx_dev_cap, bytes_rad_ind))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->walk_ang1_ind_dev,
                          &state->walk_ang1_ind_dev_cap, bytes_ang1_ind))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->walk_ang1_image_idx_dev,
                          &state->walk_ang1_image_idx_dev_cap, bytes_ang1_ind))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->walk_ang2_ind_dev,
                          &state->walk_ang2_ind_dev_cap, bytes_ang2_ind))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->walk_ang2_image_idx_dev,
                          &state->walk_ang2_image_idx_dev_cap, bytes_ang2_ind))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->walk_n_rad_pa_dev,
                          &state->walk_n_rad_pa_dev_cap, bytes_n_rad_pa))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->walk_n_ang1_pa_dev,
                          &state->walk_n_ang1_pa_dev_cap, bytes_n_ang1_pa))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->walk_n_ang2_pa_dev,
                          &state->walk_n_ang2_pa_dev_cap, bytes_n_ang2_pa))
    return -1;

  offloadMemsetAsync(state->walk_n_rad_pa_dev, 0, bytes_n_rad_pa,
                     state->stream);
  offloadMemsetAsync(state->walk_n_ang1_pa_dev, 0, bytes_n_ang1_pa,
                     state->stream);
  offloadMemsetAsync(state->walk_n_ang2_pa_dev, 0, bytes_n_ang2_pa,
                     state->stream);

  /* Overflow flag: small persistent buffer, grown once and reused. */
  if (nnp_gpu_pbuf_ensure((void **)&state->walk_overflow_dev,
                          &state->walk_overflow_cap, sizeof(int)))
    return -1;
  offloadMemsetAsync(state->walk_overflow_dev, 0, sizeof(int), state->stream);

  /* Launch */
  {
    dim3 grid(n_atoms_in_e_max, n_ele, 1);
    dim3 block(32, 1, 1);
    nnp_cell_list_walk_kernel<<<grid, block, 0, state->stream>>>(
        num_atoms_global, n_ele, n_atoms_in_e_max, max_n_neigh_rad,
        max_n_neigh_ang1, max_n_neigh_ang2, max_n_radgrp, max_n_anggrp,
        max_pair_grps, h_nbin[0], h_nbin[1], h_nbin[2], h_bin_lower[0],
        h_bin_lower[1], h_bin_lower[2], h_bin_width[0], h_bin_width[1],
        h_bin_width[2], h_bin_span[0], h_bin_span[1], h_bin_span[2], h_perd[0],
        h_perd[1], h_perd[2], h_exact_pbc_copies[0], h_exact_pbc_copies[1],
        h_exact_pbc_copies[2], state->coord_dev, state->image_table_dev,
        state->coord_scaled_dev, state->ele_ind_global_dev,
        state->cell_image_atom_dev, state->cell_image_shift_dev,
        state->cell_head_dev, state->cell_next_dev,
        state->walk_atoms_global_dev, state->walk_n_atoms_in_e_dev,
        state->pair_n_rad_dev, state->pair_n_ang1_dev, state->pair_n_ang2_dev,
        state->pair_max_cutoff_dev, state->pair_rad_groups_dev,
        state->pair_ang1_groups_dev, state->pair_ang2_groups_dev,
        state->rad_cutoff_dev, state->ang_cutoff_dev, state->walk_rad_ind_dev,
        state->walk_rad_image_idx_dev, state->walk_ang1_ind_dev,
        state->walk_ang1_image_idx_dev, state->walk_ang2_ind_dev,
        state->walk_ang2_image_idx_dev, state->walk_n_rad_pa_dev,
        state->walk_n_ang1_pa_dev, state->walk_n_ang2_pa_dev,
        state->walk_overflow_dev);
    OFFLOAD_CHECK(offloadPeekAtLastError());
  }

  /* Sync + read overflow flag back. */
  offloadStreamSynchronize(state->stream);
  offloadMemcpyDtoH(&h_overflow, state->walk_overflow_dev, sizeof(int));
  if (h_overflow) {
    fprintf(stderr, "ERROR: nnp_gpu_cell_list_walk: per-(slot, s, e) overflow "
                    "(grow max_n_neigh_rad/ang1/ang2 and retry)\n");
    return -2;
  }
  return 0;
}

/*******************************************************************************
 * \brief Device-to-host copy of the per-(slot, s, e) walk count buffers into
 *        the caller's host arrays, then synchronizes the stream so the host
 *        buffers are valid on return. Used by the Fortran host-versus-device
 *        walk comparison. Cheap (well under a megabyte for typical systems).
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
extern "C" int nnp_gpu_walk_counts_d2h_c(int n_atoms_in_e_max, int max_n_radgrp,
                                         int max_n_anggrp, int n_ele,
                                         int *h_n_rad_pa, int *h_n_ang1_pa,
                                         int *h_n_ang2_pa) {
  if (n_atoms_in_e_max <= 0 || max_n_radgrp <= 0 || max_n_anggrp <= 0 ||
      n_ele <= 0) {
    fprintf(stderr, "ERROR: nnp_gpu_walk_counts_d2h: bad sizes\n");
    return -1;
  }
  nnp_gpu_state_t *state = nnp_gpu_state_get();
  if (state == NULL)
    return -1;
  if (state->walk_n_rad_pa_dev == NULL || state->walk_n_ang1_pa_dev == NULL ||
      state->walk_n_ang2_pa_dev == NULL) {
    fprintf(stderr, "ERROR: nnp_gpu_walk_counts_d2h: walk count buffers not "
                    "allocated (call nnp_gpu_cell_list_walk_c first)\n");
    return -1;
  }
  const size_t pa_count_rad =
      (size_t)n_atoms_in_e_max * (size_t)max_n_radgrp * (size_t)n_ele;
  const size_t pa_count_ang =
      (size_t)n_atoms_in_e_max * (size_t)max_n_anggrp * (size_t)n_ele;

  offloadMemcpyAsyncDtoH(h_n_rad_pa, state->walk_n_rad_pa_dev,
                         sizeof(int) * pa_count_rad, state->stream);
  offloadMemcpyAsyncDtoH(h_n_ang1_pa, state->walk_n_ang1_pa_dev,
                         sizeof(int) * pa_count_ang, state->stream);
  offloadMemcpyAsyncDtoH(h_n_ang2_pa, state->walk_n_ang2_pa_dev,
                         sizeof(int) * pa_count_ang, state->stream);
  offloadStreamSynchronize(state->stream);
  return 0;
}

#endif /* defined(__OFFLOAD) && !defined(__NO_OFFLOAD_NNP) */
