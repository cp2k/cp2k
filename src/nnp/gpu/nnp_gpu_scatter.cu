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

#include "nnp_gpu_pack.h"

/*
 * Force assembly. The descriptor gradient the network leaves on the device
 * is multiplied by the descriptor derivatives the symmetry-function pass
 * stored per (group, element), and the products are accumulated into the
 * force array. When the caller asks for a stress tensor the same kernels
 * accumulate the analytic virial alongside the forces.
 */

/*******************************************************************************
 * \brief Descriptor-derivative to force scale factor for one descriptor row.
 *        Reproduces the host scaling applied when the network's descriptor
 *        gradient is turned into a force contribution. Scaling modes 1 and 3
 *        map the descriptor range [loc_min, loc_max] onto [scmin, scmax];
 *        mode 4 uses the descriptor standard deviation; any other mode leaves
 *        the gradient unscaled.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
static __device__ __forceinline__ double
nnp_compute_force_scale(int scaling_mode, double scmin, double scmax, int idx,
                        const double *__restrict__ scaling_loc_max,
                        const double *__restrict__ scaling_loc_min,
                        const double *__restrict__ scaling_sigma) {
  if (scaling_mode == 1 || scaling_mode == 3) {
    const double diff = scaling_loc_max[idx] - scaling_loc_min[idx];
    return (scmax - scmin) / diff;
  }
  if (scaling_mode == 4) {
    return (scmax - scmin) / scaling_sigma[idx];
  }
  return 1.0;
}

/*******************************************************************************
 * \brief Warp-reduce the 9 per-thread analytical-virial components across the
 *        32-lane scatter block and have lane 0 commit a single atomicAdd set
 *        into the per-committee accumulator virial_dev[9*c + 3*a + b]. Keeps
 *        global atomic traffic at 9 adds per (atom, sf, committee) block
 *        instead of 9 per neighbour. Assumes blockDim.x == warpSize (the
 *        scatter launches use blk(32, 1, 1)); only called from the scatter
 *        launchers, which do.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
__device__ __forceinline__ void
nnp_virial_warp_flush(const double vir[3][3], double *__restrict__ virial_dev,
                      const int c) {
#pragma unroll
  for (int a = 0; a < 3; ++a) {
#pragma unroll
    for (int b = 0; b < 3; ++b) {
      /* Butterfly XOR reduction (backend-portable) leaves the full 32-lane sum
       * in every lane; lane 0 commits the single atomicAdd. */
      double v = vir[a][b];
      for (int mask = 16; mask > 0; mask >>= 1)
        v += NNP_SHFL_XOR(v, mask);
      if (threadIdx.x == 0)
        atomicAdd(&virial_dev[9 * c + 3 * a + b], v);
    }
  }
}

/*******************************************************************************
 * \brief Radial force-scatter kernel: one block per (atom in element, symmetry
 *        function, committee member). Thread 0 adds the self contribution to
 *        the central atom; the block's threads share the per-neighbour loop,
 *        each atomic-adding one neighbour's contribution per Cartesian
 *        component into the device force accumulator. Sign convention:
 *        f[i] -= de * self, f[k_atom] += de * force. When accumulate_virial is
 *        set, each per-neighbour force contribution also accumulates
 *        -r_ij (outer) f_k into the per-committee analytical virial.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
__global__ void nnp_force_scatter_radial_kernel(
    const int n_atoms_in_e, const int n_symf, const int max_n_neigh,
    const int n_committee, const int n_input_max, const int num_atoms_global,
    const int ele_idx, const int scaling_mode, const int n_rad_offset,
    const double scmin, const double scmax,
    const int *__restrict__ symf_map_one_based,
    const int *__restrict__ atom_map, const int *__restrict__ atoms_global,
    const int *__restrict__ n_neigh_per_atom,
    const int *__restrict__ rad_ind_one_based,
    const double *__restrict__ self_partial,
    const double *__restrict__ force_rad, const double *__restrict__ dE_dG,
    const double *__restrict__ scaling_loc_max,
    const double *__restrict__ scaling_loc_min,
    const double *__restrict__ scaling_sigma, double *__restrict__ f_dev,
    /* Virial (read only when accumulate_virial != 0): rad_image_idx pairs
     * 1:1 with rad_ind_one_based; coord_dev + image_table_dev rebuild the
     * bond vector r_ij = c_i - (c_k + image), matching the SF kernel. */
    const int *__restrict__ rad_image_idx, const double *__restrict__ coord_dev,
    const double *__restrict__ image_table_dev, double *__restrict__ virial_dev,
    const int accumulate_virial) {
  const int atom_in_e = blockIdx.x;
  const int sf = blockIdx.y;
  const int c = blockIdx.z;
  if (atom_in_e >= n_atoms_in_e || sf >= n_symf || c >= n_committee)
    return;

  /* (void) suppress unused-arg warning when n_rad_offset is 0 (radial). */
  (void)n_rad_offset;

  const int i_local = atom_map[atom_in_e];
  const int i_global = atoms_global[atom_in_e];
  /* m is the descriptor row (1-based in symfgrp%symf, 0-based in
   * dE_dG / scaling tables). Radial: m = symf(sf). */
  const int m = symf_map_one_based[sf] - 1;

  /* 64-bit offsets on every per-atom-scaled buffer, matching nnp_gpu_pack.h:
   * these tables exceed 2^31 entries on large ranks. */
  const double dEdG = dE_dG[m + n_input_max * c +
                            (long long)n_input_max * n_committee * i_local];
  const int scale_idx = ele_idx * n_input_max + m;
  const double scale =
      nnp_compute_force_scale(scaling_mode, scmin, scmax, scale_idx,
                              scaling_loc_max, scaling_loc_min, scaling_sigma);
  const double de = dEdG * scale;

  /* Self contribution to atom_i: f[i] -= de * self_partial(d, sf, atom_in_e).
   * r_ii = 0, so the self term makes no virial contribution. */
  const long long self_off = 3 * (sf + (long long)n_symf * atom_in_e);
  const long long f_i_off = 3 * (i_global + (long long)num_atoms_global * c);
  if (threadIdx.x == 0) {
    atomicAdd(&f_dev[f_i_off + 0], -de * self_partial[self_off + 0]);
    atomicAdd(&f_dev[f_i_off + 1], -de * self_partial[self_off + 1]);
    atomicAdd(&f_dev[f_i_off + 2], -de * self_partial[self_off + 2]);
  }

  /* Per-thread virial accumulator vir[a][b], a = force component,
   * b = bond component. Folds onto committee_stress(b, a) via the identity
   * committee_stress(b, a) = -sum r_ij[b] * F_k[a] over neighbours. */
  double vir[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  const double ci_x = accumulate_virial ? coord_dev[3 * i_global + 0] : 0.0;
  const double ci_y = accumulate_virial ? coord_dev[3 * i_global + 1] : 0.0;
  const double ci_z = accumulate_virial ? coord_dev[3 * i_global + 2] : 0.0;

  /* Per-neighbour contribution: f[k_atom] += de * force(d, sf, j, atom_in_e).
   * One global atomicAdd per neighbour per component. */
  const int n_neigh = n_neigh_per_atom[atom_in_e];
  for (int j = threadIdx.x; j < n_neigh; j += blockDim.x) {
    const long long nb_off = j + (long long)max_n_neigh * atom_in_e;
    const int k_atom = rad_ind_one_based[nb_off] - 1;
    const long long force_off = 3 * (sf + (long long)n_symf * nb_off);
    const long long f_k_off = 3 * (k_atom + (long long)num_atoms_global * c);
    const double fk0 = +de * force_rad[force_off + 0];
    const double fk1 = +de * force_rad[force_off + 1];
    const double fk2 = +de * force_rad[force_off + 2];
    atomicAdd(&f_dev[f_k_off + 0], fk0);
    atomicAdd(&f_dev[f_k_off + 1], fk1);
    atomicAdd(&f_dev[f_k_off + 2], fk2);
    if (accumulate_virial) {
      const int img = rad_image_idx[nb_off];
      const double rb0 =
          ci_x - (coord_dev[3 * k_atom + 0] + image_table_dev[3 * img + 0]);
      const double rb1 =
          ci_y - (coord_dev[3 * k_atom + 1] + image_table_dev[3 * img + 1]);
      const double rb2 =
          ci_z - (coord_dev[3 * k_atom + 2] + image_table_dev[3 * img + 2]);
      vir[0][0] += -rb0 * fk0;
      vir[0][1] += -rb1 * fk0;
      vir[0][2] += -rb2 * fk0;
      vir[1][0] += -rb0 * fk1;
      vir[1][1] += -rb1 * fk1;
      vir[1][2] += -rb2 * fk1;
      vir[2][0] += -rb0 * fk2;
      vir[2][1] += -rb1 * fk2;
      vir[2][2] += -rb2 * fk2;
    }
  }

  if (accumulate_virial)
    nnp_virial_warp_flush(vir, virial_dev, c);
}

/*******************************************************************************
 * \brief Angular force-scatter kernel: the same block layout as the radial
 *        kernel but with two neighbour lists (ang1 and ang2). Thread 0 adds
 *        the self contribution; the threads share the ang1 loop (partials jj)
 *        and the ang2 loop (partials kk), each atomic-adding one contribution
 *        per Cartesian component. In the homonuclear case both loops read the
 *        ang1 index list. Sign convention: f[i] -= de * self,
 *        f[k_atom_j] -= de * jj, f[k_atom_k] -= de * kk.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
__global__ void nnp_force_scatter_angular_kernel(
    const int n_atoms_in_e, const int n_symf, const int max_n_ang1,
    const int max_n_kk, const int n_committee, const int n_input_max,
    const int num_atoms_global, const int ele_idx, const int scaling_mode,
    const int homo, const int n_rad_offset, const double scmin,
    const double scmax, const int *__restrict__ symf_map_one_based,
    const int *__restrict__ atom_map, const int *__restrict__ atoms_global,
    const int *__restrict__ n_ang1_per_atom,
    const int *__restrict__ n_ang2_per_atom,
    const int *__restrict__ ang1_ind_one_based,
    const int *__restrict__ ang2_ind_one_based,
    const double *__restrict__ self_partial, const double *__restrict__ jj,
    const double *__restrict__ kk, const double *__restrict__ dE_dG,
    const double *__restrict__ scaling_loc_max,
    const double *__restrict__ scaling_loc_min,
    const double *__restrict__ scaling_sigma, double *__restrict__ f_dev,
    /* Virial (read only when accumulate_virial != 0): ang1_image_idx pairs
     * with ang1_ind (bond r_ij), ang2_image_idx with ang2_ind (bond r_ik);
     * in homo mode the kk side reads ang1_image_idx to match ang1_ind. */
    const int *__restrict__ ang1_image_idx,
    const int *__restrict__ ang2_image_idx,
    const double *__restrict__ coord_dev,
    const double *__restrict__ image_table_dev, double *__restrict__ virial_dev,
    const int accumulate_virial) {
  const int atom_in_e = blockIdx.x;
  const int sf = blockIdx.y;
  const int c = blockIdx.z;
  if (atom_in_e >= n_atoms_in_e || sf >= n_symf || c >= n_committee)
    return;

  const int i_local = atom_map[atom_in_e];
  const int i_global = atoms_global[atom_in_e];
  /* Angular: m = n_rad(ele) + symf(sf), 1-based; convert to 0-based. */
  const int m = n_rad_offset + symf_map_one_based[sf] - 1;

  /* 64-bit offsets on every per-atom-scaled buffer, matching nnp_gpu_pack.h:
   * these tables exceed 2^31 entries on large ranks. */
  const double dEdG = dE_dG[m + n_input_max * c +
                            (long long)n_input_max * n_committee * i_local];
  const int scale_idx = ele_idx * n_input_max + m;
  const double scale =
      nnp_compute_force_scale(scaling_mode, scmin, scmax, scale_idx,
                              scaling_loc_max, scaling_loc_min, scaling_sigma);
  const double de = dEdG * scale;

  /* Self term: r_ii = 0, so no virial contribution. */
  const long long self_off = 3 * (sf + (long long)n_symf * atom_in_e);
  const long long f_i_off = 3 * (i_global + (long long)num_atoms_global * c);
  if (threadIdx.x == 0) {
    atomicAdd(&f_dev[f_i_off + 0], -de * self_partial[self_off + 0]);
    atomicAdd(&f_dev[f_i_off + 1], -de * self_partial[self_off + 1]);
    atomicAdd(&f_dev[f_i_off + 2], -de * self_partial[self_off + 2]);
  }

  /* Per-thread virial accumulator vir[a][b], a = force component,
   * b = bond component; folds onto committee_stress(b, a). Both the jj
   * (bond r_ij) and kk (bond r_ik) sides accumulate into it. */
  double vir[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  const double ci_x = accumulate_virial ? coord_dev[3 * i_global + 0] : 0.0;
  const double ci_y = accumulate_virial ? coord_dev[3 * i_global + 1] : 0.0;
  const double ci_z = accumulate_virial ? coord_dev[3 * i_global + 2] : 0.0;

  const int n_ang1 = n_ang1_per_atom[atom_in_e];
  const int n_ang2 = n_ang2_per_atom[atom_in_e];
  /* One global atomicAdd per neighbour per component. */
  /* ang1 neighbours j: f[k_atom] -= de * jj(d, sf, j, atom). */
  for (int j = threadIdx.x; j < n_ang1; j += blockDim.x) {
    const long long nb_off = j + (long long)max_n_ang1 * atom_in_e;
    const int k_atom = ang1_ind_one_based[nb_off] - 1;
    const long long jj_off = 3 * (sf + (long long)n_symf * nb_off);
    const long long f_k_off = 3 * (k_atom + (long long)num_atoms_global * c);
    const double fk0 = -de * jj[jj_off + 0];
    const double fk1 = -de * jj[jj_off + 1];
    const double fk2 = -de * jj[jj_off + 2];
    atomicAdd(&f_dev[f_k_off + 0], fk0);
    atomicAdd(&f_dev[f_k_off + 1], fk1);
    atomicAdd(&f_dev[f_k_off + 2], fk2);
    if (accumulate_virial) {
      const int img = ang1_image_idx[nb_off];
      const double rb0 =
          ci_x - (coord_dev[3 * k_atom + 0] + image_table_dev[3 * img + 0]);
      const double rb1 =
          ci_y - (coord_dev[3 * k_atom + 1] + image_table_dev[3 * img + 1]);
      const double rb2 =
          ci_z - (coord_dev[3 * k_atom + 2] + image_table_dev[3 * img + 2]);
      vir[0][0] += -rb0 * fk0;
      vir[0][1] += -rb1 * fk0;
      vir[0][2] += -rb2 * fk0;
      vir[1][0] += -rb0 * fk1;
      vir[1][1] += -rb1 * fk1;
      vir[1][2] += -rb2 * fk1;
      vir[2][0] += -rb0 * fk2;
      vir[2][1] += -rb1 * fk2;
      vir[2][2] += -rb2 * fk2;
    }
  }
  /* ang2 neighbours k (or ang1 in homo): f[k_atom] -= de * kk(d, sf, k, atom).
   */
  for (int kk_j = threadIdx.x; kk_j < n_ang2; kk_j += blockDim.x) {
    const int *idx_src = homo ? ang1_ind_one_based : ang2_ind_one_based;
    const long long nb_off = kk_j + (long long)max_n_kk * atom_in_e;
    const int k_atom = idx_src[nb_off] - 1;
    const long long kk_off = 3 * (sf + (long long)n_symf * nb_off);
    const long long f_k_off = 3 * (k_atom + (long long)num_atoms_global * c);
    const double fk0 = -de * kk[kk_off + 0];
    const double fk1 = -de * kk[kk_off + 1];
    const double fk2 = -de * kk[kk_off + 2];
    atomicAdd(&f_dev[f_k_off + 0], fk0);
    atomicAdd(&f_dev[f_k_off + 1], fk1);
    atomicAdd(&f_dev[f_k_off + 2], fk2);
    if (accumulate_virial) {
      const int *img_src = homo ? ang1_image_idx : ang2_image_idx;
      const int img = img_src[nb_off];
      const double rb0 =
          ci_x - (coord_dev[3 * k_atom + 0] + image_table_dev[3 * img + 0]);
      const double rb1 =
          ci_y - (coord_dev[3 * k_atom + 1] + image_table_dev[3 * img + 1]);
      const double rb2 =
          ci_z - (coord_dev[3 * k_atom + 2] + image_table_dev[3 * img + 2]);
      vir[0][0] += -rb0 * fk0;
      vir[0][1] += -rb1 * fk0;
      vir[0][2] += -rb2 * fk0;
      vir[1][0] += -rb0 * fk1;
      vir[1][1] += -rb1 * fk1;
      vir[1][2] += -rb2 * fk1;
      vir[2][0] += -rb0 * fk2;
      vir[2][1] += -rb1 * fk2;
      vir[2][2] += -rb2 * fk2;
    }
  }

  if (accumulate_virial)
    nnp_virial_warp_flush(vir, virial_dev, c);
}

/*******************************************************************************
 * \brief Seed the device force accumulator for a step. Ensures f_dev holds
 *        3 * num_atoms_global * n_committee doubles, then copies the current
 *        host force array onto the device so the scatter kernels accumulate on
 *        top of it, matching the host loop's += semantics.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
extern "C" int nnp_gpu_force_scatter_begin_c(int num_atoms_global,
                                             int n_committee,
                                             const double *h_force) {
  if (num_atoms_global <= 0 || n_committee <= 0) {
    fprintf(stderr,
            "ERROR: nnp_gpu_force_scatter_begin: bad sizes (num_atoms=%d, "
            "n_committee=%d)\n",
            num_atoms_global, n_committee);
    return -1;
  }
  nnp_gpu_state_t *state = nnp_gpu_state_get();
  if (state == NULL)
    return -1;
  const size_t bytes =
      sizeof(double) * 3 * (size_t)num_atoms_global * (size_t)n_committee;
  if (nnp_gpu_pbuf_ensure((void **)&state->f_dev, &state->f_dev_cap, bytes))
    return -1;
  /* Copy the host nnp%myforce so the scatter kernel accumulates on top,
   * matching the host loop's += semantics. */
  offloadMemcpyAsyncHtoD(state->f_dev, h_force, bytes, state->stream);
  /* Zero the per-committee analytical virial accumulator. The scatter kernels
   * atomic-add onto it when virial accumulation is requested;
   * nnp_gpu_force_scatter_virial_d2h_c reads it back. Always maintained (tiny:
   * 9 * n_committee doubles) so the buffer exists whenever a stress-requesting
   * step D2Hs it. */
  const size_t vbytes = sizeof(double) * 9 * (size_t)n_committee;
  if (nnp_gpu_pbuf_ensure((void **)&state->virial_dev, &state->virial_dev_cap,
                          vbytes))
    return -1;
  offloadMemsetAsync(state->virial_dev, 0, vbytes, state->stream);
  return 0;
}

/*******************************************************************************
 * \brief Read the accumulated device force back to the host at the end of a
 *        step. The stream synchronize here is the single per-step
 *        synchronization point: it guarantees every kernel queued this step
 *        has completed before the host reads the result.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
extern "C" int nnp_gpu_force_scatter_end_c(int num_atoms_global,
                                           int n_committee, double *h_force) {
  if (num_atoms_global <= 0 || n_committee <= 0)
    return -1;
  nnp_gpu_state_t *state = nnp_gpu_state_get();
  if (state == NULL || state->f_dev == NULL)
    return -1;
  const size_t bytes =
      sizeof(double) * 3 * (size_t)num_atoms_global * (size_t)n_committee;
  offloadMemcpyAsyncDtoH(h_force, state->f_dev, bytes, state->stream);
  offloadStreamSynchronize(state->stream);
  return 0;
}

/*******************************************************************************
 * \brief D2H the per-committee analytical virial accumulator into a host
 *        buffer laid out [9 * n_committee] (virial[9*c + 3*a + b], a = force
 *        component, b = bond component). Call immediately after
 *        nnp_gpu_force_scatter_end_c: its stream synchronize has already
 *        committed every scatter kernel's atomic adds, so a plain synchronous
 *        D2H of this tiny buffer is sufficient.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
extern "C" int nnp_gpu_force_scatter_virial_d2h_c(int n_committee,
                                                  double *h_virial) {
  if (n_committee <= 0)
    return -1;
  nnp_gpu_state_t *state = nnp_gpu_state_get();
  if (state == NULL || state->virial_dev == NULL)
    return -1;
  const size_t vbytes = sizeof(double) * 9 * (size_t)n_committee;
  offloadMemcpyDtoH(h_virial, state->virial_dev, vbytes);
  return 0;
}

/*******************************************************************************
 * \brief Radial force-scatter launcher for one (group, element) slot,
 *        companion to the radial symmetry-function launcher. Rebuilds the
 *        per-neighbour index and per-atom mapping tables (the shared index
 *        buffers were overwritten by later descriptor groups, so the
 *        rearrange is repeated here) and the 1-based symmetry-function map,
 *        then launches the radial scatter kernel to accumulate this slot's
 *        force contribution into f_dev.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
extern "C" int nnp_gpu_force_scatter_radial_c(
    const int n_atoms_in_e, const int n_symf, const int max_n_neigh_per_se,
    const int n_committee, const int n_input_max, const int num_atoms_global,
    const int ele_idx, const int scaling_mode, const double scmin,
    const double scmax, const int s_0based, const int max_n_radgrp_walk,
    const int n_atoms_in_e_max_walk, const int max_n_neigh_walk,
    const int istart_1based, const int se_idx, const int accumulate_virial) {
  if (n_atoms_in_e <= 0 || n_symf <= 0 || max_n_neigh_per_se <= 0 ||
      n_committee <= 0 || n_input_max <= 0 || num_atoms_global <= 0)
    return -1;
  if (se_idx < 0 || se_idx >= NNP_MAX_SE) {
    fprintf(stderr,
            "ERROR: nnp_gpu_force_scatter_radial: se_idx=%d out of range\n",
            se_idx);
    return -1;
  }
  nnp_gpu_state_t *state = nnp_gpu_state_get();
  if (state == NULL || state->f_dev == NULL)
    return -1;
  if (!state->scaling_uploaded)
    return -1;
  if (state->per_se_rad_self[se_idx] == NULL ||
      state->per_se_rad_force[se_idx] == NULL) {
    fprintf(stderr, "ERROR: nnp_gpu_force_scatter_radial: per_se buffers "
                    "not populated (call the symmetry-function stage first)\n");
    return -1;
  }
  if (state->walk_n_rad_pa_dev == NULL || state->walk_rad_ind_dev == NULL ||
      state->walk_atoms_global_dev == NULL || !state->se_static_uploaded) {
    fprintf(stderr, "ERROR: nnp_gpu_force_scatter_radial: walk / "
                    "se_static prerequisites not met\n");
    return -1;
  }
  /* The descriptor gradient is the network stage's output and the one input
   * this launcher does not produce itself. */
  if (state->pbuf_dEdG_out == NULL) {
    fprintf(stderr, "ERROR: nnp_gpu_force_scatter_radial: descriptor "
                    "gradient not on the device (run the network stage "
                    "first)\n");
    return -1;
  }

  const size_t neigh_count = (size_t)max_n_neigh_per_se * (size_t)n_atoms_in_e;

  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_atom_map,
                          &state->pbuf_atom_map_cap,
                          sizeof(int) * (size_t)n_atoms_in_e))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_atoms_global_rad,
                          &state->pbuf_atoms_global_rad_cap,
                          sizeof(int) * (size_t)n_atoms_in_e))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_symf_map,
                          &state->pbuf_symf_map_cap,
                          sizeof(int) * (size_t)n_symf))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_n_neigh_pa,
                          &state->pbuf_n_neigh_pa_cap,
                          sizeof(int) * (size_t)n_atoms_in_e))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->scatter_neigh_ind,
                          &state->scatter_neigh_ind_cap,
                          sizeof(int) * neigh_count))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_izeta, &state->pbuf_izeta_cap,
                          sizeof(int) * (size_t)n_symf))
    return -1;
  /* pbuf_rad_ind / pbuf_rad_image_idx are written too but unused by
   * scatter; keep their allocations alive in case the symmetry-function stage
   * has not run. */
  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_rad_ind,
                          &state->pbuf_rad_ind_cap, sizeof(int) * neigh_count))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_rad_image_idx,
                          &state->pbuf_rad_image_idx_cap,
                          sizeof(int) * neigh_count))
    return -1;

  {
    const int block = 32;
    const int j_blocks = (max_n_neigh_per_se + block - 1) / block;
    dim3 grid(n_atoms_in_e, j_blocks > 0 ? j_blocks : 1, 1);
    dim3 blk(block, 1, 1);
    nnp_pack_walk_to_pbuf_radial_kernel<<<grid, blk, 0, state->stream>>>(
        n_atoms_in_e, max_n_neigh_per_se, s_0based, ele_idx, max_n_radgrp_walk,
        n_atoms_in_e_max_walk, max_n_neigh_walk, istart_1based,
        state->walk_n_rad_pa_dev, state->walk_rad_ind_dev,
        state->walk_rad_image_idx_dev, state->walk_atoms_global_dev,
        state->pbuf_n_neigh_pa, state->pbuf_rad_ind, state->pbuf_rad_image_idx,
        state->pbuf_atoms_global_rad, state->pbuf_atom_map,
        state->scatter_neigh_ind);
    OFFLOAD_CHECK(offloadPeekAtLastError());
  }
  {
    const int block = 32;
    const int sf_blocks = (n_symf + block - 1) / block;
    dim3 grid(sf_blocks > 0 ? sf_blocks : 1, 1, 1);
    dim3 blk(block, 1, 1);
    const int se_idx_in_kind = s_0based + max_n_radgrp_walk * ele_idx;
    /* Radial scatter: pbuf_izeta = 1-based; pbuf_symf_map = 0-based scratch
     * (kernel reads pbuf_izeta and applies n_rad_offset itself). m_offset=0. */
    nnp_pack_se_static_symf_map_kernel<<<grid, blk, 0, state->stream>>>(
        se_idx_in_kind, n_symf, /*m_offset=*/0,
        state->se_static_rad_symf_packed_dev,
        state->se_static_rad_symf_offsets_dev, state->pbuf_symf_map,
        state->pbuf_izeta);
    OFFLOAD_CHECK(offloadPeekAtLastError());
  }

  {
    dim3 grid(n_atoms_in_e, n_symf, n_committee);
    dim3 blk(32, 1, 1);
    /* Scatter kernel takes 1-based symf_map (pbuf_izeta) and 1-based
     * neigh_ind (scatter_neigh_ind). */
    nnp_force_scatter_radial_kernel<<<grid, blk, 0, state->stream>>>(
        n_atoms_in_e, n_symf, max_n_neigh_per_se, n_committee, n_input_max,
        num_atoms_global, ele_idx, scaling_mode, /*n_rad_offset=*/0, scmin,
        scmax, state->pbuf_izeta, state->pbuf_atom_map,
        state->pbuf_atoms_global_rad, state->pbuf_n_neigh_pa,
        state->scatter_neigh_ind, state->per_se_rad_self[se_idx],
        state->per_se_rad_force[se_idx], state->pbuf_dEdG_out,
        state->scaling_loc_max_dev, state->scaling_loc_min_dev,
        state->scaling_sigma_dev, state->f_dev,
        /* Virial: pbuf_rad_image_idx pairs 1:1 (same out_off) with the
         * scatter_neigh_ind neighbour list the kernel reads. */
        state->pbuf_rad_image_idx, state->coord_dev, state->image_table_dev,
        state->virial_dev, accumulate_virial);
  }
  OFFLOAD_CHECK(offloadPeekAtLastError());
  return 0;
}

/*******************************************************************************
 * \brief Angular force-scatter launcher for one (group, element) slot,
 *        companion to the angular symmetry-function launcher. Rebuilds the
 *        per-atom counts, the two 1-based angular neighbour index lists and
 *        the 1-based symmetry-function map (the shared index buffers were
 *        overwritten by later descriptor groups, so the rearrange is repeated
 *        here), then launches the angular scatter kernel to accumulate this
 *        slot's force contribution into f_dev.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
extern "C" int nnp_gpu_force_scatter_angular_c(
    const int n_atoms_in_e, const int n_symf, const int max_n_ang1_per_se,
    const int max_n_kk_per_se, const int n_committee, const int n_input_max,
    const int num_atoms_global, const int ele_idx, const int n_rad_offset,
    const int scaling_mode, const int homo, const double scmin,
    const double scmax, const int s_0based, const int max_n_anggrp_walk,
    const int n_atoms_in_e_max_walk, const int max_n_neigh_ang1_walk,
    const int max_n_neigh_ang2_walk, const int istart_1based, const int se_idx,
    const int accumulate_virial) {
  if (n_atoms_in_e <= 0 || n_symf <= 0 || max_n_ang1_per_se <= 0 ||
      n_committee <= 0 || n_input_max <= 0 || num_atoms_global <= 0)
    return -1;
  if (se_idx < 0 || se_idx >= NNP_MAX_SE) {
    fprintf(stderr,
            "ERROR: nnp_gpu_force_scatter_angular: se_idx=%d out of range\n",
            se_idx);
    return -1;
  }
  nnp_gpu_state_t *state = nnp_gpu_state_get();
  if (state == NULL || state->f_dev == NULL)
    return -1;
  if (!state->scaling_uploaded)
    return -1;
  if (state->per_se_ang_self[se_idx] == NULL ||
      state->per_se_ang_jj[se_idx] == NULL ||
      state->per_se_ang_kk[se_idx] == NULL) {
    fprintf(stderr,
            "ERROR: nnp_gpu_force_scatter_angular: per_se_ang "
            "buffers not populated (call the symmetry-function stage first)\n");
    return -1;
  }
  if (state->walk_n_ang1_pa_dev == NULL || state->walk_ang1_ind_dev == NULL ||
      state->walk_atoms_global_dev == NULL || !state->se_static_uploaded) {
    fprintf(stderr, "ERROR: nnp_gpu_force_scatter_angular: walk / "
                    "se_static prerequisites not met\n");
    return -1;
  }
  /* The descriptor gradient is the network stage's output and the one input
   * this launcher does not produce itself. */
  if (state->pbuf_dEdG_out == NULL) {
    fprintf(stderr, "ERROR: nnp_gpu_force_scatter_angular: descriptor "
                    "gradient not on the device (run the network stage "
                    "first)\n");
    return -1;
  }
  /* Mirrors the bound the SF launcher enforces: the scatter reads the k-side
   * index list with the max_n_kk stride, so it must not exceed the width the
   * pack kernel wrote. */
  if (homo && max_n_kk_per_se != max_n_ang1_per_se) {
    fprintf(stderr,
            "ERROR: nnp_gpu_force_scatter_angular: homonuclear slot needs "
            "max_n_kk_per_se == max_n_ang1_per_se (%d vs %d); the k-side list "
            "is the ang1 list and is packed with that stride\n",
            max_n_kk_per_se, max_n_ang1_per_se);
    return -1;
  }

  const size_t ang1_count = (size_t)max_n_ang1_per_se * (size_t)n_atoms_in_e;
  const size_t ang2_count = (size_t)max_n_kk_per_se * (size_t)n_atoms_in_e;

  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_atom_map,
                          &state->pbuf_atom_map_cap,
                          sizeof(int) * (size_t)n_atoms_in_e))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_atoms_global_ang,
                          &state->pbuf_atoms_global_ang_cap,
                          sizeof(int) * (size_t)n_atoms_in_e))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_symf_map,
                          &state->pbuf_symf_map_cap,
                          sizeof(int) * (size_t)n_symf))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_n_ang1_pa,
                          &state->pbuf_n_ang1_pa_cap,
                          sizeof(int) * (size_t)n_atoms_in_e))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_n_ang2_pa,
                          &state->pbuf_n_ang2_pa_cap,
                          sizeof(int) * (size_t)n_atoms_in_e))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_ang1_ind,
                          &state->pbuf_ang1_ind_cap, sizeof(int) * ang1_count))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_ang1_image_idx,
                          &state->pbuf_ang1_image_idx_cap,
                          sizeof(int) * ang1_count))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_ang2_ind,
                          &state->pbuf_ang2_ind_cap,
                          sizeof(int) * (ang2_count > 0 ? ang2_count : 1)))
    return -1;
  if (nnp_gpu_pbuf_ensure((void **)&state->pbuf_ang2_image_idx,
                          &state->pbuf_ang2_image_idx_cap,
                          sizeof(int) * (ang2_count > 0 ? ang2_count : 1)))
    return -1;
  /* scatter scratch for 0-based symf_map (unused by scatter kernel; the
   * rearrange writes both 0-based and 1-based, we use only 1-based). */
  if (nnp_gpu_pbuf_ensure((void **)&state->scatter_neigh_ind,
                          &state->scatter_neigh_ind_cap,
                          sizeof(int) * (size_t)n_symf))
    return -1;

  /* Per-atom rearrange (writes pbuf_n_ang*, atom_map, atoms_global). */
  {
    const int block = 32;
    dim3 grid((n_atoms_in_e + block - 1) / block, 1, 1);
    dim3 blk(block, 1, 1);
    nnp_pack_walk_to_pbuf_angular_per_atom_kernel<<<grid, blk, 0,
                                                    state->stream>>>(
        n_atoms_in_e, homo, s_0based, ele_idx, max_n_anggrp_walk,
        n_atoms_in_e_max_walk, istart_1based, state->walk_n_ang1_pa_dev,
        state->walk_n_ang2_pa_dev, state->walk_atoms_global_dev,
        state->pbuf_n_ang1_pa, state->pbuf_n_ang2_pa,
        state->pbuf_atoms_global_ang, state->pbuf_atom_map);
    OFFLOAD_CHECK(offloadPeekAtLastError());
  }

  /* Per-(atom, j) rearrange with ind_offset = 1 (1-based for scatter). */
  {
    const int block = 32;
    const int j_blocks = (max_n_ang1_per_se + block - 1) / block;
    dim3 grid(n_atoms_in_e, j_blocks > 0 ? j_blocks : 1, 1);
    dim3 blk(block, 1, 1);
    nnp_pack_walk_to_pbuf_angular_per_neigh_kernel<<<grid, blk, 0,
                                                     state->stream>>>(
        n_atoms_in_e, max_n_ang1_per_se, s_0based, ele_idx, max_n_anggrp_walk,
        n_atoms_in_e_max_walk, max_n_neigh_ang1_walk,
        /*ind_offset=*/1, state->walk_n_ang1_pa_dev, state->walk_ang1_ind_dev,
        state->walk_ang1_image_idx_dev, state->pbuf_ang1_ind,
        state->pbuf_ang1_image_idx);
    OFFLOAD_CHECK(offloadPeekAtLastError());
  }
  if (!homo && max_n_kk_per_se > 0) {
    const int block = 32;
    const int j_blocks = (max_n_kk_per_se + block - 1) / block;
    dim3 grid(n_atoms_in_e, j_blocks > 0 ? j_blocks : 1, 1);
    dim3 blk(block, 1, 1);
    nnp_pack_walk_to_pbuf_angular_per_neigh_kernel<<<grid, blk, 0,
                                                     state->stream>>>(
        n_atoms_in_e, max_n_kk_per_se, s_0based, ele_idx, max_n_anggrp_walk,
        n_atoms_in_e_max_walk, max_n_neigh_ang2_walk,
        /*ind_offset=*/1, state->walk_n_ang2_pa_dev, state->walk_ang2_ind_dev,
        state->walk_ang2_image_idx_dev, state->pbuf_ang2_ind,
        state->pbuf_ang2_image_idx);
    OFFLOAD_CHECK(offloadPeekAtLastError());
  }

  /* Materialise the symmetry-function map: scatter wants 1-based. The kernel
   * writes both (0-based to first ptr, 1-based to second); we put 1-based
   * directly into pbuf_symf_map this time and discard 0-based to scratch. */
  {
    const int block = 32;
    const int sf_blocks = (n_symf + block - 1) / block;
    dim3 grid(sf_blocks > 0 ? sf_blocks : 1, 1, 1);
    dim3 blk(block, 1, 1);
    const int se_idx_in_kind = s_0based + max_n_anggrp_walk * ele_idx;
    /* Angular scatter: scatter_neigh_ind = 0-based scratch (discarded);
     * pbuf_symf_map = 1-based (consumed by scatter kernel which applies
     * n_rad_offset itself). m_offset=0 only affects the discarded scratch. */
    nnp_pack_se_static_symf_map_kernel<<<grid, blk, 0, state->stream>>>(
        se_idx_in_kind, n_symf, /*m_offset=*/0,
        state->se_static_ang_symf_packed_dev,
        state->se_static_ang_symf_offsets_dev, state->scatter_neigh_ind,
        state->pbuf_symf_map);
    OFFLOAD_CHECK(offloadPeekAtLastError());
  }

  {
    dim3 grid(n_atoms_in_e, n_symf, n_committee);
    dim3 blk(32, 1, 1);
    nnp_force_scatter_angular_kernel<<<grid, blk, 0, state->stream>>>(
        n_atoms_in_e, n_symf, max_n_ang1_per_se, max_n_kk_per_se, n_committee,
        n_input_max, num_atoms_global, ele_idx, scaling_mode, homo,
        n_rad_offset, scmin, scmax, state->pbuf_symf_map, state->pbuf_atom_map,
        state->pbuf_atoms_global_ang, state->pbuf_n_ang1_pa,
        state->pbuf_n_ang2_pa, state->pbuf_ang1_ind, state->pbuf_ang2_ind,
        state->per_se_ang_self[se_idx], state->per_se_ang_jj[se_idx],
        state->per_se_ang_kk[se_idx], state->pbuf_dEdG_out,
        state->scaling_loc_max_dev, state->scaling_loc_min_dev,
        state->scaling_sigma_dev, state->f_dev,
        /* Virial: pbuf_ang1/2_image_idx pair 1:1 (same out_off) with the
         * pbuf_ang1/2_ind neighbour lists the kernel reads; in homo mode the
         * kk side reads pbuf_ang1_image_idx to match its ang1_ind read. */
        state->pbuf_ang1_image_idx, state->pbuf_ang2_image_idx,
        state->coord_dev, state->image_table_dev, state->virial_dev,
        accumulate_virial);
  }
  OFFLOAD_CHECK(offloadPeekAtLastError());
  return 0;
}

#endif /* defined(__OFFLOAD) && !defined(__NO_OFFLOAD_NNP) */
