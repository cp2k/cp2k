/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#ifndef NNP_GPU_PACK_H
#define NNP_GPU_PACK_H

#include "../../offload/offload_runtime.h"
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_NNP)

#include "nnp_gpu_state.h"

/*
 * Kernels that rearrange the on-device cell-list walk output into the
 * per-(group, element) buffers the descriptor and force-assembly kernels
 * read. Both nnp_gpu_acsf.cu and nnp_gpu_scatter.cu launch them: the
 * descriptor pass fills the buffers, and the force pass rebuilds them
 * because the descriptor groups that follow overwrite the shared index
 * scratch. A kernel cannot be launched from a translation unit that does
 * not define it unless the build enables relocatable device code, so these
 * are static and each unit compiles its own copy.
 */

/*******************************************************************************
 * \brief Rearrange one (group, element) slice of the on-device walk output
 *        into the per-call pbuf_* layouts the radial symmetry-function and
 *        scatter kernels consume, removing both the device-to-host readback of
 *        the walk buffers and the host packing loop. The grid is sized so each
 *        thread handles one (atom, neighbour) pair; the j == 0 lane of each
 *        atom writes the per-atom scalar outputs.
 *
 *        Index bases. The walk stores rad_ind / rad_image_idx 0-BASED (the
 *        same convention the host uses post-conversion, where grp%rad_ind =
 *        walk + 1 is 1-based). walk_atoms_global is 1-based. The radial
 *        symmetry-function kernel expects 0-based atom and image indices, so
 *        pbuf_rad_ind / pbuf_rad_image_idx are copied through unchanged. The
 *        scatter kernel expects 1-based atom indices, so
 *        scatter_neigh_ind = walk + 1. Padding slots (j >= n_neigh) get zeros
 *        (and 0 + 1 = 1 in the 1-based buffer), harmless because downstream
 *        kernels guard on the per-atom neighbour count.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
static __global__ void nnp_pack_walk_to_pbuf_radial_kernel(
    const int n_atoms_in_e, const int max_n_neigh_per_se, const int s_0based,
    const int e_0based, const int max_n_radgrp_walk,
    const int n_atoms_in_e_max_walk, const int max_n_neigh_walk,
    const int istart_1based, const int *__restrict__ walk_n_rad_pa,
    const int *__restrict__ walk_rad_ind,
    const int *__restrict__ walk_rad_image_idx,
    const int *__restrict__ walk_atoms_global,
    int *__restrict__ pbuf_n_neigh_pa,
    int *__restrict__ pbuf_rad_ind_zero_based,
    int *__restrict__ pbuf_rad_image_idx_zero_based,
    int *__restrict__ pbuf_atoms_global_zero_based,
    int *__restrict__ pbuf_atom_map,
    int *__restrict__ scatter_neigh_ind_one_based) {
  const int atom = blockIdx.x;
  const int j = blockIdx.y * blockDim.x + threadIdx.x;
  if (atom >= n_atoms_in_e)
    return;

  /* Walk-side strides (column-major, matching state field comments). */
  const long long count_off =
      (long long)atom +
      (long long)n_atoms_in_e_max_walk *
          ((long long)s_0based +
           (long long)max_n_radgrp_walk * (long long)e_0based);
  const long long ij_base =
      (long long)max_n_neigh_walk *
      ((long long)atom +
       (long long)n_atoms_in_e_max_walk *
           ((long long)s_0based +
            (long long)max_n_radgrp_walk * (long long)e_0based));

  const int n_neigh = walk_n_rad_pa[count_off];

  /* Per-atom scalar outputs: one thread (j == 0) per atom row. */
  if (j == 0) {
    pbuf_n_neigh_pa[atom] = n_neigh;
    const int atom_global_1b =
        walk_atoms_global[atom + n_atoms_in_e_max_walk * e_0based];
    pbuf_atoms_global_zero_based[atom] = atom_global_1b - 1;
    pbuf_atom_map[atom] = atom_global_1b - istart_1based;
  }

  /* Per-(atom, j) outputs. Out-of-range slots get pad zeros so that
   * downstream kernels (which guard on n_neigh) never read garbage. */
  if (j < max_n_neigh_per_se) {
    int j_global_0b = 0;
    int img_0b = 0;
    if (j < n_neigh) {
      j_global_0b = walk_rad_ind[ij_base + j];
      img_0b = walk_rad_image_idx[ij_base + j];
    }
    const long long out_off =
        (long long)atom * (long long)max_n_neigh_per_se + (long long)j;
    /* Walk is already 0-based, so copy through directly. */
    pbuf_rad_ind_zero_based[out_off] = j_global_0b;
    pbuf_rad_image_idx_zero_based[out_off] = img_0b;
    /* Scatter wants 1-based (walk + 1). Padding (j >= n_neigh) stays 0 + 1 =
     * 1, harmless because scatter guards on n_neigh_per_atom. */
    scatter_neigh_ind_one_based[out_off] = j_global_0b + 1;
  }
}

/*******************************************************************************
 * \brief Materialise the per-(group, element) symmetry-function map into both
 *        0-based (for the symmetry-function kernel) and 1-based (for the
 *        scatter kernel) layouts in one launch, reading the once-uploaded
 *        packed table (se_symf_packed) and its prefix-sum offsets
 *        (se_symf_offsets). The packed table stores raw 1-based symf values
 *        with no descriptor-row offset. Shared by the radial and angular
 *        launchers.
 *
 *        m_offset is added to the 0-based output only. For radial groups it is
 *        0 (the descriptor row is just symf - 1). For angular groups it is
 *        n_rad(element), because angular descriptors live after all radial
 *        descriptors in an atom's G_dev row and the angular symmetry-function
 *        kernel takes no separate offset parameter, so the offset must be baked
 *        into the map. The 1-based output is never offset: both scatter kernels
 *        read it as the raw symf value and apply n_rad themselves.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
static __global__ void
nnp_pack_se_static_symf_map_kernel(const int se_idx_in_kind, const int n_symf,
                                   const int m_offset,
                                   const int *__restrict__ se_symf_packed,
                                   const int *__restrict__ se_symf_offsets,
                                   int *__restrict__ symf_map_zero_based,
                                   int *__restrict__ symf_map_one_based) {
  const int sf = blockIdx.x * blockDim.x + threadIdx.x;
  if (sf >= n_symf)
    return;
  const int off = se_symf_offsets[se_idx_in_kind];
  const int v_one = se_symf_packed[off + sf];
  symf_map_zero_based[sf] = v_one - 1 + m_offset;
  symf_map_one_based[sf] = v_one;
}

/*******************************************************************************
 * \brief Rearrange the per-atom scalar outputs (neighbour counts, 0-based
 *        center index, atom-slot map) for one angular (group, element) slice
 *        of the walk output. One block per atom row; thread 0 writes
 *        everything for that atom. homo == 1 sets pbuf_n_ang2_pa =
 *        pbuf_n_ang1_pa. walk_atoms_global is 1-based; the 0-based center index
 *        and the atom_map (atom_global - istart) are derived from it.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
static __global__ void nnp_pack_walk_to_pbuf_angular_per_atom_kernel(
    const int n_atoms_in_e, const int homo, const int s_0based,
    const int e_0based, const int max_n_anggrp_walk,
    const int n_atoms_in_e_max_walk, const int istart_1based,
    const int *__restrict__ walk_n_ang1_pa,
    const int *__restrict__ walk_n_ang2_pa,
    const int *__restrict__ walk_atoms_global, int *__restrict__ pbuf_n_ang1_pa,
    int *__restrict__ pbuf_n_ang2_pa,
    int *__restrict__ pbuf_atoms_global_zero_based,
    int *__restrict__ pbuf_atom_map) {
  const int atom = blockIdx.x * blockDim.x + threadIdx.x;
  if (atom >= n_atoms_in_e)
    return;

  const long long count_off =
      (long long)atom +
      (long long)n_atoms_in_e_max_walk *
          ((long long)s_0based +
           (long long)max_n_anggrp_walk * (long long)e_0based);

  const int n1 = walk_n_ang1_pa[count_off];
  const int n2 = homo ? n1 : walk_n_ang2_pa[count_off];
  pbuf_n_ang1_pa[atom] = n1;
  pbuf_n_ang2_pa[atom] = n2;

  const int atom_global_1b =
      walk_atoms_global[atom + n_atoms_in_e_max_walk * e_0based];
  pbuf_atoms_global_zero_based[atom] = atom_global_1b - 1;
  pbuf_atom_map[atom] = atom_global_1b - istart_1based;
}

/*******************************************************************************
 * \brief Rearrange one (group, element) slice of a single angular neighbour
 *        list (ang1 or ang2, selected by the device pointers and walk strides
 *        the caller passes) from the walk output into the per-call pbuf layout
 *        [max_n_per_se x n_atoms_in_e]. ind_offset is added to the walk index:
 *        0 for the symmetry-function kernel (0-based), 1 for the scatter kernel
 *        (1-based). image_idx is always 0-based (only the symmetry-function
 *        kernel reads it). Padding slots (j >= n_neigh) write ind = ind_offset
 *        and image = 0, inert because downstream kernels guard on the per-atom
 *        neighbour count.
 * \author Dhruv Sharma (ds2173@cam.ac.uk)
 ******************************************************************************/
static __global__ void nnp_pack_walk_to_pbuf_angular_per_neigh_kernel(
    const int n_atoms_in_e, const int max_n_per_se, const int s_0based,
    const int e_0based, const int max_n_anggrp_walk,
    const int n_atoms_in_e_max_walk, const int max_n_neigh_walk,
    const int ind_offset, const int *__restrict__ walk_n_pa,
    const int *__restrict__ walk_ind, const int *__restrict__ walk_image_idx,
    int *__restrict__ pbuf_ind, int *__restrict__ pbuf_image_idx) {
  const int atom = blockIdx.x;
  const int j = blockIdx.y * blockDim.x + threadIdx.x;
  if (atom >= n_atoms_in_e)
    return;
  if (j >= max_n_per_se)
    return;

  const long long count_off =
      (long long)atom +
      (long long)n_atoms_in_e_max_walk *
          ((long long)s_0based +
           (long long)max_n_anggrp_walk * (long long)e_0based);
  const long long ij_base =
      (long long)max_n_neigh_walk *
      ((long long)atom +
       (long long)n_atoms_in_e_max_walk *
           ((long long)s_0based +
            (long long)max_n_anggrp_walk * (long long)e_0based));

  const int n_neigh = walk_n_pa[count_off];
  int j_global_0b = 0;
  int img_0b = 0;
  if (j < n_neigh) {
    j_global_0b = walk_ind[ij_base + j];
    img_0b = walk_image_idx[ij_base + j];
  }
  const long long out_off =
      (long long)atom * (long long)max_n_per_se + (long long)j;
  pbuf_ind[out_off] = j_global_0b + ind_offset;
  pbuf_image_idx[out_off] = img_0b;
}

#endif /* defined(__OFFLOAD) && !defined(__NO_OFFLOAD_NNP) */
#endif /* NNP_GPU_PACK_H */
