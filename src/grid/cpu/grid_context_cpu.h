/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#ifndef GRID_CONTEXT_CPU_H
#define GRID_CONTEXT_CPU_H

#include "../../offload/offload_buffer.h"
#include "../common/grid_basis_set.h"

void *create_grid_context_cpu(
    const bool orthorhombic, const int ntasks, const int nlevels,
    const int natoms, const int nkinds, const int nblocks,
    const int *block_offsets, const double atom_positions[natoms][3],
    const int *const atom_kinds, const grid_basis_set **const basis_sets,
    const int *const level_list, const int *const iatom_list,
    const int *jatom_list, const int *const iset_list,
    const int *const jset_list, const int *const ipgf_list,
    const int *const jpgf_list, const int *const border_mask_list,
    const int *block_num_list, const double *const radius_list,
    const double rab_list[ntasks][3], const int npts_global[nlevels][3],
    const int npts_local[nlevels][3], const int shift_local[nlevels][3],
    const int border_width[nlevels][3], const double dh[nlevels][3][3],
    const double dh_inv[nlevels][3][3]);

void update_grid_context_cpu(
    const bool orthorhombic, const int ntasks, const int nlevels,
    const int natoms, const int nkinds, const int nblocks,
    const int *block_offsets, const double atom_positions[natoms][3],
    const int *const atom_kinds, const grid_basis_set **const basis_sets,
    const int *const level_list, const int *const iatom_list,
    const int *jatom_list, const int *const iset_list,
    const int *const jset_list, const int *const ipgf_list,
    const int *const jpgf_list, const int *const border_mask_list,
    const int *block_num_list, const double *const radius_list,
    const double rab_list[ntasks][3], const int npts_global[nlevels][3],
    const int npts_local[nlevels][3], const int shift_local[nlevels][3],
    const int border_width[nlevels][3], const double dh[nlevels][3][3],
    const double dh_inv[nlevels][3][3], void *ptr);

void initialize_grid_context_on_gpu(void *ptr, const int number_of_devices,
                                    const int *device_id);

void destroy_grid_context_cpu(void *ptr);

void apply_cutoff(void *ptr);

void update_queue_length(void *const ptr, const int queue_length);

void grid_collocate_task_list_cpu(
    void *const ptr, const bool orthorhombic, const int func, const int nlevels,
    const int npts_global[nlevels][3], const int npts_local[nlevels][3],
    const int shift_local[nlevels][3], const int border_width[nlevels][3],
    const double dh[nlevels][3][3], const double dh_inv[nlevels][3][3],
    const offload_buffer *pab_blocks, double *grid[nlevels]);
#endif
