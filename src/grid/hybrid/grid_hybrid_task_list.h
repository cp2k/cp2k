/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2021 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#ifndef GRID_HYBRID_TASK_LIST_H
#define GRID_HYBRID_TASK_LIST_H

#include "../common/grid_basis_set.h"
#include "../common/grid_buffer.h"
#include "../common/grid_constants.h"
#include <stdbool.h>

/*******************************************************************************
 * \brief Internal representation of a task list.
 ******************************************************************************/
typedef struct grid_context_ grid_hybrid_task_list;

/*******************************************************************************
 * \brief Allocates a task list for the hybrid backend.
 *        See grid_task_list.h for details.
 ******************************************************************************/
void grid_hybrid_create_task_list(
    const int ntasks, const int nlevels, const int natoms, const int nkinds,
    const int nblocks, const int block_offsets[nblocks],
    const double atom_positions[natoms][3], const int atom_kinds[natoms],
    const grid_basis_set *basis_sets[nkinds], const int level_list[ntasks],
    const int iatom_list[ntasks], const int jatom_list[ntasks],
    const int iset_list[ntasks], const int jset_list[ntasks],
    const int ipgf_list[ntasks], const int jpgf_list[ntasks],
    const int border_mask_list[ntasks], const int block_num_list[ntasks],
    const double radius_list[ntasks], const double rab_list[ntasks][3],
    grid_hybrid_task_list **task_list);

/*******************************************************************************
 * \brief Deallocates given task list, basis_sets have to be freed separately.
 ******************************************************************************/
void grid_hybrid_free_task_list(grid_hybrid_task_list *ptr);

/*******************************************************************************
 * \brief Collocate all tasks of in given list onto given grids.
 *        See grid_task_list.h for details.
 ******************************************************************************/
void grid_hybrid_collocate_task_list(
    grid_hybrid_task_list *const ptr, const bool orthorhombic,
    const enum grid_func func, const int nlevels,
    const int npts_global[nlevels][3], const int npts_local[nlevels][3],
    const int shift_local[nlevels][3], const int border_width[nlevels][3],
    const double dh[nlevels][3][3], const double dh_inv[nlevels][3][3],
    const grid_buffer *pab_blocks, double *grid[nlevels]);

void grid_hybrid_integrate_task_list(
    void *const ptr, const bool orthorhombic, const bool calculate_forces,
    const int natoms, const int nlevels, const int npts_global[nlevels][3],
    const int npts_local[nlevels][3], const int shift_local[nlevels][3],
    const int border_width[nlevels][3], const double dh[nlevels][3][3],
    const double dh_inv[nlevels][3][3], const grid_buffer *const pab_blocks,
    const double *const grid[nlevels], grid_buffer *hab_blocks,
    double forces[natoms][3], double virial[3][3]);
#endif
