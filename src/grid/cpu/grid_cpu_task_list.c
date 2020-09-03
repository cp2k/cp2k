/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../ref/grid_ref_task_list.h"
#include "grid_cpu_task_list.h"

/*******************************************************************************
 * \brief Allocates a task list which can be passed to grid_collocate_task_list.
 *        See grid_task_list.h for details.
 ******************************************************************************/
void grid_cpu_create_task_list(
    const int ntasks, const int nlevels, const int natoms, const int nkinds,
    const int nblocks, const int buffer_size, const int block_offsets[nblocks],
    const double atom_positions[natoms][3], const int atom_kinds[natoms],
    const grid_basis_set *basis_sets[nkinds], const int level_list[ntasks],
    const int iatom_list[ntasks], const int jatom_list[ntasks],
    const int iset_list[ntasks], const int jset_list[ntasks],
    const int ipgf_list[ntasks], const int jpgf_list[ntasks],
    const int border_mask_list[ntasks], const int block_num_list[ntasks],
    const double radius_list[ntasks], const double rab_list[ntasks][3],
    double **blocks_buffer, grid_cpu_task_list **task_list) {
  //
  // TODO: Replace with actual cpu backend implementation.
  //
  if (*task_list == NULL) {
    *task_list = malloc(sizeof(grid_cpu_task_list));
    (*task_list)->crutch = NULL;
  }
  grid_ref_create_task_list(
      ntasks, nlevels, natoms, nkinds, nblocks, buffer_size, block_offsets,
      atom_positions, atom_kinds, basis_sets, level_list, iatom_list,
      jatom_list, iset_list, jset_list, ipgf_list, jpgf_list, border_mask_list,
      block_num_list, radius_list, rab_list, blocks_buffer,
      &(*task_list)->crutch);
}

/*******************************************************************************
 * \brief Deallocates given task list, basis_sets have to be freed separately.
 ******************************************************************************/
void grid_cpu_free_task_list(grid_cpu_task_list *task_list) {
  //
  // TODO: Replace with actual cpu backend implementation.
  //
  grid_ref_free_task_list(task_list->crutch);
  free(task_list);
}

/*******************************************************************************
 * \brief Collocate all tasks of in given list onto given grids.
 *        See grid_task_list.h for details.
 ******************************************************************************/
void grid_cpu_collocate_task_list(
    const grid_cpu_task_list *task_list, const bool orthorhombic,
    const int func, const int nlevels, const int npts_global[nlevels][3],
    const int npts_local[nlevels][3], const int shift_local[nlevels][3],
    const int border_width[nlevels][3], const double dh[nlevels][3][3],
    const double dh_inv[nlevels][3][3], double *grid[nlevels]) {
  //
  // TODO: Replace with actual cpu backend implementation.
  //
  grid_ref_collocate_task_list(task_list->crutch, orthorhombic, func, nlevels,
                               npts_global, npts_local, shift_local,
                               border_width, dh, dh_inv, grid);
}

// EOF
