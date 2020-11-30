/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "../common/grid_library.h"
#include "grid_context_cpu.h"
#include "grid_cpu_task_list.h"

/*******************************************************************************
 * \brief Allocates a task list for the cpu backend.
 *        See grid_task_list.h for details.
 ******************************************************************************/
void grid_cpu_create_task_list(
    const int ntasks, const int nlevels, const int natoms, const int nkinds,
    const int nblocks, const int block_offsets[nblocks],
    const double atom_positions[natoms][3], const int atom_kinds[natoms],
    const grid_basis_set *basis_sets[nkinds], const int level_list[ntasks],
    const int iatom_list[ntasks], const int jatom_list[ntasks],
    const int iset_list[ntasks], const int jset_list[ntasks],
    const int ipgf_list[ntasks], const int jpgf_list[ntasks],
    const int border_mask_list[ntasks], const int block_num_list[ntasks],
    const double radius_list[ntasks], const double rab_list[ntasks][3],
    grid_cpu_task_list **task_list) {

  if (*task_list == NULL) {
    *task_list = malloc(sizeof(grid_cpu_task_list));

    (*task_list)->context = create_grid_context_cpu(
        ntasks, nlevels, natoms, nkinds, nblocks, block_offsets, atom_positions,
        atom_kinds, basis_sets, level_list, iatom_list, jatom_list, iset_list,
        jset_list, ipgf_list, jpgf_list, border_mask_list, block_num_list,
        radius_list, rab_list);
  } else {
    update_grid_context_cpu(
        ntasks, nlevels, natoms, nkinds, nblocks, block_offsets, atom_positions,
        atom_kinds, basis_sets, level_list, iatom_list, jatom_list, iset_list,
        jset_list, ipgf_list, jpgf_list, border_mask_list, block_num_list,
        radius_list, rab_list, (*task_list)->context);
  }

  const grid_library_config config = grid_library_get_config();
  if (config.apply_cutoff) {
    apply_cutoff((*task_list)->context);
  }
}

/*******************************************************************************
 * \brief Deallocates given task list, basis_sets have to be freed separately.
 ******************************************************************************/
void grid_cpu_free_task_list(grid_cpu_task_list *task_list) {
  destroy_grid_context_cpu(task_list->context);
  free(task_list);
}

/*******************************************************************************
 * \brief Collocate all tasks of in given list onto given grids.
 *        See grid_task_list.h for details.
 ******************************************************************************/
void grid_cpu_collocate_task_list(
    const grid_cpu_task_list *task_list, const bool orthorhombic,
    const enum grid_func func, const int nlevels,
    const int npts_global[nlevels][3], const int npts_local[nlevels][3],
    const int shift_local[nlevels][3], const int border_width[nlevels][3],
    const double dh[nlevels][3][3], const double dh_inv[nlevels][3][3],
    const grid_buffer *pab_blocks, double *grid[nlevels]) {

  grid_collocate_task_list_cpu(task_list->context, orthorhombic, func, nlevels,
                               npts_global, npts_local, shift_local,
                               border_width, dh, dh_inv, pab_blocks, grid);
}

// EOF
