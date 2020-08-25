/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2020  CP2K developers group                         *
 *****************************************************************************/
#ifndef GRID_CPU_TASK_LIST_H
#define GRID_CPU_TASK_LIST_H

#include <stdbool.h>

#include "../common/grid_basis_set.h"
#include "../ref/grid_ref_task_list.h"

//******************************************************************************
// \brief Internal representation of a task list.
//******************************************************************************
typedef struct {
  //
  // TODO: Replace with actual cpu backend implementation.
  //
  grid_ref_task_list *crutch;
} grid_cpu_task_list;

//******************************************************************************
// \brief Allocates a task list which can be passed to grid_collocate_task_list.
//        See grid_task_list.h for details.
//******************************************************************************
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
    double **blocks_buffer, grid_cpu_task_list **task_list);

//******************************************************************************
// \brief Deallocates given task list, basis_sets have to be freed separately.
//******************************************************************************
void grid_cpu_free_task_list(grid_cpu_task_list *task_list);

//******************************************************************************
// \brief Collocate all tasks of in given list onto given grids.
//        See grid_task_list.h for details.
//******************************************************************************
void grid_cpu_collocate_task_list(
    const grid_cpu_task_list *task_list, const bool orthorhombic,
    const int func, const int nlevels, const int npts_global[nlevels][3],
    const int npts_local[nlevels][3], const int shift_local[nlevels][3],
    const int border_width[nlevels][3], const double dh[nlevels][3][3],
    const double dh_inv[nlevels][3][3], double *grid[nlevels]);

#endif

// EOF
