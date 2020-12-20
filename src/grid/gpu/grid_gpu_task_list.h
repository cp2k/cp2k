/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/
#ifndef GRID_GPU_TASK_LIST_H
#define GRID_GPU_TASK_LIST_H

#ifdef __GRID_CUDA

#ifdef __cplusplus
extern "C" {
#endif

#include "../common/grid_basis_set.h"
#include "../common/grid_buffer.h"
#include "../common/grid_constants.h"
#include <cuda_runtime.h>
#include <stdbool.h>

/*******************************************************************************
 * \brief Internal representation of a task.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  int level;
  int iatom;
  int jatom;
  int iset;
  int jset;
  int ipgf;
  int jpgf;
  int border_mask;
  int block_num;
  double radius;
  double rab[3];
} grid_gpu_task;

/*******************************************************************************
 * \brief Internal representation of a task list.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  int ntasks;
  int nlevels;
  int natoms;
  int nkinds;
  int nblocks;
  int *tasks_per_level;
  cudaStream_t *level_streams;
  cudaStream_t main_stream;
  int lmax;
  // device pointers
  int *block_offsets_dev;
  double *atom_positions_dev;
  int *atom_kinds_dev;
  grid_basis_set *basis_sets_dev;
  grid_gpu_task *tasks_dev;
  double **grid_dev;
  size_t *grid_dev_size;
} grid_gpu_task_list;

/*******************************************************************************
 * \brief Allocates a task list for the GPU backend.
 *        See grid_task_list.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_gpu_create_task_list(
    const int ntasks, const int nlevels, const int natoms, const int nkinds,
    const int nblocks, const int block_offsets[],
    const double atom_positions[][3], const int atom_kinds[],
    const grid_basis_set *basis_sets[], const int level_list[],
    const int iatom_list[], const int jatom_list[], const int iset_list[],
    const int jset_list[], const int ipgf_list[], const int jpgf_list[],
    const int border_mask_list[], const int block_num_list[],
    const double radius_list[], const double rab_list[][3],
    grid_gpu_task_list **task_list_out);

/*******************************************************************************
 * \brief Deallocates given task list, basis_sets have to be freed separately.
 * \author Ole Schuett
 ******************************************************************************/
void grid_gpu_free_task_list(grid_gpu_task_list *task_list);

/*******************************************************************************
 * \brief Collocate all tasks of in given list onto given grids.
 *        See grid_task_list.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_gpu_collocate_task_list(
    const grid_gpu_task_list *task_list, const bool orthorhombic,
    const enum grid_func func, const int nlevels, const int npts_global[][3],
    const int npts_local[][3], const int shift_local[][3],
    const int border_width[][3], const double dh[][3][3],
    const double dh_inv[][3][3], const grid_buffer *pab_blocks, double *grid[]);

/*******************************************************************************
 * \brief Integrate all tasks of in given list onto given grids.
 *        See grid_task_list.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_gpu_integrate_task_list(
    const grid_gpu_task_list *task_list, const bool orthorhombic,
    const bool compute_tau, const int natoms, const int nlevels,
    const int npts_global[][3], const int npts_local[][3],
    const int shift_local[][3], const int border_width[][3],
    const double dh[][3][3], const double dh_inv[][3][3],
    const grid_buffer *pab_blocks, const double *grid[],
    grid_buffer *hab_blocks, double forces[][3], double virial[3][3]);

#ifdef __cplusplus
}
#endif

#endif // __GRID_CUDA
#endif

// EOF
