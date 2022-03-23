/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#ifndef GRID_GPU_TASK_LIST_H
#define GRID_GPU_TASK_LIST_H

#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_GRID)

#ifdef __cplusplus
extern "C" {
#endif

#include "../../offload/offload_buffer.h"
#include "../../offload/offload_runtime.h"
#include "../common/grid_basis_set.h"
#include "../common/grid_constants.h"
#include <stdbool.h>

/*******************************************************************************
 * \brief Internal representation of a task.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct grid_gpu_task_struct {
  bool use_orthorhombic_kernel;
  bool block_transposed;
  double radius;
  double radius2;
  double ra[3];
  double rb[3];
  double rp[3];
  double rab[3];
  double gp[3];
  double rab2;
  double zeta;
  double zetb;
  double zetp;
  double prefactor;
  double off_diag_twice;
  double dh_max;
  // angular momentum range of basis set
  int la_max_basis;
  int lb_max_basis;
  int la_min_basis;
  int lb_min_basis;
  // size of entire spherical basis
  int nsgfa;
  int nsgfb;
  // size of spherical set
  int nsgf_seta;
  int nsgf_setb;
  // start of decontracted set, ie. pab and hab
  int first_coseta;
  int first_cosetb;
  // size of decontracted set, ie. pab and hab
  int ncoseta;
  int ncosetb;
  // strides of the sphi transformation matrices
  int maxcoa;
  int maxcob;

  // offset of the pab and hab block relative to buffer pointer.
  int ab_block_offset;

  // atoms to which the forces and virial should be added
  int iatom;
  int jatom;

  // pointers basis set matrices
  const double *sphia;
  const double *sphib;

  // Stuff for the ortho kernel.
  double disr_radius;
  int cube_center_shifted[3];
  double cube_offset[3];

  // Stuff for the general kernel.
  int index_min[3];
  int index_max[3];
  int bounds_i[2];
  int bounds_j[2];
  int bounds_k[2];
} grid_gpu_task;

/*******************************************************************************
 * \brief Internal representation of a grid layout.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  int npts_global[3];
  int npts_local[3];
  int shift_local[3];
  int border_width[3];
  double dh[3][3];
  double dh_inv[3][3];
} grid_gpu_layout;

/*******************************************************************************
 * \brief Internal representation of a task list.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  bool orthorhombic;
  int ntasks;
  int nlevels;
  int natoms;
  int nkinds;
  int nblocks;
  grid_gpu_layout *layouts;
  int *tasks_per_level;
  offloadStream_t *level_streams;
  offloadStream_t main_stream;
  int lmax;
  int stats[2][20]; // [has_border_mask][lp]
  // device pointers
  double **sphis_dev;
  grid_gpu_task *tasks_dev;
} grid_gpu_task_list;

/*******************************************************************************
 * \brief Allocates a task list for the GPU backend.
 *        See grid_task_list.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_gpu_create_task_list(
    const bool orthorhombic, const int ntasks, const int nlevels,
    const int natoms, const int nkinds, const int nblocks,
    const int block_offsets[], const double atom_positions[][3],
    const int atom_kinds[], const grid_basis_set *basis_sets[],
    const int level_list[], const int iatom_list[], const int jatom_list[],
    const int iset_list[], const int jset_list[], const int ipgf_list[],
    const int jpgf_list[], const int border_mask_list[],
    const int block_num_list[], const double radius_list[],
    const double rab_list[][3], const int npts_global[][3],
    const int npts_local[][3], const int shift_local[][3],
    const int border_width[][3], const double dh[][3][3],
    const double dh_inv[][3][3], grid_gpu_task_list **task_list);

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
void grid_gpu_collocate_task_list(const grid_gpu_task_list *task_list,
                                  const enum grid_func func, const int nlevels,
                                  const offload_buffer *pab_blocks,
                                  offload_buffer *grids[]);

/*******************************************************************************
 * \brief Integrate all tasks of in given list onto given grids.
 *        See grid_task_list.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_gpu_integrate_task_list(const grid_gpu_task_list *task_list,
                                  const bool compute_tau, const int natoms,
                                  const int nlevels,
                                  const offload_buffer *pab_blocks,
                                  const offload_buffer *grids[],
                                  offload_buffer *hab_blocks,
                                  double forces[][3], double virial[3][3]);

#ifdef __cplusplus
}
#endif

#endif // defined(__OFFLOAD) && !defined(__NO_OFFLOAD_GRID)
#endif

// EOF
