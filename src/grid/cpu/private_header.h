/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#ifndef PRIVATE_HEADER_H
#define PRIVATE_HEADER_H

#include "tensor_local.h"
#include <assert.h>
#include <stdbool.h>
/* everything here is specific to the cpu and gpu backends*/
#include "../common/grid_basis_set.h"
#include "../common/grid_common.h"

enum checksum_ { task_checksum = 0x2384989, ctx_checksum = 0x2356734 };

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
  enum checksum_ checksum;
} _task;

typedef struct {
  int ntasks;  // total number of tasks
  int nlevels; // number of different grid
  int natoms;
  int nkinds;
  int nblocks;
  int nblocks_total;
  int nkinds_total;
  int nlevels_total;
  int ntasks_total;

  double *blocks_buffer;
  size_t block_buffer_size;
  size_t block_buffer_size_alloc;
  int *block_offsets;
  double *atom_positions;
  int *atom_kinds;
  grid_basis_set **basis_sets;
  _task **tasks;
  int *tasks_per_level;
  int maxco;
  bool apply_cutoff;
  bool work_on_gpu;
  int number_of_devices;
  int *device_id;
  int queue_length;
  struct collocation_integration_ **handler;
  int number_of_handler;
  tensor *grid;
  void *scratch;
  bool orthorhombic;
  enum checksum_ checksum;
} grid_context;

static inline void update_loop_index(const int xmin, const int xmax,
                                     const int global_grid_size,
                                     int *const x_offset, int *const x,
                                     int *const x1) {
  *x += xmax - xmin - 1;
  *x_offset += (xmax - xmin);

  if (*x1 == global_grid_size) {
    *x1 = -1;
  }
}

static inline int compute_next_boundaries(int *y1, int y,
                                          const int global_grid_size,
                                          const int cube_size) {
  *y1 += imin(cube_size - y, global_grid_size - *y1);
  return *y1;
}

extern void grid_transform_coef_jik_to_yxz(const double dh[3][3],
                                           const tensor *coef_xyz);
extern void grid_transform_coef_xzy_to_ikj(const double dh[3][3],
                                           const tensor *coef_xyz);
extern void compute_block_boundaries(
    const int *blockDim, const int *lb_grid, const int *grid_size,
    const int *blocked_grid_size, const int *period, const int *cube_center,
    const int *cube_size, const int *lower_boundaries_cube,
    int *lower_block_corner, int *upper_block_corner, int *pol_offsets);

extern void grid_fill_pol_dgemm(const bool transpose, const double dr,
                                const double roffset, const int pol_offset,
                                const int xmin, const int xmax, const int lp,
                                const int cmax, const double zetp,
                                double *pol_);

/* this function is not exported outside. Should move */
extern void tensor_reduction_for_collocate_integrate(
    double *scratch, const double alpha, const bool *const orthogonal,
    const struct tensor_ *Exp, const struct tensor_ *co,
    const struct tensor_ *p_alpha_beta_reduced_, struct tensor_ *cube);

extern void set_grid_parameters(
    grid_context *ctx, const int grid_level,
    const int grid_full_size[3],  /* size of the full grid */
    const int grid_local_size[3], /* size of the local grid block */
    const int shift_local[3],     /* coordinates of the lower coordinates of the
                                     local grid window */
    const int border_width[3],    /* width of the borders */
    const double
        dh[3][3], /* displacement vectors of the grid (cartesian) -> (ijk) */
    const double dh_inv[3][3], /* (ijk) -> (x,y,z) */
    double *grid_);

extern void collocate_one_grid_level_dgemm(grid_context *const ctx,
                                           const int *const, const int *const,
                                           const int func, const int level);

#endif
