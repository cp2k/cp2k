/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2021 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#ifndef CPU_PRIVATE_HEADER_H
#define CPU_PRIVATE_HEADER_H

#include "tensor_local.h"
#include <assert.h>
#include <stdbool.h>
/* everything here is specific to the cpu and gpu backends*/
#include "../common/grid_basis_set.h"
#include "../common/grid_buffer.h"
#include "../common/grid_common.h"

enum checksum_ { task_checksum = 0x2384989, ctx_checksum = 0x2356734 };

typedef struct {
  int xmin, xmax;
} Interval;

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
  double zetp;
  double zeta[2];
  double ra[3];
  double rb[3];
  double rp[3];
  int lmax[2];
  int lmin[2];
  int l1_plus_l2_;
  int offset[2];
  bool update_block_;
  double rab[3];
  double prefactor;
  enum checksum_ checksum;
} _task;

typedef struct grid_context_ {
  int ntasks;  // total number of tasks
  int nlevels; // number of different grid
  int natoms;
  int nkinds;
  int nblocks;
  int nblocks_total;
  int nkinds_total;
  int nlevels_total;
  int ntasks_total;
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

static inline void update_loop_index(const int global_grid_size, int x1,
                                     int *const x) {
  *x += global_grid_size - x1 - 1;
}

static inline Interval create_interval(const int xmin, const int xmax) {
  assert(xmax >= xmin);

  Interval t = {.xmin = xmin, .xmax = xmax};
  return t;
}

static inline bool is_point_in_interval(const int value, Interval x) {
  return (value >= x.xmin) && (value <= x.xmax);
}

static inline bool intersection_interval_is_empty(const Interval x,
                                                  const Interval y) {
  /* return true if the intersection is empty */
  if ((x.xmin > y.xmax) || (x.xmax < y.xmin))
    return true;
  else
    return false;
}

static inline Interval intersection_interval(const Interval x,
                                             const Interval y) {
  Interval z;
  z.xmin = imax(x.xmin, y.xmin);
  z.xmax = imin(x.xmax, y.xmax);
  return z;
}

static inline int compute_next_boundaries(const int y1, const int y,
                                          const int grid_size,
                                          const int cube_size) {
  return y1 + imin(cube_size - y, grid_size - y1);
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
    tensor *grid, /* tensor describing the grid */
    const bool orthorhombic,
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
                                           const int func, const int level,
                                           const grid_buffer *pab_blocks);

extern void integrate_one_grid_level_dgemm(
    grid_context *const ctx, const int level, const bool calculate_tau,
    const bool calculate_forces, const bool calculate_virial,
    const int *const shift_local, const int *const border_width,
    const grid_buffer *const pab_blocks, grid_buffer *const hab_blocks,
    tensor *forces_, tensor *virial_);

extern void compute_coefficients(grid_context *const ctx,
                                 struct collocation_integration_ *handler,
                                 const _task *previous_task, const _task *task,
                                 const grid_buffer *pab_blocks,
                                 tensor *const pab, tensor *const work,
                                 tensor *const pab_prep);

extern void extract_blocks(grid_context *const ctx, const _task *const task,
                           const grid_buffer *pab_blocks, tensor *const work,
                           tensor *const pab);

#endif
