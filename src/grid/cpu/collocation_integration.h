/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#ifndef COLLOCATION_INTEGRATION_H
#define COLLOCATION_INTEGRATION_H
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "../cpu/tensor_local.h"
#include "cpu_private_header.h"

typedef struct collocation_integration_ {
  /* number of compute device */
  int number_of_devices;
  int *device_id;

  /* number of gaussians block in each list */
  int number_of_gaussian;

  /* some scratch storage to avoid malloc / free all the time */
  tensor alpha;
  tensor pol;
  tensor coef;

  /* tensors for the grid to collocate or integrate */
  /* original grid */
  tensor grid;

  /* int period[3]; */
  /* int lb_grid[3]; */

  /* original grid decomposed in block */
  tensor blocked_grid;

  /* do we need to update the grid */
  bool grid_restored;

  /* coordinates of the blocks */
  tensor blocks_coordinates;

  double dh[3][3];
  double dh_inv[3][3];
  double dx[3];

  /* block dimensions */
  int blockDim[4];

  /* Only allocated in sequential mode */
  tensor cube;
  tensor Exp;
  size_t Exp_alloc_size;
  size_t cube_alloc_size;
  size_t coef_alloc_size;
  size_t alpha_alloc_size;
  size_t pol_alloc_size;
  size_t scratch_alloc_size;
  size_t T_alloc_size;
  size_t W_alloc_size;
  int lmax;
  /* for the spherical cutoff */
  int **map;
  void *scratch;

  bool durty;
  bool orthogonal[3];
  bool integrate;
  bool apply_cutoff;

  enum grid_func func;
  int lmin_diff[2];
  int lmax_diff[2];

  int cmax;
} collocation_integration;

extern struct collocation_integration_ *collocate_create_handle(void);
extern void collocate_synchronize(void *gaussian_handler);
extern void collocate_destroy_handle(void *gaussian_handle);
extern void calculate_collocation(void *const in);
extern void initialize_W_and_T(collocation_integration *const handler,
                               const tensor *cube, const tensor *coef);
extern void initialize_basis_vectors(collocation_integration *const handler,
                                     const double dh[3][3],
                                     const double dh_inv[3][3]);

#ifdef __cplusplus
}
#endif
#endif
