/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#ifndef COLLOCATION_INTEGRATION_H
#define COLLOCATION_INTEGRATION_H
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __GRID_CUDA
#include <cublas_v2.h>
#include <cuda.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include "../cpu/tensor_local.h"
#include "private_header.h"
#ifdef __GRID_CUDA
typedef struct pgf_list_gpu_ {
  /* number of devices */
  int number_of_devices;

  /* device_id */
  int device_id;
  /* */
  int lmax;

  /* maximum size of the batch */
  int batch_size;

  /* size of the batch */
  int list_length;

  /* number of elements occupied in the buffer */
  size_t coef_dynamic_alloc_size_gpu_;

  /*  total size of the buffer */
  size_t coef_alloc_size_gpu_;

  /* size of the previously allocated coefficent table */
  size_t coef_previous_alloc_size_;

  /* size of the previously allocated grid */
  size_t data_gpu_old_size_;

  double *coef_cpu_;
  double *coef_gpu_;

  /* Info about the cubes */
  int *coef_offset_cpu_;
  double3 *rp_cpu_;
  double *radius_cpu_, *radius_gpu_;

  int *coef_offset_gpu_;
  double3 *rp_gpu_;

  /* angular momentum */
  int *lmax_cpu_;
  int *lmax_gpu_;

  double *zeta_cpu_;
  double *zeta_gpu_;

  double *data_gpu_;

  cudaStream_t stream;
  cudaEvent_t event;
  bool job_finished;

  cublasHandle_t blas_handle;

  int3 grid_size, grid_lower_corner_position, grid_full_size, window_shift,
      window_size;

  int cmax;
  bool zeroing_grid;

  /* size of the halo when the grid is split over multiple mpi ranks */
  int *border_mask_cpu_;
  int *border_mask_gpu_;

  _task *task_list_cpu_;
  _task *task_list_gpu_;
  int3 border_width;

  struct pgf_list_gpu_ *next;
  /* if true, the grid on the gpu should be reallocated */
  bool durty;
  /* true if the buffers are used for computing already */
  bool running;
  bool apply_cutoff;
} pgf_list_gpu;
#endif

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

  int func;
  int lmin_diff[2];
  int lmax_diff[2];

  int cmax;

#ifdef __GRID_CUDA
  pgf_list_gpu *worker_list;
  int worker_list_size;
#endif

} collocation_integration;

extern struct collocation_integration_ *collocate_create_handle();
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
