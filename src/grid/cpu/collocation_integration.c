/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#ifdef __GRID_CUDA
#include <cublas_v2.h>
#include <cuda.h>
#endif

#include "../common/grid_common.h"
#include "collocation_integration.h"
#include "non_orthorombic_corrections.h"
#include "tensor_local.h"
#include "utils.h"

struct collocation_integration_ *collocate_create_handle() {
  struct collocation_integration_ *handle = NULL;
  handle = (struct collocation_integration_ *)malloc(
      sizeof(struct collocation_integration_));

  if (handle == NULL) {
    abort();
  }
  memset(handle, 0, sizeof(struct collocation_integration_));

  handle->alpha.alloc_size_ = 8192;
  handle->coef.alloc_size_ = 1024;
  handle->pol.alloc_size_ = 1024;
  /* it is a cube of size 32 x 32 x 32 */
  handle->cube.alloc_size_ = 32768;

  handle->cube_alloc_size = realloc_tensor(&handle->cube);
  handle->alpha_alloc_size = realloc_tensor(&handle->alpha);
  handle->coef_alloc_size = realloc_tensor(&handle->coef);
  handle->pol_alloc_size = realloc_tensor(&handle->pol);

  handle->scratch = memalign(4096, sizeof(double) * 32768);
  handle->scratch_alloc_size = 32768;
  handle->T_alloc_size = 8192;
  handle->W_alloc_size = 2048;
  handle->blockDim[0] = 5;
  handle->blockDim[1] = 5;
  handle->blockDim[2] = 5;
  handle->device_id = (int *)malloc(sizeof(double) * 12);
  handle->number_of_devices = 1;

  /* to suppress when we remove the spherical cutoff */
  handle->map = (int **)malloc(3 * sizeof(int *));
  handle->map[0] = (int *)malloc(sizeof(int) * 512 * 3);
  handle->map[1] = handle->map[0] + 512;
  handle->map[2] = handle->map[1] + 512;
  handle->cmax = 512 * 3;
  return handle;
}

void collocate_destroy_handle(void *gaussian_handle) {
  struct collocation_integration_ *handle =
      (struct collocation_integration_ *)gaussian_handle;
  if (handle->Exp.data)
    free(handle->Exp.data);

  if (handle->grid.data)
    free(handle->grid.data);

  free(handle->scratch);
  free(handle->pol.data);
  free(handle->cube.data);
  free(handle->blocks_coordinates.data);
  handle->alpha.data = NULL;
  handle->coef.data = NULL;
  handle->blocks_coordinates.data = NULL;
  free(handle->device_id);
  free(handle->map[0]);
  free(handle->map);
  free(handle);

  handle = NULL;
}

void initialize_W_and_T(collocation_integration *const handler,
                        const tensor *cube, const tensor *coef) {
  size_t tmp1 = compute_memory_space_tensor_3(coef->size[0] /* alpha */,
                                              coef->size[1] /* gamma */,
                                              cube->size[1] /* j */);

  size_t tmp2 = compute_memory_space_tensor_3(
      coef->size[0] /* gamma */, cube->size[1] /* j */, cube->size[2] /* i */);

  const size_t mem_alloc_size_ =
      imax(imax(tmp1 + tmp2, cube->alloc_size_), coef->alloc_size_);

  handler->T_alloc_size = tmp1;
  handler->W_alloc_size = tmp2;

  if ((mem_alloc_size_ > handler->scratch_alloc_size) ||
      (handler->scratch == NULL)) {
    handler->scratch_alloc_size = mem_alloc_size_;

    if (handler->scratch)
      free(handler->scratch);
    handler->scratch =
        memalign(64, sizeof(double) * handler->scratch_alloc_size);
    if (handler->scratch == NULL)
      abort();
  }
}

void initialize_W_and_T_integrate(collocation_integration *const handler,
                                  const int num_block, const tensor *coef,
                                  const tensor *block) {
  /* T */
  size_t tmp1 = compute_memory_space_tensor_4(num_block, block->size[0] /* k */,
                                              block->size[1] /* j */,
                                              coef->size[1] /* alpha */);

  /* W */
  size_t tmp2 = compute_memory_space_tensor_4(num_block, block->size[1] /* j */,
                                              coef->size[1] /* alpha */,
                                              coef->size[2] /* gamma */);

  const size_t mem_alloc_size_ = tmp1 + tmp2;

  handler->T_alloc_size = tmp1;
  handler->W_alloc_size = tmp2;

  if ((mem_alloc_size_ > handler->scratch_alloc_size) ||
      (handler->scratch == NULL)) {
    handler->scratch_alloc_size = mem_alloc_size_;

    if (handler->scratch)
      free(handler->scratch);
    handler->scratch =
        memalign(64, sizeof(double) * handler->scratch_alloc_size);
    if (handler->scratch == NULL)
      abort();
  }
}

void initialize_basis_vectors(collocation_integration *const handler,
                              const double dh[3][3],
                              const double dh_inv[3][3]) {
  handler->dh[0][0] = dh[0][0];
  handler->dh[0][1] = dh[0][1];
  handler->dh[0][2] = dh[0][2];
  handler->dh[1][0] = dh[1][0];
  handler->dh[1][1] = dh[1][1];
  handler->dh[1][2] = dh[1][2];
  handler->dh[2][0] = dh[2][0];
  handler->dh[2][1] = dh[2][1];
  handler->dh[2][2] = dh[2][2];

  handler->dh_inv[0][0] = dh_inv[0][0];
  handler->dh_inv[0][1] = dh_inv[0][1];
  handler->dh_inv[0][2] = dh_inv[0][2];
  handler->dh_inv[1][0] = dh_inv[1][0];
  handler->dh_inv[1][1] = dh_inv[1][1];
  handler->dh_inv[1][2] = dh_inv[1][2];
  handler->dh_inv[2][0] = dh_inv[2][0];
  handler->dh_inv[2][1] = dh_inv[2][1];
  handler->dh_inv[2][2] = dh_inv[2][2];

  /* Only used when we are in the non  orthorombic case */
  handler->dx[2] = handler->dh[0][0] * handler->dh[0][0] +
                   handler->dh[0][1] * handler->dh[0][1] +
                   handler->dh[0][2] * handler->dh[0][2];
  handler->dx[1] = handler->dh[1][0] * handler->dh[1][0] +
                   handler->dh[1][1] * handler->dh[1][1] +
                   handler->dh[1][2] * handler->dh[1][2];
  handler->dx[0] = handler->dh[2][0] * handler->dh[2][0] +
                   handler->dh[2][1] * handler->dh[2][1] +
                   handler->dh[2][2] * handler->dh[2][2];
}
