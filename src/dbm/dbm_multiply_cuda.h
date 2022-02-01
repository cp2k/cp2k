/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#ifndef DBM_MULTIPLY_CUDA_H
#define DBM_MULTIPLY_CUDA_H

#ifdef __cplusplus
extern "C" {
#endif

#include <cuda_runtime.h>

#include "dbm_multiply_internal.h"
#include "dbm_shard.h"

/*******************************************************************************
 * \brief Internal struct for storing per shard cuda objects.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  double *data; // on the device
  int data_size;
  cudaStream_t stream;
} dbm_shard_cuda_t;

/*******************************************************************************
 * \brief Internal struct for storing the cuda backend's context.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  cudaStream_t main_stream;

  int nshards;
  dbm_shard_t *shards_c_host;
  dbm_shard_cuda_t *shards_c_dev;

  dbm_pack_t pack_a_dev;
  dbm_pack_t pack_b_dev;

  int max_batch_size;
  dbm_task_t *batches_dev;
} dbm_multiply_cuda_context_t;

/*******************************************************************************
 * \brief Internal routine for intializing the cuda backend.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply_cuda_start(const int max_batch_size, const int nshards,
                             dbm_shard_t *shards_c_host,
                             dbm_multiply_cuda_context_t *ctx);

/*******************************************************************************
 * \brief Internal routine for uploading newly arrived packs onto the device.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply_cuda_upload_packs(const dbm_pack_t *pack_a,
                                    const dbm_pack_t *pack_b,
                                    dbm_multiply_cuda_context_t *ctx);

/*******************************************************************************
 * \brief Internal routine for executing the tasks in given batch on the GPU.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply_cuda_process_batch(const int ntasks, const dbm_task_t *batch,
                                     const bool transa, const bool transb,
                                     const double alpha, const int kshard,
                                     dbm_multiply_cuda_context_t *ctx);

/*******************************************************************************
 * \brief Internal routine for downloading results from the device.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply_cuda_download_results(dbm_multiply_cuda_context_t *ctx);

/*******************************************************************************
 * \brief Internal routine for shutting down the cuda backend.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply_cuda_stop(dbm_multiply_cuda_context_t *ctx);

#ifdef __cplusplus
}
#endif
#endif

// EOF
