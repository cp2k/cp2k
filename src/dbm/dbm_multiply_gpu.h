/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2024 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#ifndef DBM_MULTIPLY_GPU_H
#define DBM_MULTIPLY_GPU_H

#include "../offload/offload_runtime.h"
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_DBM)

#include "dbm_multiply_internal.h"
#include "dbm_shard.h"

/*******************************************************************************
 * \brief Internal struct for storing per shard gpu objects.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  double *data; // on the device
  int data_size;
  int data_allocated;
  offloadStream_t stream;
} dbm_shard_gpu_t;

/*******************************************************************************
 * \brief Internal struct for storing the gpu backend's context.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  offloadStream_t main_stream;

  int nshards;
  dbm_shard_t *shards_c_host;
  dbm_shard_gpu_t *shards_c_dev;

  dbm_pack_t pack_a_dev;
  dbm_pack_t pack_b_dev;

  int max_batch_size;
  dbm_task_t *batches_dev;
} dbm_multiply_gpu_context_t;

/*******************************************************************************
 * \brief Internal routine for intializing the gpu backend.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply_gpu_start(const int max_batch_size, const int nshards,
                            dbm_shard_t *shards_c_host,
                            dbm_multiply_gpu_context_t *ctx);

/*******************************************************************************
 * \brief Internal routine for uploading newly arrived packs onto the device.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply_gpu_upload_packs(const dbm_pack_t *pack_a,
                                   const dbm_pack_t *pack_b,
                                   dbm_multiply_gpu_context_t *ctx);

/*******************************************************************************
 * \brief Internal routine for executing the tasks in given batch on the GPU.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply_gpu_process_batch(const int ntasks, const dbm_task_t *batch,
                                    const int mnk_range[3][2],
                                    const double alpha, const int kshard,
                                    dbm_multiply_gpu_context_t *ctx);

/*******************************************************************************
 * \brief Internal routine for downloading results from the device.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply_gpu_download_results(dbm_multiply_gpu_context_t *ctx);

/*******************************************************************************
 * \brief Internal routine for shutting down the gpu backend.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply_gpu_stop(dbm_multiply_gpu_context_t *ctx);

#endif // defined(__OFFLOAD) && !defined(__NO_OFFLOAD_DBM)
#endif

// EOF
