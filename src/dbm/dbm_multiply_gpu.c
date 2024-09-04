/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2024 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include "../offload/offload_runtime.h"
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_DBM)

#include "../offload/offload_library.h"
#include "dbm_hyperparams.h"
#include "dbm_mempool.h"
#include "dbm_multiply_gpu.h"
#include "dbm_multiply_gpu_kernel.h"

#include <assert.h>
#include <stdio.h>

/*******************************************************************************
 * \brief Internal routine for intializing the gpu backend.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply_gpu_start(const int max_batch_size, const int nshards,
                            dbm_shard_t *shards_c_host,
                            dbm_multiply_gpu_context_t *ctx) {
  // Select GPU device.
  offload_activate_chosen_device();

  ctx->nshards = nshards;
  ctx->shards_c_host = shards_c_host;
  ctx->max_batch_size = max_batch_size;
  offloadStreamCreate(&ctx->main_stream);

  // Allocate device storage for batches.
  const size_t size = nshards * max_batch_size * sizeof(dbm_task_t);
  ctx->batches_dev = dbm_mempool_device_malloc(size);

  // Allocate and upload shards of result matrix C.
  ctx->shards_c_dev = malloc(nshards * sizeof(dbm_shard_gpu_t));
  for (int i = 0; i < nshards; i++) {
    const dbm_shard_t *shard_c_host = &ctx->shards_c_host[i];
    dbm_shard_gpu_t *shard_c_dev = &ctx->shards_c_dev[i];
    offloadStreamCreate(&shard_c_dev->stream);
    shard_c_dev->data_size = shard_c_host->data_size;
    shard_c_dev->data_allocated = shard_c_host->data_allocated;
    shard_c_dev->data =
        dbm_mempool_device_malloc(shard_c_dev->data_allocated * sizeof(double));
    offloadMemcpyAsyncHtoD(shard_c_dev->data, shard_c_host->data,
                           shard_c_dev->data_size * sizeof(double),
                           shard_c_dev->stream);
  }
}

/*******************************************************************************
 * \brief Private routine for uploading a single pack onto the device.
 * \author Ole Schuett
 ******************************************************************************/
static void upload_pack(const dbm_pack_t *pack_host, dbm_pack_t *pack_dev,
                        const offloadStream_t stream) {

  const size_t size = pack_host->data_size * sizeof(double);
  if (pack_dev->data_size < pack_host->data_size) {
    dbm_mempool_free(pack_dev->data);
    pack_dev->data = dbm_mempool_device_malloc(size);
  }
  offloadMemcpyAsyncHtoD(pack_dev->data, pack_host->data, size, stream);
}

/*******************************************************************************
 * \brief Internal routine for uploading newly arrived packs onto the device.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply_gpu_upload_packs(const dbm_pack_t *pack_a,
                                   const dbm_pack_t *pack_b,
                                   dbm_multiply_gpu_context_t *ctx) {
  // Select GPU device.
  offload_activate_chosen_device();

  // Wait for all c-streams to complete before overwriting old packs.
  offloadEvent_t event;
  offloadEventCreate(&event);
  for (int i = 0; i < ctx->nshards; i++) {
    offloadEventRecord(event, ctx->shards_c_dev[i].stream);
    offloadStreamWaitEvent(ctx->main_stream, event, 0);
  }

  upload_pack(pack_a, &ctx->pack_a_dev, ctx->main_stream);
  upload_pack(pack_b, &ctx->pack_b_dev, ctx->main_stream);

  // Have all c-streams wait until new packs are uploaded.
  offloadEventRecord(event, ctx->main_stream);
  for (int i = 0; i < ctx->nshards; i++) {
    offloadStreamWaitEvent(ctx->shards_c_dev[i].stream, event, 0);
  }
  offloadEventDestroy(event);
}

/*******************************************************************************
 * \brief Internal routine for executing the tasks in given batch on the GPU.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply_gpu_process_batch(const int ntasks, const dbm_task_t *batch,
                                    const int mnk_range[3][2],
                                    const double alpha, const int kshard,
                                    dbm_multiply_gpu_context_t *ctx) {

  if (ntasks == 0) {
    return; // Nothing to do.
  }

  // Select GPU device.
  offload_activate_chosen_device();

  const dbm_shard_t *shard_c_host = &ctx->shards_c_host[kshard];
  dbm_shard_gpu_t *shard_c_dev = &ctx->shards_c_dev[kshard];

  // Upload new batch.
  dbm_task_t *batch_dev = &ctx->batches_dev[kshard * ctx->max_batch_size];
  const size_t size = ntasks * sizeof(dbm_task_t);
  offloadMemcpyAsyncHtoD(batch_dev, batch, size, shard_c_dev->stream);
  offloadEvent_t batch_uploaded;
  offloadEventCreate(&batch_uploaded);
  offloadEventRecord(batch_uploaded, shard_c_dev->stream);

  // Reallocate shard_c_dev->data if necessary.
  if (shard_c_host->data_promised > shard_c_dev->data_allocated) {
    double *old_data_dev = shard_c_dev->data;
    shard_c_dev->data_allocated =
        ALLOCATION_FACTOR * shard_c_host->data_promised;
    shard_c_dev->data =
        dbm_mempool_device_malloc(shard_c_dev->data_allocated * sizeof(double));
    offloadMemcpyAsyncDtoD(shard_c_dev->data, old_data_dev,
                           shard_c_dev->data_size * sizeof(double),
                           shard_c_dev->stream);
    // Wait for copy to complete before freeing old buffer.
    offloadStreamSynchronize(shard_c_dev->stream);
    dbm_mempool_free(old_data_dev);
  }

  // Zero new blocks if necessary.
  if (shard_c_host->data_promised > shard_c_dev->data_size) {
    const int tail = shard_c_host->data_promised - shard_c_dev->data_size;
    offloadMemsetAsync(&shard_c_dev->data[shard_c_dev->data_size], 0,
                       tail * sizeof(double), shard_c_dev->stream);
    shard_c_dev->data_size = shard_c_host->data_promised;
  }

  // Launch kernel.
  dbm_multiply_gpu_launch_kernel(shard_c_dev->stream, mnk_range, alpha, ntasks,
                                 batch_dev, ctx->pack_a_dev.data,
                                 ctx->pack_b_dev.data, shard_c_dev->data);
  OFFLOAD_CHECK(offloadGetLastError());

  // Wait for batch to be uploaded before refilling it.
  offloadEventSynchronize(batch_uploaded);
  offloadEventDestroy(batch_uploaded);
}

/*******************************************************************************
 * \brief Internal routine for downloading results from the device.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply_gpu_download_results(dbm_multiply_gpu_context_t *ctx) {
  // Select GPU device.
  offload_activate_chosen_device();

#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < ctx->nshards; i++) {
    // Grow host buffer if necessary.
    dbm_shard_t *shard_c_host = &ctx->shards_c_host[i];
    dbm_shard_allocate_promised_blocks(shard_c_host);

    // Download results from device.
    dbm_shard_gpu_t *shard_c_dev = &ctx->shards_c_dev[i];
    assert(shard_c_host->data_size == shard_c_dev->data_size);
    const size_t size = shard_c_dev->data_size * sizeof(double);
    offloadMemcpyAsyncDtoH(shard_c_host->data, shard_c_dev->data, size,
                           shard_c_dev->stream);
  }
}

/*******************************************************************************
 * \brief Internal routine for shutting down the gpu backend.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply_gpu_stop(dbm_multiply_gpu_context_t *ctx) {
  // Select GPU device.
  offload_activate_chosen_device();

  // Wait for completion, then free gpu ressources.
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < ctx->nshards; i++) {
    dbm_shard_gpu_t *shard_c_dev = &ctx->shards_c_dev[i];
    offloadStreamSynchronize(shard_c_dev->stream);
    offloadStreamDestroy(shard_c_dev->stream);
    dbm_mempool_free(shard_c_dev->data);
  }
  free(ctx->shards_c_dev);

  dbm_mempool_free(ctx->pack_a_dev.data);
  dbm_mempool_free(ctx->pack_b_dev.data);
  dbm_mempool_free(ctx->batches_dev);
  offloadStreamDestroy(ctx->main_stream);
}

#endif // defined(__OFFLOAD) && !defined(__NO_OFFLOAD_DBM)

// EOF
