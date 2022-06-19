/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include "../offload/offload_runtime.h"

#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_DBM)

#include "../offload/offload_library.h"
#include "dbm_hyperparams.h"
#include "dbm_mempool.h"
#include "dbm_multiply_gpu.h"

#include <assert.h>
#include <stdio.h>

/*******************************************************************************
 * \brief Atomic add for doubles that also works prior to compute capability 6.
 * \author Ole Schuett
 ******************************************************************************/
__device__ static void atomicAddDouble(double *address, double val) {
  if (val == 0.0)
    return;

#if __CUDA_ARCH__ >= 600
  atomicAdd(address, val); // part of gpu library
#else
  // https://docs.nvidia.com/gpu/gpu-c-programming-guide/index.html#atomic-functions
  unsigned long long int *address_as_ull = (unsigned long long int *)address;
  unsigned long long int old = *address_as_ull, assumed;

  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,
                    __double_as_longlong(val + __longlong_as_double(assumed)));

    // Uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
  } while (assumed != old);

#endif
}

#define TILE_DIM 8

/*******************************************************************************
 * \brief A generic matrix multiplication kernel.
 * \author Ole Schuett
 ******************************************************************************/
__global__ static void process_batch_kernel(const double alpha,
                                            const dbm_task_t *batch,
                                            const double *pack_a_data,
                                            const double *pack_b_data,
                                            double *shard_c_data) {

  __shared__ double tile_a[TILE_DIM][TILE_DIM], tile_b[TILE_DIM][TILE_DIM];

  const dbm_task_t task = batch[blockIdx.x];
  const double *block_a = &pack_a_data[task.offset_a];
  const double *block_b = &pack_b_data[task.offset_b];
  double *block_c = &shard_c_data[task.offset_c];

  for (int i_tile = 0; i_tile < task.m; i_tile += TILE_DIM) {
    for (int j_tile = 0; j_tile < task.n; j_tile += TILE_DIM) {
      double result = 0.0;
      for (int l_tile = 0; l_tile < task.k; l_tile += TILE_DIM) {
        // Map indicies to threads such that memory reads are coalesced.
        const int i = i_tile + threadIdx.x;
        const int j = j_tile + threadIdx.x; // Different from j in final loop!
        const int l = l_tile + threadIdx.y;

        // Load tile_a from global into shared memory.
        const int idx_a = l * task.m + i; // transa = "N"
        const bool load_a = (l < task.k && i < task.m);
        tile_a[threadIdx.y][threadIdx.x] = (load_a) ? block_a[idx_a] : 0.0;

        // Load tile_b from global into shared memory.
        const int idx_b = l * task.n + j; // transb = "T"
        const bool load_b = (l < task.k && j < task.n);
        tile_b[threadIdx.y][threadIdx.x] = (load_b) ? block_b[idx_b] : 0.0;

        // Multiply tiles from shared memory.
        __syncthreads();
#pragma unroll
        for (int z = 0; z < TILE_DIM; z++) {
          result += tile_a[z][threadIdx.x] * tile_b[z][threadIdx.y];
        }
        __syncthreads();
      }

      // Add result tile to block_c in global memory.
      const int i = i_tile + threadIdx.x;
      const int j = j_tile + threadIdx.y; // Different from j in inner loop!
      if (i < task.m && j < task.n) {
        const int idx_c = j * task.m + i;
        // Need atomics because other thread blocks might work on same block_c.
        atomicAddDouble(&block_c[idx_c], alpha * result);
      }
    }
  }
}

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
  ctx->batches_dev = (dbm_task_t *)dbm_mempool_device_malloc(size);

  // Allocate and upload shards of result matrix C.
  ctx->shards_c_dev =
      (dbm_shard_gpu_t *)malloc(nshards * sizeof(dbm_shard_gpu_t));
  for (int i = 0; i < nshards; i++) {
    offloadStreamCreate(&ctx->shards_c_dev[i].stream);
    ctx->shards_c_dev[i].data_size = ctx->shards_c_host[i].data_size;
    ctx->shards_c_dev[i].data_allocated = ctx->shards_c_dev[i].data_size;
    const size_t size = ctx->shards_c_dev[i].data_allocated * sizeof(double);
    ctx->shards_c_dev[i].data = (double *)dbm_mempool_device_malloc(size);
    offloadMemcpyAsyncHtoD(ctx->shards_c_dev[i].data,
                           ctx->shards_c_host[i].data, size,
                           ctx->shards_c_dev[i].stream);
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
    pack_dev->data = (double *)dbm_mempool_device_malloc(size);
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

  // Reallocate shard_c_dev->data if nessecary.
  if (shard_c_host->data_promised > shard_c_dev->data_allocated) {
    double *old_data_dev = shard_c_dev->data;
    shard_c_dev->data_allocated =
        ALLOCATION_FACTOR * shard_c_host->data_promised;
    shard_c_dev->data = (double *)dbm_mempool_device_malloc(
        shard_c_dev->data_allocated * sizeof(double));
    offloadMemcpyAsyncDtoD(shard_c_dev->data, old_data_dev,
                           shard_c_dev->data_size * sizeof(double),
                           shard_c_dev->stream);
    // Wait for copy to complete before freeing old buffer.
    offloadStreamSynchronize(shard_c_dev->stream);
    dbm_mempool_free(old_data_dev);
  }

  // Zero new blocks if nessecary.
  if (shard_c_host->data_promised > shard_c_dev->data_size) {
    const int tail = shard_c_host->data_promised - shard_c_dev->data_size;
    offloadMemsetAsync(&shard_c_dev->data[shard_c_dev->data_size], 0,
                       tail * sizeof(double), shard_c_dev->stream);
    shard_c_dev->data_size = shard_c_host->data_promised;
  }

  // Launch kernel.
  const int nblocks = ntasks; // TODO tune launch parameters.
  const dim3 threads_per_block(TILE_DIM, TILE_DIM, 1);
  const size_t smem_per_block = 0;
  process_batch_kernel<<<nblocks, threads_per_block, smem_per_block,
                         shard_c_dev->stream>>>(
      alpha, batch_dev, ctx->pack_a_dev.data, ctx->pack_b_dev.data,
      shard_c_dev->data);
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
    // Grow host buffer if nessecary.
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
