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

#define NUM_THREADS 256

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

/*******************************************************************************
 * \brief Returns the larger of two given integer (missing from the C standard)
 * \author Ole Schuett
 ******************************************************************************/
__device__ constexpr static inline int imax(int x, int y) {
  return (x > y ? x : y);
}

/*******************************************************************************
 * \brief Processes a task using tile dimensions give by template parameters.
 * \author Ole Schuett
 ******************************************************************************/
template <int M, int N>
__device__ static void process_task(const int m, const int n, const int k,
                                    const double alpha, const double *block_a,
                                    const double *block_b, double *block_c,
                                    double *shared_mem) {

  constexpr int K = NUM_THREADS / imax(M, N);

  static_assert(K * imax(M, N) == NUM_THREADS, "Wasting threads.");
  static_assert(M * N <= NUM_THREADS, "Not enough threads to cover tile.");

  double *tile_a = shared_mem;
  double *tile_b = &shared_mem[K * M];

  // Layout threads to cover the entire tile.
  // Indices x,y,z mark the position within the tile.
  const int y = threadIdx.x / M;       // Can exceed N, guaranteed to reach K.
  const int x = threadIdx.x - y * M;   // Fastest index, does not exceed M.
  const int xT = threadIdx.x / N;      // Can exceed M, guaranteed to reach K.
  const int yT = threadIdx.x - xT * N; // Fastest index, does not exceed N.

  // Indices {ijl}_tile mark the beginning of the current tile.
  for (int i_tile = 0; i_tile < m; i_tile += M) {
    for (int j_tile = 0; j_tile < n; j_tile += N) {
      double result = 0.0;
      for (int l_tile = 0; l_tile < k; l_tile += K) {

        // Load tile_a from global into shared memory.
        if (x < M && y < K) {
          const int i = i_tile + x;
          const int l = l_tile + y;
          const int idx = l * m + i; // transa = "N"
          const bool load = (l < k && i < m);
          tile_a[y * M + x] = (load) ? block_a[idx] : 0.0;
        }

        // Load tile_b from global into shared memory.
        // Use transposed thread mapping to achieve coalesced memory reads.
        if (yT < N && xT < K) {
          const int j = j_tile + yT;
          const int l = l_tile + xT;
          const int idx = l * n + j; // transb = "T"
          const bool load = (l < k && j < n);
          tile_b[xT * N + yT] = (load) ? block_b[idx] : 0.0;
        }

        // Multiply tiles from shared memory.
        __syncthreads();
        if (x < M && y < N) {
#pragma unroll
          for (int z = 0; z < K; z++) {
            result += tile_a[z * M + x] * tile_b[z * N + y];
          }
        }
        __syncthreads();
      }

      // Add result tile to block_c in global memory.
      if (x < M && y < N) {
        const int i = i_tile + x;
        const int j = j_tile + y;
        if (i < m && j < n) {
          const int idx = j * m + i;
          // Need atomics because other thread blocks might work on same C.
          atomicAddDouble(&block_c[idx], alpha * result);
        }
      }
    }
  }
}

/*******************************************************************************
 * \brief A generic matrix multiplication kernel.
 * \author Ole Schuett
 ******************************************************************************/
__global__ static void process_batch_kernel(const double alpha,
                                            const dbm_task_t *batch,
                                            const double *pack_a_data,
                                            const double *pack_b_data,
                                            double *shard_c_data) {

  __shared__ double shared_mem[2 * NUM_THREADS];

  const dbm_task_t task = batch[blockIdx.x];

  const int m = task.m;
  const int n = task.n;
  const int k = task.k;

  const double *blk_a = &pack_a_data[task.offset_a];
  const double *blk_b = &pack_b_data[task.offset_b];
  double *blk_c = &shard_c_data[task.offset_c];

  if (m <= 4 && n <= 4) {
    process_task<4, 4>(m, n, k, alpha, blk_a, blk_b, blk_c, shared_mem);

  } else if (m <= 4 && n <= 8) {
    process_task<4, 8>(m, n, k, alpha, blk_a, blk_b, blk_c, shared_mem);

  } else if (m <= 4 && n <= 16) {
    process_task<4, 16>(m, n, k, alpha, blk_a, blk_b, blk_c, shared_mem);

  } else if (m <= 4 && n <= 32) {
    process_task<4, 32>(m, n, k, alpha, blk_a, blk_b, blk_c, shared_mem);

  } else if (m <= 4) {
    process_task<4, 64>(m, n, k, alpha, blk_a, blk_b, blk_c, shared_mem);

  } else if (m <= 8 && n <= 4) {
    process_task<8, 4>(m, n, k, alpha, blk_a, blk_b, blk_c, shared_mem);

  } else if (m <= 16 && n <= 4) {
    process_task<16, 4>(m, n, k, alpha, blk_a, blk_b, blk_c, shared_mem);

  } else if (m <= 32 && n <= 4) {
    process_task<32, 4>(m, n, k, alpha, blk_a, blk_b, blk_c, shared_mem);

  } else if (n <= 4) {
    process_task<64, 4>(m, n, k, alpha, blk_a, blk_b, blk_c, shared_mem);

  } else if (m <= 8 && n <= 8) {
    process_task<8, 8>(m, n, k, alpha, blk_a, blk_b, blk_c, shared_mem);

  } else if (m <= 8 && n <= 16) {
    process_task<8, 16>(m, n, k, alpha, blk_a, blk_b, blk_c, shared_mem);

  } else if (m <= 8) {
    process_task<8, 32>(m, n, k, alpha, blk_a, blk_b, blk_c, shared_mem);

  } else if (m <= 16 && n <= 8) {
    process_task<16, 8>(m, n, k, alpha, blk_a, blk_b, blk_c, shared_mem);

  } else if (n <= 8) {
    process_task<32, 8>(m, n, k, alpha, blk_a, blk_b, blk_c, shared_mem);

  } else {
    process_task<16, 16>(m, n, k, alpha, blk_a, blk_b, blk_c, shared_mem);
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
  const int threads_per_block = NUM_THREADS;
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
