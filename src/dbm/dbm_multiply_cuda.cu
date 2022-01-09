/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#ifdef __DBM_CUDA

#include <assert.h>
#include <stdio.h>

#include "../offload/offload_library.h"
#include "dbm_mempool.h"
#include "dbm_multiply_cuda.h"

/*******************************************************************************
 * \brief Check given Cuda status and upon failure abort with a nice message.
 * \author Ole Schuett
 ******************************************************************************/
#define CHECK(status)                                                          \
  if (status != cudaSuccess) {                                                 \
    fprintf(stderr, "ERROR: %s %s %d\n", cudaGetErrorString(status), __FILE__, \
            __LINE__);                                                         \
    abort();                                                                   \
  }

/*******************************************************************************
 * \brief Atomic add for doubles that also works prior to compute capability 6.
 * \author Ole Schuett
 ******************************************************************************/
__device__ static void atomicAddDouble(double *address, double val) {
  if (val == 0.0)
    return;

#if __CUDA_ARCH__ >= 600
  atomicAdd(address, val); // part of cuda library
#else
  // https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#atomic-functions
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
 * \brief Internal routine for intializing the cuda backend.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply_cuda_start(const int max_batch_size, const int nshards,
                             dbm_shard_t *shards_c_host,
                             dbm_multiply_cuda_context_t *ctx) {
  // Select GPU device.
  offload_set_device();

  ctx->nshards = nshards;
  ctx->shards_c_host = shards_c_host;
  ctx->max_batch_size = max_batch_size;
  CHECK(cudaStreamCreate(&ctx->main_stream));

  // Allocate device storage for batches.
  const size_t size = nshards * max_batch_size * sizeof(dbm_task_t);
  ctx->batches_dev = (dbm_task_t *)dbm_mempool_device_malloc(size);

  // Allocate and upload shards of result matrix C.
  ctx->shards_c_dev =
      (dbm_shard_cuda_t *)malloc(nshards * sizeof(dbm_shard_cuda_t));
  for (int i = 0; i < nshards; i++) {
    CHECK(cudaStreamCreate(&ctx->shards_c_dev[i].stream));
    ctx->shards_c_dev[i].data_size = ctx->shards_c_host[i].data_size;
    const size_t size = ctx->shards_c_dev[i].data_size * sizeof(double);
    ctx->shards_c_dev[i].data = (double *)dbm_mempool_device_malloc(size);
    CHECK(cudaMemcpyAsync(ctx->shards_c_dev[i].data, ctx->shards_c_host[i].data,
                          size, cudaMemcpyHostToDevice,
                          ctx->shards_c_dev[i].stream));
  }
}

/*******************************************************************************
 * \brief Private routine for uploading a single pack onto the device.
 * \author Ole Schuett
 ******************************************************************************/
static void upload_pack(const dbm_pack_t *pack_host, dbm_pack_t *pack_dev,
                        const cudaStream_t stream) {

  const size_t size = pack_host->data_size * sizeof(double);
  if (pack_dev->data_size < pack_host->data_size) {
    dbm_mempool_free(pack_dev->data);
    pack_dev->data = (double *)dbm_mempool_device_malloc(size);
  }
  CHECK(cudaMemcpyAsync(pack_dev->data, pack_host->data, size,
                        cudaMemcpyHostToDevice, stream));
}

/*******************************************************************************
 * \brief Internal routine for uploading newly arrived packs onto the device.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply_cuda_upload_packs(const dbm_pack_t *pack_a,
                                    const dbm_pack_t *pack_b,
                                    dbm_multiply_cuda_context_t *ctx) {
  // Select GPU device.
  offload_set_device();

  // Wait for all c-streams to complete before overwriting old packs.
  cudaEvent_t event;
  CHECK(cudaEventCreate(&event));
  for (int i = 0; i < ctx->nshards; i++) {
    CHECK(cudaEventRecord(event, ctx->shards_c_dev[i].stream))
    CHECK(cudaStreamWaitEvent(ctx->main_stream, event, 0));
  }

  upload_pack(pack_a, &ctx->pack_a_dev, ctx->main_stream);
  upload_pack(pack_b, &ctx->pack_b_dev, ctx->main_stream);

  // Have all c-streams wait until new packs are uploaded.
  CHECK(cudaEventRecord(event, ctx->main_stream))
  for (int i = 0; i < ctx->nshards; i++) {
    CHECK(cudaStreamWaitEvent(ctx->shards_c_dev[i].stream, event, 0));
  }
  CHECK(cudaEventDestroy(event));
}

/*******************************************************************************
 * \brief A very naive - but generic - matrix multiplication kernel.
 * \author Ole Schuett
 ******************************************************************************/
__global__ static void
process_batch_kernel(const bool transa, const bool transb, const double alpha,
                     const dbm_task_t *batch, const double *pack_a_data,
                     const double *pack_b_data, double *shard_c_data) {

  const dbm_task_t task = batch[blockIdx.x];
  const int lda = (transa) ? task.k : task.m;
  const int ldb = (transb) ? task.n : task.k;
  const int ldc = task.m;
  const double *data_a = &pack_a_data[task.offset_a];
  const double *data_b = &pack_b_data[task.offset_b];
  double *data_c = &shard_c_data[task.offset_c];

  for (int i = threadIdx.z; i < task.m; i += blockDim.z) {
    for (int j = threadIdx.y; j < task.n; j += blockDim.y) {
      for (int l = threadIdx.x; l < task.k; l += blockDim.x) {
        const int idx_a = (transa) ? i * lda + l : l * lda + i;
        const int idx_b = (transb) ? l * ldb + j : j * ldb + l;
        const int idx_c = j * ldc + i;
        const double val = alpha * data_a[idx_a] * data_b[idx_b];
        atomicAddDouble(&data_c[idx_c], val);
      }
    }
  }
}

/*******************************************************************************
 * \brief Internal routine for executing the tasks in given batch on the GPU.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply_cuda_process_batch(const int ntasks, const dbm_task_t *batch,
                                     const bool transa, const bool transb,
                                     const double alpha, const int kshard,
                                     dbm_multiply_cuda_context_t *ctx) {
  if (ntasks == 0) {
    return; // Nothing to do.
  }

  // Select GPU device.
  offload_set_device();

  const dbm_shard_t *shard_c_host = &ctx->shards_c_host[kshard];
  dbm_shard_cuda_t *shard_c_dev = &ctx->shards_c_dev[kshard];

  // Upload new batch.
  dbm_task_t *batch_dev = &ctx->batches_dev[kshard * ctx->max_batch_size];
  const size_t size = ntasks * sizeof(dbm_task_t);
  CHECK(cudaMemcpyAsync(batch_dev, batch, size, cudaMemcpyHostToDevice,
                        shard_c_dev->stream));
  cudaEvent_t batch_uploaded;
  CHECK(cudaEventCreate(&batch_uploaded));
  CHECK(cudaEventRecord(batch_uploaded, shard_c_dev->stream));

  // Grow shard_c_dev->data if nessecary.
  if (shard_c_dev->data_size != shard_c_host->data_promised) {
    // TODO experiment with over-allocation.
    double *old_data_dev = shard_c_dev->data;
    const size_t old_size = shard_c_dev->data_size * sizeof(double);
    shard_c_dev->data_size = shard_c_host->data_promised;
    const size_t new_size = shard_c_dev->data_size * sizeof(double);
    shard_c_dev->data = (double *)dbm_mempool_device_malloc(new_size);
    CHECK(cudaMemsetAsync(shard_c_dev->data, 0, new_size,
                          shard_c_dev->stream)); // TODO: zero only tail
    CHECK(cudaMemcpyAsync(shard_c_dev->data, old_data_dev, old_size,
                          cudaMemcpyDeviceToDevice, shard_c_dev->stream));
    // Wait for copy to complete before freeing old buffer.
    CHECK(cudaStreamSynchronize(shard_c_dev->stream));
    dbm_mempool_free(old_data_dev);
  }

  // Launch kernel.
  const int nblocks = ntasks; // TODO tune launch parameters.
  const dim3 threads_per_block(4, 4, 4);
  const size_t smem_per_block = 0;
  process_batch_kernel<<<nblocks, threads_per_block, smem_per_block,
                         shard_c_dev->stream>>>(
      transa, transb, alpha, batch_dev, ctx->pack_a_dev.data,
      ctx->pack_b_dev.data, shard_c_dev->data);
  CHECK(cudaGetLastError());

  // Wait for batch to be uploaded before refilling it.
  CHECK(cudaEventSynchronize(batch_uploaded));
  CHECK(cudaEventDestroy(batch_uploaded));
}

/*******************************************************************************
 * \brief Internal routine for downloading results from the device.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply_cuda_download_results(dbm_multiply_cuda_context_t *ctx) {
  // Select GPU device.
  offload_set_device();

#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < ctx->nshards; i++) {
    // Grow host buffer if nessecary.
    dbm_shard_t *shard_c_host = &ctx->shards_c_host[i];
    dbm_shard_allocate_promised_blocks(shard_c_host);

    // Download results from device.
    dbm_shard_cuda_t *shard_c_dev = &ctx->shards_c_dev[i];
    assert(shard_c_host->data_size == shard_c_dev->data_size);
    const size_t size = shard_c_dev->data_size * sizeof(double);
    CHECK(cudaMemcpyAsync(shard_c_host->data, shard_c_dev->data, size,
                          cudaMemcpyDeviceToHost, shard_c_dev->stream));
  }
}

/*******************************************************************************
 * \brief Internal routine for shutting down the cuda backend.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply_cuda_stop(dbm_multiply_cuda_context_t *ctx) {
  // Select GPU device.
  offload_set_device();

  // Wait for completion, then free cuda ressources.
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < ctx->nshards; i++) {
    dbm_shard_cuda_t *shard_c_dev = &ctx->shards_c_dev[i];
    CHECK(cudaStreamSynchronize(shard_c_dev->stream));
    CHECK(cudaStreamDestroy(shard_c_dev->stream));
    dbm_mempool_free(shard_c_dev->data);
  }
  free(ctx->shards_c_dev);

  dbm_mempool_free(ctx->pack_a_dev.data);
  dbm_mempool_free(ctx->pack_b_dev.data);
  dbm_mempool_free(ctx->batches_dev);
  CHECK(cudaStreamDestroy(ctx->main_stream));
}

#endif // __DBM_CUDA

// EOF
