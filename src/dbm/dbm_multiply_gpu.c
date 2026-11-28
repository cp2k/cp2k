/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include "../offload/offload_runtime.h"
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_DBM)

#include "../offload/offload_library.h"
#include "../offload/offload_mempool.h"
#include "dbm_hyperparams.h"
#include "dbm_multiply_cpu.h"
#include "dbm_multiply_gpu.h"
#include "dbm_multiply_gpu_kernel.h"

#include "../base/cp_assert.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*******************************************************************************
 * \brief Internal routine returning the maximum number of threads computing
 *        on the host while the GPU is busy (DBM_MULTIPLY_HYBRID, default: 0).
 * \author Hans Pabst
 ******************************************************************************/
int dbm_multiply_gpu_hybrid(void) {
  const char *const hybrid_env = getenv("DBM_MULTIPLY_HYBRID");
  return (NULL == hybrid_env ? 0 : atoi(hybrid_env));
}

/*******************************************************************************
 * \brief Internal routine for initializing the gpu backend.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply_gpu_start(const int max_batch_size, const int nshards,
                            dbm_shard_t *shards_c_host,
                            dbm_multiply_gpu_context_t *ctx) {
  // Select GPU device.
  offload_activate_chosen_device();

  ctx->nshards = nshards;
  ctx->max_batch_size = max_batch_size;
  ctx->shards_c_host = shards_c_host;

  // Up to |DBM_MULTIPLY_HYBRID| threads compute on the host while the GPU is
  // busy, and a negative value reports the share of the host.
  ctx->hybrid = dbm_multiply_gpu_hybrid();
  static bool hybrid_reported = false; // once per process
  if (0 > ctx->hybrid && !hybrid_reported &&
      !dbm_multiply_cpu_generated(0, 0, 0, 1.0, 0)) {
    fprintf(stderr, "INFO DBM: no generated host kernels, hence no hybrid\n");
    hybrid_reported = true;
  }
  ctx->hybrid_active = 0;
  ctx->flops[0] = ctx->flops[1] = 0;
  offloadStreamCreate(&ctx->main_stream);
  offloadEventCreate(&ctx->upload_event);

#if !defined(__OFFLOAD_UNIFIED_MEMORY)
  // Allocate device storage for batches.
  const size_t size = nshards * max_batch_size * sizeof(dbm_task_t);
  ctx->batches_dev = offload_mempool_device_malloc(size);
#endif

  // Allocate and upload shards of result matrix C.
  ctx->shards_c_dev = malloc(nshards * sizeof(dbm_shard_gpu_t));
  cp_assert(ctx->shards_c_dev != NULL || nshards == 0);
  for (int i = 0; i < nshards; i++) {
    const dbm_shard_t *const shard_c_host = &shards_c_host[i];
    dbm_shard_gpu_t *const shard_g = &ctx->shards_c_dev[i];
    shard_g->data_size = shard_c_host->data_size;
    offloadStreamCreate(&shard_g->stream);
    offloadEventCreate(&shard_g->event);
    offloadEventCreate(&shard_g->done);
    offloadEventCreate(&shard_g->download);
    shard_g->downloading = false;
    shard_g->host_data = NULL;
    shard_g->host_size = shard_g->host_allocated = 0;
    // only allocate data_size on device rather than data_allocated
    shard_g->data_allocated = shard_c_host->data_size;
    shard_g->data =
        offload_mempool_device_malloc(shard_g->data_allocated * sizeof(double));
    offloadMemcpyAsyncHtoD(shard_g->data, shard_c_host->data,
                           shard_g->data_size * sizeof(double),
                           shard_g->stream);
  }
}

/*******************************************************************************
 * \brief Private routine for uploading a single pack onto the device.
 * \author Ole Schuett
 ******************************************************************************/
static void upload_pack(const dbm_pack_t *pack_host, dbm_pack_gpu_t *pack_dev,
                        const offloadStream_t stream) {
  // Reallocate only when the device buffer is too small; data_allocated tracks
  // the device capacity independently from the host-side data_size that changes
  // each tick (previous code compared against the wrong field).
  const size_t size = pack_host->data_size * sizeof(double);
  if (pack_dev->data_allocated < pack_host->data_size) {
    offload_mempool_device_free(pack_dev->data);
    pack_dev->data = offload_mempool_device_malloc(size);
    pack_dev->data_allocated = pack_host->data_size;
  }
  offloadMemcpyAsyncHtoD(pack_dev->data, pack_host->data, size, stream);
}

/*******************************************************************************
 * \brief Internal routine for uploading newly arrived packs onto the device.
 * \author Ole Schuett and Hans Pabst
 ******************************************************************************/
bool dbm_multiply_gpu_upload_packs(const dbm_pack_t *pack_a,
                                   const dbm_pack_t *pack_b,
                                   dbm_multiply_gpu_context_t *ctx) {
  // Assume GPU device was activated earlier.
  // Wait for all c-streams to complete before overwriting old packs.
  for (int i = 0; i < ctx->nshards; i++) {
    offloadEventRecord(ctx->upload_event, ctx->shards_c_dev[i].stream);
    offloadStreamWaitEvent(ctx->main_stream, ctx->upload_event);
  }
  // Record event to check if all c-streams already completed.
  offloadEventRecord(ctx->upload_event, ctx->main_stream);
  // A growing pack frees its buffer on the host while kernels may still read
  // it (the waits above are device-side only).
  if (ctx->pack_a_dev.data_allocated < pack_a->data_size ||
      ctx->pack_b_dev.data_allocated < pack_b->data_size) {
    offloadEventSynchronize(ctx->upload_event);
  }

  bool uploaded = false;
  /*if (offloadEventQuery(ctx->upload_event))*/
  {
    upload_pack(pack_a, &ctx->pack_a_dev, ctx->main_stream);
    upload_pack(pack_b, &ctx->pack_b_dev, ctx->main_stream);

    // Have all c-streams wait until new packs are uploaded.
    offloadEventRecord(ctx->upload_event, ctx->main_stream);
    for (int i = 0; i < ctx->nshards; i++) {
      offloadStreamWaitEvent(ctx->shards_c_dev[i].stream, ctx->upload_event);
    }
    uploaded = true;
  }

  return uploaded;
}

/*******************************************************************************
 * \brief Internal routine for concluding a turn of the main loop, i.e., before
 *        the packs are reused or released and before batches are refilled.
 * \author Hans Pabst
 ******************************************************************************/
void dbm_multiply_gpu_end_turn(dbm_multiply_gpu_context_t *ctx) {
  // A turn without tasks does not otherwise wait for the asynchronous uploads.
  offloadEventSynchronize(ctx->upload_event);
#if defined(__OFFLOAD_UNIFIED_MEMORY)
  // The next turn refills the batches the previous turn's kernels read.
  for (int i = 0; i < ctx->nshards; i++) {
    offloadEventSynchronize(ctx->shards_c_dev[i].event);
    offloadEventRecord(ctx->shards_c_dev[i].event, ctx->shards_c_dev[i].stream);
  }
#endif
}

/*******************************************************************************
 * \brief Internal routine for executing the tasks in given batch on the GPU.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply_gpu_process_batch(const dbm_batch_t *batch,
                                    const double alpha, dbm_shard_t *shard_c,
                                    const int kshard, const bool finish,
                                    dbm_multiply_gpu_context_t *ctx) {
  // Assume GPU device was activated earlier.
  dbm_shard_gpu_t *const shard_g = &ctx->shards_c_dev[kshard];
  const int ntasks = batch->ntasks;
  double *old_data_dev = NULL;
  bool grown = false; // an empty shard starts without a buffer (NULL)
#if defined(__OFFLOAD_UNIFIED_MEMORY)
  const dbm_task_t *const batch_dev = batch->tasks; // read where written
#else
  dbm_task_t *const batch_dev = &ctx->batches_dev[kshard * ctx->max_batch_size];
#endif

  // Blocks promised by batches computed on the host grow the shard at finish.
  if (0 < ntasks || finish) {
    cp_assert(NULL != shard_c && NULL != shard_g);

#if !defined(__OFFLOAD_UNIFIED_MEMORY)
    if (0 < ntasks) { // Upload new batch.
      const size_t size = ntasks * sizeof(dbm_task_t);
      offloadMemcpyAsyncHtoD(batch_dev, batch->tasks, size, shard_g->stream);
    }
#endif

    // Reallocate shard_g->data if necessary.
    if (shard_c->data_promised > shard_g->data_allocated) {
      shard_g->data_allocated = DBM_ALLOCATION_FACTOR * shard_c->data_promised;
      cp_assert(shard_c->data_promised <= shard_g->data_allocated);
      old_data_dev = shard_g->data;
      shard_g->data = offload_mempool_device_malloc(shard_g->data_allocated *
                                                    sizeof(double));
      // Omit to wait for copy before freeing old buffer.
      offloadMemcpyAsyncDtoD(shard_g->data, old_data_dev,
                             shard_g->data_size * sizeof(double),
                             shard_g->stream);
      grown = true;
    }
    if (0 < ntasks || NULL != old_data_dev) {
      offloadEventRecord(shard_g->event, shard_g->stream);
    }

    // Zero the headroom once per allocation rather than new blocks per batch:
    // a memset is mostly launch latency, and kernels write promised blocks.
    if (grown) {
      const size_t headroom = shard_g->data_allocated - shard_g->data_size;
      offloadMemsetAsync(&shard_g->data[shard_g->data_size], 0,
                         headroom * sizeof(double), shard_g->stream);
    }
    if (shard_c->data_promised > shard_g->data_size) {
      shard_g->data_size = shard_c->data_promised; // within zeroed headroom
    }
  }

  if (0 < ntasks) {
    OFFLOAD_CHECK(offloadGetLastError());
    cp_assert(0 != shard_g->data_size);

    // Launch kernel.
    dbm_multiply_gpu_launch_kernel(shard_g->stream, alpha, batch, batch_dev,
                                   ctx->pack_a_dev.data, ctx->pack_b_dev.data,
                                   shard_g->data);
    OFFLOAD_CHECK(offloadGetLastError());
    offloadEventRecord(shard_g->done, shard_g->stream);
    if (0 > ctx->hybrid) {
#pragma omp atomic
      ctx->flops[0] += batch->flops;
    }
#if defined(__OFFLOAD_UNIFIED_MEMORY)
    if (!finish) { // the next batch of the shard refills what the kernel reads
      offloadStreamSynchronize(shard_g->stream);
    }
#endif
  }

  if (finish) { // Start downloading the current shard of matrix_c.
    // Growing the host buffer frees the old one, into which the previous
    // download may still write: batches computed on the host never wait.
    if (shard_g->downloading &&
        shard_c->data_allocated < shard_c->data_promised) {
      offloadEventSynchronize(shard_g->download);
      shard_g->downloading = false;
    }
    // Grow host buffer if necessary.
    dbm_shard_allocate_promised_blocks(shard_c);
    // Download results from device.
    cp_assert(shard_c->data_size == shard_g->data_size);
    offloadMemcpyAsyncDtoH(shard_c->data, shard_g->data,
                           shard_g->data_size * sizeof(double),
                           shard_g->stream);
    offloadEventRecord(shard_g->download, shard_g->stream);
    shard_g->downloading = true;
  }

  if (0 < ntasks || NULL != old_data_dev) {
    // Wait for:
    // - Batch to be uploaded (before refilling it).
    // - Safely freeing device buffer (if resized).
    offloadEventSynchronize(shard_g->event);

    if (NULL != old_data_dev) {
      offload_mempool_device_free(old_data_dev);
    }
  }
}

/*******************************************************************************
 * \brief Internal routine for executing the tasks in given batch on the host
 *        while the GPU is busy with the shard's previous batch (hybrid).
 *        Returns false if the GPU shall process the batch instead.
 * \author Hans Pabst
 ******************************************************************************/
bool dbm_multiply_gpu_process_batch_host(
    const dbm_batch_t *batch, const double alpha, const dbm_pack_t *pack_a,
    const dbm_pack_t *pack_b, dbm_shard_t *shard_c, const int kshard,
    const bool finish, const int cpu_options, dbm_multiply_gpu_context_t *ctx) {
  dbm_shard_gpu_t *const shard_g = &ctx->shards_c_dev[kshard];
  bool result = false;

  // Only generated kernels let the host rival the GPU.
  if (0 != ctx->hybrid && 0 < batch->ntasks &&
      !offloadEventQuery(shard_g->done) &&
      dbm_multiply_cpu_generated(batch->max_m, batch->max_n, batch->max_k,
                                 alpha, cpu_options)) {
    int active;
#pragma omp atomic capture
    active = ctx->hybrid_active++;
    result = (active < abs(ctx->hybrid));
    if (result) {
      // The contributions start from zero and are added to C at the end.
      if (shard_g->host_size < shard_c->data_promised) {
        if (shard_g->host_allocated < shard_c->data_promised) {
          shard_g->host_allocated =
              DBM_ALLOCATION_FACTOR * shard_c->data_promised;
          shard_g->host_data = realloc(
              shard_g->host_data, shard_g->host_allocated * sizeof(double));
          cp_assert(NULL != shard_g->host_data);
        }
        memset(&shard_g->host_data[shard_g->host_size], 0,
               (shard_c->data_promised - shard_g->host_size) * sizeof(double));
        shard_g->host_size = shard_c->data_promised;
      }
      dbm_multiply_cpu_process_tasks(batch, alpha, pack_a, pack_b,
                                     shard_g->host_data, cpu_options);
      if (0 > ctx->hybrid) {
#pragma omp atomic
        ctx->flops[1] += batch->flops;
      }
    }
#pragma omp atomic
    --ctx->hybrid_active;
    if (result && finish) { // download as if the GPU processed the batch
      dbm_batch_t empty = {0};
      dbm_multiply_gpu_process_batch(&empty, alpha, shard_c, kshard, true, ctx);
    }
  }

  return result;
}

/*******************************************************************************
 * \brief Internal routine for shutting down the gpu backend.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply_gpu_stop(dbm_multiply_gpu_context_t *ctx) {
  // Assume GPU device was activated earlier.
  // Wait for completion, then free gpu ressources.
#pragma omp parallel for DBM_OMP_SCHEDULE
  for (int i = 0; i < ctx->nshards; i++) {
    dbm_shard_gpu_t *const shard_g = &ctx->shards_c_dev[i];
    offloadStreamSynchronize(shard_g->stream);
    if (NULL != shard_g->host_data) { // add contributions computed on the host
      double *const data = ctx->shards_c_host[i].data;
      cp_assert(shard_g->host_size <= ctx->shards_c_host[i].data_size);
      for (int j = 0; j < shard_g->host_size; j++) {
        data[j] += shard_g->host_data[j];
      }
      free(shard_g->host_data);
    }
    offloadStreamDestroy(shard_g->stream);
    offloadEventDestroy(shard_g->event);
    offloadEventDestroy(shard_g->done);
    offloadEventDestroy(shard_g->download);
    offload_mempool_device_free(shard_g->data);
  }
  free(ctx->shards_c_dev);

  const int64_t flops = ctx->flops[0] + ctx->flops[1];
  if (0 > ctx->hybrid && 0 < flops) {
    fprintf(stderr, "INFO DBM: %.0f%% of %.1f GFLOP computed on the host\n",
            100.0 * ctx->flops[1] / flops, 1E-9 * flops);
  }

  offload_mempool_device_free(ctx->pack_a_dev.data);
  offload_mempool_device_free(ctx->pack_b_dev.data);
#if !defined(__OFFLOAD_UNIFIED_MEMORY)
  offload_mempool_device_free(ctx->batches_dev);
#endif
  offloadStreamDestroy(ctx->main_stream);
  offloadEventDestroy(ctx->upload_event);
}

#endif // defined(__OFFLOAD) && !defined(__NO_OFFLOAD_DBM)

// EOF
