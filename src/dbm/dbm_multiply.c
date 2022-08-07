/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../offload/offload_runtime.h"
#include "dbm_hyperparams.h"
#include "dbm_library.h"
#include "dbm_multiply.h"
#include "dbm_multiply_comm.h"
#include "dbm_multiply_cpu.h"
#include "dbm_multiply_gpu.h"
#include "dbm_multiply_internal.h"

/*******************************************************************************
 * \brief Returns the larger of two given integer (missing from the C standard)
 * \author Ole Schuett
 ******************************************************************************/
static inline int imax(int x, int y) { return (x > y ? x : y); }

/*******************************************************************************
 * \brief Private routine for computing the max filter threshold for each row.
 * \author Ole Schuett
 ******************************************************************************/
static float *compute_rows_max_eps(const bool trans, const dbm_matrix_t *matrix,
                                   const double filter_eps) {
  const int nrows = (trans) ? matrix->ncols : matrix->nrows;
  int *nblocks_per_row = calloc(nrows, sizeof(int));
  for (int ishard = 0; ishard < matrix->nshards; ishard++) {
    dbm_shard_t *shard = &matrix->shards[ishard];
    for (int iblock = 0; iblock < shard->nblocks; iblock++) {
      const dbm_block_t *blk = &shard->blocks[iblock];
      const int row = (trans) ? blk->col : blk->row;
      nblocks_per_row[row]++;
    }
  }

  dbm_mpi_sum_int(nblocks_per_row, nrows, matrix->dist->comm);

  float *row_max_eps = malloc(nrows * sizeof(float));
  for (int i = 0; i < nrows; i++) {
    const float f = ((float)filter_eps) / ((float)imax(1, nblocks_per_row[i]));
    row_max_eps[i] = f * f;
  }
  free(nblocks_per_row);
  return row_max_eps; // Ownership of row_max_eps transfers to caller.
}

/*******************************************************************************
 * \brief Private struct for storing the context of the multiplication backend.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_DBM)
  dbm_multiply_gpu_context_t gpu;
#endif
} backend_context_t;

/*******************************************************************************
 * \brief Private routine for intializing the multiplication backend.
 * \author Ole Schuett
 ******************************************************************************/
static backend_context_t *backend_start(const dbm_matrix_t *matrix_c) {
  backend_context_t *ctx = calloc(1, sizeof(backend_context_t));

#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_DBM)
  dbm_multiply_gpu_start(MAX_BATCH_SIZE, matrix_c->nshards, matrix_c->shards,
                         &ctx->gpu);
#else
  (void)matrix_c; // mark as used
#endif

  return ctx;
}

/*******************************************************************************
 * \brief Private routine for handing newly arrived packs to the backend.
 * \author Ole Schuett
 ******************************************************************************/
static void backend_upload_packs(const dbm_pack_t *pack_a,
                                 const dbm_pack_t *pack_b,
                                 backend_context_t *ctx) {

#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_DBM)
  dbm_multiply_gpu_upload_packs(pack_a, pack_b, &ctx->gpu);
#else
  (void)pack_a;   // mark as used
  (void)pack_b;
  (void)ctx;
#endif
}

/*******************************************************************************
 * \brief Private routine for sending a batch to the multiplication backend.
 * \author Ole Schuett
 ******************************************************************************/
static void backend_process_batch(const int ntasks, dbm_task_t batch[ntasks],
                                  const double alpha, const dbm_pack_t *pack_a,
                                  const dbm_pack_t *pack_b, const int kshard,
                                  dbm_shard_t *shard_c,
                                  backend_context_t *ctx) {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_DBM)
  (void)pack_a; // mark as used
  (void)pack_b;
  (void)shard_c;
  dbm_multiply_gpu_process_batch(ntasks, batch, alpha, kshard, &ctx->gpu);
#else
  (void)kshard; // mark as used
  (void)ctx;
  dbm_multiply_cpu_process_batch(ntasks, batch, alpha, pack_a, pack_b, shard_c);
#endif
}

/*******************************************************************************
 * \brief Private routine for downloading results of the multiplication backend.
 * \author Ole Schuett
 ******************************************************************************/
static void backend_download_results(backend_context_t *ctx) {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_DBM)
  dbm_multiply_gpu_download_results(&ctx->gpu);
#else
  (void)ctx; // mark as used
#endif
}

/*******************************************************************************
 * \brief Private routine for shutting down the multiplication backend.
 * \author Ole Schuett
 ******************************************************************************/
static void backend_stop(backend_context_t *ctx) {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_DBM)
  dbm_multiply_gpu_stop(&ctx->gpu);
#endif
  free(ctx);
}

/*******************************************************************************
 * \brief Private routine for multipling two packs.
 * \author Ole Schuett
 ******************************************************************************/
static void multiply_packs(const bool transa, const bool transb,
                           const double alpha, const dbm_pack_t *pack_a,
                           const dbm_pack_t *pack_b,
                           const dbm_matrix_t *matrix_a,
                           const dbm_matrix_t *matrix_b, dbm_matrix_t *matrix_c,
                           const bool retain_sparsity,
                           const float *rows_max_eps, int64_t *flop,
                           backend_context_t *ctx) {

  const float alpha2 = alpha * alpha;

  int64_t flop_sum = 0;

  const int *sum_index_sizes_a =
      (transa) ? matrix_a->row_sizes : matrix_a->col_sizes;
  const int *sum_index_sizes_b =
      (transb) ? matrix_b->col_sizes : matrix_b->row_sizes;
  const int *free_index_sizes_a =
      (transa) ? matrix_a->col_sizes : matrix_a->row_sizes;
  const int *free_index_sizes_b =
      (transb) ? matrix_b->row_sizes : matrix_b->col_sizes;

#pragma omp parallel for schedule(dynamic) reduction(+ : flop_sum)
  for (int kshard = 0; kshard < matrix_c->nshards; kshard++) {
    dbm_shard_t *shard_c = &matrix_c->shards[kshard];

    dbm_task_t batch[MAX_BATCH_SIZE];
    int ntasks = 0;

    // Essentially a merge sort (assuming blocks are pre-sorted by sum_index)
    int jblock_start = 0;
    for (int iblock = 0; iblock < pack_a->nblocks; iblock++) {
      const dbm_pack_block_t *blk_a = &pack_a->blocks[iblock];
      if (blk_a->free_index % matrix_c->nshards != kshard) {
        continue;
      }
      for (int jblock = jblock_start; jblock < pack_b->nblocks; jblock++) {
        const dbm_pack_block_t *blk_b = &pack_b->blocks[jblock];
        if (blk_a->sum_index < blk_b->sum_index) {
          break;
        }
        if (blk_a->sum_index > blk_b->sum_index) {
          jblock_start++;
          continue;
        }
        // Found block pair with blk_a->sum_index == blk_b->sum_index.

        // Check norms.
        const float result_norm = alpha2 * blk_a->norm * blk_b->norm;
        if (result_norm < rows_max_eps[blk_a->free_index]) {
          continue;
        }

        // Check block sizes.
        const int m = free_index_sizes_a[blk_a->free_index];
        const int n = free_index_sizes_b[blk_b->free_index];
        const int k = sum_index_sizes_a[blk_a->sum_index];
        assert(m == matrix_c->row_sizes[blk_a->free_index]);
        assert(n == matrix_c->col_sizes[blk_b->free_index]);
        assert(k == sum_index_sizes_b[blk_b->sum_index]);

        // Get C block.
        dbm_block_t *blk_c =
            dbm_shard_lookup(shard_c, blk_a->free_index, blk_b->free_index);
        if (blk_c == NULL && retain_sparsity) {
          continue;
        } else if (blk_c == NULL) {
          blk_c = dbm_shard_promise_new_block(shard_c, blk_a->free_index,
                                              blk_b->free_index, m * n);
        }

        // Count flops.
        assert(m * n * k > 0);
        flop_sum += 2 * m * n * k;
        dbm_library_counter_increment(m, n, k);

        // Add block multiplication to batch.
        batch[ntasks].m = m;
        batch[ntasks].n = n;
        batch[ntasks].k = k;
        batch[ntasks].offset_a = blk_a->offset;
        batch[ntasks].offset_b = blk_b->offset;
        batch[ntasks].offset_c = blk_c->offset;
        ntasks++;

        if (ntasks == MAX_BATCH_SIZE) {
          backend_process_batch(ntasks, batch, alpha, pack_a, pack_b, kshard,
                                shard_c, ctx);
          ntasks = 0;
        }
      }
    }
    backend_process_batch(ntasks, batch, alpha, pack_a, pack_b, kshard, shard_c,
                          ctx);
  }
  *flop += flop_sum;
}

/*******************************************************************************
 * \brief Performs a multiplication of two dbm_matrix_t matrices.
 *        See dbm_matrix.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply(const bool transa, const bool transb, const double alpha,
                  const dbm_matrix_t *matrix_a, const dbm_matrix_t *matrix_b,
                  const double beta, dbm_matrix_t *matrix_c,
                  const bool retain_sparsity, const double filter_eps,
                  int64_t *flop) {

  assert(omp_get_num_threads() == 1);

  // Throughout the matrix multiplication code the "sum_index" and "free_index"
  // denote the summation (aka dummy) and free index from the Einstein notation.
  const int num_sum_index_a = (transa) ? matrix_a->nrows : matrix_a->ncols;
  const int num_sum_index_b = (transb) ? matrix_b->ncols : matrix_b->nrows;
  const int num_free_index_a = (transa) ? matrix_a->ncols : matrix_a->nrows;
  const int num_free_index_b = (transb) ? matrix_b->nrows : matrix_b->ncols;

  // Sanity check matrix dimensions.
  assert(num_sum_index_a == num_sum_index_b);
  assert(num_free_index_a == matrix_c->nrows);
  assert(num_free_index_b == matrix_c->ncols);

  // Prepare matrix_c.
  dbm_scale(matrix_c, beta);

  // Start uploading matrix_c to the GPU.
  backend_context_t *ctx = backend_start(matrix_c);

  // Compute filter thresholds for each row.
  float *rows_max_eps = compute_rows_max_eps(transa, matrix_a, filter_eps);

  // Redistribute matrix_a and matrix_b across MPI ranks.
  dbm_comm_iterator_t *iter =
      dbm_comm_iterator_start(transa, transb, matrix_a, matrix_b, matrix_c);

  // Main loop.
  *flop = 0;
  dbm_pack_t *pack_a, *pack_b;
  while (dbm_comm_iterator_next(iter, &pack_a, &pack_b)) {
    backend_upload_packs(pack_a, pack_b, ctx);
    multiply_packs(transa, transb, alpha, pack_a, pack_b, matrix_a, matrix_b,
                   matrix_c, retain_sparsity, rows_max_eps, flop, ctx);
  }

  // Start downloading matrix_c from the GPU.
  backend_download_results(ctx);

  // Sanity check if matrix_c contains the correct blocks.
  for (int ishard = 0; ishard < matrix_c->nshards; ishard++) {
    dbm_shard_t *shard = &matrix_c->shards[ishard];
    for (int iblock = 0; iblock < shard->nblocks; iblock++) {
      const dbm_block_t *blk = &shard->blocks[iblock];
      const int rank = dbm_get_stored_coordinates(matrix_c, blk->row, blk->col);
      assert(rank == matrix_c->dist->my_rank);
    }
  }

  // Wait for all other MPI ranks to complete, then release ressources.
  dbm_comm_iterator_stop(iter);
  free(rows_max_eps);
  backend_stop(ctx);

  // Compute average flops per rank.
  dbm_mpi_sum_int64(flop, 1, matrix_c->dist->comm);
  *flop = (*flop + matrix_c->dist->nranks - 1) / matrix_c->dist->nranks;

  // Final filter pass.
  dbm_filter(matrix_c, filter_eps);
}

// EOF
