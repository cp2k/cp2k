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
                                  const bool transa, const bool transb,
                                  const double alpha, const dbm_pack_t *pack_a,
                                  const dbm_pack_t *pack_b, const int kshard,
                                  dbm_shard_t *shard_c,
                                  backend_context_t *ctx) {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_DBM)
  (void)pack_a; // mark as used
  (void)pack_b;
  (void)shard_c;
  dbm_multiply_gpu_process_batch(ntasks, batch, transa, transb, alpha, kshard,
                                 &ctx->gpu);
#else
  (void)kshard; // mark as used
  (void)ctx;
  dbm_multiply_cpu_process_batch(ntasks, batch, transa, transb, alpha, pack_a,
                                 pack_b, shard_c);
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
                           const float *rows_left_max_eps, int64_t *flop,
                           backend_context_t *ctx) {

  const float alpha2 = alpha * alpha;

  int64_t flop_sum = 0;
#pragma omp parallel for schedule(dynamic) reduction(+ : flop_sum)
  for (int kshard = 0; kshard < matrix_c->nshards; kshard++) {
    dbm_shard_t *shard_c = &matrix_c->shards[kshard];

    dbm_task_t batch[MAX_BATCH_SIZE];
    int ntasks = 0;

    // Essentially a merge sort (assuming blocks are pre-sorted by shared index)
    int jblock_start = 0;
    for (int iblock = 0; iblock < pack_a->nblocks; iblock++) {
      const dbm_block_t *blk_a = &pack_a->blocks[iblock];
      const int row_left = (transa) ? blk_a->col : blk_a->row;
      const int col_left = (transa) ? blk_a->row : blk_a->col;
      if (row_left % matrix_c->nshards != kshard) {
        continue;
      }
      for (int jblock = jblock_start; jblock < pack_b->nblocks; jblock++) {
        const dbm_block_t *blk_b = &pack_b->blocks[jblock];
        const int row_right = (transb) ? blk_b->col : blk_b->row;
        const int col_right = (transb) ? blk_b->row : blk_b->col;
        if (col_left < row_right) {
          break;
        }
        if (col_left > row_right) {
          jblock_start++;
          continue;
        }
        // Found block pair with col_left == row_right.

        // Check norms.
        if (alpha2 * blk_a->norm * blk_b->norm < rows_left_max_eps[row_left]) {
          continue;
        }

        // Check block sizes.
        const int row_size_a = matrix_a->row_sizes[blk_a->row];
        const int col_size_a = matrix_a->col_sizes[blk_a->col];
        const int row_size_left = (transa) ? col_size_a : row_size_a;
        const int col_size_left = (transa) ? row_size_a : col_size_a;
        const int row_size_b = matrix_b->row_sizes[blk_b->row];
        const int col_size_b = matrix_b->col_sizes[blk_b->col];
        const int row_size_right = (transb) ? col_size_b : row_size_b;
        const int col_size_right = (transb) ? row_size_b : col_size_b;
        const int row_size_c = matrix_c->row_sizes[row_left];
        const int col_size_c = matrix_c->col_sizes[col_right];
        const int m = row_size_left, n = col_size_right, k = col_size_left;
        assert(m == row_size_c);
        assert(n == col_size_c);
        assert(k == row_size_right);

        // Get C block.
        dbm_block_t *blk_c = dbm_shard_lookup(shard_c, row_left, col_right);
        if (blk_c == NULL && retain_sparsity) {
          continue;
        } else if (blk_c == NULL) {
          blk_c =
              dbm_shard_promise_new_block(shard_c, row_left, col_right, m * n);
        }

        // Count flops.
        assert(m * n * k > 0);
        flop_sum += 2 * m * n * k;

        // Invalidate norm of C block because its data is going to change.
        blk_c->norm = -1.0;

        // Add block multiplication to batch.
        batch[ntasks].m = m;
        batch[ntasks].n = n;
        batch[ntasks].k = k;
        batch[ntasks].offset_a = blk_a->offset;
        batch[ntasks].offset_b = blk_b->offset;
        batch[ntasks].offset_c = blk_c->offset;
        ntasks++;

        if (ntasks == MAX_BATCH_SIZE) {
          backend_process_batch(ntasks, batch, transa, transb, alpha, pack_a,
                                pack_b, kshard, shard_c, ctx);
          ntasks = 0;
        }
      }
    }
    backend_process_batch(ntasks, batch, transa, transb, alpha, pack_a, pack_b,
                          kshard, shard_c, ctx);
  }
  *flop += flop_sum;
}

/*******************************************************************************
 * \brief Performs a multiplication of two dbm_matrix_t matrices.
 *        See dbm_matrix.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply(const bool transa, const bool transb, const double alpha,
                  dbm_matrix_t *matrix_a, dbm_matrix_t *matrix_b,
                  const double beta, dbm_matrix_t *matrix_c,
                  const bool retain_sparsity, const double filter_eps,
                  int64_t *flop) {

  assert(omp_get_num_threads() == 1);

  // Denote left/right to matrices a/b after possible transpose.
  const int nrows_left = (transa) ? matrix_a->ncols : matrix_a->nrows;
  const int ncols_left = (transa) ? matrix_a->nrows : matrix_a->ncols;
  const int nrows_right = (transb) ? matrix_b->ncols : matrix_b->nrows;
  const int ncols_right = (transb) ? matrix_b->nrows : matrix_b->ncols;

  // Sanity check matrix dimensions.
  assert(ncols_left == nrows_right);
  assert(nrows_left == matrix_c->nrows);
  assert(ncols_right == matrix_c->ncols);

  // Prepare matrix_c.
  dbm_scale(matrix_c, beta);

  // Start uploading matrix_c to the GPU.
  backend_context_t *ctx = backend_start(matrix_c);

  // Compute norms and filter thresholds for each row.
  dbm_compute_block_norms(matrix_a);
  dbm_compute_block_norms(matrix_b);
  float *rows_left_max_eps = compute_rows_max_eps(transa, matrix_a, filter_eps);

  // Redistribute matrix_a and matrix_b across MPI ranks.
  dbm_comm_iterator_t *iter =
      dbm_comm_iterator_start(transa, transb, matrix_a, matrix_b, matrix_c);

  // Main loop.
  *flop = 0;
  dbm_pack_t *pack_a, *pack_b;
  while (dbm_comm_iterator_next(iter, &pack_a, &pack_b)) {
    backend_upload_packs(pack_a, pack_b, ctx);
    multiply_packs(transa, transb, alpha, pack_a, pack_b, matrix_a, matrix_b,
                   matrix_c, retain_sparsity, rows_left_max_eps, flop, ctx);
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
  free(rows_left_max_eps);
  backend_stop(ctx);

  // Compute average flops per rank.
  dbm_mpi_sum_int64(flop, 1, matrix_c->dist->comm);
  *flop = (*flop + matrix_c->dist->nranks - 1) / matrix_c->dist->nranks;

  // Final filter pass.
  dbm_filter(matrix_c, filter_eps);
}

// EOF
