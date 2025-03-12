/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2025 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <limits.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>

#include "../offload/offload_runtime.h"
#include "dbm_hyperparams.h"
#include "dbm_internal.h"
#include "dbm_library.h"
#include "dbm_multiply.h"
#include "dbm_multiply_comm.h"
#include "dbm_multiply_cpu.h"
#include "dbm_multiply_gpu.h"

#if defined(__LIBXSMM)
#include <libxsmm.h>
#endif

#if !defined(DBM_VALIDATE_AGAINST_LIBXSMM) && 0
#define DBM_VALIDATE_AGAINST_LIBXSMM
#endif

/*******************************************************************************
 * \brief Private routine for computing the max filter threshold for each row.
 * \author Ole Schuett
 ******************************************************************************/
static float *compute_rows_max_eps(const bool trans, const dbm_matrix_t *matrix,
                                   const double filter_eps) {
  const int nrows = (trans) ? matrix->ncols : matrix->nrows;
  int *nblocks_per_row = calloc(nrows, sizeof(int));
  float *row_max_eps = malloc(nrows * sizeof(float));
  assert(row_max_eps != NULL);

#pragma omp parallel
  {
#pragma omp for
    for (int ishard = 0; ishard < dbm_get_num_shards(matrix); ishard++) {
      dbm_shard_t *shard = &matrix->shards[ishard];
      for (int iblock = 0; iblock < shard->nblocks; iblock++) {
        const dbm_block_t *blk = &shard->blocks[iblock];
        const int row = (trans) ? blk->col : blk->row;
#pragma omp atomic
        nblocks_per_row[row]++;
      }
    }
#pragma omp single
    dbm_mpi_sum_int(nblocks_per_row, nrows, matrix->dist->comm);
#pragma omp barrier
#pragma omp for
    for (int i = 0; i < nrows; i++) {
      const float f =
          ((float)filter_eps) / ((float)imax(1, nblocks_per_row[i]));
      row_max_eps[i] = f * f;
    }
  } // end of omp parallel region

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
 * \brief Private routine for initializing the multiplication backend.
 * \author Ole Schuett
 ******************************************************************************/
static backend_context_t *backend_start(const dbm_matrix_t *matrix_c) {
  backend_context_t *ctx = calloc(1, sizeof(backend_context_t));

#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_DBM)
  dbm_multiply_gpu_start(DBM_MAX_BATCH_SIZE, dbm_get_num_shards(matrix_c),
                         matrix_c->shards, &ctx->gpu);
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
  (void)pack_a; // mark as used
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
  dbm_multiply_gpu_process_batch(ntasks, batch, alpha, kshard, &ctx->gpu);
#if defined(DBM_VALIDATE_AGAINST_LIBXSMM) && defined(__LIBXSMM)
  dbm_shard_gpu_t *const shard_g = &ctx->gpu.shards_c_dev[kshard];
  dbm_shard_t shard_r;
  dbm_shard_allocate_promised_blocks(shard_c);
  /* start transferring GPU result to host */
  assert(shard_c->data_size == shard_g->data_size);
  dbm_shard_init(&shard_r);
  dbm_shard_copy(&shard_r, shard_c);
  offloadMemcpyAsyncDtoH(shard_c->data, shard_g->data,
                         shard_c->data_size * sizeof(double), shard_g->stream);
  dbm_multiply_cpu_process_batch(ntasks, batch, alpha, pack_a, pack_b,
                                 &shard_r);
  /* finish transferring GPU result to host */
  offloadStreamSynchronize(shard_g->stream);
  libxsmm_matdiff_info diff;
  libxsmm_matdiff_clear(&diff);
  for (int itask = 0; itask < ntasks; ++itask) {
    const dbm_task_t task = batch[itask];
    const double *const tst = &shard_c->data[task.offset_c];
    const double *const ref = &shard_r.data[task.offset_c];
    libxsmm_matdiff_info d;
    if (EXIT_SUCCESS == libxsmm_matdiff(&d, LIBXSMM_DATATYPE(double), task.m,
                                        task.n, ref, tst, NULL /*ldref*/,
                                        NULL /*ldtst*/)) {
      libxsmm_matdiff_reduce(&diff, &d);
    }
  }
  const char *const maxeps_env = getenv("DBM_MULTIPLY_MAXEPS");
  const double maxeps = (NULL == maxeps_env ? 1E-13 : fabs(atof(maxeps_env)));
  const double epsilon = libxsmm_matdiff_epsilon(&diff);
  if (maxeps < epsilon) {
    fprintf(stderr, "INFO ACC/LIBDBM: diff=%g\n", epsilon);
  }
  dbm_shard_release(&shard_r);
#else
  (void)pack_a;
  (void)pack_b;
  (void)shard_c; // mark as used
#endif
#else
  (void)kshard;
  (void)ctx; // mark as used
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

  const int nshard_rows = matrix_c->dist->rows.nshards;
  const int nshard_cols = matrix_c->dist->cols.nshards;
  int *shard_row_start = calloc(nshard_rows, sizeof(int));
  int *shard_col_start = calloc(nshard_cols, sizeof(int));
  assert(NULL != shard_row_start && NULL != shard_col_start);

  const int *sum_index_sizes_a =
      (transa) ? matrix_a->row_sizes : matrix_a->col_sizes;
  const int *sum_index_sizes_b =
      (transb) ? matrix_b->col_sizes : matrix_b->row_sizes;
  const int *free_index_sizes_a =
      (transa) ? matrix_a->col_sizes : matrix_a->row_sizes;
  const int *free_index_sizes_b =
      (transb) ? matrix_b->row_sizes : matrix_b->col_sizes;

#pragma omp parallel reduction(+ : flop_sum)
  {
    // Thread-private array covering given work in piece-wise fashion.
    dbm_task_t *batch =
#if (201811 /*v5.0*/ <= _OPENMP)
        omp_alloc(sizeof(dbm_task_t) * DBM_MAX_BATCH_SIZE, omp_null_allocator);
#else
        malloc(sizeof(dbm_task_t) * DBM_MAX_BATCH_SIZE);
#endif
    assert(NULL != batch);

    // Blocks are ordered first by shard. Creating lookup tables of boundaries.
#pragma omp for nowait
    for (int iblock = 1; iblock < pack_a->nblocks; iblock++) {
      const int shard_row = pack_a->blocks[iblock].free_index % nshard_rows;
      const int prev_shard_row =
          pack_a->blocks[iblock - 1].free_index % nshard_rows;
      if (prev_shard_row != shard_row) {
        shard_row_start[shard_row] = iblock;
      }
    }
#pragma omp for
    for (int jblock = 1; jblock < pack_b->nblocks; jblock++) {
      const int shard_col = pack_b->blocks[jblock].free_index % nshard_cols;
      const int prev_shard_col =
          pack_b->blocks[jblock - 1].free_index % nshard_cols;
      if (prev_shard_col != shard_col) {
        shard_col_start[shard_col] = jblock;
      }
    }

#pragma omp for collapse(2) DBM_OMP_SCHEDULE
    for (int shard_row = 0; shard_row < nshard_rows; shard_row++) {
      for (int shard_col = 0; shard_col < nshard_cols; shard_col++) {
        const int ishard = shard_row * nshard_cols + shard_col;
        dbm_shard_t *shard_c = &matrix_c->shards[ishard];
        int ntasks = 0;

        // Use a merge-join to find pairs of blocks with matching sum indices.
        // This utilizes that blocks within a shard are ordered by sum_index.
        const int iblock_start = shard_row_start[shard_row];
        int jblock_start = shard_col_start[shard_col];
        for (int iblock = iblock_start; iblock < pack_a->nblocks; iblock++) {
          const dbm_pack_block_t *blk_a = &pack_a->blocks[iblock];
          if (blk_a->free_index % nshard_rows != shard_row) {
            break;
          }
          for (int jblock = jblock_start; jblock < pack_b->nblocks; jblock++) {
            const dbm_pack_block_t *blk_b = &pack_b->blocks[jblock];
            if (blk_b->free_index % nshard_cols != shard_col) {
              break;
            }
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
            const int row = blk_a->free_index, col = blk_b->free_index;
            dbm_block_t *blk_c = dbm_shard_lookup(shard_c, row, col);
            if (blk_c == NULL && retain_sparsity) {
              continue;
            } else if (blk_c == NULL) {
              assert(dbm_get_shard_index(matrix_c, row, col) == ishard);
              assert(dbm_get_stored_coordinates(matrix_c, row, col) ==
                     matrix_c->dist->my_rank);
              blk_c = dbm_shard_promise_new_block(shard_c, row, col, m * n);
            }

            // Count flops.
            const int64_t task_flops = 2LL * m * n * k;
            if (task_flops == 0) {
              continue;
            }
            flop_sum += task_flops;
            dbm_library_counter_increment(m, n, k);

            // Add block multiplication to batch.
            batch[ntasks].m = m;
            batch[ntasks].n = n;
            batch[ntasks].k = k;
            batch[ntasks].offset_a = blk_a->offset;
            batch[ntasks].offset_b = blk_b->offset;
            batch[ntasks].offset_c = blk_c->offset;
            ++ntasks;

            if (ntasks == DBM_MAX_BATCH_SIZE) {
              backend_process_batch(ntasks, batch, alpha, pack_a, pack_b,
                                    ishard, shard_c, ctx);
              ntasks = 0;
            }
          }
        }
        backend_process_batch(ntasks, batch, alpha, pack_a, pack_b, ishard,
                              shard_c, ctx);
      }
    }

#if (201811 /*v5.0*/ <= _OPENMP)
    omp_free(batch, omp_null_allocator);
#else
    free(batch);
#endif
  }

  free(shard_row_start);
  free(shard_col_start);

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

  // Wait for all other MPI ranks to complete, then release ressources.
  dbm_comm_iterator_stop(iter);
  free(rows_max_eps);
  backend_stop(ctx);

  // Final filter pass.
  dbm_filter(matrix_c, filter_eps);
}

// EOF
