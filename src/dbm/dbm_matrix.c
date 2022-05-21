/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "dbm_hyperparams.h"
#include "dbm_matrix.h"
#include "dbm_mpi.h"

/*******************************************************************************
 * \brief Creates a new matrix.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_create(dbm_matrix_t **matrix_out, dbm_distribution_t *dist,
                const char name[], const int nrows, const int ncols,
                const int row_sizes[nrows], const int col_sizes[ncols]) {
  assert(omp_get_num_threads() == 1);

  dbm_matrix_t *matrix = calloc(1, sizeof(dbm_matrix_t));

  assert(dist->rows.length == nrows);
  assert(dist->cols.length == ncols);
  dbm_distribution_hold(dist);
  matrix->dist = dist;

  size_t size = (strlen(name) + 1) * sizeof(char);
  matrix->name = malloc(size);
  memcpy(matrix->name, name, size);

  matrix->nrows = nrows;
  matrix->ncols = ncols;

  size = nrows * sizeof(int);
  matrix->row_sizes = malloc(size);
  memcpy(matrix->row_sizes, row_sizes, size);

  size = ncols * sizeof(int);
  matrix->col_sizes = malloc(size);
  memcpy(matrix->col_sizes, col_sizes, size);

  matrix->nshards = SHARDS_PER_THREAD * omp_get_max_threads();
  matrix->shards = malloc(matrix->nshards * sizeof(dbm_shard_t));
#pragma omp parallel for
  for (int ishard = 0; ishard < matrix->nshards; ishard++) {
    dbm_shard_init(&matrix->shards[ishard]);
  }

  assert(*matrix_out == NULL);
  *matrix_out = matrix;
}

/*******************************************************************************
 * \brief Releases a matrix and all its ressources.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_release(dbm_matrix_t *matrix) {
  assert(omp_get_num_threads() == 1);
  dbm_distribution_release(matrix->dist);
  free(matrix->name);
  free(matrix->row_sizes);
  free(matrix->col_sizes);
  for (int ishard = 0; ishard < matrix->nshards; ishard++) {
    dbm_shard_release(&matrix->shards[ishard]);
  }
  free(matrix->shards);
  free(matrix);
}

/*******************************************************************************
 * \brief Copies content of matrix_b into matrix_a.
 *        Matrices must have the same row/col block sizes and distribution.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_copy(dbm_matrix_t *matrix_a, const dbm_matrix_t *matrix_b) {
  assert(omp_get_num_threads() == 1);

  assert(matrix_b->nrows == matrix_a->nrows);
  for (int i = 0; i < matrix_b->nrows; i++) {
    assert(matrix_b->row_sizes[i] == matrix_a->row_sizes[i]);
  }
  assert(matrix_b->ncols == matrix_a->ncols);
  for (int i = 0; i < matrix_b->ncols; i++) {
    assert(matrix_b->col_sizes[i] == matrix_a->col_sizes[i]);
  }

  assert(matrix_a->nshards == matrix_b->nshards);
  assert(matrix_a->dist == matrix_b->dist);

#pragma omp parallel for schedule(dynamic)
  for (int ishard = 0; ishard < matrix_a->nshards; ishard++) {
    dbm_shard_t *shard_a = &matrix_a->shards[ishard];
    const dbm_shard_t *shard_b = &matrix_b->shards[ishard];

    free(shard_a->blocks);
    shard_a->nblocks = shard_b->nblocks;
    shard_a->nblocks_allocated = shard_b->nblocks_allocated;
    shard_a->blocks = malloc(shard_b->nblocks_allocated * sizeof(dbm_block_t));
    memcpy(shard_a->blocks, shard_b->blocks,
           shard_b->nblocks * sizeof(dbm_block_t));

    free(shard_a->hashtable);
    shard_a->hashtable_size = shard_b->hashtable_size;
    shard_a->hashtable = malloc(shard_b->hashtable_size * sizeof(int));
    memcpy(shard_a->hashtable, shard_b->hashtable,
           shard_b->hashtable_size * sizeof(int));

    free(shard_a->data);
    shard_a->data_allocated = shard_b->data_allocated;
    shard_a->data = malloc(shard_b->data_allocated * sizeof(double));
    shard_a->data_size = shard_b->data_size;
    memcpy(shard_a->data, shard_b->data, shard_b->data_size * sizeof(double));
  }
}

/*******************************************************************************
 * \brief Copies content of matrix_b into matrix_a.
 *        Matrices may have different distributions.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_redistribute(const dbm_matrix_t *matrix, dbm_matrix_t *redist) {
  assert(omp_get_num_threads() == 1);
  assert(matrix->nrows == redist->nrows);
  for (int i = 0; i < matrix->nrows; i++) {
    assert(matrix->row_sizes[i] == redist->row_sizes[i]);
  }
  assert(matrix->ncols == redist->ncols);
  for (int i = 0; i < matrix->ncols; i++) {
    assert(matrix->col_sizes[i] == redist->col_sizes[i]);
  }

  assert(dbm_mpi_comms_are_similar(matrix->dist->comm, redist->dist->comm));
  const dbm_mpi_comm_t comm = redist->dist->comm;
  const int nranks = dbm_mpi_comm_size(comm);

  // 1st pass: Compute send_count.
  int send_count[nranks];
  memset(send_count, 0, nranks * sizeof(int));
  for (int ishard = 0; ishard < matrix->nshards; ishard++) {
    dbm_shard_t *shard = &matrix->shards[ishard];
    for (int iblock = 0; iblock < shard->nblocks; iblock++) {
      const dbm_block_t *blk = &shard->blocks[iblock];
      const int row_size = matrix->row_sizes[blk->row];
      const int col_size = matrix->col_sizes[blk->col];
      const int block_size = row_size * col_size;
      const int rank = dbm_get_stored_coordinates(redist, blk->row, blk->col);
      assert(0 <= rank && rank < nranks);
      send_count[rank] += 2 + block_size;
    }
  }

  // 1st communication: Exchange counts.
  int recv_count[nranks];
  dbm_mpi_alltoall_int(send_count, 1, recv_count, 1, comm);

  // Compute displacements and allocate data buffers.
  int send_displ[nranks + 1], recv_displ[nranks + 1];
  send_displ[0] = recv_displ[0] = 0;
  for (int irank = 1; irank < nranks + 1; irank++) {
    send_displ[irank] = send_displ[irank - 1] + send_count[irank - 1];
    recv_displ[irank] = recv_displ[irank - 1] + recv_count[irank - 1];
  }
  const int total_send_count = send_displ[nranks];
  const int total_recv_count = recv_displ[nranks];
  double *data_send = malloc(total_send_count * sizeof(double));
  double *data_recv = malloc(total_recv_count * sizeof(double));

  // 2nd pass: Fill send_data.
  int send_data_positions[nranks];
  memcpy(send_data_positions, send_displ, nranks * sizeof(int));
  for (int ishard = 0; ishard < matrix->nshards; ishard++) {
    dbm_shard_t *shard = &matrix->shards[ishard];
    for (int iblock = 0; iblock < shard->nblocks; iblock++) {
      const dbm_block_t *blk = &shard->blocks[iblock];
      const double *blk_data = &shard->data[blk->offset];
      const int row_size = matrix->row_sizes[blk->row];
      const int col_size = matrix->col_sizes[blk->col];
      const int block_size = row_size * col_size;
      const int rank = dbm_get_stored_coordinates(redist, blk->row, blk->col);
      const int pos = send_data_positions[rank];
      data_send[pos + 0] = blk->row; // send integers as doubles
      data_send[pos + 1] = blk->col;
      memcpy(&data_send[pos + 2], blk_data, block_size * sizeof(double));
      send_data_positions[rank] += 2 + block_size;
    }
  }
  for (int irank = 0; irank < nranks; irank++) {
    assert(send_data_positions[irank] == send_displ[irank + 1]);
  }

  // 2nd communication: Exchange data.
  dbm_mpi_alltoallv_double(data_send, send_count, send_displ, data_recv,
                           recv_count, recv_displ, comm);
  free(data_send);

  // 3rd pass: Unpack data.
  dbm_clear(redist);
  int recv_data_pos = 0;
  while (recv_data_pos < total_recv_count) {
    const int row = (int)data_recv[recv_data_pos + 0];
    const int col = (int)data_recv[recv_data_pos + 1];
    assert(data_recv[recv_data_pos + 0] - (double)row == 0.0);
    assert(data_recv[recv_data_pos + 1] - (double)col == 0.0);
    dbm_put_block(redist, row, col, false, &data_recv[recv_data_pos + 2]);
    const int row_size = matrix->row_sizes[row];
    const int col_size = matrix->col_sizes[col];
    const int block_size = row_size * col_size;
    recv_data_pos += 2 + block_size;
  }
  assert(recv_data_pos == total_recv_count);
  free(data_recv);
}

/*******************************************************************************
 * \brief Looks up a block from given matrics. This routine is thread-safe.
 *        If the block is not found then a null pointer is returned.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_get_block_p(dbm_matrix_t *matrix, const int row, const int col,
                     double **block, int *row_size, int *col_size) {
  assert(0 <= row && row < matrix->nrows);
  assert(0 <= col && col < matrix->ncols);
  assert(dbm_get_stored_coordinates(matrix, row, col) == matrix->dist->my_rank);
  *row_size = matrix->row_sizes[row];
  *col_size = matrix->col_sizes[col];
  *block = NULL;

  const int ishard = row % matrix->nshards;
  const dbm_shard_t *shard = &matrix->shards[ishard];
  dbm_block_t *blk = dbm_shard_lookup(shard, row, col);
  if (blk != NULL) {
    blk->norm = -1.0; // Invalidate norm because caller might modify block data.
    *block = &shard->data[blk->offset];
  }
}

/*******************************************************************************
 * \brief Adds a block to given matrix. This routine is thread-safe.
 *        If block already exist then it gets overwritten (or summed).
 * \author Ole Schuett
 ******************************************************************************/
void dbm_put_block(dbm_matrix_t *matrix, const int row, const int col,
                   const bool summation, const double *block) {
  assert(0 <= row && row < matrix->nrows);
  assert(0 <= col && col < matrix->ncols);
  assert(dbm_get_stored_coordinates(matrix, row, col) == matrix->dist->my_rank);
  const int row_size = matrix->row_sizes[row];
  const int col_size = matrix->col_sizes[col];
  const int block_size = row_size * col_size;

  const int ishard = row % matrix->nshards;
  dbm_shard_t *shard = &matrix->shards[ishard];
  omp_set_lock(&shard->lock);
  dbm_block_t *blk =
      dbm_shard_get_or_allocate_block(shard, row, col, block_size);
  double *blk_data = &shard->data[blk->offset];
  if (summation) {
    for (int i = 0; i < block_size; i++) {
      blk_data[i] += block[i];
    }
  } else {
    memcpy(blk_data, block, block_size * sizeof(double));
  }
  blk->norm = -1.0; // Invalidate norm because block data has changed.
  omp_unset_lock(&shard->lock);
}

/*******************************************************************************
 * \brief Remove all blocks from matrix, but does not release underlying memory.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_clear(dbm_matrix_t *matrix) {
  assert(omp_get_num_threads() == 1);

#pragma omp parallel for schedule(dynamic)
  for (int ishard = 0; ishard < matrix->nshards; ishard++) {
    dbm_shard_t *shard = &matrix->shards[ishard];
    shard->nblocks = 0;
    shard->data_size = 0;
    shard->data_promised = 0;
    // Does not deallocate memory, hence data_allocated remains unchanged.
    memset(shard->hashtable, 0, shard->hashtable_size * sizeof(int));
  }
}

/*******************************************************************************
 * \brief Internal routine for re-computing any invalid block norms.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_compute_block_norms(dbm_matrix_t *matrix) {
  assert(omp_get_num_threads() == 1);

#pragma omp parallel for schedule(dynamic)
  for (int ishard = 0; ishard < matrix->nshards; ishard++) {
    dbm_shard_t *shard = &matrix->shards[ishard];
    for (int iblock = 0; iblock < shard->nblocks; iblock++) {
      dbm_block_t *blk = &shard->blocks[iblock];
      if (blk->norm < 0.0) { // negative norms are invalid
        const double *blk_data = &shard->data[blk->offset];
        const int row_size = matrix->row_sizes[blk->row];
        const int col_size = matrix->col_sizes[blk->col];
        double norm = 0.0; // Compute as double ...
        for (int i = 0; i < row_size * col_size; i++) {
          norm += blk_data[i] * blk_data[i];
        }
        blk->norm = (float)norm; // ...store as float.
      }
    }
  }
}

/*******************************************************************************
 * \brief Removes all blocks from the matrix whose norm is below the threshold.
 *        Blocks of size zero are always kept.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_filter(dbm_matrix_t *matrix, const double eps) {
  assert(omp_get_num_threads() == 1);

  if (eps == 0.0) {
    return;
  }

  dbm_compute_block_norms(matrix);

#pragma omp parallel for schedule(dynamic)
  for (int ishard = 0; ishard < matrix->nshards; ishard++) {
    dbm_shard_t *shard = &matrix->shards[ishard];
    const int old_nblocks = shard->nblocks;
    shard->nblocks = 0;
    shard->data_promised = 0;
    memset(shard->hashtable, 0, shard->hashtable_size * sizeof(int));

    for (int iblock = 0; iblock < old_nblocks; iblock++) {
      const dbm_block_t old_blk = shard->blocks[iblock];
      const int row_size = matrix->row_sizes[old_blk.row];
      const int col_size = matrix->col_sizes[old_blk.col];
      const int block_size = row_size * col_size;
      // For historic reasons zero-sized blocks are never filtered.
      if (sqrt(old_blk.norm) < eps && block_size > 0) {
        continue; // filter the block
      }
      // Re-create block.
      dbm_block_t *new_blk = dbm_shard_promise_new_block(
          shard, old_blk.row, old_blk.col, block_size);
      new_blk->norm = old_blk.norm;
      assert(new_blk->offset <= old_blk.offset);
      if (new_blk->offset != old_blk.offset) {
        // Using memmove instead of memcpy because it handles overlap correctly.
        const double *old_blk_data = &shard->data[old_blk.offset];
        double *new_blk_data = &shard->data[new_blk->offset];
        memmove(new_blk_data, old_blk_data, block_size * sizeof(double));
      }
    }
    shard->data_size = shard->data_promised;
    // TODO: Could call realloc to release excess memory.
  }
}

/*******************************************************************************
 * \brief Adds list of blocks efficiently. The blocks will be filled with zeros.
 *        This routine must always be called within an OpenMP parallel region.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_reserve_blocks(dbm_matrix_t *matrix, const int nblocks,
                        const int rows[], const int cols[]) {
  assert(omp_get_num_threads() == omp_get_max_threads() &&
         "Please call dbm_reserve_blocks within an OpenMP parallel region.");
  const int my_rank = matrix->dist->my_rank;
  for (int i = 0; i < nblocks; i++) {
    const int ishard = rows[i] % matrix->nshards;
    dbm_shard_t *shard = &matrix->shards[ishard];
    omp_set_lock(&shard->lock);
    assert(0 <= rows[i] && rows[i] < matrix->nrows);
    assert(0 <= cols[i] && cols[i] < matrix->ncols);
    assert(dbm_get_stored_coordinates(matrix, rows[i], cols[i]) == my_rank);
    const int row_size = matrix->row_sizes[rows[i]];
    const int col_size = matrix->col_sizes[cols[i]];
    const int block_size = row_size * col_size;
    dbm_shard_get_or_promise_block(shard, rows[i], cols[i], block_size);
    omp_unset_lock(&shard->lock);
  }
#pragma omp barrier

#pragma omp for schedule(dynamic)
  for (int ishard = 0; ishard < matrix->nshards; ishard++) {
    dbm_shard_t *shard = &matrix->shards[ishard];
    dbm_shard_allocate_promised_blocks(shard);
  }
}

/*******************************************************************************
 * \brief Multiplies all entries in the given matrix by the given factor alpha.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_scale(dbm_matrix_t *matrix, const double alpha) {
  assert(omp_get_num_threads() == 1);
  if (alpha == 1.0) {
    return;
  }
  if (alpha == 0.0) {
    dbm_zero(matrix);
    return;
  }

#pragma omp parallel for schedule(dynamic)
  for (int ishard = 0; ishard < matrix->nshards; ishard++) {
    dbm_shard_t *shard = &matrix->shards[ishard];
    for (int i = 0; i < shard->data_size; i++) {
      shard->data[i] *= alpha;
    }
    const double alpha2 = alpha * alpha;
    for (int iblock = 0; iblock < shard->nblocks; iblock++) {
      dbm_block_t *blk = &shard->blocks[iblock];
      blk->norm *= alpha2;
    }
  }
}

/*******************************************************************************
 * \brief Sets all blocks in the given matrix to zero.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_zero(dbm_matrix_t *matrix) {
  assert(omp_get_num_threads() == 1);

#pragma omp parallel for schedule(dynamic)
  for (int ishard = 0; ishard < matrix->nshards; ishard++) {
    dbm_shard_t *shard = &matrix->shards[ishard];
    memset(shard->data, 0, shard->data_size * sizeof(double));
    for (int iblock = 0; iblock < shard->nblocks; iblock++) {
      shard->blocks[iblock].norm = 0.0;
    }
  }
}

/*******************************************************************************
 * \brief Adds matrix_b to matrix_a.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_add(dbm_matrix_t *matrix_a, const dbm_matrix_t *matrix_b) {
  assert(omp_get_num_threads() == 1);
  assert(matrix_a->nshards == matrix_b->nshards);
  assert(matrix_a->dist == matrix_b->dist);

#pragma omp parallel for schedule(dynamic)
  for (int ishard = 0; ishard < matrix_b->nshards; ishard++) {
    dbm_shard_t *shard_a = &matrix_a->shards[ishard];
    const dbm_shard_t *shard_b = &matrix_b->shards[ishard];
    for (int iblock = 0; iblock < shard_b->nblocks; iblock++) {
      const dbm_block_t blk_b = shard_b->blocks[iblock];

      const int row_size = matrix_b->row_sizes[blk_b.row];
      const int col_size = matrix_b->col_sizes[blk_b.col];
      assert(row_size == matrix_a->row_sizes[blk_b.row]);
      assert(col_size == matrix_a->col_sizes[blk_b.col]);
      const int block_size = row_size * col_size;
      dbm_block_t *blk_a = dbm_shard_get_or_allocate_block(
          shard_a, blk_b.row, blk_b.col, block_size);
      double *data_a = &shard_a->data[blk_a->offset];
      const double *data_b = &shard_b->data[blk_b.offset];
      for (int i = 0; i < block_size; i++) {
        data_a[i] += data_b[i];
      }
      blk_a->norm = -1.0; // Invalidate norm because block data has changed.
    }
  }
}

/*******************************************************************************
 * \brief Creates an iterator for the blocks of the given matrix.
 *        The iteration order is not stable.
 *        This routine must always be called within an OpenMP parallel region.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_iterator_start(dbm_iterator_t **iter_out, const dbm_matrix_t *matrix) {
  assert(omp_get_num_threads() == omp_get_max_threads() &&
         "Please call dbm_iterator_start within an OpenMP parallel region.");
  dbm_iterator_t *iter = malloc(sizeof(dbm_iterator_t));
  iter->matrix = matrix;
  iter->next_block = 0;
  iter->next_shard = omp_get_thread_num();
  while (iter->next_shard < matrix->nshards &&
         matrix->shards[iter->next_shard].nblocks == 0) {
    iter->next_shard += omp_get_num_threads();
  }
  assert(*iter_out == NULL);
  *iter_out = iter;
}

/*******************************************************************************
 * \brief Returns number of blocks the iterator will provide to calling thread.
 * \author Ole Schuett
 ******************************************************************************/
int dbm_iterator_num_blocks(const dbm_iterator_t *iter) {
  int num_blocks = 0;
  for (int ishard = omp_get_thread_num(); ishard < iter->matrix->nshards;
       ishard += omp_get_num_threads()) {
    num_blocks += iter->matrix->shards[ishard].nblocks;
  }
  return num_blocks;
}

/*******************************************************************************
 * \brief Tests whether the given iterator has any block left.
 * \author Ole Schuett
 ******************************************************************************/
bool dbm_iterator_blocks_left(const dbm_iterator_t *iter) {
  return iter->next_shard < iter->matrix->nshards;
}

/*******************************************************************************
 * \brief Returns the next block from the given iterator.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_iterator_next_block(dbm_iterator_t *iter, int *row, int *col,
                             double **block, int *row_size, int *col_size) {
  const dbm_matrix_t *matrix = iter->matrix;
  assert(iter->next_shard < matrix->nshards);
  const dbm_shard_t *shard = &matrix->shards[iter->next_shard];
  assert(iter->next_block < shard->nblocks);
  dbm_block_t *blk = &shard->blocks[iter->next_block];

  *row = blk->row;
  *col = blk->col;
  *row_size = matrix->row_sizes[blk->row];
  *col_size = matrix->col_sizes[blk->col];
  *block = &shard->data[blk->offset];
  blk->norm = -1.0; // Invalidate norm because caller might modify block data.

  iter->next_block++;
  if (iter->next_block >= shard->nblocks) {
    // Advance to the next non-empty shard...
    iter->next_shard += omp_get_num_threads();
    while (iter->next_shard < matrix->nshards &&
           matrix->shards[iter->next_shard].nblocks == 0) {
      iter->next_shard += omp_get_num_threads();
    }
    iter->next_block = 0; // ...and continue with its first block.
  }
}

/*******************************************************************************
 * \brief Releases the given iterator.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_iterator_stop(dbm_iterator_t *iter) { free(iter); }

/*******************************************************************************
 * \brief Computes a checksum of the given matrix.
 * \author Ole Schuett
 ******************************************************************************/
double dbm_checksum(const dbm_matrix_t *matrix) {
  double checksum = 0.0;
  for (int ishard = 0; ishard < matrix->nshards; ishard++) {
    const dbm_shard_t *shard = &matrix->shards[ishard];
    for (int i = 0; i < shard->data_size; i++) {
      checksum += shard->data[i] * shard->data[i];
    }
  }
  dbm_mpi_sum_double(&checksum, 1, matrix->dist->comm);
  return checksum;
}

/*******************************************************************************
 * \brief Returns the absolute value of the larges element of the entire matrix.
 * \author Ole Schuett
 ******************************************************************************/
double dbm_maxabs(const dbm_matrix_t *matrix) {
  double maxabs = 0.0;
  for (int ishard = 0; ishard < matrix->nshards; ishard++) {
    dbm_shard_t *shard = &matrix->shards[ishard];
    for (int i = 0; i < shard->data_size; i++) {
      maxabs = fmax(maxabs, fabs(shard->data[i]));
    }
  }
  dbm_mpi_max_double(&maxabs, 1, matrix->dist->comm);
  return maxabs;
}

/*******************************************************************************
 * \brief Returns the name of the matrix of the given matrix.
 * \author Ole Schuett
 ******************************************************************************/
const char *dbm_get_name(const dbm_matrix_t *matrix) { return matrix->name; }

/*******************************************************************************
 * \brief Returns the number of local Non-Zero Elements of the given matrix.
 * \author Ole Schuett
 ******************************************************************************/
int dbm_get_nze(const dbm_matrix_t *matrix) {
  int nze = 0;
  for (int ishard = 0; ishard < matrix->nshards; ishard++) {
    nze += matrix->shards[ishard].data_size;
  }
  return nze;
}

/*******************************************************************************
 * \brief Returns the number of local blocks of the given matrix.
 * \author Ole Schuett
 ******************************************************************************/
int dbm_get_num_blocks(const dbm_matrix_t *matrix) {
  int nblocks = 0;
  for (int ishard = 0; ishard < matrix->nshards; ishard++) {
    nblocks += matrix->shards[ishard].nblocks;
  }
  return nblocks;
}

/*******************************************************************************
 * \brief Returns the row block sizes of the given matrix.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_get_row_sizes(const dbm_matrix_t *matrix, int *nrows,
                       const int **row_sizes) {
  *nrows = matrix->nrows;
  *row_sizes = matrix->row_sizes;
}

/*******************************************************************************
 * \brief Returns the column block sizes of the given matrix.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_get_col_sizes(const dbm_matrix_t *matrix, int *ncols,
                       const int **col_sizes) {
  *ncols = matrix->ncols;
  *col_sizes = matrix->col_sizes;
}

/*******************************************************************************
 * \brief Returns the local row block sizes of the given matrix.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_get_local_rows(const dbm_matrix_t *matrix, int *nlocal_rows,
                        const int **local_rows) {
  *nlocal_rows = matrix->dist->rows.nlocals;
  *local_rows = matrix->dist->rows.local_indicies;
}

/*******************************************************************************
 * \brief Returns the local column block sizes of the given matrix.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_get_local_cols(const dbm_matrix_t *matrix, int *nlocal_cols,
                        const int **local_cols) {
  *nlocal_cols = matrix->dist->cols.nlocals;
  *local_cols = matrix->dist->cols.local_indicies;
}

/*******************************************************************************
 * \brief Returns the MPI rank on which the given block should be stored.
 * \author Ole Schuett
 ******************************************************************************/
int dbm_get_stored_coordinates(const dbm_matrix_t *matrix, const int row,
                               const int col) {
  return dbm_distribution_stored_coords(matrix->dist, row, col);
}

/*******************************************************************************
 * \brief Returns the distribution of the given matrix.
 * \author Ole Schuett
 ******************************************************************************/
const dbm_distribution_t *dbm_get_distribution(const dbm_matrix_t *matrix) {
  return matrix->dist;
}

// EOF
