/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#ifndef DBM_MATRIX_H
#define DBM_MATRIX_H

#include <stdbool.h>

#include "dbm_distribution.h"
#include "dbm_shard.h"

/*******************************************************************************
 * \brief Internal struct for storing a matrix.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  dbm_distribution_t *dist;
  char *name;
  int nrows;
  int ncols;
  int *row_sizes;
  int *col_sizes;

  int nshards;
  dbm_shard_t *shards;
} dbm_matrix_t;

/*******************************************************************************
 * \brief Internal struct for storing a block iterator.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  const dbm_matrix_t *matrix;
  int next_block;
  int next_shard;
} dbm_iterator_t;

/*******************************************************************************
 * \brief Creates a new matrix.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_create(dbm_matrix_t **matrix_out, dbm_distribution_t *dist,
                const char name[], const int nrows, const int ncols,
                const int row_sizes[nrows], const int col_sizes[ncols]);

/*******************************************************************************
 * \brief Releases a matrix and all its ressources.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_release(dbm_matrix_t *matrix);

/*******************************************************************************
 * \brief Copies content of matrix_b into matrix_a.
 *        Matrices must have the same row/col block sizes and distribution.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_copy(dbm_matrix_t *matrix_a, const dbm_matrix_t *matrix_b);

/*******************************************************************************
 * \brief Copies content of matrix_b into matrix_a.
 *        Matrices may have different distributions.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_redistribute(const dbm_matrix_t *matrix, dbm_matrix_t *redist);

/*******************************************************************************
 * \brief Looks up a block from given matrics. This routine is thread-safe.
 *        If the block is not found then a null pointer is returned.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_get_block_p(dbm_matrix_t *matrix, const int row, const int col,
                     double **block, int *row_size, int *col_size);

/*******************************************************************************
 * \brief Adds a block to given matrix. This routine is thread-safe.
 *        If block already exist then it gets overwritten (or summed).
 * \author Ole Schuett
 ******************************************************************************/
void dbm_put_block(dbm_matrix_t *matrix, const int row, const int col,
                   const bool summation, const double *block);

/*******************************************************************************
 * \brief Remove all blocks from matrix, but does not release underlying memory.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_clear(dbm_matrix_t *matrix);

/*******************************************************************************
 * \brief Internal routine for re-computing any invalid block norms.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_compute_block_norms(dbm_matrix_t *matrix);

/*******************************************************************************
 * \brief Removes all blocks from the matrix whose norm is below the threshold.
 *        Blocks of size zero are always kept.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_filter(dbm_matrix_t *matrix, const double eps);

/*******************************************************************************
 * \brief Adds list of blocks efficiently. The blocks will be filled with zeros.
 *        This routine must always be called within an OpenMP parallel region.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_reserve_blocks(dbm_matrix_t *matrix, const int nblocks,
                        const int rows[], const int cols[]);

/*******************************************************************************
 * \brief Multiplies all entries in the given matrix by the given factor alpha.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_scale(dbm_matrix_t *matrix, const double alpha);

/*******************************************************************************
 * \brief Sets all blocks in the given matrix to zero.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_zero(dbm_matrix_t *matrix);

/*******************************************************************************
 * \brief Adds matrix_b to matrix_a.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_add(dbm_matrix_t *matrix_a, const dbm_matrix_t *matrix_b);

/*******************************************************************************
 * \brief Creates an iterator for the blocks of the given matrix.
 *        The iteration order is not stable.
 *        This routine must always be called within an OpenMP parallel region.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_iterator_start(dbm_iterator_t **iter_out, const dbm_matrix_t *matrix);

/*******************************************************************************
 * \brief Returns number of blocks the iterator will provide to calling thread.
 * \author Ole Schuett
 ******************************************************************************/
int dbm_iterator_num_blocks(const dbm_iterator_t *iter);

/*******************************************************************************
 * \brief Tests whether the given iterator has any block left.
 * \author Ole Schuett
 ******************************************************************************/
bool dbm_iterator_blocks_left(const dbm_iterator_t *iter);

/*******************************************************************************
 * \brief Returns the next block from the given iterator.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_iterator_next_block(dbm_iterator_t *iter, int *row, int *col,
                             double **block, int *row_size, int *col_size);

/*******************************************************************************
 * \brief Releases the given iterator.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_iterator_stop(dbm_iterator_t *iter);

/*******************************************************************************
 * \brief Computes a checksum of the given matrix.
 * \author Ole Schuett
 ******************************************************************************/
double dbm_checksum(const dbm_matrix_t *matrix);

/*******************************************************************************
 * \brief Returns the absolute value of the larges element of the entire matrix.
 * \author Ole Schuett
 ******************************************************************************/
double dbm_maxabs(const dbm_matrix_t *matrix);

/*******************************************************************************
 * \brief Returns the name of the matrix of the given matrix.
 * \author Ole Schuett
 ******************************************************************************/
const char *dbm_get_name(const dbm_matrix_t *matrix);

/*******************************************************************************
 * \brief Returns the number of local Non-Zero Elements of the given matrix.
 * \author Ole Schuett
 ******************************************************************************/
int dbm_get_nze(const dbm_matrix_t *matrix);

/*******************************************************************************
 * \brief Returns the number of local blocks of the given matrix.
 * \author Ole Schuett
 ******************************************************************************/
int dbm_get_num_blocks(const dbm_matrix_t *matrix);

/*******************************************************************************
 * \brief Returns the row block sizes of the given matrix.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_get_row_sizes(const dbm_matrix_t *matrix, int *nrows,
                       const int **row_sizes);

/*******************************************************************************
 * \brief Returns the column block sizes of the given matrix.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_get_col_sizes(const dbm_matrix_t *matrix, int *ncols,
                       const int **col_sizes);

/*******************************************************************************
 * \brief Returns the local row block sizes of the given matrix.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_get_local_rows(const dbm_matrix_t *matrix, int *nlocal_rows,
                        const int **local_rows);

/*******************************************************************************
 * \brief Returns the local column block sizes of the given matrix.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_get_local_cols(const dbm_matrix_t *matrix, int *nlocal_cols,
                        const int **local_cols);

/*******************************************************************************
 * \brief Returns the MPI rank on which the given block should be stored.
 * \author Ole Schuett
 ******************************************************************************/
int dbm_get_stored_coordinates(const dbm_matrix_t *matrix, const int row,
                               const int col);

/*******************************************************************************
 * \brief Returns the distribution of the given matrix.
 * \author Ole Schuett
 ******************************************************************************/
const dbm_distribution_t *dbm_get_distribution(const dbm_matrix_t *matrix);

#endif

// EOF
