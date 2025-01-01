/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2025 CP2K developers group <https://cp2k.org>              */
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

#include "dbm_distribution.h"
#include "dbm_hyperparams.h"

/*******************************************************************************
 * \brief Private routine for creating a new one dimensional distribution.
 * \author Ole Schuett
 ******************************************************************************/
static void dbm_dist_1d_new(dbm_dist_1d_t *dist, const int length,
                            const int coords[length], const dbm_mpi_comm_t comm,
                            const int nshards) {
  dist->comm = comm;
  dist->nshards = nshards;
  dist->my_rank = dbm_mpi_comm_rank(comm);
  dist->nranks = dbm_mpi_comm_size(comm);
  dist->length = length;
  dist->index2coord = malloc(length * sizeof(int));
  assert(dist->index2coord != NULL);
  memcpy(dist->index2coord, coords, length * sizeof(int));

  // Check that cart coordinates and ranks are equivalent.
  int cart_dims[1], cart_periods[1], cart_coords[1];
  dbm_mpi_cart_get(comm, 1, cart_dims, cart_periods, cart_coords);
  assert(dist->nranks == cart_dims[0]);
  assert(dist->my_rank == cart_coords[0]);

  // Count local rows/columns.
  for (int i = 0; i < length; i++) {
    assert(0 <= coords[i] && coords[i] < dist->nranks);
    if (coords[i] == dist->my_rank) {
      dist->nlocals++;
    }
  }

  // Store local rows/columns.
  dist->local_indicies = malloc(dist->nlocals * sizeof(int));
  assert(dist->local_indicies != NULL);
  int j = 0;
  for (int i = 0; i < length; i++) {
    if (coords[i] == dist->my_rank) {
      dist->local_indicies[j++] = i;
    }
  }
  assert(j == dist->nlocals);
}

/*******************************************************************************
 * \brief Private routine for releasing a one dimensional distribution.
 * \author Ole Schuett
 ******************************************************************************/
static void dbm_dist_1d_free(dbm_dist_1d_t *dist) {
  free(dist->index2coord);
  free(dist->local_indicies);
  dbm_mpi_comm_free(&dist->comm);
}

/*******************************************************************************
 * \brief Returns the larger of two given integer (missing from the C standard)
 * \author Ole Schuett
 ******************************************************************************/
static inline int imax(int x, int y) { return (x > y ? x : y); }

/*******************************************************************************
 * \brief Private routine for finding the optimal number of shard rows.
 * \author Ole Schuett
 ******************************************************************************/
static int find_best_nrow_shards(const int nshards, const int nrows,
                                 const int ncols) {
  const double target = (double)imax(nrows, 1) / (double)imax(ncols, 1);
  int best_nrow_shards = nshards;
  double best_error = fabs(log(target / (double)nshards));

  for (int nrow_shards = 1; nrow_shards <= nshards; nrow_shards++) {
    const int ncol_shards = nshards / nrow_shards;
    if (nrow_shards * ncol_shards != nshards)
      continue; // Not a factor of nshards.
    const double ratio = (double)nrow_shards / (double)ncol_shards;
    const double error = fabs(log(target / ratio));
    if (error < best_error) {
      best_error = error;
      best_nrow_shards = nrow_shards;
    }
  }
  return best_nrow_shards;
}

/*******************************************************************************
 * \brief Creates a new two dimensional distribution.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_distribution_new(dbm_distribution_t **dist_out, const int fortran_comm,
                          const int nrows, const int ncols,
                          const int row_dist[nrows],
                          const int col_dist[ncols]) {
  assert(omp_get_num_threads() == 1);
  dbm_distribution_t *dist = calloc(1, sizeof(dbm_distribution_t));
  dist->ref_count = 1;

  dist->comm = dbm_mpi_comm_f2c(fortran_comm);
  dist->my_rank = dbm_mpi_comm_rank(dist->comm);
  dist->nranks = dbm_mpi_comm_size(dist->comm);

  const int row_dim_remains[2] = {1, 0};
  const dbm_mpi_comm_t row_comm = dbm_mpi_cart_sub(dist->comm, row_dim_remains);

  const int col_dim_remains[2] = {0, 1};
  const dbm_mpi_comm_t col_comm = dbm_mpi_cart_sub(dist->comm, col_dim_remains);

  const int nshards = SHARDS_PER_THREAD * omp_get_max_threads();
  const int nrow_shards = find_best_nrow_shards(nshards, nrows, ncols);
  const int ncol_shards = nshards / nrow_shards;

  dbm_dist_1d_new(&dist->rows, nrows, row_dist, row_comm, nrow_shards);
  dbm_dist_1d_new(&dist->cols, ncols, col_dist, col_comm, ncol_shards);

  assert(*dist_out == NULL);
  *dist_out = dist;
}

/*******************************************************************************
 * \brief Increases the reference counter of the given distribution.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_distribution_hold(dbm_distribution_t *dist) {
  assert(dist->ref_count > 0);
  dist->ref_count++;
}

/*******************************************************************************
 * \brief Decreases the reference counter of the given distribution.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_distribution_release(dbm_distribution_t *dist) {
  assert(dist->ref_count > 0);
  dist->ref_count--;
  if (dist->ref_count == 0) {
    dbm_dist_1d_free(&dist->rows);
    dbm_dist_1d_free(&dist->cols);
    free(dist);
  }
}

/*******************************************************************************
 * \brief Returns the rows of the given distribution.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_distribution_row_dist(const dbm_distribution_t *dist, int *nrows,
                               const int **row_dist) {
  assert(dist->ref_count > 0);
  *nrows = dist->rows.length;
  *row_dist = dist->rows.index2coord;
}

/*******************************************************************************
 * \brief Returns the columns of the given distribution.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_distribution_col_dist(const dbm_distribution_t *dist, int *ncols,
                               const int **col_dist) {
  assert(dist->ref_count > 0);
  *ncols = dist->cols.length;
  *col_dist = dist->cols.index2coord;
}

/*******************************************************************************
 * \brief Returns the MPI rank on which the given block should be stored.
 * \author Ole Schuett
 ******************************************************************************/
int dbm_distribution_stored_coords(const dbm_distribution_t *dist,
                                   const int row, const int col) {
  assert(dist->ref_count > 0);
  assert(0 <= row && row < dist->rows.length);
  assert(0 <= col && col < dist->cols.length);
  int coords[2] = {dist->rows.index2coord[row], dist->cols.index2coord[col]};
  return dbm_mpi_cart_rank(dist->comm, coords);
}

// EOF
