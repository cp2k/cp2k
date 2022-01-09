/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#ifndef DBM_DISTRIBUTION_H
#define DBM_DISTRIBUTION_H

#include "dbm_mpi.h"

/*******************************************************************************
 * \brief Internal struct for storing a one dimensional distribution.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  int length;       // global number of rows/cols
  int *index2coord; // maps row/col indicies to cart coordinate
  int nlocals;
  int *local_indicies; // list of row/col indicies that reside locally.
  dbm_mpi_comm_t comm; // 1D communicator
  int nranks;
  int my_rank;
} dbm_dist_1d_t;

/*******************************************************************************
 * \brief Internal struct for storing a two dimensional distribution.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  int ref_count;
  dbm_dist_1d_t rows;
  dbm_dist_1d_t cols;
  dbm_mpi_comm_t comm;
  int nranks;
  int my_rank;
} dbm_distribution_t;

/*******************************************************************************
 * \brief Creates a new two dimensional distribution.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_distribution_new(dbm_distribution_t **dist_out, const int fortran_comm,
                          const int nrows, const int ncols,
                          const int row_dist[nrows], const int col_dist[ncols]);

/*******************************************************************************
 * \brief Increases the reference counter of the given distribution.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_distribution_hold(dbm_distribution_t *dist);

/*******************************************************************************
 * \brief Decreases the reference counter of the given distribution.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_distribution_release(dbm_distribution_t *dist);

/*******************************************************************************
 * \brief Returns the rows of the given distribution.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_distribution_row_dist(const dbm_distribution_t *dist, int *nrows,
                               const int **row_dist);

/*******************************************************************************
 * \brief Returns the columns of the given distribution.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_distribution_col_dist(const dbm_distribution_t *dist, int *ncols,
                               const int **col_dist);

/*******************************************************************************
 * \brief Returns the MPI rank on which the given block should be stored.
 * \author Ole Schuett
 ******************************************************************************/
int dbm_distribution_stored_coords(const dbm_distribution_t *dist,
                                   const int row, const int col);

#endif

// EOF
