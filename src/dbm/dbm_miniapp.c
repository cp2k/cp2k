/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../offload/offload_library.h"
#include "dbm_library.h"
#include "dbm_matrix.h"
#include "dbm_mpi.h"

/*******************************************************************************
 * \brief Private routine for creating a distribution and an empty matrix.
 * \author Ole Schuett
 ******************************************************************************/
static dbm_matrix_t *create_some_matrix(const dbm_mpi_comm_t comm) {
  const int nrow = 200;
  const int ncol = 200;

  int cart_dims[2], cart_periods[2], cart_coords[2];
  dbm_mpi_cart_get(comm, 2, cart_dims, cart_periods, cart_coords);

  // Create distribution.
  int row_dist[nrow];
  int col_dist[ncol];
  for (int i = 0; i < nrow; i++) {
    row_dist[i] = i % cart_dims[0];
  }
  for (int i = 0; i < ncol; i++) {
    col_dist[i] = i % cart_dims[1];
  }
  const int fortran_comm = dbm_mpi_comm_c2f(comm);
  dbm_distribution_t *dist = NULL;
  dbm_distribution_new(&dist, fortran_comm, nrow, ncol, row_dist, col_dist);

  // Create matrix.
  int row_sizes[nrow];
  int col_sizes[ncol];
  for (int i = 0; i < nrow; i++) {
    row_sizes[i] = 18;
  }
  for (int i = 0; i < ncol; i++) {
    col_sizes[i] = 18;
  }
  dbm_matrix_t *matrix = NULL;
  dbm_create(&matrix, dist, "some name", nrow, ncol, row_sizes, col_sizes);
  dbm_distribution_release(dist);
  return matrix;
}

/*******************************************************************************
 * \brief Private routine for reserving all blocks of the given matrix.
 * \author Ole Schuett
 ******************************************************************************/
static void reserve_all_blocks(dbm_matrix_t *matrix) {
  int nrows, ncols;
  const int *row_sizes, *col_sizes;
  dbm_get_row_sizes(matrix, &nrows, &row_sizes);
  dbm_get_col_sizes(matrix, &ncols, &col_sizes);

#pragma omp parallel
  {
    const int nblocks = nrows * ncols;
    int reserve_row[nblocks / omp_get_num_threads() + 1];
    int reserve_col[nblocks / omp_get_num_threads() + 1];
    int nblocks_reserve = 0;
#pragma omp for
    for (int i = 0; i < nblocks; i++) {
      const int row = i % nrows;
      const int col = i % ncols;
      if (dbm_get_stored_coordinates(matrix, row, col) ==
          matrix->dist->my_rank) {
        reserve_row[nblocks_reserve] = row;
        reserve_col[nblocks_reserve] = col;
        nblocks_reserve++;
      }
    }
    dbm_reserve_blocks(matrix, nblocks_reserve, reserve_row, reserve_col);
  }
}

/*******************************************************************************
 * \brief Private routine for setting all blocks to 1.0.
 * \author Ole Schuett
 ******************************************************************************/
static void set_all_blocks(dbm_matrix_t *matrix) {
#pragma omp parallel
  {
    dbm_iterator_t *iter = NULL;
    dbm_iterator_start(&iter, matrix);
    while (dbm_iterator_blocks_left(iter)) {
      int row, col, row_size, col_size;
      double *block;
      dbm_iterator_next_block(iter, &row, &col, &block, &row_size, &col_size);
      const int block_size = row_size * col_size;
      for (int i = 0; i < block_size; i++) {
        block[i] = 1.0;
      }
    }
    dbm_iterator_stop(iter);
  }
}

/*******************************************************************************
 * \brief Stand-alone miniapp for smoke-testing and benchmarking dbm_multiply.
 * \author Ole Schuett
 ******************************************************************************/
int main(int argc, char *argv[]) {
  dbm_mpi_init(&argc, &argv);

  const dbm_mpi_comm_t world_comm = dbm_mpi_get_comm_world();
  const int nranks = dbm_mpi_comm_size(world_comm);
  const int my_rank = dbm_mpi_comm_rank(world_comm);

  if (my_rank == 0) {
    printf("MPI ranks:      %i\n", nranks);
    printf("OpenMP threads: %i\n", omp_get_max_threads());
  }

  if (offload_get_device_count() > 0) {
    offload_set_chosen_device(my_rank % offload_get_device_count());
  }

  // Create 2D cart.
  const int dims[2] = {nranks, 1};
  const int periods[2] = {true, true};
  dbm_mpi_comm_t comm =
      dbm_mpi_cart_create(world_comm, 2, dims, periods, false);

  dbm_library_init();

  // Benchmark reserve blocks.
  const double time_start_reserve = omp_get_wtime();
  for (int icycle = 0; icycle < 50; icycle++) {
    dbm_matrix_t *matrix = create_some_matrix(comm);
    reserve_all_blocks(matrix);
    dbm_release(matrix);
  }
  const double time_end_reserve = omp_get_wtime();
  if (my_rank == 0) {
    const double t = time_end_reserve - time_start_reserve;
    printf("reserve blocks: %.3f seconds\n", t);
  }

  // Benchmark matrix multiply.
  dbm_matrix_t *matrix_a = create_some_matrix(comm);
  dbm_matrix_t *matrix_b = create_some_matrix(comm);
  dbm_matrix_t *matrix_c = create_some_matrix(comm);
  reserve_all_blocks(matrix_a);
  reserve_all_blocks(matrix_b);
  set_all_blocks(matrix_a);
  set_all_blocks(matrix_b);
  int64_t flop;
  const double time_start_multiply = omp_get_wtime();
  dbm_multiply(false, false, 1.0, matrix_a, matrix_b, 1.0, matrix_c, false,
               1e-8, &flop);
  const double time_end_multiply = omp_get_wtime();
  dbm_release(matrix_a);
  dbm_release(matrix_b);
  dbm_release(matrix_c);

  if (my_rank == 0) {
    const double t = time_end_multiply - time_start_multiply;
    printf("matrix multiply: %.3f s, %.1f MFLOP/s \n", t, 1e-6 * flop / t);
    printf("done :-)\n");
  }

  dbm_library_finalize();
  dbm_mpi_comm_free(&comm);
  dbm_mpi_finalize();
  return 0;
}

// EOF
