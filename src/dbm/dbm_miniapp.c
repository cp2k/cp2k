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
 * \brief Wrapper for printf, passed to dbm_library_print_stats.
 * \author Ole Schuett
 ******************************************************************************/
static void print_func(char *message, int output_unit) {
  (void)output_unit; // mark used
  printf("%s", message);
}

/*******************************************************************************
 * \brief Returns the smaller of two given integer (missing from the C standard)
 * \author Ole Schuett
 ******************************************************************************/
static inline int imin(int x, int y) { return (x < y ? x : y); }

/*******************************************************************************
 * \brief Private routine for creating a distribution and an empty matrix.
 * \author Ole Schuett
 ******************************************************************************/
static dbm_matrix_t *create_some_matrix(const int row_size, const int col_size,
                                        const dbm_mpi_comm_t comm) {
  const int N = 8000;
  const int nrow = imin(500, N / row_size);
  const int ncol = imin(500, N / col_size);

  int cart_dims[2], cart_periods[2], cart_coords[2];
  dbm_mpi_cart_get(comm, 2, cart_dims, cart_periods, cart_coords);

  // Create distribution.
  int *row_dist = malloc(nrow * sizeof(int));
  int *col_dist = malloc(ncol * sizeof(int));
  for (int i = 0; i < nrow; i++) {
    row_dist[i] = i % cart_dims[0];
  }
  for (int i = 0; i < ncol; i++) {
    col_dist[i] = i % cart_dims[1];
  }
  const int fortran_comm = dbm_mpi_comm_c2f(comm);
  dbm_distribution_t *dist = NULL;
  dbm_distribution_new(&dist, fortran_comm, nrow, ncol, row_dist, col_dist);
  free(row_dist);
  free(col_dist);

  // Create matrix.
  int *row_sizes = malloc(nrow * sizeof(int));
  int *col_sizes = malloc(ncol * sizeof(int));
  for (int i = 0; i < nrow; i++) {
    row_sizes[i] = row_size;
  }
  for (int i = 0; i < ncol; i++) {
    col_sizes[i] = col_size;
  }
  dbm_matrix_t *matrix = NULL;
  dbm_create(&matrix, dist, "some name", nrow, ncol, row_sizes, col_sizes);
  dbm_distribution_release(dist);
  free(row_sizes);
  free(col_sizes);
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
    int nblocks = 0;
#pragma omp for collapse(2)
    for (int row = 0; row < nrows; row++) {
      for (int col = 0; col < ncols; col++) {
        if (dbm_get_stored_coordinates(matrix, row, col) ==
            matrix->dist->my_rank) {
          nblocks++;
        }
      }
    }
    int *reserve_row = malloc(nblocks * sizeof(int));
    int *reserve_col = malloc(nblocks * sizeof(int));
    int iblock = 0;
#pragma omp for collapse(2)
    for (int row = 0; row < nrows; row++) {
      for (int col = 0; col < ncols; col++) {
        if (dbm_get_stored_coordinates(matrix, row, col) ==
            matrix->dist->my_rank) {
          reserve_row[iblock] = row;
          reserve_col[iblock] = col;
          iblock++;
        }
      }
    }
    assert(iblock == nblocks);
    dbm_reserve_blocks(matrix, nblocks, reserve_row, reserve_col);
    free(reserve_row);
    free(reserve_col);
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
 * \brief Run a benchmark of dbm_multiply with given block sizes.
 * \author Ole Schuett
 ******************************************************************************/
void bechmark_multiply(const int m, const int n, const int k,
                       const dbm_mpi_comm_t comm) {
  dbm_matrix_t *matrix_a = create_some_matrix(m, k, comm);
  dbm_matrix_t *matrix_b = create_some_matrix(k, n, comm);
  dbm_matrix_t *matrix_c = create_some_matrix(m, n, comm);
  reserve_all_blocks(matrix_a);
  reserve_all_blocks(matrix_b);
  set_all_blocks(matrix_a);
  set_all_blocks(matrix_b);

  int64_t flop = 0;
  const double time_start_multiply = omp_get_wtime();
  dbm_multiply(false, false, 1.0, matrix_a, matrix_b, 1.0, matrix_c, false,
               1e-8, &flop);
  const double time_end_multiply = omp_get_wtime();

  dbm_release(matrix_a);
  dbm_release(matrix_b);
  dbm_release(matrix_c);

  if (dbm_mpi_comm_rank(comm) == 0) {
    const double duration = time_end_multiply - time_start_multiply;
    printf("multiply  %3i  x  %3i  x  %3i : %6.3f s  =>  %5.1f GFLOP/s \n", m,
           n, k, duration, 1e-9 * flop / duration);
  }
}

/*******************************************************************************
 * \brief Stand-alone miniapp for smoke-testing and benchmarking dbm_multiply.
 * \author Ole Schuett
 ******************************************************************************/
int main(int argc, char *argv[]) {
  dbm_mpi_init(&argc, &argv);
  dbm_library_init();

  const dbm_mpi_comm_t world_comm = dbm_mpi_get_comm_world();
  const int nranks = dbm_mpi_comm_size(world_comm);
  const int my_rank = dbm_mpi_comm_rank(world_comm);

  if (my_rank == 0) {
    printf("MPI-ranks: %i  OpenMP-threads: %i\n\n", nranks,
           omp_get_max_threads());
  }

  if (offload_get_device_count() > 0) {
    offload_set_chosen_device(my_rank % offload_get_device_count());
  }

  // Create 2D cart.
  const int dims[2] = {nranks, 1};
  const int periods[2] = {true, true};
  dbm_mpi_comm_t comm =
      dbm_mpi_cart_create(world_comm, 2, dims, periods, false);

  bechmark_multiply(4, 4, 4, comm);

  bechmark_multiply(128, 4, 4, comm);
  bechmark_multiply(4, 128, 4, comm);
  bechmark_multiply(4, 4, 128, comm);
  bechmark_multiply(4, 128, 128, comm);
  bechmark_multiply(128, 4, 128, comm);
  bechmark_multiply(128, 128, 4, comm);
  bechmark_multiply(128, 128, 128, comm);

  dbm_library_print_stats(dbm_mpi_comm_c2f(comm), &print_func, 0);
  dbm_library_finalize();
  dbm_mpi_comm_free(&comm);
  dbm_mpi_finalize();
  return 0;
}

// EOF
