/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2025 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(__LIBXSMM)
#include <libxsmm.h>
#endif

#include "../offload/offload_library.h"
#include "dbm_library.h"
#include "dbm_matrix.h"
#include "dbm_mpi.h"

/*******************************************************************************
 * \brief Wrapper for printf, passed to dbm_library_print_stats.
 * \author Ole Schuett
 ******************************************************************************/
static void print_func(char *message, int output_unit) {
  if (output_unit == 0) { // i.e. my_rank == 0
    printf("%s", message);
  }
}

/*******************************************************************************
 * \brief Returns the smaller of the two integers (missing from the C standard).
 * \author Ole Schuett
 ******************************************************************************/
static inline int imin(int x, int y) { return (x < y ? x : y); }

/*******************************************************************************
 * \brief Private routine for creating a distribution.
 * \author Ole Schuett
 ******************************************************************************/
static dbm_distribution_t *create_dist(const int nrows, const int ncols,
                                       const dbm_mpi_comm_t comm) {
  int cart_dims[2], cart_periods[2], cart_coords[2];
  dbm_mpi_cart_get(comm, 2, cart_dims, cart_periods, cart_coords);

  // Create distribution.
  assert(0 < nrows && 0 < ncols);
  int *row_dist = malloc(nrows * sizeof(int));
  int *col_dist = malloc(ncols * sizeof(int));
  assert(row_dist != NULL && col_dist != NULL);
  for (int i = 0; i < nrows; i++) {
    row_dist[i] = i % cart_dims[0];
  }
  for (int i = 0; i < ncols; i++) {
    col_dist[i] = i % cart_dims[1];
  }
  const int fortran_comm = dbm_mpi_comm_c2f(comm);
  dbm_distribution_t *dist = NULL;
  dbm_distribution_new(&dist, fortran_comm, nrows, ncols, row_dist, col_dist);
  free(row_dist);
  free(col_dist);
  return dist;
}

/*******************************************************************************
 * \brief Private routine for creating a distribution and an empty matrix.
 * \author Ole Schuett
 ******************************************************************************/
static dbm_matrix_t *
create_some_matrix(const int nrows, const int ncols, const int nrows_min,
                   const int nrows_max, const int ncols_min,
                   const int ncols_max, const dbm_mpi_comm_t comm) {
  // Create distribution.
  dbm_distribution_t *dist = create_dist(nrows, ncols, comm);

  // Create matrix.
  assert(0 < nrows && 0 < ncols);
  int *row_sizes = malloc(nrows * sizeof(int));
  int *col_sizes = malloc(ncols * sizeof(int));
  assert(row_sizes != NULL && col_sizes != NULL);
  assert(0 < nrows_min && nrows_min <= nrows_max);
  assert(0 < ncols_min && ncols_min <= ncols_max);
  if (nrows_min != nrows_max) {
    const int row_size = nrows_max - nrows_min + 1;
    for (int i = 0; i < nrows; i++) {
      row_sizes[i] = rand() % row_size + 1;
    }
  } else {
    for (int i = 0; i < nrows; i++) {
      row_sizes[i] = nrows_max;
    }
  }
  if (ncols_min != ncols_max) {
    const int col_size = ncols_max - ncols_min + 1;
    for (int i = 0; i < ncols; i++) {
      col_sizes[i] = rand() % col_size + 1;
    }
  } else {
    for (int i = 0; i < ncols; i++) {
      col_sizes[i] = ncols_max;
    }
  }
  dbm_matrix_t *matrix = NULL;
  dbm_create(&matrix, dist, "some name", nrows, ncols, row_sizes, col_sizes);
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
          ++nblocks;
        }
      }
    }
    assert(0 < nblocks);
    int *reserve_row = malloc(nblocks * sizeof(int));
    int *reserve_col = malloc(nblocks * sizeof(int));
    assert(reserve_row != NULL && reserve_col != NULL);
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
 * \brief Private routine for setting all blocks.
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
        block[i] = 1.0 / (i + 1);
      }
    }
    dbm_iterator_stop(iter);
  }
}
/*******************************************************************************
 * \brief Run a benchmark of dbm_multiply with given block sizes.
 * \author Ole Schuett
 ******************************************************************************/
void benchmark_multiply(const int M, const int N, const int K, const int m,
                        const int n, const int k, const dbm_mpi_comm_t comm) {
  dbm_matrix_t *matrix_a = create_some_matrix(M, K, 1, m, k, k, comm);
  dbm_matrix_t *matrix_b = create_some_matrix(K, N, k, k, 1, n, comm);
  dbm_distribution_t *dist_c = create_dist(M, N, comm);
  dbm_matrix_t *matrix_c = NULL, *matrix_d = NULL;
  dbm_create(&matrix_c, dist_c, "result", M, N, matrix_a->row_sizes,
             matrix_b->col_sizes);
  dbm_distribution_release(dist_c);

  reserve_all_blocks(matrix_a);
  reserve_all_blocks(matrix_b);
  set_all_blocks(matrix_a);
  set_all_blocks(matrix_b);

  dbm_distribution_t *const dist_shared = matrix_c->dist;
  dbm_create(&matrix_d, dist_shared, matrix_c->name, matrix_c->nrows,
             matrix_c->ncols, matrix_c->row_sizes, matrix_c->col_sizes);
  dbm_copy(matrix_d, matrix_c);

  int64_t flop = 0;
  const double time_start_multiply = omp_get_wtime();
  dbm_multiply(false, false, 1.0, matrix_a, matrix_b, 1.0, matrix_c, false,
               1e-8, &flop);
  const double time_end_multiply = omp_get_wtime();

  // Calculate result on the host for validation
  dbm_multiply(false, false, 1.0, matrix_a, matrix_b, 1.0, matrix_d, false,
               1e-8, NULL);

  if (dbm_mpi_comm_rank(comm) == 0) {
    printf("%5i x %5i x %5i  with  %3i x %3i x %3i blocks: ", M, N, K, m, n, k);
  }

  // Validate result
  const double maxeps = 1E-5, epsilon = dbm_maxeps(matrix_d, matrix_c);
  if (maxeps >= epsilon) {
    dbm_mpi_sum_int64(&flop, 1, comm);
    if (dbm_mpi_comm_rank(comm) == 0) {
      const double duration = time_end_multiply - time_start_multiply;
      printf("%6.3f s =>  %6.1f GFLOP/s\n", duration, 1e-9 * flop / duration);
      fflush(stdout);
    }
  } else {
    printf("ERROR\n");
    fprintf(stderr, "Failed validation (epsilon=%f).\n", epsilon);
    exit(1);
  }

  dbm_release(matrix_a);
  dbm_release(matrix_b);
  dbm_release(matrix_c);
  dbm_release(matrix_d);
}

/*******************************************************************************
 * \brief Stand-alone miniapp for smoke-testing and benchmarking dbm_multiply.
 * \author Ole Schuett
 ******************************************************************************/
int main(int argc, char *argv[]) {
  int result = EXIT_SUCCESS;

  srand(25071975); // seed rng

  dbm_mpi_init(&argc, &argv);
  dbm_library_init();

  const dbm_mpi_comm_t world_comm = dbm_mpi_get_comm_world();
  const int nranks = dbm_mpi_comm_size(world_comm);
  const int my_rank = dbm_mpi_comm_rank(world_comm);

  if (offload_get_device_count() > 0) {
    offload_set_chosen_device(my_rank % offload_get_device_count());
  }

  // Create 2D cart.
  int dims[2] = {0, 0};
  dbm_mpi_dims_create(nranks, 2, dims);
  const int periods[2] = {true, true};
  dbm_mpi_comm_t comm =
      dbm_mpi_cart_create(world_comm, 2, dims, periods, false);

  if (my_rank == 0) {
    printf("OpenMP-threads: %i  GPUs: %i", omp_get_max_threads(),
           imin(offload_get_device_count(), nranks));
#if defined(__LIBXSMM)
    printf("  Libxsmm: %s", LIBXSMM_VERSION);
#else
    printf("  Libxsmm: n/a");
#endif
#if defined(__parallel)
    printf("  MPI-ranks: %i  MPI-cart: %i x %i", nranks, dims[0], dims[1]);
#else
    printf("  MPI: n/a");
#endif
    printf("\n\n");
    fflush(stdout);
  }

  if (1 >= argc) {
    benchmark_multiply(16384, 128, 128, 4, 4, 4, comm);
    benchmark_multiply(128, 16384, 128, 4, 4, 4, comm);
    benchmark_multiply(128, 128, 16384, 4, 4, 4, comm);
    benchmark_multiply(645, 645, 645, 4, 4, 4, comm);
    if (my_rank == 0)
      printf("\n");

    benchmark_multiply(60, 500, 500, 128, 4, 4, comm);
    benchmark_multiply(500, 60, 500, 4, 128, 4, comm);
    benchmark_multiply(500, 500, 60, 4, 4, 128, comm);
    if (my_rank == 0)
      printf("\n");

    benchmark_multiply(500, 60, 60, 4, 128, 128, comm);
    benchmark_multiply(60, 500, 60, 128, 4, 128, comm);
    benchmark_multiply(60, 60, 500, 128, 128, 4, comm);
    if (my_rank == 0)
      printf("\n");

    benchmark_multiply(350, 350, 350, 23, 23, 23, comm);
    benchmark_multiply(250, 250, 250, 32, 32, 32, comm);
    benchmark_multiply(60, 60, 60, 128, 128, 128, comm);
  } else { /* read triplet(s) from file or one triplet from command line */
    FILE *const file = fopen(argv[1], "r"); /* try 1st arg as filename */
    char buffer[1024];
    const char delims[] = "x,;:|/\t ";
    int mnk[] = {0, 0, 0}, i = 1, j = 0;
    while (i < argc &&
           (NULL == file || NULL != fgets(buffer, sizeof(buffer), file))) {
      const char *arg = strtok(NULL != file ? buffer : argv[i], delims);
      for (; NULL != arg && j < 3; arg = strtok(NULL, delims), ++j) {
        mnk[j] = atoi(arg);
      }
      if (NULL != file) {
        j = 0;
      } else if (++i < argc) {
        continue;
      }
      if (0 < mnk[0]) { /* valid MxNxK? */
        const int m = mnk[0];
        const int n = (0 < mnk[1] ? mnk[1] : m);
        const int k = (0 < mnk[2] ? mnk[2] : m);
        int M = (NULL == arg ? 0 : atoi(arg)), N, K;
        if (0 < M) {
          arg = strtok(NULL, delims);
          N = (NULL == arg ? 1 : atoi(arg));
          arg = strtok(NULL, delims);
          K = (NULL == arg ? 1 : atoi(arg));
        } else { /* default */
          M = N = K = 128;
        }
        benchmark_multiply(M, N, K, m, n, k, comm);
        mnk[0] = mnk[1] = mnk[2] = 0;
      } else {
        fprintf(stderr, "ERROR: invalid argument(s)\n");
        result = EXIT_FAILURE;
        i = argc; /* break */
      }
    }
    if (NULL != file) {
      fclose(file);
    }
  }

  if (EXIT_SUCCESS == result) {
    dbm_library_print_stats(dbm_mpi_comm_c2f(comm), &print_func, my_rank);
  }
  dbm_library_finalize();
  dbm_mpi_comm_free(&comm);
  dbm_mpi_finalize();
  return result;
}

// EOF
