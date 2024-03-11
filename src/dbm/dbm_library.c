/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2024 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <inttypes.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dbm_library.h"
#include "dbm_mempool.h"
#include "dbm_mpi.h"

#define DBM_NUM_COUNTERS 64

static int64_t **per_thread_counters = NULL;
static bool library_initialized = false;
static int max_threads = 0;

#if !defined(_OPENMP)
#error "OpenMP is required. Please add -fopenmp to your C compiler flags."
#endif

/*******************************************************************************
 * \brief Initializes the DBM library.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_library_init(void) {
  assert(omp_get_num_threads() == 1);

  if (library_initialized) {
    fprintf(stderr, "DBM library was already initialized.\n");
    abort();
  }

  max_threads = omp_get_max_threads();
  per_thread_counters = malloc(max_threads * sizeof(int64_t *));

  // Using parallel regions to ensure memory is allocated near a thread's core.
#pragma omp parallel default(none) shared(per_thread_counters)                 \
    num_threads(max_threads)
  {
    const int ithread = omp_get_thread_num();
    const size_t counters_size = DBM_NUM_COUNTERS * sizeof(int64_t);
    per_thread_counters[ithread] = malloc(counters_size);
    memset(per_thread_counters[ithread], 0, counters_size);
  }

  library_initialized = true;
}

/*******************************************************************************
 * \brief Finalizes the DBM library.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_library_finalize(void) {
  assert(omp_get_num_threads() == 1);

  if (!library_initialized) {
    fprintf(stderr, "Error: DBM library is not initialized.\n");
    abort();
  }

  for (int i = 0; i < max_threads; i++) {
    free(per_thread_counters[i]);
  }
  free(per_thread_counters);
  per_thread_counters = NULL;

  dbm_mempool_clear();
  library_initialized = false;
}

/*******************************************************************************
 * \brief Computes min(3, floor(log10(x))).
 * \author Ole Schuett
 ******************************************************************************/
static int floorlog10(const int x) {
  if (x >= 1000) {
    return 3;
  }
  if (x >= 100) {
    return 2;
  }
  if (x >= 10) {
    return 1;
  }
  return 0;
}

/*******************************************************************************
 * \brief Add given block multiplication to stats. This routine is thread-safe.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_library_counter_increment(const int m, const int n, const int k) {
  const int ithread = omp_get_thread_num();
  assert(ithread < max_threads);
  const int idx = 16 * floorlog10(m) + 4 * floorlog10(n) + floorlog10(k);
  per_thread_counters[ithread][idx]++;
}

/*******************************************************************************
 * \brief Comperator passed to qsort to compare two counters.
 * \author Ole Schuett
 ******************************************************************************/
static int compare_counters(const void *a, const void *b) {
  return *(const int64_t *)b - *(const int64_t *)a;
}

/*******************************************************************************
 * \brief Prints statistics gathered by the DBM library.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_library_print_stats(const int fortran_comm,
                             void (*print_func)(char *, int),
                             const int output_unit) {
  assert(omp_get_num_threads() == 1);

  if (!library_initialized) {
    fprintf(stderr, "Error: DBM library is not initialized.\n");
    abort();
  }

  const dbm_mpi_comm_t comm = dbm_mpi_comm_f2c(fortran_comm);
  // Sum all counters across threads and mpi ranks.
  int64_t counters[DBM_NUM_COUNTERS][2];
  memset(counters, 0, DBM_NUM_COUNTERS * 2 * sizeof(int64_t));
  double total = 0.0;
  for (int i = 0; i < DBM_NUM_COUNTERS; i++) {
    counters[i][1] = i; // needed as inverse index after qsort
    for (int j = 0; j < max_threads; j++) {
      counters[i][0] += per_thread_counters[j][i];
    }
    dbm_mpi_sum_int64(&counters[i][0], 1, comm);
    total += counters[i][0];
  }

  // Sort counters.
  qsort(counters, DBM_NUM_COUNTERS, 2 * sizeof(int64_t), &compare_counters);

  // Print counters.
  print_func("\n", output_unit);
  print_func(" ----------------------------------------------------------------"
             "---------------\n",
             output_unit);
  print_func(" -                                                               "
             "              -\n",
             output_unit);
  print_func(" -                                DBM STATISTICS                 "
             "              -\n",
             output_unit);
  print_func(" -                                                               "
             "              -\n",
             output_unit);
  print_func(" ----------------------------------------------------------------"
             "---------------\n",
             output_unit);
  print_func("    M  x    N  x    K                                          "
             "COUNT     PERCENT\n",
             output_unit);

  const char *labels[] = {"?", "??", "???", ">999"};
  for (int i = 0; i < DBM_NUM_COUNTERS; i++) {
    if (counters[i][0] == 0) {
      continue; // skip empty counters
    }
    const double percent = 100.0 * counters[i][0] / total;
    const int idx = counters[i][1];
    const int m = (idx % 64) / 16;
    const int n = (idx % 16) / 4;
    const int k = (idx % 4) / 1;
    char buffer[100];
    snprintf(buffer, sizeof(buffer),
             " %4s  x %4s  x %4s %46" PRId64 " %10.2f%%\n", labels[m],
             labels[n], labels[k], counters[i][0], percent);
    print_func(buffer, output_unit);
  }

  print_func(" ----------------------------------------------------------------"
             "---------------\n",
             output_unit);
}

// EOF
