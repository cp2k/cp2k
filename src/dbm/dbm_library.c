/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#include "dbm_library.h"
#include "../mpiwrap/cp_mpi.h"
#include "../offload/offload_library.h"
#include "../offload/offload_mempool.h"

#include <assert.h>
#include <inttypes.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DBM_LIBRARY_PRINT(FN, MSG, OUTPUT_UNIT)                                \
  ((FN)(MSG, (int)strlen(MSG), OUTPUT_UNIT))

static int64_t **per_thread_counters = NULL;
static double phase_seconds[DBM_NUM_PHASES] = {0};
static int64_t phase_bytes[DBM_NUM_PHASES] = {0};
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

  offload_init();

  max_threads = omp_get_max_threads();
  memset(phase_seconds, 0, sizeof(phase_seconds));
  memset(phase_bytes, 0, sizeof(phase_bytes));
  per_thread_counters = malloc(max_threads * sizeof(int64_t *));
  assert(per_thread_counters != NULL);

  // Using parallel regions to ensure memory is allocated near a thread's core.
#pragma omp parallel default(none) shared(per_thread_counters)                 \
    num_threads(max_threads)
  {
    const int ithread = omp_get_thread_num();
    const size_t counters_size = DBM_NUM_COUNTERS * sizeof(int64_t);
    per_thread_counters[ithread] = malloc(counters_size);
    assert(per_thread_counters[ithread] != NULL);
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

  offload_mempool_clear();
  library_initialized = false;
}

/*******************************************************************************
 * \brief Add the calling thread's counts to the stats. This routine is
 *        thread-safe.
 * \author Hans Pabst
 ******************************************************************************/
void dbm_library_counters_add(const int64_t counters[DBM_NUM_COUNTERS]) {
  const int ithread = omp_get_thread_num();
  assert(ithread < max_threads);
  for (int i = 0; i < DBM_NUM_COUNTERS; i++) {
    per_thread_counters[ithread][i] += counters[i];
  }
}

/*******************************************************************************
 * \brief Add the durations (seconds) of a multiplication's phases to the
 *        stats. Called once per dbm_multiply outside of parallel regions.
 * \author Hans Pabst
 ******************************************************************************/
void dbm_library_phases_add(const double seconds[DBM_NUM_PHASES],
                            const int64_t bytes[DBM_NUM_PHASES]) {
  assert(omp_get_num_threads() == 1);
  for (int i = 0; i < DBM_NUM_PHASES; i++) {
    phase_seconds[i] += seconds[i];
    phase_bytes[i] += bytes[i];
  }
}

/*******************************************************************************
 * \brief Comperator passed to qsort to compare two counters.
 * \author Ole Schuett
 ******************************************************************************/
static int compare_counters(const void *a, const void *b) {
  // Descending order: the difference of two counts may overflow int.
  const int64_t x = *(const int64_t *)a, y = *(const int64_t *)b;
  return (x < y) - (x > y);
}

/*******************************************************************************
 * \brief Prints statistics gathered by the DBM library.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_library_print_stats(const int fortran_comm,
                             void (*print_func)(const char *, int, int),
                             const int output_unit) {
  assert(omp_get_num_threads() == 1);

  if (!library_initialized) {
    fprintf(stderr, "Error: DBM library is not initialized.\n");
    abort();
  }

  const cp_mpi_comm_t comm = cp_mpi_comm_f2c(fortran_comm);
  // Sum all counters across threads and mpi ranks.
  int64_t counters[DBM_NUM_COUNTERS][2] = {{0}};
  double total = 0.0;
  for (int i = 0; i < DBM_NUM_COUNTERS; i++) {
    counters[i][1] = i; // needed as inverse index after qsort
    for (int j = 0; j < max_threads; j++) {
      counters[i][0] += per_thread_counters[j][i];
    }
    cp_mpi_sum_int64(&counters[i][0], 1, comm);
    total += counters[i][0];
  }

  // Phases per rank: minimum, average, and maximum over ranks (last: total).
  double phase_min[DBM_NUM_PHASES + 1], phase_avg[DBM_NUM_PHASES + 1];
  double phase_max[DBM_NUM_PHASES + 1];
  phase_max[DBM_NUM_PHASES] = 0.0;
  for (int i = 0; i < DBM_NUM_PHASES; i++) {
    phase_max[i] = phase_seconds[i];
    phase_max[DBM_NUM_PHASES] += phase_seconds[i];
  }
  for (int i = 0; i <= DBM_NUM_PHASES; i++) {
    phase_min[i] = -phase_max[i]; // minimum as maximum of the negated values
    phase_avg[i] = phase_max[i];
  }
  cp_mpi_max_double(phase_max, DBM_NUM_PHASES + 1, comm);
  cp_mpi_max_double(phase_min, DBM_NUM_PHASES + 1, comm);
  cp_mpi_sum_double(phase_avg, DBM_NUM_PHASES + 1, comm);
  int64_t bytes[DBM_NUM_PHASES];
  memcpy(bytes, phase_bytes, sizeof(bytes));
  cp_mpi_sum_int64(bytes, DBM_NUM_PHASES, comm);
  const int nranks = cp_mpi_comm_size(comm);

  // Sort counters.
  qsort(counters, DBM_NUM_COUNTERS, 2 * sizeof(int64_t), &compare_counters);

  // Determine if anything needs to be printed.
  bool print = false;
  for (int i = 0; i < DBM_NUM_COUNTERS && !print; i++) {
    if (counters[i][0] != 0) {
      print = true;
    }
  }
  if (!print) {
    return; // nothing to be printed
  }

  // Print counters.
  DBM_LIBRARY_PRINT(print_func, "\n", output_unit);
  DBM_LIBRARY_PRINT(
      print_func,
      " ----------------------------------------------------------------"
      "---------------\n",
      output_unit);
  DBM_LIBRARY_PRINT(
      print_func,
      " -                                                               "
      "              -\n",
      output_unit);
  DBM_LIBRARY_PRINT(
      print_func,
      " -                                DBM STATISTICS                 "
      "              -\n",
      output_unit);
  DBM_LIBRARY_PRINT(
      print_func,
      " -                                                               "
      "              -\n",
      output_unit);
  DBM_LIBRARY_PRINT(
      print_func,
      " ----------------------------------------------------------------"
      "---------------\n",
      output_unit);
  DBM_LIBRARY_PRINT(
      print_func,
      "    M  x    N  x    K                                          "
      "COUNT     PERCENT\n",
      output_unit);

  const char *labels[] = {"?", "??", "???", ">999"};
  char buffer[100];
  for (int i = 0; i < DBM_NUM_COUNTERS; i++) {
    if (counters[i][0] == 0) {
      continue; // skip empty counters
    }
    const double percent = 100.0 * counters[i][0] / total;
    const int idx = counters[i][1];
    const int m = (idx % 64) / 16;
    const int n = (idx % 16) / 4;
    const int k = (idx % 4) / 1;
    snprintf(buffer, sizeof(buffer),
             " %4s  x %4s  x %4s %46" PRId64 " %10.2f%%\n", labels[m],
             labels[n], labels[k], counters[i][0], percent);
    DBM_LIBRARY_PRINT(print_func, buffer, output_unit);
  }

  // Print phases, i.e., where a rank spends the time of dbm_multiply.
  if (0.0 < phase_max[DBM_NUM_PHASES]) {
    const char *const phases[] = {"setup",    "pack",   "alltoall",
                                  "sort",     "shift",  "upload",
                                  "multiply", "finish", "total"};
    DBM_LIBRARY_PRINT(
        print_func,
        " ----------------------------------------------------------------"
        "---------------\n",
        output_unit);
    snprintf(buffer, sizeof(buffer), "    %-24s %12s %12s %12s %12s\n",
             "PHASE PER RANK", "MIN [s]", "AVG [s]", "MAX [s]", "GB/s");
    DBM_LIBRARY_PRINT(print_func, buffer, output_unit);
    for (int i = 0; i <= DBM_NUM_PHASES; i++) {
      // Bytes sent over the time all ranks spent in the phase.
      const bool rate =
          (i < DBM_NUM_PHASES && 0 < bytes[i] && 0.0 < phase_avg[i]);
      const int n = snprintf(
          buffer, sizeof(buffer), "    %-24s %12.3f %12.3f %12.3f", phases[i],
          -phase_min[i], phase_avg[i] / nranks, phase_max[i]);
      if (rate) {
        snprintf(buffer + n, sizeof(buffer) - n, " %12.1f\n",
                 1E-9 * bytes[i] / phase_avg[i]);
      } else {
        snprintf(buffer + n, sizeof(buffer) - n, "\n");
      }
      DBM_LIBRARY_PRINT(print_func, buffer, output_unit);
    }
  }

  DBM_LIBRARY_PRINT(
      print_func,
      " ----------------------------------------------------------------"
      "---------------\n",
      output_unit);
}

// EOF
