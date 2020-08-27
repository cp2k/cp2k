/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2020  CP2K developers group                         *
 *****************************************************************************/

#include <assert.h>
#include <omp.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "grid_constants.h"
#include "grid_library.h"

static grid_library_stats **per_thread_stats = NULL;
static bool library_initialized = false;
static grid_library_config config = {.backend = GRID_BACKEND_AUTO,
                                     .validate = false};

//******************************************************************************
// \brief Initializes the grid library.
// \author Ole Schuett
//******************************************************************************
void grid_library_init() {
  if (library_initialized) {
    printf("Error: Grid library was already initialized.\n");
    abort();
  }

  per_thread_stats =
      malloc(sizeof(grid_library_stats *) * omp_get_max_threads());

// Using parallel regions to ensure memory is allocated near a thread's core.
#pragma omp parallel default(none) shared(per_thread_stats)
  {
    const int ithread = omp_get_thread_num();
    per_thread_stats[ithread] = malloc(sizeof(grid_library_stats));
    memset(per_thread_stats[ithread], 0, sizeof(grid_library_stats));
  }

  library_initialized = true;
}

//******************************************************************************
// \brief Finalizes the grid library.
// \author Ole Schuett
//******************************************************************************
void grid_library_finalize() {
  if (!library_initialized) {
    printf("Error: Grid library is not initialized.\n");
    abort();
  }

  for (int i = 0; i < omp_get_max_threads(); i++) {
    free(per_thread_stats[i]);
  }
  free(per_thread_stats);
  per_thread_stats = NULL;
  library_initialized = false;
}

//******************************************************************************
// \brief Configures the grid library.
// \author Ole Schuett
//******************************************************************************
void grid_library_set_config(const int backend, const bool validate) {
  config.backend = backend;
  config.validate = validate;
}

//******************************************************************************
// \brief Returns the library config.
// \author Ole Schuett
//******************************************************************************
grid_library_config grid_library_get_config() { return config; }

//******************************************************************************
// \brief Internal helper for summing two sets of counters.
// \author Ole Schuett
//******************************************************************************
static void sum_stats(const grid_library_stats increment,
                      grid_library_stats *accumulator) {
  accumulator->ref_collocate_ortho += increment.ref_collocate_ortho;
  accumulator->ref_collocate_general += increment.ref_collocate_general;
}

//******************************************************************************
// \brief Increment global counters by given values.
// \author Ole Schuett
//******************************************************************************
void grid_library_gather_stats(const grid_library_stats increment) {
  if (!library_initialized) {
    printf("Error: Grid library is not initialized.\n");
    abort();
  }
  sum_stats(increment, per_thread_stats[omp_get_thread_num()]);
}

//******************************************************************************
// \brief Prints statistics gathered by the grid library.
// \author Ole Schuett
//******************************************************************************
void grid_library_print_stats(void (*mpi_sum_func)(long *),
                              void (*print_func)(char *)) {
  if (!library_initialized) {
    printf("Error: Grid library is not initialized.\n");
    abort();
  }
  print_func("\n");
  print_func(" ----------------------------------------------------------------"
             "---------------\n");
  print_func(" -                                                               "
             "              -\n");
  print_func(" -                                GRID STATISTICS                "
             "              -\n");
  print_func(" -                                                               "
             "              -\n");
  print_func(" ----------------------------------------------------------------"
             "---------------\n");
  print_func(" COUNTER                                                         "
             "          VALUE\n");

  grid_library_stats totals;
  memset(&totals, 0, sizeof(grid_library_stats));

  for (int i = 0; i < omp_get_max_threads(); i++) {
    sum_stats(*per_thread_stats[i], &totals);
  }

  char buffer[100];
  mpi_sum_func(&totals.ref_collocate_ortho);
  snprintf(buffer, sizeof(buffer), " %-58s %20li\n", "ref_collocate_ortho",
           totals.ref_collocate_ortho);
  print_func(buffer);

  mpi_sum_func(&totals.ref_collocate_general);
  snprintf(buffer, sizeof(buffer), " %-58s %20li\n", "ref_collocate_general",
           totals.ref_collocate_general);
  print_func(buffer);

  print_func(" ----------------------------------------------------------------"
             "---------------\n");
}

// EOF
