/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <omp.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "grid_common.h"
#include "grid_constants.h"
#include "grid_library.h"

static grid_library_globals **per_thread_globals = NULL;
static bool library_initialized = false;
static grid_library_config config = {.backend = GRID_BACKEND_AUTO,
                                     .device_id = -1,
                                     .validate = false,
                                     .apply_cutoff = false,
                                     .queue_length = 8192};

/*******************************************************************************
 * \brief Initializes the grid library.
 * \author Ole Schuett
 ******************************************************************************/
void grid_library_init() {
  if (library_initialized) {
    printf("Error: Grid library was already initialized.\n");
    abort();
  }

  per_thread_globals =
      malloc(sizeof(grid_library_globals *) * omp_get_max_threads());

// Using parallel regions to ensure memory is allocated near a thread's core.
#pragma omp parallel default(none) shared(per_thread_globals)
  {
    const int ithread = omp_get_thread_num();
    per_thread_globals[ithread] = malloc(sizeof(grid_library_globals));
    memset(per_thread_globals[ithread], 0, sizeof(grid_library_globals));
  }

  library_initialized = true;
}

/*******************************************************************************
 * \brief Finalizes the grid library.
 * \author Ole Schuett
 ******************************************************************************/
void grid_library_finalize() {
  if (!library_initialized) {
    printf("Error: Grid library is not initialized.\n");
    abort();
  }

  for (int i = 0; i < omp_get_max_threads(); i++) {
    grid_sphere_cache_free(&per_thread_globals[i]->sphere_cache);
    free(per_thread_globals[i]);
  }
  free(per_thread_globals);
  per_thread_globals = NULL;
  library_initialized = false;
}

/*******************************************************************************
 * \brief Returns a pointer to the thread local sphere cache.
 * \author Ole Schuett
 ******************************************************************************/
grid_sphere_cache *grid_library_get_sphere_cache() {
  return &per_thread_globals[omp_get_thread_num()]->sphere_cache;
}

/*******************************************************************************
 * \brief Configures the grid library.
 * \author Ole Schuett
 ******************************************************************************/
void grid_library_set_config(const int backend, const int device_id,
                             const bool validate, const bool apply_cutoff) {
  config.backend = backend;
  config.device_id = device_id;
  config.validate = validate;
  config.apply_cutoff = apply_cutoff;
}

/*******************************************************************************
 * \brief Returns the library config.
 * \author Ole Schuett
 ******************************************************************************/
grid_library_config grid_library_get_config() { return config; }

/*******************************************************************************
 * \brief Increment specified counter, see grid_library.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_library_increment_counter(int lp, int kern, int op) {
  lp = imin(lp, 19);
  per_thread_globals[omp_get_thread_num()]->counters[lp][kern][op]++;
}

/*******************************************************************************
 * \brief Prints statistics gathered by the grid library.
 * \author Ole Schuett
 ******************************************************************************/
void grid_library_print_stats(void (*mpi_sum_func)(long *, int),
                              const int mpi_comm,
                              void (*print_func)(char *, int),
                              const int output_unit) {
  if (!library_initialized) {
    printf("Error: Grid library is not initialized.\n");
    abort();
  }

  // Sum all counters across threads and mpi ranks.
  long counters[20][2][2] = {0};
  double total = 0.0;
  for (int lp = 0; lp < 20; lp++) {
    for (int kern = 0; kern < 2; kern++) {
      for (int op = 0; op < 2; op++) {
        for (int i = 0; i < omp_get_max_threads(); i++) {
          counters[lp][kern][op] +=
              per_thread_globals[i]->counters[lp][kern][op];
        }
        mpi_sum_func(&counters[lp][kern][op], mpi_comm);
        total += counters[lp][kern][op];
      }
    }
  }

  // Print counters.
  print_func("\n", output_unit);
  print_func(" ----------------------------------------------------------------"
             "---------------\n",
             output_unit);
  print_func(" -                                                               "
             "              -\n",
             output_unit);
  print_func(" -                                GRID STATISTICS                "
             "              -\n",
             output_unit);
  print_func(" -                                                               "
             "              -\n",
             output_unit);
  print_func(" ----------------------------------------------------------------"
             "---------------\n",
             output_unit);
  print_func(" LP    KERNEL    OPERATION                                     "
             "COUNT     PERCENT\n",
             output_unit);

  for (int lp = 0; lp < 20; lp++) {
    for (int kern = 0; kern < 2; kern++) {
      for (int op = 0; op < 2; op++) {
        if (counters[lp][kern][op] == 0)
          continue; // skip empty counters
        const char *op_str = (op == 1) ? "collocate" : "integrate";
        const char *kern_str = (kern == 1) ? "ortho" : "general";
        char buffer[100];
        double percent = 100.0 * counters[lp][kern][op] / total;
        snprintf(buffer, sizeof(buffer), " %-5i %-9s %-12s %38li %10.2f%%\n",
                 lp, kern_str, op_str, counters[lp][kern][op], percent);
        print_func(buffer, output_unit);
      }
    }
  }

  print_func(" ----------------------------------------------------------------"
             "---------------\n",
             output_unit);
}

// EOF
