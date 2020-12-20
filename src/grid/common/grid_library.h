/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/
#ifndef GRID_LIBRARY_H
#define GRID_LIBRARY_H

#ifdef __cplusplus
extern "C" {
#endif

#include "grid_constants.h"
#include "grid_sphere_cache.h"
#include <stdbool.h>

/*******************************************************************************
 * \brief Initializes the grid library.
 * \author Ole Schuett
 ******************************************************************************/
void grid_library_init();

/*******************************************************************************
 * \brief Finalizes the grid library.
 * \author Ole Schuett
 ******************************************************************************/
void grid_library_finalize();

/*******************************************************************************
 * \brief Configuration of the grid library.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  enum grid_backend
      backend;       // Selectes the backend to be used by the grid library.
  int device_id;     // gpu id
  bool validate;     // When true the reference backend runs in shadow mode.
  bool apply_cutoff; // only important for the dgemm and gpu backends
  int queue_length;  // Length of the queue for the gpu backend
} grid_library_config;

/*******************************************************************************
 * \brief Configures the grid library.
 * \author Ole Schuett
 ******************************************************************************/
void grid_library_set_config(const enum grid_backend backend,
                             const int device_id, const bool validate,
                             const bool apply_cutoff, const int queue_length);

/*******************************************************************************
 * \brief Returns the library config.
 * \author Ole Schuett
 ******************************************************************************/
grid_library_config grid_library_get_config();

/*******************************************************************************
 * \brief Prints statistics gathered by the grid library.
 * \author Ole Schuett
 ******************************************************************************/
void grid_library_print_stats(void (*mpi_sum_func)(long *, int), int mpi_comm,
                              void (*print_func)(char *, int), int output_unit);

/*******************************************************************************
 * \brief All exiting counters. When adding a counter also update functions
 *        internal_add_stats() and grid_library_print_counters().
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  grid_sphere_cache sphere_cache;
  long counters[20][2][2]; // [lp][kernel][operation]
} grid_library_globals;

/*******************************************************************************
 * \brief Returns a pointer to the thread local sphere cache.
 * \author Ole Schuett
 ******************************************************************************/
grid_sphere_cache *grid_library_get_sphere_cache();

/*******************************************************************************
 * \brief Increment specified counter.
 * \param lp    Total angular momentum: la_max + lb_max.
 * \param kern  Kernel: 0=ortho, 1=general.
 * \param op    Operation: 0=integrate, 1=collocate.
 * \author Ole Schuett
 ******************************************************************************/
void grid_library_increment_counter(int collocate, int ortho, int lp);

#ifdef __cplusplus
}
#endif

#endif // GRID_LIBRARY_H

// EOF
