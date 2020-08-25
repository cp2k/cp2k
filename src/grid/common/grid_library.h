/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2020  CP2K developers group                         *
 *****************************************************************************/
#ifndef GRID_LIBRARY_H
#define GRID_LIBRARY_H

#include <stdbool.h>

//******************************************************************************
// \brief Initializes the grid library.
// \author Ole Schuett
//******************************************************************************
void grid_library_init();

//******************************************************************************
// \brief Finalizes the grid library.
// \author Ole Schuett
//******************************************************************************
void grid_library_finalize();

//******************************************************************************
// \brief Configuration of the grid library.
// \author Ole Schuett
//******************************************************************************
typedef struct {
  int backend;   // Selectes the backend to be used by the grid library.
  bool validate; // When true the reference backend runs in shadow mode.
} grid_library_config;

//******************************************************************************
// \brief Configures the grid library.
// \author Ole Schuett
//******************************************************************************
void grid_library_set_config(int backend, bool validate);

//******************************************************************************
// \brief Returns the library config.
// \author Ole Schuett
//******************************************************************************
grid_library_config grid_library_get_config();

//******************************************************************************
// \brief Prints statistics gathered by the grid library.
// \author Ole Schuett
//******************************************************************************
void grid_library_print_stats(void (*mpi_sum_func)(long *),
                              void (*print_func)(char *));

//******************************************************************************
// \brief All exiting counters. When adding a counter also update functions
//        internal_add_stats() and grid_library_print_counters().
// \author Ole Schuett
//******************************************************************************
typedef struct {
  long ref_collocate_ortho;
  long ref_collocate_general;
} grid_library_stats;

//******************************************************************************
// \brief Increment global counters by given values.
// \author Ole Schuett
//******************************************************************************
void grid_library_gather_stats(grid_library_stats increment);

#endif // GRID_LIBRARY_H

// EOF
