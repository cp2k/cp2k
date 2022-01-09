/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "dbm_library.h"
#include "dbm_mempool.h"

static bool library_initialized = false;

#if !defined(_OPENMP)
#error "OpenMP is required. Please add -fopenmp to your C compiler flags."
#endif

/*******************************************************************************
 * \brief Initializes the DBM library.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_library_init(void) {
  if (library_initialized) {
    fprintf(stderr, "DBM library was already initialized.\n");
    abort();
  }

  // Nothing to do, yet.
  library_initialized = true;
}

/*******************************************************************************
 * \brief Finalizes the DBM library.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_library_finalize(void) {
  if (!library_initialized) {
    fprintf(stderr, "Error: DBM library is not initialized.\n");
    abort();
  }

  dbm_mempool_clear();
  library_initialized = false;
}

// EOF
