/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#ifndef DBM_LIBRARY_H
#define DBM_LIBRARY_H

#include "dbm_multiply.h"

/*******************************************************************************
 * \brief Initializes the DBM library.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_library_init(void);

/*******************************************************************************
 * \brief Finalizes the DBM library.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_library_finalize(void);

/*******************************************************************************
 * \brief Add given block multiplication to the stats.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_library_counter_increment(const int m, const int n, const int k);

/*******************************************************************************
 * \brief Prints statistics gathered by the DBM library.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_library_print_stats(const int fortran_comm,
                             void (*print_func)(char *, int),
                             const int output_unit);

#endif

// EOF
