/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2025 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#ifndef GRID_REPLAY_H
#define GRID_REPLAY_H

#include <stdbool.h>

/*******************************************************************************
 * \brief Reads a .task file, collocates/integrates it, and compares results.
 *
 * \param filename          Name of the task file.
 * \param cycles            Number of times the task should be collocated.
 * \param collocate         When true collocate is called otherwise integrate.
 * \param batch             When false grid_ref_collocate_pgf_product is called.
 *                          When true grid_collocate_task_list is called.
 * \param cycles_per_block  Number of cycles per matrix block decontraction.
 * \param tolerance         Tolerance for comparing floating point results.
 * \returns                 Returns true iff the test passed.
 *
 * \author Ole Schuett
 ******************************************************************************/
bool grid_replay(const char *filename, const int cycles, const bool collocate,
                 const bool batch, const int cycles_per_block,
                 const double tolerance);

#endif

// EOF
