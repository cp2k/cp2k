/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2025 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#ifndef GRID_REF_PREPARE_PAB_H
#define GRID_REF_PREPARE_PAB_H

#include "../common/grid_constants.h"

/*******************************************************************************
 * \brief Returns block size changes due to transformation grid_prepare_pab.
 * \author Ole Schuett
 ******************************************************************************/
void grid_ref_prepare_get_ldiffs(const enum grid_func func, int *la_min_diff,
                                 int *la_max_diff, int *lb_min_diff,
                                 int *lb_max_diff);

/*******************************************************************************
 * \brief Selects and transforms a sub-block of the given density matrix block.
 *
 * \param func          Transformation function to apply, one of GRID_FUNC_*.
 * \param o{1,2}        Offsets of the sub-block within the matrix block.
 * \param l{a,b}_max    Max angular momentum to collocate for give atom.
 * \param l{a,b}_min    Lowest angular momentum to collocate for give atom.
 * \param zet_{a,b}     Gaussian's exponent of given atom.
 * \param n{1,2}        Dimensions of input matrix block.
 * \param pab           Input matrix block.
 * \param n{1,2}_prep   Dimensions of the transformed matrix sub-block.
 * \param pab_prep      Resulting transformed matrix sub-block.
 *
 * \author Ole Schuett
 ******************************************************************************/
void grid_ref_prepare_pab(const enum grid_func func, const int o1, const int o2,
                          const int la_max, const int la_min, const int lb_max,
                          const int lb_min, const double zeta,
                          const double zetb, const int n1, const int n2,
                          const double pab[n2][n1], const int n1_prep,
                          const int n2_prep, double pab_prep[n2_prep][n1_prep]);

#endif

// EOF
