/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2021 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#ifndef GRID_REF_INTEGRATE_H
#define GRID_REF_INTEGRATE_H

#include <stdbool.h>

/*******************************************************************************
 * \brief Integrates a single task. A task consists of a pair of atoms each
 *        with a position, Gaussian exponent, and a range of angular momentum.
 *        This function then integrates all combinations of spherical harmonics.
 *        Arguments are identical with grid_collocate_pgf_product except for:
 *
 * \param grid          Input grid array.
 * \param hab           Output Hamiltonian matrix block.
 *
 * \param pab           Optional input density matrix block.
 * \param forces        Optional output forces, requires pab.
 * \param virials       Optional output virials, requires pab.
 * \param hdab          Optional output derivative d(hab)/da, requires pab.
 * \param a_hdab        Optional output virial of hab, requires pab.
 *
 * \author Ole Schuett
 ******************************************************************************/
void grid_ref_integrate_pgf_product(
    const bool orthorhombic, const bool compute_tau, const int border_mask,
    const int la_max, const int la_min, const int lb_max, const int lb_min,
    const double zeta, const double zetb, const double dh[3][3],
    const double dh_inv[3][3], const double ra[3], const double rab[3],
    const int npts_global[3], const int npts_local[3], const int shift_local[3],
    const int border_width[3], const double radius, const int o1, const int o2,
    const int n1, const int n2, const double *grid, double hab[n2][n1],
    const double pab[n2][n1], double forces[2][3], double virials[2][3][3],
    double hdab[n2][n1][3], double hadb[n2][n1][3],
    double a_hdab[n2][n1][3][3]);

#endif
// EOF
