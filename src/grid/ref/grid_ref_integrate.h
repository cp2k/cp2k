/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#ifndef GRID_REF_INTEGRATE_H
#define GRID_REF_INTEGRATE_H

#include <stdbool.h>

/*******************************************************************************
 * \brief Integrates a single task. A task consists of a pair of atoms each
 *        with a position, Gaussian exponent, and a range of angular momentum.
 *        This function then integrates all combinations of spherical harmonics.
 *        The arguments are mostly identical with grid_collocate_pgf_product.
 *
 * \param grid          The input grid array to integrate from.
 * \param vab           The output atom-pair's potential matrix block V_{ab}.
 *
 * \author Ole Schuett
 ******************************************************************************/
void grid_ref_integrate_vab(const bool orthorhombic, const int border_mask,
                            const int la_max, const int la_min,
                            const int lb_max, const int lb_min,
                            const double zeta, const double zetb,
                            const double rscale, const double dh[3][3],
                            const double dh_inv[3][3], const double ra[3],
                            const double rab[3], const int npts_global[3],
                            const int npts_local[3], const int shift_local[3],
                            const int border_width[3], const double radius,
                            const double *grid, double *vab);

#endif
// EOF
