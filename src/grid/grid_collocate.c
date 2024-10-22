/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2024 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include "grid_collocate.h"
#include "cpu/grid_cpu_collocate.h"

/*******************************************************************************
 * \brief Public entry point. A thin wrapper to prevent circular dependencies.
 * \author Frederick Stein
 ******************************************************************************/
void grid_collocate_pgf_product(
    const grid_multigrid *multigrid, const int ilevel, const int border_mask,
    const enum grid_func func, const int la_max, const int la_min,
    const int lb_max, const int lb_min, const double zeta, const double zetb,
    const double rscale, const double ra[3], const double rab[3],
    const double radius, const int o1, const int o2, const int n1, const int n2,
    const double pab[n2][n1], double *grid) {

  grid_cpu_collocate_pgf_product(multigrid->orthorhombic, border_mask, func,
                                 la_max, la_min, lb_max, lb_min, zeta, zetb,
                                 rscale, *(multigrid->dh + (ilevel - 1)),
                                 *(multigrid->dh_inv + (ilevel - 1)), ra, rab,
                                 *(multigrid->npts_global + (ilevel - 1)),
                                 *(multigrid->npts_local + (ilevel - 1)),
                                 *(multigrid->shift_local + (ilevel - 1)),
                                 *(multigrid->border_width + (ilevel - 1)),
                                 radius, o1, o2, n1, n2, pab, grid);
}

// EOF
