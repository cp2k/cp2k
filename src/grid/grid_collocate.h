/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2024 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#ifndef GRID_COLLOCATE_H
#define GRID_COLLOCATE_H

#include "common/grid_constants.h"
#include "grid_multigrid.h"

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
    const double pab[n2][n1], double *grid);

#endif

// EOF
