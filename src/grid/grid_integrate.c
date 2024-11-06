/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2024 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include "grid_integrate.h"
#include "cpu/grid_cpu_integrate.h"

#include <assert.h>
#include <stdio.h>

/*******************************************************************************
 * \brief Integrates a single task. Thin wrapper to prevent circular
 *dependencies. \author Frederick Stein
 ******************************************************************************/
void grid_integrate_pgf_product(
    const grid_multigrid *multigrid, const int ilevel, const bool compute_tau,
    const int border_mask, const int la_max, const int la_min, const int lb_max,
    const int lb_min, const double zeta, const double zetb, const double ra[3],
    const double rab[3], const double radius, const int o1, const int o2,
    const int n1, const int n2, double hab[n2][n1], const double pab[n2][n1],
    double forces[2][3], double virials[2][3][3], double hdab[n2][n1][3],
    double hadb[n2][n1][3], double a_hdab[n2][n1][3][3]) {

  assert(multigrid->cpu->singlegrids != NULL);
  assert(multigrid->cpu->singlegrids[ilevel - 1] != NULL);
  grid_cpu_integrate_pgf_product_multigrid(
      multigrid->cpu, ilevel, compute_tau, border_mask, la_max, la_min, lb_max,
      lb_min, zeta, zetb, ra, rab, radius, o1, o2, n1, n2,
      offload_get_buffer_host_pointer(multigrid->grids[ilevel - 1]), hab, pab,
      forces, virials, hdab, hadb, a_hdab);
}

// EOF
