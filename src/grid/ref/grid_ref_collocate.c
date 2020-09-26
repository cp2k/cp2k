/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define GRID_DO_COLLOCATE 1
#include "../common/grid_common.h"
#include "grid_ref_collint.h"
#include "grid_ref_collocate.h"
#include "grid_ref_prepare_pab.h"

/*******************************************************************************
 * \brief Collocates a single product of primitiv Gaussians.
 *        See grid_collocate.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_ref_collocate_pgf_product(
    const bool orthorhombic, const int border_mask, const int func,
    const int la_max, const int la_min, const int lb_max, const int lb_min,
    const double zeta, const double zetb, const double rscale,
    const double dh[3][3], const double dh_inv[3][3], const double ra[3],
    const double rab[3], const int npts_global[3], const int npts_local[3],
    const int shift_local[3], const int border_width[3], const double radius,
    const int o1, const int o2, const int n1, const int n2,
    const double pab[n2][n1], double *grid) {

  int la_min_diff, la_max_diff, lb_min_diff, lb_max_diff;
  grid_ref_prepare_get_ldiffs(func, &la_min_diff, &la_max_diff, &lb_min_diff,
                              &lb_max_diff);

  const int la_min_cab = imax(la_min + la_min_diff, 0);
  const int lb_min_cab = imax(lb_min + lb_min_diff, 0);
  const int la_max_cab = la_max + la_max_diff;
  const int lb_max_cab = lb_max + lb_max_diff;
  const int n1_cab = ncoset[la_max_cab];
  const int n2_cab = ncoset[lb_max_cab];

  const size_t cab_size = n2_cab * n1_cab;
  double cab[cab_size];
  memset(cab, 0, cab_size * sizeof(double));

  grid_ref_prepare_pab(func, o1, o2, la_max, la_min, lb_max, lb_min, zeta, zetb,
                       n1, n2, pab, n1_cab, n2_cab, (double(*)[n1_cab])cab);
  cab_to_grid(orthorhombic, border_mask, la_max_cab, la_min_cab, lb_max_cab,
              lb_min_cab, zeta, zetb, rscale, dh, dh_inv, ra, rab, npts_global,
              npts_local, shift_local, border_width, radius, cab, grid);
}

// EOF
