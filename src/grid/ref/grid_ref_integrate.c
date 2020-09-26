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

#define GRID_DO_COLLOCATE 0
#include "../common/grid_common.h"
#include "grid_ref_collint.h"
#include "grid_ref_integrate.h"

/*******************************************************************************
 * \brief Integrates a single task. See grid_ref_integrate.h for details.
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
                            const double *grid, double *vab) {

  cab_to_grid(orthorhombic, border_mask, la_max, la_min, lb_max, lb_min, zeta,
              zetb, rscale, dh, dh_inv, ra, rab, npts_global, npts_local,
              shift_local, border_width, radius, vab, grid);
}

// EOF
