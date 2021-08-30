/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2021 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#ifndef GRID_COLLOCATE_DGEMM_H
#define GRID_COLLOCATE_DGEMM_H

#include "../common/grid_constants.h"
#ifdef __cplusplus
extern "C" {
#endif

void grid_collocate_pgf_product_cpu_dgemm(
    const bool orthorhombic, const int border_mask, const enum grid_func func,
    const int la_max, const int la_min, const int lb_max, const int lb_min,
    const double zeta, const double zetb, const double rscale,
    const double dh[3][3], const double dh_inv[3][3], const double ra[3],
    const double rab[3], const int npts_global[3], const int npts_local[3],
    const int shift_local[3], const int border_width[3], const double radius,
    const int o1, const int o2, const int n1, const int n2,
    const double pab_[n2][n1], double *const grid);
#ifdef __cplusplus
}
#endif

#endif
// EOF
