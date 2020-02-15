/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2020  CP2K developers group                         *
 *****************************************************************************/

#ifndef GRID_COLLOCATE_REPLAY_H
#define GRID_COLLOCATE_REPLAY_H

#include <stdbool.h>

void grid_collocate_record(const bool orthorhombic,
                           const bool use_subpatch,
                           const int subpatch,
                           const int border,
                           const int func,
                           const int la_max,
                           const int la_min,
                           const int lb_max,
                           const int lb_min,
                           const double zeta,
                           const double zetb,
                           const double rscale,
                           const double dh[3][3],
                           const double dh_inv[3][3],
                           const double ra[3],
                           const double rab[3],
                           const int npts_global[3],
                           const int npts_local[3],
                           const int shift_local[3],
                           const double radius,
                           const int o1,
                           const int o2,
                           const int n1,
                           const int n2,
                           const double pab[n2][n1],
                           const double* grid);

double grid_collocate_replay(const char* filename, int cycles);

#endif

//EOF
