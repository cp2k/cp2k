/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2020  CP2K developers group                         *
 *****************************************************************************/
#ifndef GRID_PREPARE_PAB_H
#define GRID_PREPARE_PAB_H

#include "../common/grid_constants.h"
#include "tensor_local.h"
void grid_prepare_get_ldiffs_dgemm(const int func, int *const lmin_diff,
                                   int *const lmax_diff);

void grid_prepare_pab_dgemm(const int func, const int *const offset,
                            const int *const lmax, const int *const lmin,
                            const double *const zeta, tensor *const pab,
                            tensor *const pab_prep);

#endif

// EOF
