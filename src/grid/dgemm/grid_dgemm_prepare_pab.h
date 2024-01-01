/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2024 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#ifndef GRID_DGEMM_PREPARE_PAB_H
#define GRID_DGEMM_PREPARE_PAB_H

#include "../common/grid_constants.h"
#include "grid_dgemm_tensor_local.h"
void grid_prepare_get_ldiffs_dgemm(const enum grid_func func,
                                   int *const lmin_diff, int *const lmax_diff);

void grid_prepare_pab_dgemm(const enum grid_func func, const int *const offset,
                            const int *const lmax, const int *const lmin,
                            const double *const zeta, tensor *const pab,
                            tensor *const pab_prep);

#endif
