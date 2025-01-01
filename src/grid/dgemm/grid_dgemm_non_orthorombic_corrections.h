/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2025 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#ifndef GRID_DGEMM_NON_ORTHOROMBIC_CORRECTIONS_H
#define GRID_DGEMM_NON_ORTHOROMBIC_CORRECTIONS_H

#include <stdbool.h>

#include "grid_dgemm_tensor_local.h"
extern void calculate_non_orthorombic_corrections_tensor(
    const double mu_mean, const double *r_ab, const double basis[3][3],
    const int *const xmin, const int *const xmax, bool plane[3],
    tensor *const Exp);
extern void calculate_non_orthorombic_corrections_tensor_blocked(
    const double mu_mean, const double *r_ab, const double basis[3][3],
    const int *const lower_block, const int *const upper_block,
    const int *const block_size, const int *const offset, const int *const xmin,
    const int *const xmax, bool *plane, tensor *const Exp);

extern void apply_non_orthorombic_corrections(const bool plane[3],
                                              const tensor *const Exp,
                                              tensor *const cube);
extern void
apply_non_orthorombic_corrections_xy(const int x, const int y,
                                     const struct tensor_ *const Exp,
                                     struct tensor_ *const m);
extern void
apply_non_orthorombic_corrections_xz(const int x, const int z,
                                     const struct tensor_ *const Exp,
                                     struct tensor_ *const m);
extern void
apply_non_orthorombic_corrections_yz(const int y, const int z,
                                     const struct tensor_ *const Exp,
                                     struct tensor_ *const m);
extern void
apply_non_orthorombic_corrections_xy_blocked(const struct tensor_ *const Exp,
                                             struct tensor_ *const m);
extern void
apply_non_orthorombic_corrections_xz_blocked(const struct tensor_ *const Exp,
                                             struct tensor_ *const m);
extern void
apply_non_orthorombic_corrections_yz_blocked(const struct tensor_ *const Exp,
                                             struct tensor_ *const m);
extern void apply_non_orthorombic_corrections_xz_yz_blocked(
    const struct tensor_ *const Exp_xz, const struct tensor_ *const Exp_yz,
    struct tensor_ *const m);
#endif
