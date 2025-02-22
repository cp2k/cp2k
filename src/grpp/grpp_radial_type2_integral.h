/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2025 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: MIT                                              */
/*----------------------------------------------------------------------------*/

/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2023 Alexander Oleynichenko
 */

#ifndef LIBGRPP_RADIAL_TYPE2_INTEGRAL_H
#define LIBGRPP_RADIAL_TYPE2_INTEGRAL_H

#include "libgrpp_types.h"

typedef struct {
  int lambda1_max;
  int lambda2_max;
  int n_max;
  double *radial_integrals;
} radial_type2_table_t;

radial_type2_table_t *libgrpp_tabulate_radial_type2_integrals(
    int lambda1_max, int lambda2_max, int n_max, double CA_2, double CB_2,
    libgrpp_potential_t *potential, libgrpp_shell_t *bra, libgrpp_shell_t *ket);

double libgrpp_get_radial_type2_integral(radial_type2_table_t *table,
                                         int lambda1, int lambda2, int n);

void libgrpp_delete_radial_type2_integrals(radial_type2_table_t *table);

#endif // LIBGRPP_RADIAL_TYPE2_INTEGRAL_H
