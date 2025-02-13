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

#ifndef LIBGRPP_RADIAL_TYPE1_INTEGRAL_H
#define LIBGRPP_RADIAL_TYPE1_INTEGRAL_H

#include "libgrpp.h"

typedef struct {
  int lambda_max;
  int n_max;
  double *radial_integrals;
} radial_type1_table_t;

radial_type1_table_t *libgrpp_tabulate_radial_type1_integrals(
    int lambda_max, int n_max, double CA_2, double CB_2, double alpha_A,
    double alpha_B, double k, double prefactor,
    double (*potential)(double r, void *params), void *potential_params);

void libgrpp_delete_radial_type1_integrals(radial_type1_table_t *table);

double libgrpp_get_radial_type1_integral(radial_type1_table_t *table,
                                         int lambda, int n);

#endif // LIBGRPP_RADIAL_TYPE1_INTEGRAL_H
