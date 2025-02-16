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

#ifndef LIBGRPP_PARAMETERS_H
#define LIBGRPP_PARAMETERS_H

enum { LIBGRPP_CART_ORDER_DIRAC, LIBGRPP_CART_ORDER_TURBOMOLE };

typedef struct {
  double radial_tolerance;
  double angular_screening_tolerance;
  double modified_bessel_tolerance;
  int (*cartesian_generator)(int L, int *cart_list);
} libgrpp_parameters_t;

extern libgrpp_parameters_t libgrpp_params;

void libgrpp_set_default_parameters();

void libgrpp_set_radial_tolerance(double tolerance);

void libgrpp_set_angular_screening_tolerance(double tolerance);

void libgrpp_set_modified_bessel_tolerance(double tolerance);

void libgrpp_set_cartesian_order(int order);

void libgrpp_set_cartesian_generator(
    int (*cartesian_generator)(int L, int *cart_list));

#endif // LIBGRPP_PARAMETERS_H
