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

#ifndef LIBGRPP_SPHERICAL_HARMONICS_H
#define LIBGRPP_SPHERICAL_HARMONICS_H

/*
 * Tables with pretabulated expansion coefficients
 */
typedef struct {
  int L;
  int n_cart_comb;
  int *cartesian_comb;
  double *coeffs;
} rsh_coef_table_t;

double libgrpp_spherical_to_cartesian_coef(int l, int m, int lx, int ly);

double libgrpp_evaluate_real_spherical_harmonic(int l, int m, const double *k);

void libgrpp_evaluate_real_spherical_harmonics_array(int l, const double *k,
                                                     double *rsh_array);

void libgrpp_create_real_spherical_harmonic_coeffs_tables(int Lmax);

rsh_coef_table_t *libgrpp_get_real_spherical_harmonic_table(int L);

#endif // LIBGRPP_SPHERICAL_HARMONICS_H
