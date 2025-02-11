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

double spherical_to_cartesian_coef(int l, int m, int lx, int ly, int lz);

double evaluate_real_spherical_harmonic(int l, int m, const double *k);

void evaluate_real_spherical_harmonics_array(int l, const double *k,
                                             double *rsh_array);

void create_real_spherical_harmonic_coeffs_tables(int Lmax);

rsh_coef_table_t *get_real_spherical_harmonic_table(int L);

#endif // LIBGRPP_SPHERICAL_HARMONICS_H
