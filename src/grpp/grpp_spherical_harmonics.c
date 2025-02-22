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

/*
 * Constructs tables with the expansion coefficients of real spherical harmonics
 * in the basis of (cartesian) unitary spherical polynomials.
 *
 * For more details about the algorithm used, see:
 * R. Flores-Moreno et al, J. Comput. Chem. 27, 1009 (2006),
 * doi: 10.1002/jcc.20410
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440
#endif

#include "grpp_binomial.h"
#include "grpp_factorial.h"
#include "grpp_spherical_harmonics.h"
#include "libgrpp.h"
/*
 * Tables with pretabulated expansion coefficients
 */

static rsh_coef_table_t **rsh_coef_tables = NULL;
static int rsh_tables_lmax = -1;

/*
 * Function pre-definitions
 */
rsh_coef_table_t *libgrpp_tabulate_real_spherical_harmonic_coeffs(int L);

static int *generate_cartesian_combinations(int L, int *num);

/**
 * Constructs the set of tables with C_{l,m}^{lx,ly,lz} coefficients
 * (up to maximum angular momentum Lmax).
 * (pretabulation step)
 */
void libgrpp_create_real_spherical_harmonic_coeffs_tables(int Lmax) {
  if (Lmax <= rsh_tables_lmax) {
    // nothing to do
  } else {
    // expand tables: realloc memory and add tables for the highest L values
    rsh_coef_tables = (rsh_coef_table_t **)realloc(
        rsh_coef_tables, (Lmax + 1) * sizeof(rsh_coef_table_t *));

    for (int L = rsh_tables_lmax + 1; L <= Lmax; L++) {
      rsh_coef_tables[L] = libgrpp_tabulate_real_spherical_harmonic_coeffs(L);
    }
    rsh_tables_lmax = Lmax;
  }
}

/**
 * Calculates all C_{l,m}^{lx,ly,lz} coefficients for the given angular momentum
 * L.
 */
rsh_coef_table_t *libgrpp_tabulate_real_spherical_harmonic_coeffs(int L) {
  int ncart = (L + 1) * (L + 2) / 2;

  rsh_coef_table_t *coef_table =
      (rsh_coef_table_t *)calloc(1, sizeof(rsh_coef_table_t));
  coef_table->n_cart_comb = ncart;
  coef_table->cartesian_comb = generate_cartesian_combinations(L, &ncart);
  coef_table->coeffs = (double *)calloc((2 * L + 1) * ncart, sizeof(double));

  for (int m = -L; m <= L; m++) {
    for (int icomb = 0; icomb < ncart; icomb++) {
      int lx = coef_table->cartesian_comb[3 * icomb];
      int ly = coef_table->cartesian_comb[3 * icomb + 1];
      // int lz = coef_table->cartesian_comb[3 * icomb + 2];
      double u_lm_lx_ly_lz = libgrpp_spherical_to_cartesian_coef(L, m, lx, ly);
      int index = (m + L) * ncart + icomb;
      coef_table->coeffs[index] = u_lm_lx_ly_lz;
    }
  }

  return coef_table;
}

/**
 * Access to the table for the angular momentum value L.
 */
rsh_coef_table_t *libgrpp_get_real_spherical_harmonic_table(int L) {
  if (L > rsh_tables_lmax) {
    printf("get_real_spherical_harmonic_table(): %d > Lmax\n", L);
    return NULL;
  }

  return rsh_coef_tables[L];
}

/**
 * For the given real spherical harmonic (RSH) S_lm calculates the coefficient
 * C_{l,m}^{lx,ly,lz} before the unitary spherical polynomial (USP) in its
 * expansion.
 *
 * The formula is taken from:
 * R. Flores-Moreno et al, J. Comput. Chem. 27, 1009 (2006)
 * doi: 10.1002/jcc.20410
 * (formula 32)
 */
double libgrpp_spherical_to_cartesian_coef(int l, int m, int lx, int ly) {
  int j = lx + ly - abs(m);
  if (j % 2 != 0) {
    return 0.0;
  }
  j /= 2;

  if (!((m > 0 && (abs(m) - lx) % 2 == 0) || (m == 0 && lx % 2 == 0) ||
        (m < 0 && (abs(m) - lx) % 2 != 0))) {
    return 0.0;
  }

  double prefactor =
      sqrt((2 * l + 1) / (2 * M_PI) * libgrpp_factorial(l - abs(m)) /
           libgrpp_factorial(l + abs(m)));
  prefactor /= pow(2, l) * libgrpp_factorial(l);

  double u_lm_lx_ly_lz = 0.0;
  for (int i = j; i <= (l - abs(m)) / 2; i++) {
    // any term that implies the factorial of a negative number is neglected
    if (2 * l - 2 * i < 0) {
      u_lm_lx_ly_lz = 0.0;
      break;
    }
    if (l - abs(m) - 2 * i < 0) {
      u_lm_lx_ly_lz = 0.0;
      break;
    }

    double factor_1 =
        libgrpp_binomial(l, i) * libgrpp_binomial(i, j) * pow(-1, i) *
        libgrpp_factorial_ratio(2 * l - 2 * i, l - abs(m) - 2 * i);

    double sum = 0.0;
    for (int k = 0; k <= j; k++) {
      sum += libgrpp_binomial(j, k) * libgrpp_binomial(abs(m), lx - 2 * k) *
             pow(-1, (abs(m) - lx + 2 * k) / 2);
    }

    u_lm_lx_ly_lz += factor_1 * sum;
  }

  u_lm_lx_ly_lz *= prefactor;
  if (m == 0 && (lx % 2 == 0)) {
    u_lm_lx_ly_lz *= M_SQRT1_2; // x 1/sqrt(2)
  }

  return u_lm_lx_ly_lz;
}

/**
 * Calculates value of the real spherical harmonic S_lm at the point k/|k| of
 * the unit sphere.
 */
double libgrpp_evaluate_real_spherical_harmonic(const int l, const int m,
                                                const double *k) {
  double unitary_kx;
  double unitary_ky;
  double unitary_kz;
  double kx_powers[200];
  double ky_powers[200];
  double kz_powers[200];

  rsh_coef_table_t *rsh_coef_l = libgrpp_get_real_spherical_harmonic_table(l);

  double length_k = sqrt(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]);

  if (length_k > LIBGRPP_ZERO_THRESH) {
    unitary_kx = k[0] / length_k;
    unitary_ky = k[1] / length_k;
    unitary_kz = k[2] / length_k;
  } else {
    unitary_kx = 0.0;
    unitary_ky = 0.0;
    unitary_kz = 0.0;
  }

  kx_powers[0] = 1.0;
  ky_powers[0] = 1.0;
  kz_powers[0] = 1.0;

  for (int i = 1; i <= l; i++) {
    kx_powers[i] = kx_powers[i - 1] * unitary_kx;
    ky_powers[i] = ky_powers[i - 1] * unitary_ky;
    kz_powers[i] = kz_powers[i - 1] * unitary_kz;
  }

  double value = 0.0;

  int ncart = rsh_coef_l->n_cart_comb;
  for (int icomb = 0; icomb < ncart; icomb++) {
    int r = rsh_coef_l->cartesian_comb[3 * icomb];
    int s = rsh_coef_l->cartesian_comb[3 * icomb + 1];
    int t = rsh_coef_l->cartesian_comb[3 * icomb + 2];
    double y_lm_rst = rsh_coef_l->coeffs[(m + l) * ncart + icomb];

    value += y_lm_rst * kx_powers[r] * ky_powers[s] * kz_powers[t];
  }

  return value;
}

/**
 * Calculates values of the real spherical harmonic S_lm at the point k/|k| of
 * the unit sphere for all m = -l, ..., +l
 */
void libgrpp_evaluate_real_spherical_harmonics_array(const int l,
                                                     const double *k,
                                                     double *rsh_array) {
  double unitary_kx;
  double unitary_ky;
  double unitary_kz;
  double kx_powers[200];
  double ky_powers[200];
  double kz_powers[200];

  rsh_coef_table_t *rsh_coef_l = libgrpp_get_real_spherical_harmonic_table(l);

  double length_k = sqrt(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]);

  if (length_k > LIBGRPP_ZERO_THRESH) {
    double inv_length = 1.0 / length_k;
    unitary_kx = k[0] * inv_length;
    unitary_ky = k[1] * inv_length;
    unitary_kz = k[2] * inv_length;
  } else {
    unitary_kx = 0.0;
    unitary_ky = 0.0;
    unitary_kz = 0.0;
  }

  kx_powers[0] = 1.0;
  ky_powers[0] = 1.0;
  kz_powers[0] = 1.0;

  for (int i = 1; i <= l; i++) {
    kx_powers[i] = kx_powers[i - 1] * unitary_kx;
    ky_powers[i] = ky_powers[i - 1] * unitary_ky;
    kz_powers[i] = kz_powers[i - 1] * unitary_kz;
  }

  memset(rsh_array, 0, (2 * l + 1) * sizeof(double));

  int ncart = rsh_coef_l->n_cart_comb;
  int *rst_array = rsh_coef_l->cartesian_comb;

  for (int icomb = 0; icomb < ncart; icomb++) {
    int r = rst_array[3 * icomb];
    int s = rst_array[3 * icomb + 1];
    int t = rst_array[3 * icomb + 2];

    double k_xyz = kx_powers[r] * ky_powers[s] * kz_powers[t];

    for (int m = -l; m <= l; m++) {
      double y_lm_rst = rsh_coef_l->coeffs[(m + l) * ncart + icomb];
      rsh_array[m + l] += y_lm_rst * k_xyz;
    }
  }
}

static int *generate_cartesian_combinations(int L, int *num) {
  *num = (L + 1) * (L + 2) / 2;

  int *combinations = (int *)calloc(*num, 3 * sizeof(int));

  int n = 0;
  for (int i = 0; i <= L; i++) {
    for (int j = 0; j <= L; j++) {
      for (int k = 0; k <= L; k++) {
        if (i + j + k == L) {
          combinations[3 * n + 0] = i;
          combinations[3 * n + 1] = j;
          combinations[3 * n + 2] = k;
          n++;
        }
      }
    }
  }

  return combinations;
}
