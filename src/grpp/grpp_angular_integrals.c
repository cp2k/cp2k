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
 * This file contains subroutines for evaluation of angular integrals of the 1st
 * and 2nd type. It also contains functions for construction of matrices of the
 * angular momentum operator in the bases of either real or complex spherical
 * harmonics.
 *
 * For more details on the angular integrals used in RPP integration, see:
 *
 * L. E. McMurchie, E. R. Davidson. Calculation of integrals over ab initio
 * pseudopotentials. J. Comput. Phys. 44(2), 289 (1981).
 * doi: 10.1016/0021-9991(81)90053-x
 */

#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "grpp_angular_integrals.h"
#include "grpp_factorial.h"
#include "grpp_spherical_harmonics.h"
#include "libgrpp.h"

static double integrate_unitary_sphere_polynomial(int i, int j, int k);

/**
 * Type 1 angular integral.
 * (see MrMurchie & Davidson, formula (28))
 */
double libgrpp_angular_type1_integral(int lambda, int II, int JJ, int KK,
                                      double *k) {
  double sum = 0.0;

  for (int mu = -lambda; mu <= lambda; mu++) {
    double sum2 = 0.0;
    for (int r = 0; r <= lambda; r++) {
      for (int s = 0; s <= lambda; s++) {
        for (int t = 0; t <= lambda; t++) {
          if (r + s + t == lambda) {
            double y_lm_rst =
                libgrpp_spherical_to_cartesian_coef(lambda, mu, r, s);
            double usp_int =
                integrate_unitary_sphere_polynomial(II + r, JJ + s, KK + t);
            sum2 += y_lm_rst * usp_int;
          }
        }
      }
    }
    sum += sum2 * libgrpp_evaluate_real_spherical_harmonic(lambda, mu, k);
  }

  return sum;
}

/**
 * Type 2 angular integral.
 * (see MrMurchie & Davidson, formula (29))
 */
double libgrpp_angular_type2_integral(const int lambda, const int L,
                                      const int m, const int a, const int b,
                                      const int c, const double *rsh_values) {
  double sum = 0.0;

  rsh_coef_table_t *rsh_coef_lambda =
      libgrpp_get_real_spherical_harmonic_table(lambda);
  rsh_coef_table_t *rsh_coef_L = libgrpp_get_real_spherical_harmonic_table(L);

  int ncomb_rst = rsh_coef_lambda->n_cart_comb;
  int ncomb_uvw = rsh_coef_L->n_cart_comb;

  for (int mu = -lambda; mu <= lambda; mu++) {
    double sum2 = 0.0;
    double rsh_value_k = rsh_values[mu + lambda];
    if (fabs(rsh_value_k) < LIBGRPP_ZERO_THRESH) {
      continue;
    }

    for (int icomb_uvw = 0; icomb_uvw < ncomb_uvw; icomb_uvw++) {

      int u = rsh_coef_L->cartesian_comb[3 * icomb_uvw];
      int v = rsh_coef_L->cartesian_comb[3 * icomb_uvw + 1];
      int w = rsh_coef_L->cartesian_comb[3 * icomb_uvw + 2];
      double y_lm_uvw = rsh_coef_L->coeffs[(m + L) * ncomb_uvw + icomb_uvw];
      if (fabs(y_lm_uvw) < LIBGRPP_ZERO_THRESH) {
        continue;
      }

      for (int icomb_rst = 0; icomb_rst < ncomb_rst; icomb_rst++) {

        int r = rsh_coef_lambda->cartesian_comb[3 * icomb_rst];
        int s = rsh_coef_lambda->cartesian_comb[3 * icomb_rst + 1];
        int t = rsh_coef_lambda->cartesian_comb[3 * icomb_rst + 2];
        double y_lam_mu_rst =
            rsh_coef_lambda->coeffs[(mu + lambda) * ncomb_rst + icomb_rst];
        if (fabs(y_lam_mu_rst) < LIBGRPP_ZERO_THRESH) {
          continue;
        }

        double usp_int = integrate_unitary_sphere_polynomial(
            a + r + u, b + s + v, c + t + w);

        sum2 += y_lam_mu_rst * y_lm_uvw * usp_int;
      }
    }

    sum += sum2 * rsh_value_k;
  }

  return sum;
}

/**
 * Integral of the unitary sphere polynomial over full solid angle.
 * (see MrMurchie & Davidson, formula (30))
 */
static double integrate_unitary_sphere_polynomial(int i, int j, int k) {
  if ((i % 2 == 0) && (j % 2 == 0) && (k % 2 == 0)) {
    double dfac_i = (double)libgrpp_double_factorial(i - 1);
    double dfac_j = (double)libgrpp_double_factorial(j - 1);
    double dfac_k = (double)libgrpp_double_factorial(k - 1);
    double dfac_ijk = (double)libgrpp_double_factorial(i + j + k + 1);
    return 4 * M_PI * dfac_i * dfac_j * dfac_k / dfac_ijk;
  } else {
    return 0.0;
  }
}
