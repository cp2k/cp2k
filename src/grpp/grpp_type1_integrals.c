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

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "grpp_angular_integrals.h"
#include "grpp_binomial.h"
#include "grpp_norm_gaussian.h"
#include "grpp_radial_type1_integral.h"
#include "grpp_type1_mcmurchie_davidson.h"
#include "grpp_utils.h"
#include "libgrpp.h"
#include "libgrpp_types.h"
/* for the old (numerical) version:
void evaluate_type1_integral_primitive_gaussians(double *A, int n_cart_A, int
*cart_list_A, double alpha_A, double *B, int n_cart_B, int *cart_list_B, double
alpha_B, double *C, libgrpp_potential_t *potential, double *matrix);
*/

extern void libgrpp_delete_radial_type1_integrals(radial_type1_table_t *table);

void libgrpp_evaluate_radially_local_potential_integral_primitive_gaussians(
    double *A, int n_cart_A, int *cart_list_A, double alpha_A, double *B,
    int n_cart_B, int *cart_list_B, double alpha_B, double *C,
    double (*potential)(double r, void *params), void *potential_params,
    double *matrix);

static double evaluate_pseudopotential(double r, void *params);

/**
 * Evaluation of type 1 RPP integrals (scalar-relativistic radially local RPP).
 */
void libgrpp_type1_integrals(libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B,
                             double *rpp_origin, libgrpp_potential_t *potential,
                             double *matrix) {
  assert(libgrpp_is_initialized());

  int size_A = libgrpp_get_shell_size(shell_A);
  int size_B = libgrpp_get_shell_size(shell_B);
  memset(matrix, 0, size_A * size_B * sizeof(double));

  if (potential == NULL) {
    return;
  }

  /*
   * RPP terms with n = 1, 2 are evaluated in a completely analytic manner
   * using the Obara-Saika-like recurrence relations
   */

  double *buf = calloc(size_A * size_B, sizeof(double));

  for (int k = 0; k < potential->num_primitives; k++) {
    double pot_coef = potential->coeffs[k];
    double pot_alpha = potential->alpha[k];
    int pot_n = potential->powers[k];

    libgrpp_type1_integrals_mcmurchie_davidson_1978(
        shell_A, shell_B, rpp_origin, pot_alpha, pot_n, buf);

    libgrpp_daxpy(size_A * size_B, pot_coef, buf, matrix);
  }

  free(buf);

  /*
   * old (numerical) version
   */

  /*for (int i = 0; i < shell_A->num_primitives; i++) {
      for (int j = 0; j < shell_B->num_primitives; j++) {
          double coef_A_i = shell_A->coeffs[i];
          double coef_B_j = shell_B->coeffs[j];

          if (fabs(coef_A_i * coef_B_j) < 1e-15) {
              continue;
          }

          evaluate_type1_integral_primitive_gaussians(
                  shell_A->origin, size_A, shell_A->cart_list,
  shell_A->alpha[i], shell_B->origin, size_B, shell_B->cart_list,
  shell_B->alpha[j], rpp_origin, potential, buf
          );

          libgrpp_daxpy(size_A * size_B, coef_A_i * coef_B_j, buf, matrix);
      }
  }*/
}

/**
 * Evaluation of type 1 RPP integrals (scalar-relativistic radially local RPP)
 * for the pair of shells constructed from primitive Gaussians.
 */
void evaluate_type1_integral_primitive_gaussians(
    double *A, int n_cart_A, int *cart_list_A, double alpha_A, double *B,
    int n_cart_B, int *cart_list_B, double alpha_B, double *C,
    libgrpp_potential_t *potential, double *matrix) {
  libgrpp_potential_t *potential_shrinked = libgrpp_shrink_potential(potential);

  libgrpp_evaluate_radially_local_potential_integral_primitive_gaussians(
      A, n_cart_A, cart_list_A, alpha_A, B, n_cart_B, cart_list_B, alpha_B, C,
      evaluate_pseudopotential, potential_shrinked, matrix);

  libgrpp_delete_potential(potential_shrinked);
}

static double evaluate_pseudopotential(double r, void *params) {
  libgrpp_potential_t *potential = (libgrpp_potential_t *)params;

  double u = libgrpp_potential_value(potential, r);

  return u;
}

/**
 * Evaluation of AO integrals for an arbitrary radially-local operator
 * for the pair of shells constructed from primitive Gaussians.
 */
void libgrpp_evaluate_radially_local_potential_integral_primitive_gaussians(
    double *A, int n_cart_A, int *cart_list_A, double alpha_A, double *B,
    int n_cart_B, int *cart_list_B, double alpha_B, double *C,
    double (*potential)(double r, void *params), void *potential_params,
    double *matrix) {
  assert(n_cart_A > 0);
  assert(n_cart_B > 0);

  memset(matrix, 0, n_cart_A * n_cart_B * sizeof(double));

  double CA_x = C[0] - A[0];
  double CA_y = C[1] - A[1];
  double CA_z = C[2] - A[2];
  double CB_x = C[0] - B[0];
  double CB_y = C[1] - B[1];
  double CB_z = C[2] - B[2];
  double CA_2 = CA_x * CA_x + CA_y * CA_y + CA_z * CA_z;
  double CB_2 = CB_x * CB_x + CB_y * CB_y + CB_z * CB_z;

  double kx = -2.0 * (alpha_A * CA_x + alpha_B * CB_x);
  double ky = -2.0 * (alpha_A * CA_y + alpha_B * CB_y);
  double kz = -2.0 * (alpha_A * CA_z + alpha_B * CB_z);
  double k = sqrt(kx * kx + ky * ky + kz * kz);
  double kvec[3];
  kvec[0] = kx;
  kvec[1] = ky;
  kvec[2] = kz;

  int L_A = cart_list_A[0] + cart_list_A[1] + cart_list_A[2];
  int L_B = cart_list_B[0] + cart_list_B[1] + cart_list_B[2];

  double N_A = libgrpp_gaussian_norm_factor(L_A, 0, 0, alpha_A);
  double N_B = libgrpp_gaussian_norm_factor(L_B, 0, 0, alpha_B);
  double D_ABC = 4 * M_PI * N_A * N_B;

  int lambda_max = L_A + L_B;
  int n_max = lambda_max;
  // create_real_spherical_harmonic_coeffs_tables(lambda_max);

  /*
   * pre-compute type 1 radial integrals
   */
  radial_type1_table_t *radial_table = libgrpp_tabulate_radial_type1_integrals(
      lambda_max, n_max, CA_2, CB_2, alpha_A, alpha_B, k, D_ABC, potential,
      potential_params);

  /*
   * main loop
   * over shell pairs
   */
  for (int icart = 0; icart < n_cart_A; icart++) {
    for (int jcart = 0; jcart < n_cart_B; jcart++) {

      double chi_AB = 0.0;

      int n_A = cart_list_A[3 * icart + 0];
      int l_A = cart_list_A[3 * icart + 1];
      int m_A = cart_list_A[3 * icart + 2];
      int n_B = cart_list_B[3 * jcart + 0];
      int l_B = cart_list_B[3 * jcart + 1];
      int m_B = cart_list_B[3 * jcart + 2];

      for (int a = 0; a <= n_A; a++) {
        double C_nA_a = libgrpp_binomial(n_A, a);
        double pow_CA_x = pow(CA_x, n_A - a);

        for (int b = 0; b <= l_A; b++) {
          double C_lA_b = libgrpp_binomial(l_A, b);
          double pow_CA_y = pow(CA_y, l_A - b);

          for (int c = 0; c <= m_A; c++) {
            double C_mA_c = libgrpp_binomial(m_A, c);
            double pow_CA_z = pow(CA_z, m_A - c);

            for (int d = 0; d <= n_B; d++) {
              double C_nB_d = libgrpp_binomial(n_B, d);
              double pow_CB_x = pow(CB_x, n_B - d);

              for (int e = 0; e <= l_B; e++) {
                double C_lB_e = libgrpp_binomial(l_B, e);
                double pow_CB_y = pow(CB_y, l_B - e);

                for (int f = 0; f <= m_B; f++) {
                  double C_mB_f = libgrpp_binomial(m_B, f);
                  double pow_CB_z = pow(CB_z, m_B - f);

                  double factor = C_nA_a * C_lA_b * C_mA_c * C_nB_d * C_lB_e *
                                  C_mB_f * pow_CA_x * pow_CA_y * pow_CA_z *
                                  pow_CB_x * pow_CB_y * pow_CB_z;

                  if (fabs(factor) < 1e-13) {
                    continue;
                  }

                  int N = a + b + c + d + e + f;
                  double sum_omega_Q = 0.0;
                  for (int lambda = 0; lambda <= lambda_max; lambda++) {

                    double Q = libgrpp_get_radial_type1_integral(radial_table,
                                                                 lambda, N);
                    if (fabs(Q) < 1e-16) {
                      continue;
                    }

                    double omega = libgrpp_angular_type1_integral(
                        lambda, a + d, b + e, c + f, kvec);

                    sum_omega_Q += omega * Q;
                  }

                  chi_AB += factor * sum_omega_Q;
                }
              }
            }
          }
        }
      }

      matrix[icart * n_cart_B + jcart] = chi_AB;
    }
  }

  libgrpp_delete_radial_type1_integrals(radial_table);
}
