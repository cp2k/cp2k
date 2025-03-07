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
#include <string.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "grpp_angular_integrals.h"
#include "grpp_binomial.h"
#include "grpp_radial_type2_integral.h"
#include "grpp_spherical_harmonics.h"
#include "grpp_utils.h"
#include "libgrpp.h"

#define LMAX (2 * LIBGRPP_MAX_BASIS_L + LIBGRPP_MAX_RPP_L)

static double type2_angular_sum(int L, int lambda_1, int a, int b, int c,
                                int lambda_2, int d, int e, int f,
                                double *rsh_values_kA, double *rsh_values_kB);

/**
 * Evaluation of type 2 RPP integrals (scalar-relativistic semilocal RPP with
 * L-projectors).
 */
void libgrpp_type2_integrals(libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B,
                             double *rpp_origin, libgrpp_potential_t *potential,
                             double *matrix) {
  assert(libgrpp_is_initialized());

  int size_A = libgrpp_get_shell_size(shell_A);
  int size_B = libgrpp_get_shell_size(shell_B);

  memset(matrix, 0, size_A * size_B * sizeof(double));

  int L = potential->L;
  int L_A =
      shell_A->cart_list[0] + shell_A->cart_list[1] + shell_A->cart_list[2];
  int L_B =
      shell_B->cart_list[0] + shell_B->cart_list[1] + shell_B->cart_list[2];

  double *A = shell_A->origin;
  double *B = shell_B->origin;
  double *C = rpp_origin;

  double CA_x = C[0] - A[0];
  double CA_y = C[1] - A[1];
  double CA_z = C[2] - A[2];
  double CB_x = C[0] - B[0];
  double CB_y = C[1] - B[1];
  double CB_z = C[2] - B[2];
  double CA_2 = CA_x * CA_x + CA_y * CA_y + CA_z * CA_z;
  double CB_2 = CB_x * CB_x + CB_y * CB_y + CB_z * CB_z;

  double alpha_A = shell_A->alpha[0];
  double alpha_B = shell_B->alpha[0];
  double kA_x = -2.0 * (alpha_A * CA_x);
  double kA_y = -2.0 * (alpha_A * CA_y);
  double kA_z = -2.0 * (alpha_A * CA_z);
  double kB_x = -2.0 * (alpha_B * CB_x);
  double kB_y = -2.0 * (alpha_B * CB_y);
  double kB_z = -2.0 * (alpha_B * CB_z);
  double kA_vec[3];
  kA_vec[0] = kA_x;
  kA_vec[1] = kA_y;
  kA_vec[2] = kA_z;
  double kB_vec[3];
  kB_vec[0] = kB_x;
  kB_vec[1] = kB_y;
  kB_vec[2] = kB_z;

  int lambda1_max = L + L_A;
  int lambda2_max = L + L_B;
  int N_max = L_A + L_B;

  /*
   * for further evaluation of angular integrals
   */
  int lmax = int_max3(lambda1_max, lambda2_max, L);
  // create_real_spherical_harmonic_coeffs_tables(lmax);

  /*
   * pre-compute type 2 radial integrals
   */
  radial_type2_table_t *radial_table = libgrpp_tabulate_radial_type2_integrals(
      lambda1_max, lambda2_max, N_max, CA_2, CB_2, potential, shell_A, shell_B);

  /*
   * pre-calculate values of real spherical harmonics for different L
   */
  double rsh_values_kA[3 * LMAX][6 * LMAX];
  double rsh_values_kB[3 * LMAX][6 * LMAX];

  for (int lambda = 0; lambda <= lmax; lambda++) {
    libgrpp_evaluate_real_spherical_harmonics_array(lambda, kA_vec,
                                                    rsh_values_kA[lambda]);
    libgrpp_evaluate_real_spherical_harmonics_array(lambda, kB_vec,
                                                    rsh_values_kB[lambda]);
  }

  /*
   * main loop
   * over shell pairs
   */
  for (int icart = 0; icart < size_A; icart++) {
    for (int jcart = 0; jcart < size_B; jcart++) {

      double gamma_AB = 0.0;

      int n_A = shell_A->cart_list[3 * icart + 0];
      int l_A = shell_A->cart_list[3 * icart + 1];
      int m_A = shell_A->cart_list[3 * icart + 2];
      int n_B = shell_B->cart_list[3 * jcart + 0];
      int l_B = shell_B->cart_list[3 * jcart + 1];
      int m_B = shell_B->cart_list[3 * jcart + 2];

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

                  int lambda1_lower = int_max2(L - a - b - c, 0);
                  int lambda2_lower = int_max2(L - d - e - f, 0);
                  int lambda1_upper = L + a + b + c;
                  int lambda2_upper = L + d + e + f;

                  for (int lambda_1 = lambda1_lower; lambda_1 <= lambda1_upper;
                       lambda_1++) {
                    if ((L + a + b + c - lambda_1) % 2 != 0) {
                      continue;
                    }

                    for (int lambda_2 = lambda2_lower;
                         lambda_2 <= lambda2_upper; lambda_2++) {
                      if ((L + d + e + f - lambda_2) % 2 != 0) {
                        continue;
                      }

                      double QN = libgrpp_get_radial_type2_integral(
                          radial_table, lambda_1, lambda_2, N);
                      if (fabs(QN) < 1e-16) {
                        continue;
                      }

                      double sum_angular = type2_angular_sum(
                          L, lambda_1, a, b, c, lambda_2, d, e, f,
                          rsh_values_kA[lambda_1], rsh_values_kB[lambda_2]);

                      sum_omega_Q += QN * sum_angular;
                    } // loop over lambda_2
                  } // loop over lambda_1

                  gamma_AB += factor * sum_omega_Q;
                }
              }
            }
          }
        }
      }

      gamma_AB *= 16 * M_PI * M_PI;

      matrix[icart * size_B + jcart] = gamma_AB;
    }
  }

  libgrpp_delete_radial_type2_integrals(radial_table);
}

/*
 * Sum of products of type 2 angular integrals
 * (McMurchie, Davidson, 1981, formulas (23) and (24))
 */
static double type2_angular_sum(int L, int lambda_1, int a, int b, int c,
                                int lambda_2, int d, int e, int f,
                                double *rsh_values_kA, double *rsh_values_kB) {
  double sum_angular = 0.0;

  /*
   * contract tensors with angular integrals
   */
  for (int m = -L; m <= L; m++) {
    double omega_1 =
        libgrpp_angular_type2_integral(lambda_1, L, m, a, b, c, rsh_values_kA);
    if (fabs(omega_1) < 1e-16) {
      continue;
    }

    double omega_2 =
        libgrpp_angular_type2_integral(lambda_2, L, m, d, e, f, rsh_values_kB);

    sum_angular += omega_1 * omega_2;
  }

  return sum_angular;
}
