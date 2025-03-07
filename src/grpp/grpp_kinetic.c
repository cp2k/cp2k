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

/**
 * Calculation of kinetic-energy integrals.
 *
 * For details, see:
 * T. Helgaker, P. Jorgensen, J. Olsen, Molecular Electronic-Structure Theory,
 * John Wiley & Sons Ltd, 2000.
 * Chapter 9.3.4, "Momentum and kinetic-energy integrals"
 */
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "grpp_kinetic.h"

#include "grpp_norm_gaussian.h"
#include "grpp_utils.h"
#include "libgrpp.h"

static void kinetic_energy_integrals_shell_pair_obara_saika(
    libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B, double alpha_A,
    double alpha_B, double *kinetic_matrix);

void libgrpp_kinetic_energy_integrals(libgrpp_shell_t *shell_A,
                                      libgrpp_shell_t *shell_B,
                                      double *kinetic_matrix) {
  int size_A = libgrpp_get_shell_size(shell_A);
  int size_B = libgrpp_get_shell_size(shell_B);

  double *buf = calloc(size_A * size_B, sizeof(double));

  memset(kinetic_matrix, 0, size_A * size_B * sizeof(double));

  // loop over primitives in contractions
  for (int i = 0; i < shell_A->num_primitives; i++) {
    for (int j = 0; j < shell_B->num_primitives; j++) {
      double alpha_i = shell_A->alpha[i];
      double alpha_j = shell_B->alpha[j];
      double coef_A_i = shell_A->coeffs[i];
      double coef_B_j = shell_B->coeffs[j];

      kinetic_energy_integrals_shell_pair_obara_saika(shell_A, shell_B, alpha_i,
                                                      alpha_j, buf);

      libgrpp_daxpy(size_A * size_B, coef_A_i * coef_B_j, buf, kinetic_matrix);
    }
  }

  free(buf);
}

static void kinetic_energy_integrals_shell_pair_obara_saika(
    libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B, double alpha_A,
    double alpha_B, double *kinetic_matrix) {
  int size_A = libgrpp_get_shell_size(shell_A);
  int size_B = libgrpp_get_shell_size(shell_B);
  int L_A = shell_A->L;
  int L_B = shell_B->L;
  double N_A = libgrpp_gaussian_norm_factor(L_A, 0, 0, alpha_A);
  double N_B = libgrpp_gaussian_norm_factor(L_B, 0, 0, alpha_B);

  double p = alpha_A + alpha_B;
  double mu = alpha_A * alpha_B / (alpha_A + alpha_B);
  double *A = shell_A->origin;
  double *B = shell_B->origin;

  // calculate S_ij
  double S[3][LIBGRPP_MAX_BASIS_L + 2][LIBGRPP_MAX_BASIS_L + 2];

  for (int coord = 0; coord < 3; coord++) {
    double P = (alpha_A * A[coord] + alpha_B * B[coord]) / p;

    double X_AB = A[coord] - B[coord];
    double X_PA = P - A[coord];
    double X_PB = P - B[coord];
    double pfac = 1.0 / (2.0 * p);

    for (int i = 0; i <= L_A + 2; i++) {
      for (int j = 0; j <= L_B + 2; j++) {
        double S_ij = 0.0;

        if (i + j == 0) {
          S[coord][0][0] = sqrt(M_PI / p) * exp(-mu * X_AB * X_AB);
          continue;
        }

        if (i == 0) { // upward by j
          S_ij += X_PB * S[coord][i][j - 1];
          if (j - 1 > 0) {
            S_ij += (j - 1) * pfac * S[coord][i][j - 2];
          }
        } else { // upward by i
          S_ij += X_PA * S[coord][i - 1][j];
          if (i - 1 > 0) {
            S_ij += (i - 1) * pfac * S[coord][i - 2][j];
          }
          if (j > 0) {
            S_ij += j * pfac * S[coord][i - 1][j - 1];
          }
        }

        S[coord][i][j] = S_ij;
      }
    }
  }

  // calculate D^2_ij

  double D2[3][LIBGRPP_MAX_BASIS_L][LIBGRPP_MAX_BASIS_L];

  for (int coord = 0; coord < 3; coord++) {
    for (int i = 0; i <= L_A; i++) {
      for (int j = 0; j <= L_B; j++) {

        double D2_ij = 0.0;
        D2_ij += 4.0 * alpha_A * alpha_A * S[coord][i + 2][j];
        D2_ij -= 2.0 * alpha_A * (2 * i + 1) * S[coord][i][j];
        if (i >= 2) {
          D2_ij += i * (i - 1) * S[coord][i - 2][j];
        }

        D2[coord][i][j] = D2_ij;
      }
    }
  }

  // loop over cartesian functions inside the shells
  for (int m = 0; m < size_A; m++) {
    for (int n = 0; n < size_B; n++) {
      int n_A = shell_A->cart_list[3 * m + 0];
      int l_A = shell_A->cart_list[3 * m + 1];
      int m_A = shell_A->cart_list[3 * m + 2];
      int n_B = shell_B->cart_list[3 * n + 0];
      int l_B = shell_B->cart_list[3 * n + 1];
      int m_B = shell_B->cart_list[3 * n + 2];

      kinetic_matrix[m * size_B + n] =
          -0.5 * N_A * N_B *
          (D2[0][n_A][n_B] * S[1][l_A][l_B] * S[2][m_A][m_B] +
           S[0][n_A][n_B] * D2[1][l_A][l_B] * S[2][m_A][m_B] +
           S[0][n_A][n_B] * S[1][l_A][l_B] * D2[2][m_A][m_B]);
    }
  }
}
