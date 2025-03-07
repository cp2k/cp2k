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
 * Implementation of the McMurchie-Davidson (MMD) recurrence relations
 * for radially local RPP integrals.
 *
 * Operators to be integrated are:
 * n = 2: e^{-ar^2}
 * n = 1: e^{-ar^2}/r^1
 * n = 0: e^{-ar^2}/r^2
 *
 * General description of the MMD scheme is given in:
 * (1) T. Helgaker, P. Jorgensen, J. Olsen, Molecular Electronic-Structure
 * Theory, John Wiley & Sons Ltd, 2000. (2) L. E. McMurchie, E. R. Davidson,
 *     One- and two-electron integrals over cartesian gaussian functions.
 *     J. Comput. Phys. 26(2), 218 (1978).
 *     doi: 10.1016/0021-9991(78)90092-X
 *
 * More on the evaluation of integrals over the 1/r^2 operator:
 * (1) J. 0. Jensen, A. H. Cameri, C. P. Vlahacos, D. Zeroka, H. F. Hameka, C.
 * N. Merrow, Evaluation of one-electron integrals for arbitrary operators V(r)
 * over cartesian Gaussians: Application to inverse-square distance and Yukawa
 * operators. J. Comput. Chem. 14(8), 986 (1993). doi: 10.1002/jcc.540140814 (2)
 * B. Gao, A. J. Thorvaldsen, K. Ruud, GEN1INT: A unified procedure for the
 * evaluation of one-electron integrals over Gaussian basis functions and their
 * geometric derivatives. Int. J. Quantum Chem. 111(4), 858 (2011).
 *     doi: 10.1002/qua.22886
 */
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "grpp_type1_mcmurchie_davidson.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "grpp_norm_gaussian.h"
#include "grpp_specfunc.h"
#include "grpp_utils.h"
#include "libgrpp.h"

#define LMAX LIBGRPP_MAX_BASIS_L

/*
 * temporary data for the McMurchie-Davidson algorithm
 */
struct mmd_data {
  double A[3];
  double B[3];
  double C[3];
  double Q[3];
  double q;
  double K_abc[3];
  double R_QC_2;
  double boys[100];
  double E[3][LMAX][LMAX][2 * LMAX];
  double R[2 * LMAX][2 * LMAX][2 * LMAX][2 * LMAX];
};

/*
 * functions used below in the file
 */

static void evaluate_rpp_type1_mmd_n2_primitive_shell_pair(
    libgrpp_shell_t *shell_A, double alpha_A, libgrpp_shell_t *shell_B,
    double alpha_B, double *rpp_origin, double rpp_alpha, double *rpp_matrix);

void libgrpp_evaluate_rpp_type1_mmd_n1_primitive_shell_pair(
    libgrpp_shell_t *shell_A, double alpha_A, libgrpp_shell_t *shell_B,
    double alpha_B, double *rpp_origin, double rpp_alpha, double *rpp_matrix);

static void evaluate_rpp_type1_mmd_n0_primitive_shell_pair(
    libgrpp_shell_t *shell_A, double alpha_A, libgrpp_shell_t *shell_B,
    double alpha_B, double *rpp_origin, double rpp_alpha, double *rpp_matrix);

static void setup_E_array(struct mmd_data *data, int L_A, int L_B);

static void setup_R_array(struct mmd_data *data, int L_A, int L_B);

static void setup_G_array(struct mmd_data *data, int L_A, int L_B);

/**
 * General interface for the McMurchie-Davidson algorithm for integrals
 * over the radially local RPP operator.
 */
void libgrpp_type1_integrals_mcmurchie_davidson_1978(
    libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B, double *origin_C,
    double alpha_C, int ecp_power, double *rpp_matrix) {
  assert((ecp_power == 0) || (ecp_power == 1) || (ecp_power == 2));

  int size_A = libgrpp_get_shell_size(shell_A);
  int size_B = libgrpp_get_shell_size(shell_B);

  double *buf = calloc(size_A * size_B, sizeof(double));

  memset(rpp_matrix, 0, size_A * size_B * sizeof(double));

  /*
   * loop over primitives in contractions
   */
  for (int i = 0; i < shell_A->num_primitives; i++) {
    double coef_A_i = shell_A->coeffs[i];
    if (fabs(coef_A_i) < LIBGRPP_ZERO_THRESH) {
      continue;
    }

    for (int j = 0; j < shell_B->num_primitives; j++) {
      double coef_B_j = shell_B->coeffs[j];
      if (fabs(coef_B_j) < LIBGRPP_ZERO_THRESH) {
        continue;
      }

      if (ecp_power == 2) {
        evaluate_rpp_type1_mmd_n2_primitive_shell_pair(
            shell_A, shell_A->alpha[i], shell_B, shell_B->alpha[j], origin_C,
            alpha_C, buf);
      } else if (ecp_power == 1) {
        libgrpp_evaluate_rpp_type1_mmd_n1_primitive_shell_pair(
            shell_A, shell_A->alpha[i], shell_B, shell_B->alpha[j], origin_C,
            alpha_C, buf);
      } else if (ecp_power == 0) {
        evaluate_rpp_type1_mmd_n0_primitive_shell_pair(
            shell_A, shell_A->alpha[i], shell_B, shell_B->alpha[j], origin_C,
            alpha_C, buf);
      }

      libgrpp_daxpy(size_A * size_B, coef_A_i * coef_B_j, buf, rpp_matrix);
    }
  }

  free(buf);
}

/**
 * Integrals of the operator: e^{-ar^2}/r^2
 */
static void evaluate_rpp_type1_mmd_n0_primitive_shell_pair(
    libgrpp_shell_t *shell_A, double alpha_A, libgrpp_shell_t *shell_B,
    double alpha_B, double *rpp_origin, double rpp_alpha, double *rpp_matrix) {
  double a = alpha_A;
  double b = alpha_B;
  double c = rpp_alpha;

  double *A = shell_A->origin;
  double *B = shell_B->origin;
  double *C = rpp_origin;

  double q = a + b + c;
  double mu_AB = a * b / q;
  double mu_AC = a * c / q;
  double mu_BC = b * c / q;

  struct mmd_data data;

  for (int i = 0; i < 3; i++) {
    data.A[i] = A[i];
    data.B[i] = B[i];
    data.C[i] = C[i];
    data.Q[i] = (a * A[i] + b * B[i] + c * C[i]) / q;

    double X_AB = A[i] - B[i];
    double X_AC = A[i] - C[i];
    double X_BC = B[i] - C[i];

    double K_ab = exp(-mu_AB * X_AB * X_AB);
    double K_ac = exp(-mu_AC * X_AC * X_AC);
    double K_bc = exp(-mu_BC * X_BC * X_BC);

    data.K_abc[i] = K_ab * K_ac * K_bc;
  }

  data.q = q;
  data.R_QC_2 = distance_squared(data.Q, data.C);
  int L_A = shell_A->L;
  int L_B = shell_B->L;
  int Nmax = L_A + L_B;
  libgrpp_gfun_values(q * data.R_QC_2, Nmax, data.boys);

  /*
   * setup E array
   */
  setup_E_array(&data, L_A, L_B);

  /*
   * setup R array
   */
  setup_G_array(&data, L_A, L_B);

  /*
   * loop over cartesian functions inside the shells
   */
  int size_A = libgrpp_get_shell_size(shell_A);
  int size_B = libgrpp_get_shell_size(shell_B);
  double N_A = libgrpp_gaussian_norm_factor(L_A, 0, 0, alpha_A);
  double N_B = libgrpp_gaussian_norm_factor(L_B, 0, 0, alpha_B);

  for (int m = 0; m < size_A; m++) {
    for (int n = 0; n < size_B; n++) {
      int n_A = shell_A->cart_list[3 * m + 0];
      int l_A = shell_A->cart_list[3 * m + 1];
      int m_A = shell_A->cart_list[3 * m + 2];
      int n_B = shell_B->cart_list[3 * n + 0];
      int l_B = shell_B->cart_list[3 * n + 1];
      int m_B = shell_B->cart_list[3 * n + 2];

      double s = 0.0;

      for (int t = 0; t <= n_A + n_B; t++) {
        double Ex = data.E[0][n_A][n_B][t];

        for (int u = 0; u <= l_A + l_B; u++) {
          double Ey = data.E[1][l_A][l_B][u];

          for (int v = 0; v <= m_A + m_B; v++) {
            double Ez = data.E[2][m_A][m_B][v];

            double R_tuv = data.R[t][u][v][0];

            s += Ex * Ey * Ez * R_tuv;
          }
        }
      }

      s *= N_A * N_B * 2.0 * pow(M_PI, 1.5) / sqrt(q);
      rpp_matrix[m * size_B + n] = s;
    }
  }
}

/**
 * Integrals of the operator: e^{-ar^2}/r
 */
void libgrpp_evaluate_rpp_type1_mmd_n1_primitive_shell_pair(
    libgrpp_shell_t *shell_A, double alpha_A, libgrpp_shell_t *shell_B,
    double alpha_B, double *rpp_origin, double rpp_alpha, double *rpp_matrix) {
  double a = alpha_A;
  double b = alpha_B;
  double c = rpp_alpha;

  double *A = shell_A->origin;
  double *B = shell_B->origin;
  double *C = rpp_origin;

  double q = a + b + c;
  double mu_AB = a * b / q;
  double mu_AC = a * c / q;
  double mu_BC = b * c / q;

  struct mmd_data data;

  for (int i = 0; i < 3; i++) {
    data.A[i] = A[i];
    data.B[i] = B[i];
    data.C[i] = C[i];
    data.Q[i] = (a * A[i] + b * B[i] + c * C[i]) / q;

    double X_AB = A[i] - B[i];
    double X_AC = A[i] - C[i];
    double X_BC = B[i] - C[i];

    double K_ab = exp(-mu_AB * X_AB * X_AB);
    double K_ac = exp(-mu_AC * X_AC * X_AC);
    double K_bc = exp(-mu_BC * X_BC * X_BC);

    data.K_abc[i] = K_ab * K_ac * K_bc;
  }

  data.q = q;
  data.R_QC_2 = distance_squared(data.Q, data.C);
  int L_A = shell_A->L;
  int L_B = shell_B->L;
  int Nmax = L_A + L_B;

  libgrpp_boys_values(q * data.R_QC_2, Nmax, data.boys);

  /*
   * setup E array
   */
  setup_E_array(&data, L_A, L_B);

  /*
   * setup R array
   */
  setup_R_array(&data, L_A, L_B);

  /*
   * loop over cartesian functions inside the shells
   */
  int size_A = libgrpp_get_shell_size(shell_A);
  int size_B = libgrpp_get_shell_size(shell_B);
  double N_A = libgrpp_gaussian_norm_factor(L_A, 0, 0, alpha_A);
  double N_B = libgrpp_gaussian_norm_factor(L_B, 0, 0, alpha_B);

  for (int m = 0; m < size_A; m++) {
    for (int n = 0; n < size_B; n++) {
      int n_A = shell_A->cart_list[3 * m + 0];
      int l_A = shell_A->cart_list[3 * m + 1];
      int m_A = shell_A->cart_list[3 * m + 2];
      int n_B = shell_B->cart_list[3 * n + 0];
      int l_B = shell_B->cart_list[3 * n + 1];
      int m_B = shell_B->cart_list[3 * n + 2];

      double s = 0.0;

      for (int t = 0; t <= n_A + n_B; t++) {
        double Ex = data.E[0][n_A][n_B][t];

        for (int u = 0; u <= l_A + l_B; u++) {
          double Ey = data.E[1][l_A][l_B][u];

          for (int v = 0; v <= m_A + m_B; v++) {
            double Ez = data.E[2][m_A][m_B][v];

            double R_tuv = data.R[t][u][v][0];

            s += Ex * Ey * Ez * R_tuv;
          }
        }
      }

      s *= N_A * N_B * 2.0 * M_PI / q;
      rpp_matrix[m * size_B + n] = s;
    }
  }
}

/**
 * Integrals of the operator: e^{-ar^2}
 */
static void evaluate_rpp_type1_mmd_n2_primitive_shell_pair(
    libgrpp_shell_t *shell_A, double alpha_A, libgrpp_shell_t *shell_B,
    double alpha_B, double *rpp_origin, double rpp_alpha, double *rpp_matrix) {
  double a = alpha_A;
  double b = alpha_B;
  double c = rpp_alpha;

  double *A = shell_A->origin;
  double *B = shell_B->origin;
  double *C = rpp_origin;

  double q = a + b + c;
  double mu_AB = a * b / q;
  double mu_AC = a * c / q;
  double mu_BC = b * c / q;

  struct mmd_data data;

  for (int i = 0; i < 3; i++) {
    data.A[i] = A[i];
    data.B[i] = B[i];
    data.C[i] = C[i];
    data.Q[i] = (a * A[i] + b * B[i] + c * C[i]) / q;

    double X_AB = A[i] - B[i];
    double X_AC = A[i] - C[i];
    double X_BC = B[i] - C[i];

    double K_ab = exp(-mu_AB * X_AB * X_AB);
    double K_ac = exp(-mu_AC * X_AC * X_AC);
    double K_bc = exp(-mu_BC * X_BC * X_BC);

    data.K_abc[i] = K_ab * K_ac * K_bc;
  }

  data.q = q;
  data.R_QC_2 = distance_squared(data.Q, data.C);

  int L_A = shell_A->L;
  int L_B = shell_B->L;

  /*
   * setup E array
   */
  setup_E_array(&data, L_A, L_B);

  /*
   * loop over cartesian functions inside the shells
   */
  int size_A = libgrpp_get_shell_size(shell_A);
  int size_B = libgrpp_get_shell_size(shell_B);
  double N_A = libgrpp_gaussian_norm_factor(L_A, 0, 0, alpha_A);
  double N_B = libgrpp_gaussian_norm_factor(L_B, 0, 0, alpha_B);

  for (int m = 0; m < size_A; m++) {
    for (int n = 0; n < size_B; n++) {
      int n_A = shell_A->cart_list[3 * m + 0];
      int l_A = shell_A->cart_list[3 * m + 1];
      int m_A = shell_A->cart_list[3 * m + 2];
      int n_B = shell_B->cart_list[3 * n + 0];
      int l_B = shell_B->cart_list[3 * n + 1];
      int m_B = shell_B->cart_list[3 * n + 2];

      double E_ij_0 = data.E[0][n_A][n_B][0];
      double E_kl_0 = data.E[1][l_A][l_B][0];
      double E_mn_0 = data.E[2][m_A][m_B][0];

      double s = N_A * N_B * E_ij_0 * E_kl_0 * E_mn_0 * pow(M_PI / q, 1.5);

      rpp_matrix[m * size_B + n] = s;
    }
  }
}

/**
 * Calculation of the R_{tuv}^N auxiliary integrals for the 1/r operator
 */
static void setup_R_array(struct mmd_data *data, int L_A, int L_B) {
  double q = data->q;
  double X_QC = data->Q[0] - data->C[0];
  double Y_QC = data->Q[1] - data->C[1];
  double Z_QC = data->Q[2] - data->C[2];

  int nmax = L_A + L_B;
  int L_SUM = L_A + L_B;

  for (int t = 0; t <= L_SUM; t++) {
    for (int u = 0; u <= L_SUM; u++) {
      for (int v = 0; v <= L_SUM; v++) {

        if (t + u + v > L_SUM) {
          continue;
        }

        for (int n = 0; n <= nmax - (t + u + v); n++) {

          double val = 0.0;

          if (t + u + v == 0) {
            val = pow(-2.0 * q, n) * data->boys[n];
          } else if (t + u == 0) {
            if (v > 1) {
              val += (v - 1) * data->R[t][u][v - 2][n + 1];
            }
            val += Z_QC * data->R[t][u][v - 1][n + 1];
          } else if (t == 0) {
            if (u > 1) {
              val += (u - 1) * data->R[t][u - 2][v][n + 1];
            }
            val += Y_QC * data->R[t][u - 1][v][n + 1];
          } else {
            if (t > 1) {
              val += (t - 1) * data->R[t - 2][u][v][n + 1];
            }
            val += X_QC * data->R[t - 1][u][v][n + 1];
          }

          data->R[t][u][v][n] = val;
        }
      }
    }
  }
}

/**
 * Calculation of the R_{tuv}^N auxiliary integrals for the 1/r^2 operator
 */
static void setup_G_array(struct mmd_data *data, int L_A, int L_B) {
  double q = data->q;
  double X_QC = data->Q[0] - data->C[0];
  double Y_QC = data->Q[1] - data->C[1];
  double Z_QC = data->Q[2] - data->C[2];

  int nmax = L_A + L_B;
  int L_SUM = L_A + L_B;

  for (int t = 0; t <= L_SUM; t++) {
    for (int u = 0; u <= L_SUM; u++) {
      for (int v = 0; v <= L_SUM; v++) {

        if (t + u + v > L_SUM) {
          continue;
        }

        for (int n = 0; n <= nmax - (t + u + v); n++) {

          double val = 0.0;

          if (t + u + v == 0) {
            val = pow(2.0 * q, n) * data->boys[n];
          } else if (t + u == 0) {
            if (v > 1) {
              val += (v - 1) * (data->R[t][u][v - 2][n + 1] -
                                2.0 * q * data->R[t][u][v - 2][n]);
            }
            val += Z_QC * (data->R[t][u][v - 1][n + 1] -
                           2.0 * q * data->R[t][u][v - 1][n]);
          } else if (t == 0) {
            if (u > 1) {
              val += (u - 1) * (data->R[t][u - 2][v][n + 1] -
                                2.0 * q * data->R[t][u - 2][v][n]);
            }
            val += Y_QC * (data->R[t][u - 1][v][n + 1] -
                           2.0 * q * data->R[t][u - 1][v][n]);
          } else {
            if (t > 1) {
              val += (t - 1) * (data->R[t - 2][u][v][n + 1] -
                                2.0 * q * data->R[t - 2][u][v][n]);
            }
            val += X_QC * (data->R[t - 1][u][v][n + 1] -
                           2.0 * q * data->R[t - 1][u][v][n]);
          }

          data->R[t][u][v][n] = val;
        }
      }
    }
  }
}

/**
 * Calculates E^{ij}_t coefficients in the MMD scheme.
 */
static void setup_E_array(struct mmd_data *data, int L_A, int L_B) {
  double q = data->q;

  for (int coord = 0; coord < 3; coord++) {
    for (int i = 0; i <= L_A; i++) {
      for (int j = 0; j <= L_B; j++) {

        for (int t = 0; t <= i + j; t++) {

          if (t == 0) {
            if (i == 0 && j == 0) {
              data->E[coord][0][0][0] = data->K_abc[coord];
              continue;
            } else if (i == 1 && j == 0) {
              double X_QA = data->Q[coord] - data->A[coord];
              data->E[coord][1][0][0] = X_QA * data->E[coord][0][0][0];
              continue;
            } else if (i == 0 && j == 1) {
              double X_QB = data->Q[coord] - data->B[coord];
              data->E[coord][0][1][0] = X_QB * data->E[coord][0][0][0];
              continue;
            } else if (i == 0) {
              double X_QB = data->Q[coord] - data->B[coord];
              data->E[coord][i][j][0] = X_QB * data->E[coord][i][j - 1][0] +
                                        data->E[coord][i][j - 1][1];
              continue;
            } else {
              double X_QA = data->Q[coord] - data->A[coord];
              data->E[coord][i][j][0] = X_QA * data->E[coord][i - 1][j][0] +
                                        data->E[coord][i - 1][j][1];
              continue;
            }
          } else {
            double E_ijt = 0.0;
            double factor = 1.0 / (2.0 * q * t);

            if (i > 0) {
              E_ijt += factor * i * data->E[coord][i - 1][j][t - 1];
            }
            if (j > 0) {
              E_ijt += factor * j * data->E[coord][i][j - 1][t - 1];
            }

            data->E[coord][i][j][t] = E_ijt;
          }
        }
      }
    }
  }
}
