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
 * Integration of the non-local terms in the GRPP operator.
 * These integrals are reduced to the type 1 integrals and overlap integrals.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "grpp_factorial.h"
#include "grpp_lmatrix.h"
#include "grpp_overlap.h"
#include "grpp_spherical_harmonics.h"
#include "grpp_utils.h"
#include "libgrpp.h"

/*
 * pre-definitions of function used below
 */

void libgrpp_outercore_potential_integrals_part_1(
    libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B, double *C,
    libgrpp_potential_t *oc_potential, libgrpp_shell_t *oc_shell,
    double *arep_matrix, double *so_x_matrix, double *so_y_matrix,
    double *so_z_matrix);

void libgrpp_outercore_potential_integrals_part_2(
    libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B, double *C,
    libgrpp_potential_t *oc_potential_1, libgrpp_shell_t *oc_shell_1,
    libgrpp_potential_t *oc_potential_2, libgrpp_shell_t *oc_shell_2,
    double *arep_matrix, double *so_x_matrix, double *so_y_matrix,
    double *so_z_matrix);

static double calculate_delta_integral(libgrpp_potential_t *oc_pot_1,
                                       libgrpp_shell_t *oc_shell_1,
                                       libgrpp_potential_t *oc_pot_2,
                                       libgrpp_shell_t *oc_shell_2);

static void transform_to_sph_basis_ket(int dim_bra, int dim_ket_cart,
                                       int dim_ket_sph, double *A_in,
                                       double *A_out, double *S_lm_coef);

static void transform_to_sph_basis_bra(int dim_bra_cart, int dim_bra_sph,
                                       int dim_ket, double *A_in, double *A_out,
                                       double *S_lm_coef);

static double ang_norm_factor(int lx, int ly, int lz);

static double analytic_one_center_rpp_integral_contracted(
    libgrpp_shell_t *bra, libgrpp_shell_t *ket, libgrpp_potential_t *pot);

static double analytic_one_center_rpp_integral_primitive(int L, double alpha1,
                                                         double alpha2, int n,
                                                         double zeta);

static double radial_gto_norm_factor(int L, double alpha);

/**
 * Calculates non-local contributions to the scalar-relativistic ECP and
 * effective spin-orbit interaction matrices from the outercore (OC) potentials.
 */
void libgrpp_outercore_potential_integrals(
    libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B, double *rpp_origin,
    int num_oc_shells, libgrpp_potential_t **oc_potentials,
    libgrpp_shell_t **oc_shells, double *arep, double *esop_x, double *esop_y,
    double *esop_z) {
  assert(libgrpp_is_initialized());

  // clear output matrices
  int size_A = libgrpp_get_shell_size(shell_A);
  int size_B = libgrpp_get_shell_size(shell_B);
  memset(arep, 0, size_A * size_B * sizeof(double));
  memset(esop_x, 0, size_A * size_B * sizeof(double));
  memset(esop_y, 0, size_A * size_B * sizeof(double));
  memset(esop_z, 0, size_A * size_B * sizeof(double));

  // bind outercore shells of the GRPP to the RPP center
  for (int ioc = 0; ioc < num_oc_shells; ioc++) {
    oc_shells[ioc]->origin[0] = rpp_origin[0];
    oc_shells[ioc]->origin[1] = rpp_origin[1];
    oc_shells[ioc]->origin[2] = rpp_origin[2];
  }

  // 1. the U * |nlj><nlj| + |nlj><nlj| * U part
  for (int ioc = 0; ioc < num_oc_shells; ioc++) {
    libgrpp_potential_t *pot = oc_potentials[ioc];
    libgrpp_shell_t *nlj = oc_shells[ioc];

    libgrpp_outercore_potential_integrals_part_1(
        shell_A, shell_B, rpp_origin, pot, nlj, arep, esop_x, esop_y, esop_z);
  }

  // 2. the |nlj><nlj| U |n'lj><n'lj| part
  for (int ioc = 0; ioc < num_oc_shells; ioc++) {
    for (int joc = 0; joc < num_oc_shells; joc++) {

      libgrpp_potential_t *pot_1 = oc_potentials[ioc];
      libgrpp_potential_t *pot_2 = oc_potentials[joc];
      libgrpp_shell_t *nlj_1 = oc_shells[ioc];
      libgrpp_shell_t *nlj_2 = oc_shells[joc];

      libgrpp_outercore_potential_integrals_part_2(
          shell_A, shell_B, rpp_origin, pot_1, nlj_1, pot_2, nlj_2, arep,
          esop_x, esop_y, esop_z);
    }
  }
}

/**
 * integration of the non-local outercore potential:
 * the U*|nlj><nlj| + |nlj><nlj|*U part
 */
void libgrpp_outercore_potential_integrals_part_1(
    libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B, double *C,
    libgrpp_potential_t *oc_potential, libgrpp_shell_t *oc_shell,
    double *arep_matrix, double *so_x_matrix, double *so_y_matrix,
    double *so_z_matrix) {
  int L = oc_shell->L;
  double J = oc_potential->J / 2.0;
  int size_A = libgrpp_get_shell_size(shell_A);
  int size_B = libgrpp_get_shell_size(shell_B);
  int size_nlj = libgrpp_get_shell_size(oc_shell);

  /*
   * foolproof: bind outercore shells to the ECP center
   */
  oc_shell->origin[0] = C[0];
  oc_shell->origin[1] = C[1];
  oc_shell->origin[2] = C[2];

  /*
   * Transformation coefficients: Cartesian -> Real spherical
   * Rows: real spherical harmonics S_lm
   * Columns: cartesians x^r y^s z^t
   */
  double *S_lm_coef = (double *)calloc((2 * L + 1) * size_nlj, sizeof(double));
  for (int m = -L; m <= +L; m++) {
    for (int icart = 0; icart < size_nlj; icart++) {
      int r = oc_shell->cart_list[3 * icart + 0];
      int s = oc_shell->cart_list[3 * icart + 1];
      int t = oc_shell->cart_list[3 * icart + 2];
      double u = libgrpp_spherical_to_cartesian_coef(L, m, r, s) /
                 ang_norm_factor(r, s, t);
      S_lm_coef[size_nlj * (m + L) + icart] = u;
    }
  }

  /*
   * Overlap integrals of the A and B shells with the outercore shell |nlj>:
   * <A|nljm> and <nljm|B>
   *
   * Integrals are evaluated first for the cartesian components of |nlj>, and
   * then transformed to the basis of real spherical components |nljm>.
   * Resulting integrals are stored in the 'S_a_nljm' and 'S_nljm_b' arrays.
   */
  double *S_nljm_b_cart =
      alloc_zeros_1d(size_nlj * size_B); // in cart_list basis
  double *S_a_nljm_cart = alloc_zeros_1d(size_A * size_nlj);
  double *S_nljm_b = alloc_zeros_1d((2 * L + 1) * size_B); // in spherical basis
  double *S_a_nljm = alloc_zeros_1d(size_A * (2 * L + 1));
  libgrpp_overlap_integrals(oc_shell, shell_B, S_nljm_b_cart); // <nljm|B>
  libgrpp_overlap_integrals(shell_A, oc_shell, S_a_nljm_cart); // <A|nljm>
  transform_to_sph_basis_ket(size_A, size_nlj, 2 * L + 1, S_a_nljm_cart,
                             S_a_nljm, S_lm_coef);
  transform_to_sph_basis_bra(size_nlj, 2 * L + 1, size_B, S_nljm_b_cart,
                             S_nljm_b, S_lm_coef);
  free(S_nljm_b_cart);
  free(S_a_nljm_cart);

  /*
   * ECP type-2 (semilocal) integrals of the A and B shells with the outercore
   * shell |nlj>: <A|U(r)P_L|nljm> and <nljm|U(r)P_L|B>
   *
   * Integrals are evaluated first for the cartesian components of |nlj>, and
   * then transformed to the basis of real spherical components |nljm>.
   * Resulting integrals are stored in the 'U_a_nljm' and 'U_nljm_b' arrays.
   */
  double *U_nljm_b_cart =
      alloc_zeros_1d(size_nlj * size_B); // in cart_list basis
  double *U_a_nljm_cart = alloc_zeros_1d(size_A * size_nlj);
  double *U_nljm_b = alloc_zeros_1d((2 * L + 1) * size_B); // in spherical basis
  double *U_a_nljm = alloc_zeros_1d(size_A * (2 * L + 1));
  libgrpp_type1_integrals(oc_shell, shell_B, C, oc_potential,
                          U_nljm_b_cart); // <nljm|U(r)P_L|B>
  libgrpp_type1_integrals(shell_A, oc_shell, C, oc_potential,
                          U_a_nljm_cart); // <A|U(r)P_L|nljm>
  transform_to_sph_basis_ket(size_A, size_nlj, 2 * L + 1, U_a_nljm_cart,
                             U_a_nljm, S_lm_coef);
  transform_to_sph_basis_bra(size_nlj, 2 * L + 1, size_B, U_nljm_b_cart,
                             U_nljm_b, S_lm_coef);
  free(U_nljm_b_cart);
  free(U_a_nljm_cart);

  /*
   * Construct outercore AREP matrix elements
   * < a | P_nlj U(r) P_L + U(r) P_nlj P_L | b >
   */
  double arep_factor =
      (J < L) ? (L / (2.0 * L + 1)) : ((L + 1) / (2.0 * L + 1));
  double *buf = alloc_zeros_1d(size_A * size_B);
  libgrpp_multiply_matrices(size_A, size_B, 2 * L + 1, S_a_nljm, U_nljm_b, buf);
  libgrpp_multiply_matrices(size_A, size_B, 2 * L + 1, U_a_nljm, S_nljm_b, buf);
  libgrpp_daxpy(size_A * size_B, arep_factor, buf, arep_matrix);
  free(buf);

  /*
   * Construct outercore effective SO potential matrix elements
   * < a | P_nlj U(r) L P_L + U(r) P_nlj L P_L | b >
   */
  double **L_matrices = alloc_zeros_2d(3, (2 * L + 1) * (2 * L + 1));
  libgrpp_construct_angular_momentum_matrices_rsh(L, L_matrices[0],
                                                  L_matrices[1], L_matrices[2]);

  double **so_buf = alloc_zeros_2d(3, size_A * size_B);
  buf = alloc_zeros_1d((2 * L + 1) * int_max2(size_A, size_B));

  for (int icoord = 0; icoord < 3; icoord++) {
    // U(r) P_nlj
    memset(buf, 0, (2 * L + 1) * size_B * sizeof(double));
    libgrpp_multiply_matrices(2 * L + 1, size_B, 2 * L + 1, L_matrices[icoord],
                              S_nljm_b, buf);
    libgrpp_multiply_matrices(size_A, size_B, 2 * L + 1, U_a_nljm, buf,
                              so_buf[icoord]);

    // P_nlj U(r)
    memset(buf, 0, (2 * L + 1) * size_A * sizeof(double));
    libgrpp_multiply_matrices(size_A, 2 * L + 1, 2 * L + 1, S_a_nljm,
                              L_matrices[icoord], buf);
    libgrpp_multiply_matrices(size_A, size_B, 2 * L + 1, buf, U_nljm_b,
                              so_buf[icoord]);
  }

  free(buf);

  double esop_factor = (J < L) ? (-2.0 / (2 * L + 1)) : (+2.0 / (2 * L + 1));
  libgrpp_daxpy(size_A * size_B, esop_factor, so_buf[0], so_x_matrix);
  libgrpp_daxpy(size_A * size_B, esop_factor, so_buf[1], so_y_matrix);
  libgrpp_daxpy(size_A * size_B, esop_factor, so_buf[2], so_z_matrix);

  /*
   * Cleanup
   */
  free(S_lm_coef);
  free(S_a_nljm);
  free(S_nljm_b);
  free(U_a_nljm);
  free(U_nljm_b);
  free_2d(L_matrices, 3);
  free_2d(so_buf, 3);
}

/**
 * integration of the non-local outercore potential:
 * the |nlj><nlj| U |n'lj><n'lj| part
 */
void libgrpp_outercore_potential_integrals_part_2(
    libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B, double *C,
    libgrpp_potential_t *oc_potential_1, libgrpp_shell_t *oc_shell_1,
    libgrpp_potential_t *oc_potential_2, libgrpp_shell_t *oc_shell_2,
    double *arep_matrix, double *so_x_matrix, double *so_y_matrix,
    double *so_z_matrix) {
  int size_A = libgrpp_get_shell_size(shell_A);
  int size_B = libgrpp_get_shell_size(shell_B);

  if (oc_potential_1->L != oc_potential_2->L) {
    return;
  }
  if (oc_potential_1->J != oc_potential_2->J) {
    return;
  }

  /*
   * foolproof: bind outercore shells to the ECP center
   */
  oc_shell_1->origin[0] = C[0];
  oc_shell_1->origin[1] = C[1];
  oc_shell_1->origin[2] = C[2];
  oc_shell_2->origin[0] = C[0];
  oc_shell_2->origin[1] = C[1];
  oc_shell_2->origin[2] = C[2];

  /*
   * just to be consistent with the MOLGEP code by A. Titov, N. Mosyagin, A.
   * Petrov: off-diagonal elements <n'lj|U|n''lj> = 0
   */
  /*if (!ecps_are_equal(oc_potential_1, oc_potential_2)) {
      return;
  }*/

  int L = oc_potential_1->L;
  double J = oc_potential_1->J / 2.0;
  int size_nlj = libgrpp_get_shell_size(oc_shell_1);

  double delta = calculate_delta_integral(oc_potential_1, oc_shell_1,
                                          oc_potential_2, oc_shell_2);

  /*
   * Transformation coefficients: Cartesian -> Real spherical
   * Rows: real spherical harmonics S_lm
   * Columns: cartesians x^r y^s z^t
   */
  double *S_lm_coef = (double *)calloc((2 * L + 1) * size_nlj, sizeof(double));
  for (int m = -L; m <= +L; m++) {
    for (int icart = 0; icart < size_nlj; icart++) {
      int r = oc_shell_1->cart_list[3 * icart + 0];
      int s = oc_shell_1->cart_list[3 * icart + 1];
      int t = oc_shell_1->cart_list[3 * icart + 2];
      double u = libgrpp_spherical_to_cartesian_coef(L, m, r, s) /
                 ang_norm_factor(r, s, t);
      S_lm_coef[size_nlj * (m + L) + icart] = u;
    }
  }

  /*
   * Overlap integrals of the A and B shells with the outercore shells |nlj> and
   * |n'lj>: <A|nljm> and <n'ljm|B>
   *
   * Integrals are evaluated first for the cartesian components of |nlj>, and
   * then transformed to the basis of real spherical components |nljm>.
   * Resulting integrals are stored in the 'S_a_nljm1' and 'S_nljm2_b' arrays.
   */
  double *S_nljm2_b_cart =
      alloc_zeros_1d(size_nlj * size_B); // in cart_list basis
  double *S_a_nljm1_cart = alloc_zeros_1d(size_A * size_nlj);
  double *S_nljm2_b =
      alloc_zeros_1d((2 * L + 1) * size_B); // in spherical basis
  double *S_a_nljm1 = alloc_zeros_1d(size_A * (2 * L + 1));
  libgrpp_overlap_integrals(oc_shell_2, shell_B, S_nljm2_b_cart); // <nljm|B>
  libgrpp_overlap_integrals(shell_A, oc_shell_1, S_a_nljm1_cart); // <A|nljm>
  transform_to_sph_basis_ket(size_A, size_nlj, 2 * L + 1, S_a_nljm1_cart,
                             S_a_nljm1, S_lm_coef);
  transform_to_sph_basis_bra(size_nlj, 2 * L + 1, size_B, S_nljm2_b_cart,
                             S_nljm2_b, S_lm_coef);
  free(S_nljm2_b_cart);
  free(S_a_nljm1_cart);

  /*
   * Construct outercore AREP matrix elements
   * < a | P_nlj U(r) P_n'lj | b >
   */
  double arep_factor =
      (J < L) ? (L / (2.0 * L + 1)) : ((L + 1) / (2.0 * L + 1));
  double *buf = alloc_zeros_1d(size_A * size_B);
  libgrpp_multiply_matrices(size_A, size_B, 2 * L + 1, S_a_nljm1, S_nljm2_b,
                            buf);
  libgrpp_daxpy(size_A * size_B, (-1.0) * delta * arep_factor, buf,
                arep_matrix);
  free(buf);

  /*
   * Construct outercore effective SO potential matrix elements
   * < a | P_nlj U(r) L P_L + U(r) P_nlj L P_L | b >
   */
  double **L_matrices = alloc_zeros_2d(3, (2 * L + 1) * (2 * L + 1));
  libgrpp_construct_angular_momentum_matrices_rsh(L, L_matrices[0],
                                                  L_matrices[1], L_matrices[2]);

  double **so_buf = alloc_zeros_2d(3, size_A * size_B);
  buf = alloc_zeros_1d((2 * L + 1) * size_B);
  for (int icoord = 0; icoord < 3; icoord++) {
    memset(buf, 0, (2 * L + 1) * size_B * sizeof(double));
    libgrpp_multiply_matrices(2 * L + 1, size_B, 2 * L + 1, L_matrices[icoord],
                              S_nljm2_b, buf);
    libgrpp_multiply_matrices(size_A, size_B, 2 * L + 1, S_a_nljm1, buf,
                              so_buf[icoord]);
  }
  free(buf);

  double esop_factor = (J < L) ? (-2.0 / (2 * L + 1)) : (+2.0 / (2 * L + 1));
  libgrpp_daxpy(size_A * size_B, (-1.0) * delta * esop_factor, so_buf[0],
                so_x_matrix);
  libgrpp_daxpy(size_A * size_B, (-1.0) * delta * esop_factor, so_buf[1],
                so_y_matrix);
  libgrpp_daxpy(size_A * size_B, (-1.0) * delta * esop_factor, so_buf[2],
                so_z_matrix);

  free_2d(so_buf, 3);
  free_2d(L_matrices, 3);

  free(S_nljm2_b);
  free(S_a_nljm1);
  free(S_lm_coef);
}

/**
 * Calculation of radial "delta" integrals.
 * Is performed analytically.
 */
static double calculate_delta_integral(libgrpp_potential_t *oc_pot_1,
                                       libgrpp_shell_t *oc_shell_1,
                                       libgrpp_potential_t *oc_pot_2,
                                       libgrpp_shell_t *oc_shell_2) {
  // both shells must have equal L,J quantum numbers, otherwise
  // the < nlj | U | n'l'j' > integral is strictly zero
  if (oc_pot_1->L != oc_pot_2->L || oc_pot_1->J != oc_pot_2->J) {
    return 0.0;
  }

  double U1 = analytic_one_center_rpp_integral_contracted(oc_shell_1,
                                                          oc_shell_2, oc_pot_1);
  double U2 = analytic_one_center_rpp_integral_contracted(oc_shell_1,
                                                          oc_shell_2, oc_pot_2);

  return 0.5 * (U1 + U2);
}

/**
 * analytic formula for one-center RPP integral between two contracted gaussian
 * functions.
 */
static double analytic_one_center_rpp_integral_contracted(
    libgrpp_shell_t *bra, libgrpp_shell_t *ket, libgrpp_potential_t *pot) {
  double sum = 0.0;
  int L = pot->L;

  for (int i = 0; i < bra->num_primitives; i++) {
    double coef_i = bra->coeffs[i];
    double alpha_i = bra->alpha[i];
    double N_i = radial_gto_norm_factor(L, alpha_i);

    for (int j = 0; j < ket->num_primitives; j++) {
      double coef_j = ket->coeffs[j];
      double alpha_j = ket->alpha[j];
      double N_j = radial_gto_norm_factor(L, alpha_j);
      double factor = N_i * N_j * coef_i * coef_j;

      for (int k = 0; k < pot->num_primitives; k++) {
        double coef_k = pot->coeffs[k];
        double zeta = pot->alpha[k];
        int n_rpp = pot->powers[k];

        double u_ijk = analytic_one_center_rpp_integral_primitive(
            L, alpha_i, alpha_j, n_rpp, zeta);

        sum += factor * coef_k * u_ijk;
      }
    }
  }

  return sum;
}

/**
 * analytic formula for one-center RPP integral between two gaussian primitives.
 * normalization factors are omitted here.
 */
static double analytic_one_center_rpp_integral_primitive(int L, double alpha1,
                                                         double alpha2, int n,
                                                         double zeta) {
  double a = alpha1 + alpha2 + zeta;

  if (n % 2 == 0) { // even n
    int k = L + n / 2;
    return libgrpp_double_factorial(2 * k - 1) / (pow(2.0, k + 1) * pow(a, k)) *
           sqrt(M_PI / a);
  } else { // odd n
    int k = L + (n - 1) / 2;
    return libgrpp_factorial(k) / (2.0 * pow(a, k + 1));
  }
}

/**
 * calculate normalization factor for the radial Gaussian-type orbital:
 * G(r) = N * r^L * exp(-alpha * r^2)
 */
static double radial_gto_norm_factor(int L, double alpha) {
  // pre-tabulated factors for calculation of normalization constants
  // (for each value of L)
  static const double factors[] = {
      2.5264751109842587,    2.9173221708553032,    2.6093322745198853,
      1.9724697960897537,    1.3149798640598356,    7.9296269381073192e-1,
      4.3985656185609934e-1, 2.2714095183849672e-1, 1.1017954545099481e-1,
      5.0553842554329785e-2, 2.2063505731056757e-2, 9.2011179391124215e-3,
      3.6804471756449694e-3, 1.4166047783978804e-3, 5.2611380677564405e-4,
      1.8898565833279173e-4, 6.5796360823633550e-5, 2.2243229718298637e-5,
      7.3135288801774484e-6, 2.3422037547660024e-6};

  return factors[L] * pow(alpha, 0.75 + L / 2.0);
}

/**
 * Transforms matrix from the basis of unitary sphere polynomials
 * to the basis of real spherical harmonics S_lm
 * (separately for 'bra' and 'ket' vectors)
 */
static void transform_to_sph_basis_ket(int dim_bra, int dim_ket_cart,
                                       int dim_ket_sph, double *A_in,
                                       double *A_out, double *S_lm_coef) {
  for (int i = 0; i < dim_bra; i++) {
    for (int j = 0; j < dim_ket_sph; j++) {
      double s = 0.0;

      for (int icart = 0; icart < dim_ket_cart; icart++) {
        double u_nljm = S_lm_coef[dim_ket_cart * j + icart];
        s += u_nljm * A_in[i * dim_ket_cart + icart];
      }

      A_out[i * dim_ket_sph + j] = s;
    }
  }
}

static void transform_to_sph_basis_bra(int dim_bra_cart, int dim_bra_sph,
                                       int dim_ket, double *A_in, double *A_out,
                                       double *S_lm_coef) {
  for (int j = 0; j < dim_ket; j++) {
    for (int i = 0; i < dim_bra_sph; i++) {

      double s = 0.0;

      for (int icart = 0; icart < dim_bra_cart; icart++) {
        double u_nljm = S_lm_coef[dim_bra_cart * i + icart];
        s += u_nljm * A_in[icart * dim_ket + j];
      }

      A_out[i * dim_ket + j] = s;
    }
  }
}

static double ang_norm_factor(int lx, int ly, int lz) {
  int L = lx + ly + lz;

  return 1.0 / (2.0 * sqrt(M_PI)) *
         sqrt(libgrpp_double_factorial(2 * L + 1)
              /*(double_factorial(2 * lx - 1) * double_factorial(2 * ly - 1) *
                 double_factorial(2 * lz - 1))*/
         );
}
