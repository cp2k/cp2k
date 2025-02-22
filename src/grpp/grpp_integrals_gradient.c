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

#include "libgrpp.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "grpp_diff_gaussian.h"
#include "grpp_utils.h"

void grpp_gradient_diff_bra_contribution(
    libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B,
    libgrpp_grpp_t *grpp_operator, double *grpp_origin, double **grad_arep,
    double **grad_so_x, double **grad_so_y, double **grad_so_z, double factor);

void grpp_gradient_diff_bra_grpp_integrals(
    libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B,
    libgrpp_grpp_t *grpp_operator, double *grpp_origin,
    double **arep_matrix_down, double **so_x_matrix_down,
    double **so_y_matrix_down, double **so_z_matrix_down,
    double **arep_matrix_up, double **so_x_matrix_up, double **so_y_matrix_up,
    double **so_z_matrix_up, int *cart_size_down, int *cart_size_up);

void grpp_gradient_diff_ket_contribution(
    libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B,
    libgrpp_grpp_t *grpp_operator, double *grpp_origin, double **grad_arep,
    double **grad_so_x, double **grad_so_y, double **grad_so_z, double factor);

void grpp_gradient_diff_ket_grpp_integrals(
    libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B,
    libgrpp_grpp_t *grpp_operator, double *grpp_origin,
    double **arep_matrix_down, double **so_x_matrix_down,
    double **so_y_matrix_down, double **so_z_matrix_down,
    double **arep_matrix_up, double **so_x_matrix_up, double **so_y_matrix_up,
    double **so_z_matrix_up, int *cart_size_down, int *cart_size_up);

void grpp_gradient_contribution(libgrpp_shell_t *shell_A,
                                libgrpp_shell_t *shell_B,
                                libgrpp_grpp_t *grpp_operator,
                                double *grpp_origin, double **grad_arep,
                                double **grad_so_x, double **grad_so_y,
                                double **grad_so_z, int diff_bra,
                                double factor);

void grpp_gradient_diff_gaussian(
    libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B,
    libgrpp_grpp_t *grpp_operator, double *grpp_origin,
    double **arep_matrix_down, double **so_x_matrix_down,
    double **so_y_matrix_down, double **so_z_matrix_down,
    double **arep_matrix_up, double **so_x_matrix_up, double **so_y_matrix_up,
    double **so_z_matrix_up, int *cart_size_down, int *cart_size_up,
    int diff_bra);

extern int libgrpp_nlm_to_linear(int *nlm);

double **libgrpp_alloc_gradients(libgrpp_shell_t *bra, libgrpp_shell_t *ket);

void libgrpp_dealloc_gradients(double **grad);

/**
 * Analytic calculation of gradients of LOCAL potential integrals for a given
 * shell pair with respect to the point 'point_3d'.
 */
void libgrpp_type1_integrals_gradient(libgrpp_shell_t *shell_A,
                                      libgrpp_shell_t *shell_B,
                                      double *grpp_origin,
                                      libgrpp_potential_t *potential,
                                      double *point_3d, double **grad_arep) {
  libgrpp_grpp_t *grpp_operator = libgrpp_new_grpp();
  libgrpp_grpp_set_local_potential(grpp_operator, potential);

  /*
   * these arrays are not actually used.
   * they are needed only in order to use the
   * libgrpp_full_grpp_integrals_gradient() routine.
   */
  double **stub_grad_so_x = libgrpp_alloc_gradients(shell_A, shell_B);
  double **stub_grad_so_y = libgrpp_alloc_gradients(shell_A, shell_B);
  double **stub_grad_so_z = libgrpp_alloc_gradients(shell_A, shell_B);

  libgrpp_full_grpp_integrals_gradient(
      shell_A, shell_B, grpp_operator, grpp_origin, point_3d, grad_arep,
      stub_grad_so_x, stub_grad_so_y, stub_grad_so_z);

  libgrpp_dealloc_gradients(stub_grad_so_x);
  libgrpp_dealloc_gradients(stub_grad_so_y);
  libgrpp_dealloc_gradients(stub_grad_so_z);

  grpp_operator->U_L = NULL;
  libgrpp_delete_grpp(grpp_operator);
}

/**
 * Analytic calculation of gradients of SEMI-LOCAL potential integrals for a
 * given shell pair with respect to the point 'point_3d'.
 */
void libgrpp_type2_integrals_gradient(libgrpp_shell_t *shell_A,
                                      libgrpp_shell_t *shell_B,
                                      double *grpp_origin,
                                      libgrpp_potential_t *potential,
                                      double *point_3d, double **grad_arep) {
  libgrpp_grpp_t *grpp_operator = libgrpp_new_grpp();
  libgrpp_grpp_add_averaged_potential(grpp_operator, potential);

  /*
   * these arrays are not actually used.
   * they are needed only in order to use the
   * libgrpp_full_grpp_integrals_gradient() routine.
   */
  double **stub_grad_so_x = libgrpp_alloc_gradients(shell_A, shell_B);
  double **stub_grad_so_y = libgrpp_alloc_gradients(shell_A, shell_B);
  double **stub_grad_so_z = libgrpp_alloc_gradients(shell_A, shell_B);

  libgrpp_full_grpp_integrals_gradient(
      shell_A, shell_B, grpp_operator, grpp_origin, point_3d, grad_arep,
      stub_grad_so_x, stub_grad_so_y, stub_grad_so_z);

  libgrpp_dealloc_gradients(stub_grad_so_x);
  libgrpp_dealloc_gradients(stub_grad_so_y);
  libgrpp_dealloc_gradients(stub_grad_so_z);

  grpp_operator->n_arep = 0;
  libgrpp_delete_grpp(grpp_operator);
}

/**
 * Analytic calculation of gradients of integrals over the effective spin-orbit
 * operator (potential) for a given shell pair (with respect to the point
 * 'point_3d').
 */
void libgrpp_spin_orbit_integrals_gradient(
    libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B, double *grpp_origin,
    libgrpp_potential_t *potential, double *point_3d, double **grad_so_x,
    double **grad_so_y, double **grad_so_z) {
  libgrpp_grpp_t *grpp_operator = libgrpp_new_grpp();
  libgrpp_grpp_add_spin_orbit_potential(grpp_operator, potential);

  /*
   * this array is not actually used and is needed only in order
   * to use the libgrpp_full_grpp_integrals_gradient() routine.
   */
  double **stub_grad_arep = libgrpp_alloc_gradients(shell_A, shell_B);

  libgrpp_full_grpp_integrals_gradient(shell_A, shell_B, grpp_operator,
                                       grpp_origin, point_3d, stub_grad_arep,
                                       grad_so_x, grad_so_y, grad_so_z);

  /*
   * inside the libgrpp_full_grpp_integrals_gradient() function
   * the SO potential was scaled by 2/(2L+1). Thus the result has to be
   * re-scaled by (2L+1)/2 to get rid of any problems with pre-factor
   */
  int L = potential->L;
  int buf_size = shell_A->cart_size * shell_B->cart_size;
  for (int icoord = 0; icoord < 3; icoord++) {
    for (int i = 0; i < buf_size; i++) {
      grad_so_x[icoord][i] *= (2.0 * L + 1.0) / 2.0;
      grad_so_y[icoord][i] *= (2.0 * L + 1.0) / 2.0;
      grad_so_z[icoord][i] *= (2.0 * L + 1.0) / 2.0;
    }
  }

  libgrpp_dealloc_gradients(stub_grad_arep);

  grpp_operator->n_esop = 0;
  libgrpp_delete_grpp(grpp_operator);
}

void libgrpp_outercore_potential_integrals_gradient(
    libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B, double *rpp_origin,
    int num_oc_shells, libgrpp_potential_t **oc_potentials,
    libgrpp_shell_t **oc_shells, double *point_3d, double **grad_arep,
    double **grad_so_x, double **grad_so_y, double **grad_so_z) {
  libgrpp_grpp_t *grpp_operator = libgrpp_new_grpp();
  for (int ioc = 0; ioc < num_oc_shells; ioc++) {
    libgrpp_grpp_add_outercore_potential(grpp_operator, oc_potentials[ioc],
                                         oc_shells[ioc]);
  }

  libgrpp_full_grpp_integrals_gradient(shell_A, shell_B, grpp_operator,
                                       rpp_origin, point_3d, grad_arep,
                                       grad_so_x, grad_so_y, grad_so_z);

  grpp_operator->n_oc_shells = 0;
  libgrpp_delete_grpp(grpp_operator);
}

/**
 * Analytic calculation of gradients of GRPP integrals for a given shell pair
 * with respect to the point 'point_3d'.
 * (for the full GRPP operator which includes local, semi-local and non-local
 * parts)
 */
void libgrpp_full_grpp_integrals_gradient(
    libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B,
    libgrpp_grpp_t *grpp_operator, double *grpp_origin, double *point_3d,
    double **grad_arep, double **grad_so_x, double **grad_so_y,
    double **grad_so_z) {
  int cart_size_A = shell_A->cart_size;
  int cart_size_B = shell_B->cart_size;
  int buf_size = cart_size_A * cart_size_B;

  /*
   * initialization: set all gradients to zero
   */
  for (int icoord = 0; icoord < 3; icoord++) {
    memset(grad_arep[icoord], 0, sizeof(double) * buf_size);
    memset(grad_so_x[icoord], 0, sizeof(double) * buf_size);
    memset(grad_so_y[icoord], 0, sizeof(double) * buf_size);
    memset(grad_so_z[icoord], 0, sizeof(double) * buf_size);
  }

  /*
   * d<AAA>/d... = 0
   */
  if (points_are_equal(shell_A->origin, grpp_origin) &&
      points_are_equal(shell_B->origin, grpp_origin)) {
    return;
  }

  /*
   * d<ACB>/dD = 0
   */
  if (!points_are_equal(shell_A->origin, point_3d) &&
      !points_are_equal(shell_B->origin, point_3d) &&
      !points_are_equal(grpp_origin, point_3d)) {
    return;
  }

  double *A = shell_A->origin;
  double *B = shell_B->origin;
  double *C = grpp_origin;
  double *D = point_3d;

  const int diff_bra = 1;
  const int diff_ket = 0;

  /*
   * Type ACB
   */
  if (!points_are_equal(A, C) && !points_are_equal(C, B) &&
      !points_are_equal(A, B)) {
    if (points_are_equal(A, D)) {
      grpp_gradient_contribution(shell_A, shell_B, grpp_operator, grpp_origin,
                                 grad_arep, grad_so_x, grad_so_y, grad_so_z,
                                 diff_bra, +1.0);
    }
    if (points_are_equal(B, D)) {
      grpp_gradient_contribution(shell_A, shell_B, grpp_operator, grpp_origin,
                                 grad_arep, grad_so_x, grad_so_y, grad_so_z,
                                 diff_ket, +1.0);
    }
    if (points_are_equal(C, D)) {
      grpp_gradient_contribution(shell_A, shell_B, grpp_operator, grpp_origin,
                                 grad_arep, grad_so_x, grad_so_y, grad_so_z,
                                 diff_bra, -1.0);
      grpp_gradient_contribution(shell_A, shell_B, grpp_operator, grpp_origin,
                                 grad_arep, grad_so_x, grad_so_y, grad_so_z,
                                 diff_ket, -1.0);
    }
  }

  /*
   * Type ACA
   */
  if (points_are_equal(A, B) && !points_are_equal(A, C)) {
    if (points_are_equal(A, D)) {
      grpp_gradient_diff_bra_contribution(shell_A, shell_B, grpp_operator,
                                          grpp_origin, grad_arep, grad_so_x,
                                          grad_so_y, grad_so_z, +1.0);
      grpp_gradient_diff_ket_contribution(shell_A, shell_B, grpp_operator,
                                          grpp_origin, grad_arep, grad_so_x,
                                          grad_so_y, grad_so_z, +1.0);
    } else {
      grpp_gradient_diff_bra_contribution(shell_A, shell_B, grpp_operator,
                                          grpp_origin, grad_arep, grad_so_x,
                                          grad_so_y, grad_so_z, -1.0);
      grpp_gradient_diff_ket_contribution(shell_A, shell_B, grpp_operator,
                                          grpp_origin, grad_arep, grad_so_x,
                                          grad_so_y, grad_so_z, -1.0);
    }
  }

  /*
   * Type ACC
   */
  if (!points_are_equal(A, C) && points_are_equal(C, B)) {
    if (points_are_equal(A, D)) {
      grpp_gradient_contribution(shell_A, shell_B, grpp_operator, grpp_origin,
                                 grad_arep, grad_so_x, grad_so_y, grad_so_z,
                                 diff_bra, +1.0);
    } else {
      grpp_gradient_contribution(shell_A, shell_B, grpp_operator, grpp_origin,
                                 grad_arep, grad_so_x, grad_so_y, grad_so_z,
                                 diff_bra, -1.0);
    }
  }

  /*
   * Type CCB
   */
  if (points_are_equal(A, C) && !points_are_equal(C, B)) {
    if (points_are_equal(B, D)) {
      grpp_gradient_contribution(shell_A, shell_B, grpp_operator, grpp_origin,
                                 grad_arep, grad_so_x, grad_so_y, grad_so_z,
                                 diff_ket, +1.0);
    } else {
      grpp_gradient_contribution(shell_A, shell_B, grpp_operator, grpp_origin,
                                 grad_arep, grad_so_x, grad_so_y, grad_so_z,
                                 diff_ket, -1.0);
    }
  }
}

/**
 * Calculates contribution to gradients arising from the < df/dA | V | g > term:
 *
 * grad += factor * < df/dA | V | g >
 *
 * (bra basis function is differentiated).
 */
void grpp_gradient_diff_bra_contribution(
    libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B,
    libgrpp_grpp_t *grpp_operator, double *grpp_origin, double **grad_arep,
    double **grad_so_x, double **grad_so_y, double **grad_so_z, double factor) {
  /*
   * calculate integrals < df/dA | V | B >
   */
  double *arep_matrix_down = NULL;
  double *so_x_matrix_down = NULL;
  double *so_y_matrix_down = NULL;
  double *so_z_matrix_down = NULL;
  double *arep_matrix_up = NULL;
  double *so_x_matrix_up = NULL;
  double *so_y_matrix_up = NULL;
  double *so_z_matrix_up = NULL;
  int cart_size_down = 0;
  int cart_size_up = 0;

  grpp_gradient_diff_bra_grpp_integrals(
      shell_A, shell_B, grpp_operator, grpp_origin, &arep_matrix_down,
      &so_x_matrix_down, &so_y_matrix_down, &so_z_matrix_down, &arep_matrix_up,
      &so_x_matrix_up, &so_y_matrix_up, &so_z_matrix_up, &cart_size_down,
      &cart_size_up);

  /*
   * construct contributions to gradients:
   * d<A|V|B>/dA += < df/dA | V | B >
   */
  for (int icoord = 0; icoord < 3; icoord++) {
    for (int i = 0; i < shell_A->cart_size; i++) {
      for (int j = 0; j < shell_B->cart_size; j++) {

        int bra_nlm[3];
        bra_nlm[0] = shell_A->cart_list[3 * i + 0];
        bra_nlm[1] = shell_A->cart_list[3 * i + 1];
        bra_nlm[2] = shell_A->cart_list[3 * i + 2];

        int ket_nlm[3];
        ket_nlm[0] = shell_B->cart_list[3 * j + 0];
        ket_nlm[1] = shell_B->cart_list[3 * j + 1];
        ket_nlm[2] = shell_B->cart_list[3 * j + 2];

        int index = i * shell_B->cart_size + j;

        /*
         * contribution from the L-1 gaussian
         */
        if (shell_A->L > 0) {
          bra_nlm[icoord] -= 1;
          int bra_index = libgrpp_nlm_to_linear(bra_nlm);
          int ket_index = libgrpp_nlm_to_linear(ket_nlm);
          bra_nlm[icoord] += 1;

          grad_arep[icoord][index] -=
              factor * bra_nlm[icoord] *
              arep_matrix_down[shell_B->cart_size * bra_index + ket_index];
          grad_so_x[icoord][index] -=
              factor * bra_nlm[icoord] *
              so_x_matrix_down[shell_B->cart_size * bra_index + ket_index];
          grad_so_y[icoord][index] -=
              factor * bra_nlm[icoord] *
              so_y_matrix_down[shell_B->cart_size * bra_index + ket_index];
          grad_so_z[icoord][index] -=
              factor * bra_nlm[icoord] *
              so_z_matrix_down[shell_B->cart_size * bra_index + ket_index];
        }

        /*
         * contribution from the L+1 gaussian
         */
        bra_nlm[icoord] += 1;
        int bra_index = libgrpp_nlm_to_linear(bra_nlm);
        int ket_index = libgrpp_nlm_to_linear(ket_nlm);
        bra_nlm[icoord] -= 1;

        grad_arep[icoord][index] +=
            factor * arep_matrix_up[shell_B->cart_size * bra_index + ket_index];
        grad_so_x[icoord][index] +=
            factor * so_x_matrix_up[shell_B->cart_size * bra_index + ket_index];
        grad_so_y[icoord][index] +=
            factor * so_y_matrix_up[shell_B->cart_size * bra_index + ket_index];
        grad_so_z[icoord][index] +=
            factor * so_z_matrix_up[shell_B->cart_size * bra_index + ket_index];
      }
    }
  }

  if (arep_matrix_down) {
    free(arep_matrix_down);
    free(so_x_matrix_down);
    free(so_y_matrix_down);
    free(so_z_matrix_down);
  }
  free(arep_matrix_up);
  free(so_x_matrix_up);
  free(so_y_matrix_up);
  free(so_z_matrix_up);
}

/**
 * To assemble the contribution < df/dA | V | g > to gradients, one have to
 * differentiate a Gaussian function. Such a differentiation yields two
 * Gaussians with angular momenta L-1 ("down") and L+1 ("up"): dG/dA -> G(L-1)
 * and G(L+1)
 *
 * This function constructs overlap matrices with these "downgraded" and
 * "upgraded" Gaussian functions: < G(L-1) | V | G' > and < G(L+1) | V | G' >
 *
 */
void grpp_gradient_diff_bra_grpp_integrals(
    libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B,
    libgrpp_grpp_t *grpp_operator, double *grpp_origin,
    double **arep_matrix_down, double **so_x_matrix_down,
    double **so_y_matrix_down, double **so_z_matrix_down,
    double **arep_matrix_up, double **so_x_matrix_up, double **so_y_matrix_up,
    double **so_z_matrix_up, int *cart_size_down, int *cart_size_up) {
  /*
   * differentiation of contracted Gaussian functions
   */
  libgrpp_shell_t *shell_A_down = NULL;
  libgrpp_shell_t *shell_A_up = NULL;
  libgrpp_differentiate_shell(shell_A, &shell_A_down, &shell_A_up);

  *cart_size_down = 0;
  if (shell_A_down != NULL) {
    *cart_size_down = shell_A_down->cart_size;
  }
  *cart_size_up = shell_A_up->cart_size;

  /*
   * matrix < L-1 | V | L >
   */
  if (shell_A_down != NULL) {
    size_t mat_size_down = shell_A_down->cart_size * shell_B->cart_size;
    *arep_matrix_down = (double *)calloc(mat_size_down, sizeof(double));
    *so_x_matrix_down = (double *)calloc(mat_size_down, sizeof(double));
    *so_y_matrix_down = (double *)calloc(mat_size_down, sizeof(double));
    *so_z_matrix_down = (double *)calloc(mat_size_down, sizeof(double));

    libgrpp_full_grpp_integrals(
        shell_A_down, shell_B, grpp_operator, grpp_origin, *arep_matrix_down,
        *so_x_matrix_down, *so_y_matrix_down, *so_z_matrix_down);
  } else {
    *arep_matrix_down = NULL;
    *so_x_matrix_down = NULL;
    *so_y_matrix_down = NULL;
    *so_z_matrix_down = NULL;
  }

  /*
   * matrix < L+1 | V | L >
   */
  size_t mat_size_up = shell_A_up->cart_size * shell_B->cart_size;
  *arep_matrix_up = (double *)calloc(mat_size_up, sizeof(double));
  *so_x_matrix_up = (double *)calloc(mat_size_up, sizeof(double));
  *so_y_matrix_up = (double *)calloc(mat_size_up, sizeof(double));
  *so_z_matrix_up = (double *)calloc(mat_size_up, sizeof(double));

  libgrpp_full_grpp_integrals(shell_A_up, shell_B, grpp_operator, grpp_origin,
                              *arep_matrix_up, *so_x_matrix_up, *so_y_matrix_up,
                              *so_z_matrix_up);

  /*
   * clean up
   */
  if (shell_A_down) {
    libgrpp_delete_shell(shell_A_down);
  }
  libgrpp_delete_shell(shell_A_up);
}

/**
 * Calculates contribution to gradients arising from the < df/dA | V | g > term:
 *
 * grad += factor * < f | V | dg/dA >
 *
 * (bra basis function is differentiated).
 */
void grpp_gradient_diff_ket_contribution(
    libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B,
    libgrpp_grpp_t *grpp_operator, double *grpp_origin, double **grad_arep,
    double **grad_so_x, double **grad_so_y, double **grad_so_z, double factor) {
  /*
   * calculate integrals < df/dA | V | B >
   */
  double *arep_matrix_down = NULL;
  double *so_x_matrix_down = NULL;
  double *so_y_matrix_down = NULL;
  double *so_z_matrix_down = NULL;
  double *arep_matrix_up = NULL;
  double *so_x_matrix_up = NULL;
  double *so_y_matrix_up = NULL;
  double *so_z_matrix_up = NULL;
  int cart_size_down = 0;
  int cart_size_up = 0;

  grpp_gradient_diff_ket_grpp_integrals(
      shell_A, shell_B, grpp_operator, grpp_origin, &arep_matrix_down,
      &so_x_matrix_down, &so_y_matrix_down, &so_z_matrix_down, &arep_matrix_up,
      &so_x_matrix_up, &so_y_matrix_up, &so_z_matrix_up, &cart_size_down,
      &cart_size_up);

  /*
   * construct contributions to gradients:
   * d<A|B>/dA += < df/dA | V | B >
   */
  for (int icoord = 0; icoord < 3; icoord++) {
    for (int i = 0; i < shell_A->cart_size; i++) {
      for (int j = 0; j < shell_B->cart_size; j++) {

        int bra_nlm[3];
        bra_nlm[0] = shell_A->cart_list[3 * i + 0];
        bra_nlm[1] = shell_A->cart_list[3 * i + 1];
        bra_nlm[2] = shell_A->cart_list[3 * i + 2];

        int ket_nlm[3];
        ket_nlm[0] = shell_B->cart_list[3 * j + 0];
        ket_nlm[1] = shell_B->cart_list[3 * j + 1];
        ket_nlm[2] = shell_B->cart_list[3 * j + 2];

        int index = i * shell_B->cart_size + j;

        /*
         * contribution from the L-1 gaussian
         */
        if (shell_B->L > 0) {
          ket_nlm[icoord] -= 1;
          int bra_index = libgrpp_nlm_to_linear(bra_nlm);
          int ket_index = libgrpp_nlm_to_linear(ket_nlm);
          ket_nlm[icoord] += 1;

          grad_arep[icoord][index] -=
              factor * ket_nlm[icoord] *
              arep_matrix_down[cart_size_down * bra_index + ket_index];
          grad_so_x[icoord][index] -=
              factor * ket_nlm[icoord] *
              so_x_matrix_down[cart_size_down * bra_index + ket_index];
          grad_so_y[icoord][index] -=
              factor * ket_nlm[icoord] *
              so_y_matrix_down[cart_size_down * bra_index + ket_index];
          grad_so_z[icoord][index] -=
              factor * ket_nlm[icoord] *
              so_z_matrix_down[cart_size_down * bra_index + ket_index];
        }

        /*
         * contribution from the L+1 gaussian
         */
        ket_nlm[icoord] += 1;
        int bra_index = libgrpp_nlm_to_linear(bra_nlm);
        int ket_index = libgrpp_nlm_to_linear(ket_nlm);
        ket_nlm[icoord] -= 1;

        grad_arep[icoord][index] +=
            factor * arep_matrix_up[cart_size_up * bra_index + ket_index];
        grad_so_x[icoord][index] +=
            factor * so_x_matrix_up[cart_size_up * bra_index + ket_index];
        grad_so_y[icoord][index] +=
            factor * so_y_matrix_up[cart_size_up * bra_index + ket_index];
        grad_so_z[icoord][index] +=
            factor * so_z_matrix_up[cart_size_up * bra_index + ket_index];
      }
    }
  }

  if (arep_matrix_down) {
    free(arep_matrix_down);
    free(so_x_matrix_down);
    free(so_y_matrix_down);
    free(so_z_matrix_down);
  }
  free(arep_matrix_up);
  free(so_x_matrix_up);
  free(so_y_matrix_up);
  free(so_z_matrix_up);
}

/**
 * To assemble the contribution < df/dA | V | g > to gradients, one have to
 * differentiate Gaussian function. Such a differentiation yields two Gaussians
 * with angular momenta L-1 ("down") and L+1 ("up"): dG/dA -> G(L-1) and G(L+1)
 *
 * This function constructs matrices with these "downgraded" and "upgraded"
 * Gaussian functions:
 * < G(L-1) | V | G' > and < G(L+1) | V | G' >
 *
 */
void grpp_gradient_diff_ket_grpp_integrals(
    libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B,
    libgrpp_grpp_t *grpp_operator, double *grpp_origin,
    double **arep_matrix_down, double **so_x_matrix_down,
    double **so_y_matrix_down, double **so_z_matrix_down,
    double **arep_matrix_up, double **so_x_matrix_up, double **so_y_matrix_up,
    double **so_z_matrix_up, int *cart_size_down, int *cart_size_up) {
  /*
   * differentiation of contracted Gaussian functions
   */
  libgrpp_shell_t *shell_B_down = NULL;
  libgrpp_shell_t *shell_B_up = NULL;
  libgrpp_differentiate_shell(shell_B, &shell_B_down, &shell_B_up);

  *cart_size_down = 0;
  if (shell_B_down != NULL) {
    *cart_size_down = shell_B_down->cart_size;
  }
  *cart_size_up = shell_B_up->cart_size;

  /*
   * matrix < L-1 | L>
   */
  if (shell_B_down != NULL) {
    size_t mat_size_down = shell_A->cart_size * shell_B_down->cart_size;
    *arep_matrix_down = (double *)calloc(mat_size_down, sizeof(double));
    *so_x_matrix_down = (double *)calloc(mat_size_down, sizeof(double));
    *so_y_matrix_down = (double *)calloc(mat_size_down, sizeof(double));
    *so_z_matrix_down = (double *)calloc(mat_size_down, sizeof(double));

    libgrpp_full_grpp_integrals(
        // evaluate_grpp_integrals_shell_pair(
        shell_A, shell_B_down, grpp_operator, grpp_origin, *arep_matrix_down,
        *so_x_matrix_down, *so_y_matrix_down, *so_z_matrix_down);
  } else {
    *arep_matrix_down = NULL;
    *so_x_matrix_down = NULL;
    *so_y_matrix_down = NULL;
    *so_z_matrix_down = NULL;
  }

  /*
   * matrix < L+1 | L>
   */
  size_t mat_size_up = shell_A->cart_size * shell_B_up->cart_size;
  *arep_matrix_up = (double *)calloc(mat_size_up, sizeof(double));
  *so_x_matrix_up = (double *)calloc(mat_size_up, sizeof(double));
  *so_y_matrix_up = (double *)calloc(mat_size_up, sizeof(double));
  *so_z_matrix_up = (double *)calloc(mat_size_up, sizeof(double));

  libgrpp_full_grpp_integrals(
      // evaluate_grpp_integrals_shell_pair(
      shell_A, shell_B_up, grpp_operator, grpp_origin, *arep_matrix_up,
      *so_x_matrix_up, *so_y_matrix_up, *so_z_matrix_up);

  /*
   * clean up
   */
  if (shell_B_down) {
    libgrpp_delete_shell(shell_B_down);
  }
  libgrpp_delete_shell(shell_B_up);
}

void grpp_gradient_contribution(libgrpp_shell_t *shell_A,
                                libgrpp_shell_t *shell_B,
                                libgrpp_grpp_t *grpp_operator,
                                double *grpp_origin, double **grad_arep,
                                double **grad_so_x, double **grad_so_y,
                                double **grad_so_z, int diff_bra,
                                double factor) {
  //  int diff_ket = 0;

  if (diff_bra == 0) {
    diff_bra = 0;
    //    diff_ket = 1;
  } else {
    diff_bra = 1;
    //    diff_ket = 0;
  }

  /*
   * calculate overlap integrals < df/dA | V | B >
   */
  double *arep_matrix_down = NULL;
  double *so_x_matrix_down = NULL;
  double *so_y_matrix_down = NULL;
  double *so_z_matrix_down = NULL;
  double *arep_matrix_up = NULL;
  double *so_x_matrix_up = NULL;
  double *so_y_matrix_up = NULL;
  double *so_z_matrix_up = NULL;
  int cart_size_down = 0;
  int cart_size_up = 0;

  grpp_gradient_diff_gaussian(
      shell_A, shell_B, grpp_operator, grpp_origin, &arep_matrix_down,
      &so_x_matrix_down, &so_y_matrix_down, &so_z_matrix_down, &arep_matrix_up,
      &so_x_matrix_up, &so_y_matrix_up, &so_z_matrix_up, &cart_size_down,
      &cart_size_up, diff_bra);

  /*
   * construct contributions to gradients:
   * d<A|U|B>/dA += < df/dA | U | B >
   */
  for (int icoord = 0; icoord < 3; icoord++) {
    for (int i = 0; i < shell_A->cart_size; i++) {
      for (int j = 0; j < shell_B->cart_size; j++) {

        int bra_nlm[3];
        bra_nlm[0] = shell_A->cart_list[3 * i + 0];
        bra_nlm[1] = shell_A->cart_list[3 * i + 1];
        bra_nlm[2] = shell_A->cart_list[3 * i + 2];

        int ket_nlm[3];
        ket_nlm[0] = shell_B->cart_list[3 * j + 0];
        ket_nlm[1] = shell_B->cart_list[3 * j + 1];
        ket_nlm[2] = shell_B->cart_list[3 * j + 2];

        int index = i * shell_B->cart_size + j;

        int *diff_nlm = diff_bra ? bra_nlm : ket_nlm;

        /*
         * contribution from the L-1 gaussian
         */
        if (cart_size_down > 0) {
          diff_nlm[icoord] -= 1;
          int bra_index = libgrpp_nlm_to_linear(bra_nlm);
          int ket_index = libgrpp_nlm_to_linear(ket_nlm);
          diff_nlm[icoord] += 1;

          int n = diff_nlm[icoord];
          int row_len = diff_bra ? shell_B->cart_size : cart_size_down;
          int index_down = row_len * bra_index + ket_index;

          grad_arep[icoord][index] -= factor * n * arep_matrix_down[index_down];
          grad_so_x[icoord][index] -= factor * n * so_x_matrix_down[index_down];
          grad_so_y[icoord][index] -= factor * n * so_y_matrix_down[index_down];
          grad_so_z[icoord][index] -= factor * n * so_z_matrix_down[index_down];
        }

        /*
         * contribution from the L+1 gaussian
         */
        diff_nlm[icoord] += 1;
        int bra_index = libgrpp_nlm_to_linear(bra_nlm);
        int ket_index = libgrpp_nlm_to_linear(ket_nlm);
        diff_nlm[icoord] -= 1;

        int row_len = diff_bra ? shell_B->cart_size : cart_size_up;
        int index_up = row_len * bra_index + ket_index;

        grad_arep[icoord][index] += factor * arep_matrix_up[index_up];
        grad_so_x[icoord][index] += factor * so_x_matrix_up[index_up];
        grad_so_y[icoord][index] += factor * so_y_matrix_up[index_up];
        grad_so_z[icoord][index] += factor * so_z_matrix_up[index_up];
      }
    }
  }

  if (arep_matrix_down) {
    free(arep_matrix_down);
    free(so_x_matrix_down);
    free(so_y_matrix_down);
    free(so_z_matrix_down);
  }
  free(arep_matrix_up);
  free(so_x_matrix_up);
  free(so_y_matrix_up);
  free(so_z_matrix_up);
}

/**
 * To assemble the contribution < df/dA | V | g > to gradients, one have to
 * differentiate Gaussian function. Such a differentiation yields two Gaussians
 * with angular momenta L-1 ("down") and L+1 ("up"): dG/dA -> G(L-1) and G(L+1)
 *
 * This function constructs matrices with these "downgraded" and "upgraded"
 * Gaussian functions:
 * < G(L-1) | V | G' > and < G(L+1) | V | G' >
 *
 */
void grpp_gradient_diff_gaussian(
    libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B,
    libgrpp_grpp_t *grpp_operator, double *grpp_origin,
    double **arep_matrix_down, double **so_x_matrix_down,
    double **so_y_matrix_down, double **so_z_matrix_down,
    double **arep_matrix_up, double **so_x_matrix_up, double **so_y_matrix_up,
    double **so_z_matrix_up, int *cart_size_down, int *cart_size_up,
    int diff_bra) {
  // int diff_ket = 0;

  if (diff_bra == 0) {
    diff_bra = 0;
    //    diff_ket = 1;
  } else {
    diff_bra = 1;
    //    diff_ket = 0;
  }

  /*
   * which shell should be differentiated, bra or ket
   */
  libgrpp_shell_t *const_shell = NULL;
  libgrpp_shell_t *diff_shell = NULL;
  if (diff_bra) {
    diff_shell = shell_A;
    const_shell = shell_B;
  } else {
    diff_shell = shell_B;
    const_shell = shell_A;
  }

  /*
   * differentiation of contracted Gaussian functions
   */
  libgrpp_shell_t *shell_down = NULL;
  libgrpp_shell_t *shell_up = NULL;
  libgrpp_differentiate_shell(diff_shell, &shell_down, &shell_up);

  *cart_size_down = 0;
  if (shell_down != NULL) {
    *cart_size_down = shell_down->cart_size;
  }
  *cart_size_up = shell_up->cart_size;

  /*
   * GRPP matrix:
   * < L-1 | U | L > or < L | U | L-1 >
   */
  if (shell_down != NULL) {
    size_t mat_size_down = const_shell->cart_size * shell_down->cart_size;
    *arep_matrix_down = (double *)calloc(mat_size_down, sizeof(double));
    *so_x_matrix_down = (double *)calloc(mat_size_down, sizeof(double));
    *so_y_matrix_down = (double *)calloc(mat_size_down, sizeof(double));
    *so_z_matrix_down = (double *)calloc(mat_size_down, sizeof(double));

    libgrpp_full_grpp_integrals(
        diff_bra ? shell_down : shell_A, diff_bra ? shell_B : shell_down,
        grpp_operator, grpp_origin, *arep_matrix_down, *so_x_matrix_down,
        *so_y_matrix_down, *so_z_matrix_down);
  } else {
    *arep_matrix_down = NULL;
    *so_x_matrix_down = NULL;
    *so_y_matrix_down = NULL;
    *so_z_matrix_down = NULL;
  }

  /*
   * GRPP matrix:
   * < L+1 | U | L > or < L | U | L+1 >
   */
  size_t mat_size_up = const_shell->cart_size * shell_up->cart_size;
  *arep_matrix_up = (double *)calloc(mat_size_up, sizeof(double));
  *so_x_matrix_up = (double *)calloc(mat_size_up, sizeof(double));
  *so_y_matrix_up = (double *)calloc(mat_size_up, sizeof(double));
  *so_z_matrix_up = (double *)calloc(mat_size_up, sizeof(double));

  libgrpp_full_grpp_integrals(diff_bra ? shell_up : shell_A,
                              diff_bra ? shell_B : shell_up, grpp_operator,
                              grpp_origin, *arep_matrix_up, *so_x_matrix_up,
                              *so_y_matrix_up, *so_z_matrix_up);

  /*
   * clean up
   */
  if (shell_down) {
    libgrpp_delete_shell(shell_down);
  }
  libgrpp_delete_shell(shell_up);
}

/**
 * Allocates memory for gradients for a given shell pair.
 */
double **libgrpp_alloc_gradients(libgrpp_shell_t *bra, libgrpp_shell_t *ket) {
  size_t size = bra->cart_size * ket->cart_size;

  double **grad = (double **)calloc(3, sizeof(double *));
  for (int icoord = 0; icoord < 3; icoord++) {
    grad[icoord] = (double *)calloc(size, sizeof(double));
  }

  return grad;
}

/**
 * Deallocates arrays containing gradients of AO integrals.
 */
void libgrpp_dealloc_gradients(double **grad) {
  free(grad[0]);
  free(grad[1]);
  free(grad[2]);
  free(grad);
}
