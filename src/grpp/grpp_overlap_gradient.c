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

#include <stdlib.h>
#include <string.h>

#include "grpp_overlap_gradient.h"
#include "libgrpp.h"

#include "grpp_diff_gaussian.h"
#include "grpp_utils.h"

static void overlap_gradient_diff_bra_contribution(libgrpp_shell_t *shell_A,
                                                   libgrpp_shell_t *shell_B,
                                                   double **grad,
                                                   double factor);

static void overlap_gradient_diff_bra_overlap_integrals(
    libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B, double **overlap_down,
    double **overlap_up, int *cart_size_down, int *cart_size_up);

int nlm_to_linear(int *nlm);

/**
 * Analytic calculation of gradients of overlap integrals for a given shell pair
 * with respect to the point 'point_3d'.
 */
void libgrpp_overlap_integrals_gradient(libgrpp_shell_t *shell_A,
                                        libgrpp_shell_t *shell_B,
                                        double *point_3d, double **grad) {
  int cart_size_A = libgrpp_get_shell_size(shell_A);
  int cart_size_B = libgrpp_get_shell_size(shell_B);
  int buf_size = cart_size_A * cart_size_B;

  /*
   * initializations: set gradients to zero
   */
  for (int icoord = 0; icoord < 3; icoord++) {
    memset(grad[icoord], 0, sizeof(double) * buf_size);
  }

  /*
   * integrals are zero:
   * (1) for 1-center integrals <A|A> (due to the translational invariance)
   * (2) d<A|B> / dC = 0 (integral is constant for the given 'point_3d')
   */
  if (points_are_equal(shell_A->origin, shell_B->origin)) {
    return;
  }

  /*
   * construct gradients:
   * d<A|B>/dA = + < df/dA | B >
   * d<A|B>/dB = - < df/dA | B >
   *
   * note that due to the property of translational invariance,
   * d<A|B>/dB = - d<A|B>/dA
   */
  if (points_are_equal(shell_A->origin, point_3d)) {
    overlap_gradient_diff_bra_contribution(shell_A, shell_B, grad, +1.0);
  }
  if (points_are_equal(shell_B->origin, point_3d)) {
    overlap_gradient_diff_bra_contribution(shell_A, shell_B, grad, -1.0);
  }
}

/**
 * Calculates contribution to gradients arising from the < df/dA | g > term:
 *
 * grad += factor * < df/dA | g >
 *
 * (bra basis function is differentiated).
 */
static void overlap_gradient_diff_bra_contribution(libgrpp_shell_t *shell_A,
                                                   libgrpp_shell_t *shell_B,
                                                   double **grad,
                                                   double factor) {
  /*
   * calculate overlap integrals < df/dA | B >
   */
  double *overlap_down = NULL;
  double *overlap_up = NULL;
  int cart_size_down = 0;
  int cart_size_up = 0;

  overlap_gradient_diff_bra_overlap_integrals(shell_A, shell_B, &overlap_down,
                                              &overlap_up, &cart_size_down,
                                              &cart_size_up);

  /*
   * construct contributions to gradients:
   * d<A|B>/dA += < df/dA | B >
   */
  for (int icoord = 0; icoord < 3; icoord++) {
    for (int i = 0; i < shell_A->cart_size; i++) {
      for (int j = 0; j < shell_B->cart_size; j++) {

        int *bra_nlm = shell_A->cart_list + 3 * i;
        int *ket_nlm = shell_B->cart_list + 3 * j;
        int index = i * shell_B->cart_size + j;

        /*
         * contribution from the L-1 gaussian
         */
        if (shell_A->L > 0) {
          bra_nlm[icoord] -= 1;
          int bra_index = libgrpp_nlm_to_linear(bra_nlm);
          int ket_index = libgrpp_nlm_to_linear(ket_nlm);
          bra_nlm[icoord] += 1;

          grad[icoord][index] -=
              factor * bra_nlm[icoord] *
              overlap_down[shell_B->cart_size * bra_index + ket_index];
        }

        /*
         * contribution from the L+1 gaussian
         */
        bra_nlm[icoord] += 1;
        int bra_index = libgrpp_nlm_to_linear(bra_nlm);
        int ket_index = libgrpp_nlm_to_linear(ket_nlm);
        bra_nlm[icoord] -= 1;

        grad[icoord][index] +=
            factor * overlap_up[shell_B->cart_size * bra_index + ket_index];
      }
    }
  }

  if (overlap_down) {
    free(overlap_down);
  }
  free(overlap_up);
}

/**
 * To assemble the contribution < df/dA | g > to gradients, one have to
 * differentiate Gaussian function. Such a differentiation yields two Gaussians
 * with angular momenta L-1 ("down") and L+1 ("up"): dG/dA -> G(L-1) and G(L+1)
 *
 * This function constructs overlap matrices with these "downgraded" and
 * "upgraded" Gaussian functions: < G(L-1) | G' > and < G(L+1) | G' >
 *
 */
static void overlap_gradient_diff_bra_overlap_integrals(
    libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B, double **overlap_down,
    double **overlap_up, int *cart_size_down, int *cart_size_up) {
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
   * overlap matrix:
   * < L-1 | L>
   */
  if (shell_A_down != NULL) {
    *overlap_down = (double *)calloc(
        shell_A_down->cart_size * shell_B->cart_size, sizeof(double));
    libgrpp_overlap_integrals(shell_A_down, shell_B, *overlap_down);
  } else {
    *overlap_down = NULL;
  }

  /*
   * overlap matrix:
   * < L+1 | L>
   */
  *overlap_up = (double *)calloc(shell_A_up->cart_size * shell_B->cart_size,
                                 sizeof(double));
  libgrpp_overlap_integrals(shell_A_up, shell_B, *overlap_up);

  /*
   * clean up
   */
  if (shell_A_down) {
    libgrpp_delete_shell(shell_A_down);
  }
  libgrpp_delete_shell(shell_A_up);
}

/**
 * calculates sequential ("linear") index of the (n,l,m) primitive in the
 * cartesian shell
 */
int libgrpp_nlm_to_linear(int *nlm) {
  int n = nlm[0];
  int l = nlm[1];
  int m = nlm[2];

  int L = n + l + m;
  int cart_size = (L + 1) * (L + 2) / 2;
  int *cart_list = libgrpp_generate_shell_cartesians(L);

  int index = 0;
  for (index = 0; index < cart_size; index++) {
    if (cart_list[3 * index + 0] == n && cart_list[3 * index + 1] == l &&
        cart_list[3 * index + 2] == m) {
      break;
    }
  }

  free(cart_list);

  return index;
}
