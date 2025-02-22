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
 * Differentiation of contracted Gaussian functions.
 * Derivatives are then used to calculate analytic gradients of 1-el integrals.
 */
#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "grpp_diff_gaussian.h"

#include "libgrpp.h"

static double norm_factor(double alpha, int L);

/**
 * Performs differentiation of a contracted Gaussian.
 *
 * Note that the "2 alpha" factors are absorbed into coefficients, while the 'n'
 * factor is not. The latter must be accounted for explicitly at the stage of
 * gradient construction. For more details, see: T. Helgaker, P. Jorgensen, J.
 * Olsen, Molecular Electronic-Structure Theory, John Wiley & Sons Ltd, 2000.
 * Chapter 9.2.2, "Recurrence relations for Cartesian Gaussians"
 *
 */
void libgrpp_differentiate_shell(libgrpp_shell_t *shell,
                                 libgrpp_shell_t **shell_minus,
                                 libgrpp_shell_t **shell_plus) {
  // downwards
  if (shell->L > 0) {
    *shell_minus =
        libgrpp_new_shell(shell->origin, shell->L - 1, shell->num_primitives,
                          shell->coeffs, shell->alpha);

    for (int i = 0; i < shell->num_primitives; i++) {

      double alpha = shell->alpha[i];
      double L = shell->L;

      (*shell_minus)->coeffs[i] *=
          norm_factor(alpha, L) / norm_factor(alpha, L - 1);
    }
  } else {
    *shell_minus = NULL;
  }

  // upwards
  *shell_plus =
      libgrpp_new_shell(shell->origin, shell->L + 1, shell->num_primitives,
                        shell->coeffs, shell->alpha);
  for (int i = 0; i < shell->num_primitives; i++) {

    double alpha = shell->alpha[i];
    double L = shell->L;

    (*shell_plus)->coeffs[i] *=
        2.0 * alpha * norm_factor(alpha, L) / norm_factor(alpha, L + 1);
  }
}

/**
 * Calculates normalization factor for the primitive Gaussian
 * with the exponential parameter 'alpha' and angular momentum L.
 */
static double norm_factor(double alpha, int L) {
  return pow(2 * alpha / M_PI, 0.75) * pow(4 * alpha, 0.5 * L);
}
