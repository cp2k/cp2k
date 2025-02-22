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
 * representation of atom-centered shell of contracted Gaussian functions
 */
#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433
#endif

#include "libgrpp.h"

#include "grpp_norm_gaussian.h"

/**
 * constructs new object representing a shell; returns pointer to it.
 */
libgrpp_shell_t *libgrpp_new_shell(double *origin, int L, int num_primitives,
                                   double *coeffs, double *alpha) {
  libgrpp_shell_t *shell = (libgrpp_shell_t *)malloc(sizeof(libgrpp_shell_t));

  shell->L = L;
  shell->origin[0] = origin[0];
  shell->origin[1] = origin[1];
  shell->origin[2] = origin[2];
  shell->cart_size = (L + 1) * (L + 2) / 2;
  shell->cart_list = libgrpp_generate_shell_cartesians(L);

  shell->num_primitives = num_primitives;
  shell->coeffs = (double *)calloc(num_primitives, sizeof(double));
  shell->alpha = (double *)calloc(num_primitives, sizeof(double));
  for (int i = 0; i < num_primitives; i++) {
    shell->coeffs[i] = coeffs[i];
    shell->alpha[i] = alpha[i];
  }

  return shell;
}

/**
 * creates deep copy of the 'libgrpp_shell_t' object
 */
libgrpp_shell_t *libgrpp_shell_deep_copy(libgrpp_shell_t *src_shell) {
  libgrpp_shell_t *new_shell = libgrpp_new_shell(
      src_shell->origin, src_shell->L, src_shell->num_primitives,
      src_shell->coeffs, src_shell->alpha);

  return new_shell;
}

/**
 * removes primitive gaussians (from the contracted function)
 * with zero coefficients
 */
void libgrpp_shell_shrink(libgrpp_shell_t *shell) {
  int nprim = 0;

  for (int i = 0; i < shell->num_primitives; i++) {
    if (fabs(shell->coeffs[i]) > LIBGRPP_ZERO_THRESH) {
      shell->coeffs[nprim] = shell->coeffs[i];
      shell->alpha[nprim] = shell->alpha[i];
      nprim++;
    }
  }

  shell->num_primitives = nprim;
}

/**
 * multiplies coefficients of the primitive gaussians by their normalization
 * factors
 */
void libgrpp_shell_mult_normcoef(libgrpp_shell_t *shell) {
  for (int i = 0; i < shell->num_primitives; i++) {
    double norm_factor =
        libgrpp_gaussian_norm_factor(shell->L, 0, 0, shell->alpha[i]);
    shell->coeffs[i] *= norm_factor;
  }
}

/**
 * returns number of Cartesian primitives encapsulated inside the shell
 */
int libgrpp_get_shell_size(libgrpp_shell_t *shell) { return shell->cart_size; }

/**
 * destructor for the shell object
 */
void libgrpp_delete_shell(libgrpp_shell_t *shell) {
  free(shell->cart_list);
  free(shell->coeffs);
  free(shell->alpha);
  free(shell);
}

int *libgrpp_generate_shell_cartesians(int L) {
  int ncart = (L + 1) * (L + 2) / 2;

  int *cart_list = (int *)calloc(3 * ncart, sizeof(int));
  libgrpp_params.cartesian_generator(L, cart_list);

  return cart_list;
}
