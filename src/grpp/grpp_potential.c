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
 * representation of (generalized) effective core potentials
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433
#endif

#include "libgrpp.h"

/**
 * constructor for the pseudopotential
 */
libgrpp_potential_t *libgrpp_new_potential(int L, int J, int num_primitives,
                                           int *powers, double *coeffs,
                                           double *alpha) {
  libgrpp_potential_t *pot =
      (libgrpp_potential_t *)calloc(1, sizeof(libgrpp_potential_t));

  pot->L = L;
  pot->J = J;
  pot->powers = (int *)calloc(num_primitives, sizeof(int));
  pot->coeffs = (double *)calloc(num_primitives, sizeof(double));
  pot->alpha = (double *)calloc(num_primitives, sizeof(double));

  pot->num_primitives = 0;
  for (int i = 0; i < num_primitives; i++) {
    if (fabs(coeffs[i]) < LIBGRPP_ZERO_THRESH) {
      continue;
    }

    pot->coeffs[pot->num_primitives] = coeffs[i];
    pot->powers[pot->num_primitives] = powers[i];
    pot->alpha[pot->num_primitives] = alpha[i];

    pot->num_primitives++;
  }

  return pot;
}

/*
 * destructor for the pseudopotential
 */
void libgrpp_delete_potential(libgrpp_potential_t *potential) {
  if (potential == NULL) {
    return;
  }

  free(potential->powers);
  free(potential->coeffs);
  free(potential->alpha);
  free(potential);
}

/*
 * calculates value of the pseudopotential at the point 'r'
 *
 * TODO: remove the invocation of 'pow()'
 */
double libgrpp_potential_value(libgrpp_potential_t *potential, double r) {
  double val = 0.0;
  double r_2 = r * r;

  for (int i = 0; i < potential->num_primitives; i++) {
    int n = potential->powers[i];
    val +=
        potential->coeffs[i] * pow(r, n - 2) * exp(-potential->alpha[i] * r_2);
  }

  return val;
}

/*
 * removes redundant (zero) primitives from the RPP.
 * argument remains constant.
 */
libgrpp_potential_t *
libgrpp_shrink_potential(libgrpp_potential_t *src_potential) {
  int n = src_potential->num_primitives;
  int *new_powers = calloc(n, sizeof(int));
  double *new_coeffs = calloc(n, sizeof(double));
  double *new_alpha = calloc(n, sizeof(double));

  int n_nonzero_primitives = 0;
  for (int i = 0; i < n; i++) {
    if (fabs(src_potential->coeffs[i]) > LIBGRPP_ZERO_THRESH) {
      new_powers[n_nonzero_primitives] = src_potential->powers[i];
      new_coeffs[n_nonzero_primitives] = src_potential->coeffs[i];
      new_alpha[n_nonzero_primitives] = src_potential->alpha[i];
      n_nonzero_primitives++;
    }
  }

  libgrpp_potential_t *new_pot = libgrpp_new_potential(
      src_potential->L, src_potential->J, n_nonzero_primitives, new_powers,
      new_coeffs, new_alpha);

  free(new_powers);
  free(new_coeffs);
  free(new_alpha);

  return new_pot;
}
