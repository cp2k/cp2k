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

#include <stdlib.h>

#define LIBGRPP_MAX_NUMBER_POTENTIALS 100

libgrpp_grpp_t *libgrpp_new_grpp() {
  return (libgrpp_grpp_t *)calloc(1, sizeof(libgrpp_grpp_t));
}

void libgrpp_grpp_set_local_potential(libgrpp_grpp_t *grpp,
                                      libgrpp_potential_t *pot) {
  if (grpp->U_L != NULL) {
    libgrpp_delete_potential(grpp->U_L);
  }

  grpp->U_L = pot;
}

void libgrpp_grpp_add_averaged_potential(libgrpp_grpp_t *grpp,
                                         libgrpp_potential_t *pot) {
  if (grpp->n_arep == 0) {
    grpp->U_arep = (libgrpp_potential_t **)calloc(
        LIBGRPP_MAX_NUMBER_POTENTIALS, sizeof(libgrpp_potential_t *));
  }

  grpp->U_arep[grpp->n_arep++] = pot;
}

void libgrpp_grpp_add_spin_orbit_potential(libgrpp_grpp_t *grpp,
                                           libgrpp_potential_t *pot) {
  if (grpp->n_esop == 0) {
    grpp->U_esop = (libgrpp_potential_t **)calloc(
        LIBGRPP_MAX_NUMBER_POTENTIALS, sizeof(libgrpp_potential_t *));
  }

  grpp->U_esop[grpp->n_esop++] = pot;
}

void libgrpp_grpp_add_outercore_potential(libgrpp_grpp_t *grpp,
                                          libgrpp_potential_t *pot,
                                          libgrpp_shell_t *oc_shell) {
  if (grpp->n_oc_shells == 0) {
    grpp->U_oc = (libgrpp_potential_t **)calloc(LIBGRPP_MAX_NUMBER_POTENTIALS,
                                                sizeof(libgrpp_potential_t *));
    grpp->oc_shells = (libgrpp_shell_t **)calloc(LIBGRPP_MAX_NUMBER_POTENTIALS,
                                                 sizeof(libgrpp_shell_t *));
  }

  grpp->U_oc[grpp->n_oc_shells] = pot;
  grpp->oc_shells[grpp->n_oc_shells] = oc_shell;
  grpp->n_oc_shells++;
}

void libgrpp_delete_grpp(libgrpp_grpp_t *grpp) {
  if (grpp == NULL) {
    return;
  }

  /*
   * scalar-relativistic part
   */
  if (grpp->U_L != NULL) {
    libgrpp_delete_potential(grpp->U_L);
  }

  for (int i = 0; i < grpp->n_arep; i++) {
    if (grpp->U_arep[i] != NULL) {
      libgrpp_delete_potential(grpp->U_arep[i]);
    }
  }
  free(grpp->U_arep);

  /*
   * effective spin-orbit operator
   */
  for (int i = 0; i < grpp->n_esop; i++) {
    if (grpp->U_esop[i] != NULL) {
      libgrpp_delete_potential(grpp->U_esop[i]);
    }
  }
  free(grpp->U_esop);

  /*
   * outercore shells and potentials
   */
  for (int i = 0; i < grpp->n_oc_shells; i++) {
    if (grpp->U_oc[i] != NULL) {
      libgrpp_delete_potential(grpp->U_oc[i]);
    }
    if (grpp->oc_shells[i] != NULL) {
      libgrpp_delete_shell(grpp->oc_shells[i]);
    }
  }
  free(grpp->U_oc);
  free(grpp->oc_shells);

  free(grpp);
}
