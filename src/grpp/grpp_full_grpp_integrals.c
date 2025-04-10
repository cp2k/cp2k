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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "grpp_utils.h"

/**
 * Evaluates integrals over the full GRPP operator consisting of three parts:
 * - scalar relativistic (local)
 * - scalar relativistic (semi-local)
 * - effective spin-orbit (semi-local)
 * - outercore potentials (non-local)
 *
 * See libgrpp.h for the definition of the libgrpp_grpp_t structure.
 */
void libgrpp_full_grpp_integrals(libgrpp_shell_t *shell_A,
                                 libgrpp_shell_t *shell_B,
                                 libgrpp_grpp_t *grpp_operator,
                                 double *grpp_origin, double *arep_matrix,
                                 double *so_x_matrix, double *so_y_matrix,
                                 double *so_z_matrix) {
  assert(libgrpp_is_initialized());

  size_t size = shell_A->cart_size * shell_B->cart_size;
  double *buf_arep = (double *)calloc(size, sizeof(double));
  double *buf_so_x = (double *)calloc(size, sizeof(double));
  double *buf_so_y = (double *)calloc(size, sizeof(double));
  double *buf_so_z = (double *)calloc(size, sizeof(double));

  memset(arep_matrix, 0, sizeof(double) * size);
  memset(so_x_matrix, 0, sizeof(double) * size);
  memset(so_y_matrix, 0, sizeof(double) * size);
  memset(so_z_matrix, 0, sizeof(double) * size);

  /*
   * radially-local ("type-1") integrals
   */
  libgrpp_type1_integrals(shell_A, shell_B, grpp_origin, grpp_operator->U_L,
                          buf_arep);
  libgrpp_daxpy(size, 1.0, buf_arep, arep_matrix);

  /*
   * semilocal AREP ("type-2") integrals
   */
  for (int L = 0; L < grpp_operator->n_arep; L++) {
    libgrpp_type2_integrals(shell_A, shell_B, grpp_origin,
                            grpp_operator->U_arep[L], buf_arep);
    libgrpp_daxpy(size, 1.0, buf_arep, arep_matrix);
  }

  /*
   * semilocal SO ("type-3") integrals
   */
  for (int i_so = 0; i_so < grpp_operator->n_esop; i_so++) {
    libgrpp_potential_t *so_potential = grpp_operator->U_esop[i_so];
    libgrpp_spin_orbit_integrals(shell_A, shell_B, grpp_origin, so_potential,
                                 buf_so_x, buf_so_y, buf_so_z);

    int L = so_potential->L;
    libgrpp_daxpy(size, 2.0 / (2 * L + 1), buf_so_x, so_x_matrix);
    libgrpp_daxpy(size, 2.0 / (2 * L + 1), buf_so_y, so_y_matrix);
    libgrpp_daxpy(size, 2.0 / (2 * L + 1), buf_so_z, so_z_matrix);
  }

  /*
   * integrals over outercore non-local potentials,
   * the part specific for GRPP.
   *
   * note that proper pre-factors for the SO part are calculated inside
   * the libgrpp_outercore_potential_integrals() procedure.
   */
  libgrpp_outercore_potential_integrals(
      shell_A, shell_B, grpp_origin, grpp_operator->n_oc_shells,
      grpp_operator->U_oc, grpp_operator->oc_shells, buf_arep, buf_so_x,
      buf_so_y, buf_so_z);

  libgrpp_daxpy(size, 1.0, buf_arep, arep_matrix);
  libgrpp_daxpy(size, 1.0, buf_so_x, so_x_matrix);
  libgrpp_daxpy(size, 1.0, buf_so_y, so_y_matrix);
  libgrpp_daxpy(size, 1.0, buf_so_z, so_z_matrix);

  /*
   * cleanup
   */
  free(buf_arep);
  free(buf_so_x);
  free(buf_so_y);
  free(buf_so_z);
}
