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
 * Calculation of nuclear attraction integrals.
 *
 * For the point charge nuclear model the recursive Obara-Saika scheme is used
 * to calculate nuclear attraction integrals. For details, see
 * T. Helgaker, P. Jorgensen, J. Olsen, Molecular Electronic-Structure Theory,
 * John Wiley & Sons Ltd, 2000.
 * Chapter 9.10.1, "The Obara-Saika scheme for one-electron Coulomb integrals"
 *
 * For the other three models,
 * - uniformly charged ball
 * - Gaussian nucleus
 * - Fermi nucleus,
 * the scheme is actually the same as for the type 1 (radially-local) ECP
 * integrals. Electrostatic potential V(r) induced by the finite nuclear charge
 * distribution is integrated numerically on a radial grid.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "grpp_norm_gaussian.h"
#include "grpp_nuclear_models.h"
#include "grpp_utils.h"
#include "libgrpp.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void libgrpp_evaluate_radially_local_potential_integral_primitive_gaussians(
    double *A, int n_cart_A, int *cart_list_A, double alpha_A, double *B,
    int n_cart_B, int *cart_list_B, double alpha_B, double *C,
    double (*potential)(double r, void *params), void *potential_params,
    double *matrix);

void libgrpp_evaluate_rpp_type1_mmd_n1_primitive_shell_pair(
    libgrpp_shell_t *shell_A, double alpha_A, libgrpp_shell_t *shell_B,
    double alpha_B, double *rpp_origin, double rpp_alpha, double *rpp_matrix);

static double wrapper_coulomb_potential_point(double r, void *params);

static double wrapper_coulomb_potential_ball(double r, void *params);

static double wrapper_coulomb_potential_gaussian(double r, void *params);

static double wrapper_coulomb_potential_fermi(double r, void *params);

static double wrapper_coulomb_potential_fermi_bubble(double r, void *params);

/**
 * Calculates nuclear attraction integral between two shells
 * represented by contracted Gaussian functions.
 *
 * nuclear model should be one of:
 *   LIBGRPP_NUCLEAR_MODEL_POINT_CHARGE
 *   LIBGRPP_NUCLEAR_MODEL_CHARGED_BALL
 *   LIBGRPP_NUCLEAR_MODEL_GAUSSIAN
 *   LIBGRPP_NUCLEAR_MODEL_FERMI
 *   LIBGRPP_NUCLEAR_MODEL_FERMI_BUBBLE
 *   LIBGRPP_NUCLEAR_MODEL_POINT_CHARGE_NUMERICAL
 */
void libgrpp_nuclear_attraction_integrals(libgrpp_shell_t *shell_A,
                                          libgrpp_shell_t *shell_B,
                                          double *charge_origin, int charge,
                                          int nuclear_model,
                                          double *model_params,
                                          double *coulomb_matrix) {
  assert(libgrpp_is_initialized());

  int size_A = libgrpp_get_shell_size(shell_A);
  int size_B = libgrpp_get_shell_size(shell_B);

  double *buf = calloc(size_A * size_B, sizeof(double));

  memset(coulomb_matrix, 0, size_A * size_B * sizeof(double));

  // loop over primitives in contractions
  for (int i = 0; i < shell_A->num_primitives; i++) {
    double coef_A_i = shell_A->coeffs[i];

    for (int j = 0; j < shell_B->num_primitives; j++) {
      double coef_B_j = shell_B->coeffs[j];

      if (fabs(coef_A_i * coef_B_j) < 1e-13) {
        continue;
      }

      if (nuclear_model == LIBGRPP_NUCLEAR_MODEL_POINT_CHARGE) {

        // use code for RPP type-1 integrals with RPP exponent = 0.0
        libgrpp_evaluate_rpp_type1_mmd_n1_primitive_shell_pair(
            shell_A, shell_A->alpha[i], shell_B, shell_B->alpha[j],
            charge_origin, 0.0, buf);

        for (int k = 0; k < size_A * size_B; k++) {
          buf[k] *= (-1) * charge;
        }
      } else if (nuclear_model == LIBGRPP_NUCLEAR_MODEL_CHARGED_BALL ||
                 nuclear_model == LIBGRPP_NUCLEAR_MODEL_GAUSSIAN ||
                 nuclear_model == LIBGRPP_NUCLEAR_MODEL_FERMI ||
                 nuclear_model == LIBGRPP_NUCLEAR_MODEL_FERMI_BUBBLE ||
                 nuclear_model ==
                     LIBGRPP_NUCLEAR_MODEL_POINT_CHARGE_NUMERICAL) {

        double params[10];
        params[0] = charge;

        /*
         * choose nuclear model
         */
        double (*electrostatic_potential_fun)(double, void *) = NULL;

        if (nuclear_model == LIBGRPP_NUCLEAR_MODEL_POINT_CHARGE_NUMERICAL) {
          // printf("charge distribution: point\n");
          electrostatic_potential_fun = wrapper_coulomb_potential_point;
        } else if (nuclear_model == LIBGRPP_NUCLEAR_MODEL_CHARGED_BALL) {
          params[1] = model_params[0]; // R_rms
          electrostatic_potential_fun = wrapper_coulomb_potential_ball;
        } else if (nuclear_model == LIBGRPP_NUCLEAR_MODEL_GAUSSIAN) {
          params[1] = model_params[0]; // R_rms
          electrostatic_potential_fun = wrapper_coulomb_potential_gaussian;
        } else if (nuclear_model == LIBGRPP_NUCLEAR_MODEL_FERMI) {
          params[1] = model_params[0]; // c
          params[2] = model_params[1]; // a
          electrostatic_potential_fun = wrapper_coulomb_potential_fermi;
        } else {
          params[1] = model_params[0]; // c
          params[2] = model_params[1]; // a
          params[3] = model_params[2]; // k
          electrostatic_potential_fun = wrapper_coulomb_potential_fermi_bubble;
        }

        /*
         * calculate integrals for the shell pair
         */
        libgrpp_evaluate_radially_local_potential_integral_primitive_gaussians(
            shell_A->origin, size_A, shell_A->cart_list, shell_A->alpha[i],
            shell_B->origin, size_B, shell_B->cart_list, shell_B->alpha[j],
            charge_origin, electrostatic_potential_fun, params, buf);
      } else {
        printf("LIBGRPP: unknown finite nuclear charge distribution model!\n");
        exit(0);
      }

      libgrpp_daxpy(size_A * size_B, coef_A_i * coef_B_j, buf, coulomb_matrix);
    }
  }

  free(buf);
}

/**
 * Calculates nuclear attraction integral between two shells
 * represented by contracted Gaussian functions for the electrostatic potential
 * generated by the point charge.
 */
void libgrpp_nuclear_attraction_integrals_point_charge(libgrpp_shell_t *shell_A,
                                                       libgrpp_shell_t *shell_B,
                                                       double *charge_origin,
                                                       int charge,
                                                       double *coulomb_matrix) {
  libgrpp_nuclear_attraction_integrals(shell_A, shell_B, charge_origin, charge,
                                       LIBGRPP_NUCLEAR_MODEL_POINT_CHARGE, NULL,
                                       coulomb_matrix);
}

/**
 * Calculates nuclear attraction integral between two shells
 * represented by contracted Gaussian functions for the electrostatic potential
 * generated by the charged ball.
 *
 * r_rms stands for the root mean square radius (in bohrs)
 */
void libgrpp_nuclear_attraction_integrals_charged_ball(libgrpp_shell_t *shell_A,
                                                       libgrpp_shell_t *shell_B,
                                                       double *charge_origin,
                                                       int charge, double r_rms,
                                                       double *coulomb_matrix) {
  double params[10];
  params[0] = charge;
  params[1] = r_rms;

  libgrpp_nuclear_attraction_integrals(shell_A, shell_B, charge_origin, charge,
                                       LIBGRPP_NUCLEAR_MODEL_CHARGED_BALL,
                                       params, coulomb_matrix);
}

/**
 * Calculates nuclear attraction integral between two shells
 * represented by contracted Gaussian functions for the electrostatic potential
 * generated by the Gaussian distribution.
 *
 * r_rms stands for the root mean square radius (in bohrs)
 */
void libgrpp_nuclear_attraction_integrals_gaussian_model(
    libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B, double *charge_origin,
    int charge, double r_rms, double *coulomb_matrix) {
  double params[10];
  params[0] = charge;
  params[1] = r_rms;

  libgrpp_nuclear_attraction_integrals(shell_A, shell_B, charge_origin, charge,
                                       LIBGRPP_NUCLEAR_MODEL_GAUSSIAN, params,
                                       coulomb_matrix);
}

/**
 * Calculates nuclear attraction integral between two shells
 * represented by contracted Gaussian functions for the electrostatic potential
 * generated by the Fermi distribution.
 *
 * Model parameters 'c' and 'a' must be given in bohrs.
 */
void libgrpp_nuclear_attraction_integrals_fermi_model(
    libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B, double *charge_origin,
    int charge, double fermi_param_c, double fermi_param_a,
    double *coulomb_matrix) {
  double params[10];
  params[0] = charge;
  params[1] = fermi_param_c;
  params[2] = fermi_param_a;

  libgrpp_nuclear_attraction_integrals(shell_A, shell_B, charge_origin, charge,
                                       LIBGRPP_NUCLEAR_MODEL_FERMI, params,
                                       coulomb_matrix);
}

/**
 * Calculates nuclear attraction integral between two shells
 * represented by contracted Gaussian functions for the electrostatic potential
 * generated by the "Fermi + bubble" distribution.
 *
 * Model parameters 'c' and 'a' must be given in bohrs.
 * The 'k' constant is dimensionless.
 */
void libgrpp_nuclear_attraction_integrals_fermi_bubble_model(
    libgrpp_shell_t *shell_A, libgrpp_shell_t *shell_B, double *charge_origin,
    int charge, double param_c, double param_a, double param_k,
    double *coulomb_matrix) {
  double params[10];
  params[0] = charge;
  params[1] = param_c;
  params[2] = param_a;
  params[3] = param_k;

  libgrpp_nuclear_attraction_integrals(shell_A, shell_B, charge_origin, charge,
                                       LIBGRPP_NUCLEAR_MODEL_FERMI_BUBBLE,
                                       params, coulomb_matrix);
}

/**
 * wrappers for charge distribution functions.
 * are used to provide a unified interface to radially-local potentials.
 * the 'params' argument is unpacked, then the specific routines are invoked.
 */

static double wrapper_coulomb_potential_point(double r, void *params) {
  double Z = ((double *)params)[0];

  return libgrpp_coulomb_potential_point(r, Z);
}

static double wrapper_coulomb_potential_ball(double r, void *params) {
  double Z = ((double *)params)[0];
  double R_rms = ((double *)params)[1];

  return libgrpp_coulomb_potential_ball(r, Z, R_rms);
}

static double wrapper_coulomb_potential_gaussian(double r, void *params) {
  double Z = ((double *)params)[0];
  double R_rms = ((double *)params)[1];

  return libgrpp_coulomb_potential_gaussian(r, Z, R_rms);
}

static double wrapper_coulomb_potential_fermi(double r, void *params) {
  double Z = ((double *)params)[0];
  double c = ((double *)params)[1];
  double a = ((double *)params)[2];

  return libgrpp_coulomb_potential_fermi(r, Z, c, a);
}

static double wrapper_coulomb_potential_fermi_bubble(double r, void *params) {
  double Z = ((double *)params)[0];
  double c = ((double *)params)[1];
  double a = ((double *)params)[2];
  double k = ((double *)params)[3];

  return libgrpp_coulomb_potential_fermi_bubble(r, Z, c, a, k);
}
