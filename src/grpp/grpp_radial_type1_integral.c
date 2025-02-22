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
 * Evaluation of type 1 radial integrals.
 *
 * The procedure in general follows that described in:
 * R. Flores-Moreno et al. Half-numerical evaluation of pseudopotential
 integrals.
 * J. Comp. Chem. 27, 1009 (2006)
 * (see formulas (12) and (13) for radial integrals invoking contracted Gaussian
 functions and RPPs).

 * In contrast to type 2 integrals, the special case of type 1 integrals Bessel
 functions
 * cannot be factorized, and one cannot use contracted Gaussians directly
 * (and we have to use primitive Gaussians instead).
 * However, the RPP radial function still can be used as a whole in the
 integrand.
 *
 * The Log3 integration scheme used here is detailed in:
 * C.-K. Skylaris et al. An efficient method for calculating effective core
 potential integrals
 * which involve projection operators.
 * Chem. Phys. Lett. 296, 445 (1998)
 */
#include <math.h>
#include <stdlib.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "grpp_radial_type1_integral.h"
#include "libgrpp.h"

#include "grpp_specfunc.h"
#include "grpp_utils.h"

#define MIN_GRID 2047
#define MAX_GRID 10000

typedef struct {
  double k;
  double alpha_A;
  double alpha_B;
  double CA_2;
  double CB_2;
  double prefactor;
  double (*potential)(double r, void *params);
  void *potential_params;
} radial_type1_params_t;

typedef struct {
  int nr;
  int n_max;
  int lambda_max;
  double *r;
  double *w;
  double *pot_values;
  double *gto_values;
  double **r_N;
  double **mod_bessel;
  radial_type1_params_t *params;
} radial_type1_grid_t;

static radial_type1_grid_t *
create_radial_type1_grid(int lambda_max, int n_max,
                         radial_type1_params_t *params);

static void expand_radial_type1_grid(radial_type1_grid_t *grid, int nr);

static void delete_radial_type1_grid(radial_type1_grid_t *grid);

static double calculate_radial_type1_integral(radial_type1_grid_t *grid, int n,
                                              int lambda, double tolerance,
                                              int *converged);

radial_type1_table_t *libgrpp_tabulate_radial_type1_integrals(
    int lambda_max, int n_max, double CA_2, double CB_2, double alpha_A,
    double alpha_B, double k, double prefactor,
    double (*potential)(double r, void *params), void *potential_params) {
  radial_type1_table_t *table;
  double const tolerance = libgrpp_params.radial_tolerance;

  table = (radial_type1_table_t *)calloc(1, sizeof(radial_type1_table_t));
  table->lambda_max = lambda_max;
  table->n_max = n_max;
  table->radial_integrals =
      (double *)calloc((lambda_max + 1) * (n_max + 1), sizeof(double));

  radial_type1_params_t params;
  params.CA_2 = CA_2;
  params.CB_2 = CB_2;
  params.alpha_A = alpha_A;
  params.alpha_B = alpha_B;
  params.k = k;
  params.prefactor = prefactor;
  params.potential = potential;
  params.potential_params = potential_params;

  radial_type1_grid_t *grid =
      create_radial_type1_grid(lambda_max, n_max, &params);

  for (int lambda = 0; lambda <= lambda_max; lambda++) {
    for (int n = 0; n <= n_max; n++) {

      int converged;
      double Q = calculate_radial_type1_integral(grid, n, lambda, tolerance,
                                                 &converged);

      table->radial_integrals[lambda * (lambda_max + 1) + n] = Q;
    }
  }

  delete_radial_type1_grid(grid);

  return table;
}

void libgrpp_delete_radial_type1_integrals(radial_type1_table_t *table) {
  free(table->radial_integrals);
  free(table);
}

double libgrpp_get_radial_type1_integral(radial_type1_table_t *table,
                                         int lambda, int n) {
  int lambda_max = table->lambda_max;
  return table->radial_integrals[lambda * (lambda_max + 1) + n];
}

static double radial_type1_integrand_fun(double r,
                                         radial_type1_params_t *params) {
  double alpha_A = params->alpha_A;
  double alpha_B = params->alpha_B;
  double k = params->k;
  double CA_2 = params->CA_2;
  double CB_2 = params->CB_2;
  double prefactor = params->prefactor;

  double power = k * r - (alpha_A + alpha_B) * r * r - alpha_A * CA_2 -
                 alpha_B * CB_2; // + N * log(r);

  return prefactor * exp(power);
}

static radial_type1_grid_t *
create_radial_type1_grid(int lambda_max, int n_max,
                         radial_type1_params_t *params) {
  radial_type1_grid_t *grid =
      (radial_type1_grid_t *)calloc(1, sizeof(radial_type1_grid_t));

  grid->nr = MIN_GRID;
  grid->n_max = n_max;
  grid->lambda_max = lambda_max;
  grid->params = params;

  grid->r = (double *)calloc(MAX_GRID, sizeof(double));
  grid->w = (double *)calloc(MAX_GRID, sizeof(double));
  grid->pot_values = (double *)calloc(MAX_GRID, sizeof(double));
  grid->gto_values = (double *)calloc(MAX_GRID, sizeof(double));
  grid->r_N = alloc_zeros_2d(n_max + 1, MAX_GRID);
  grid->mod_bessel = alloc_zeros_2d(lambda_max + 1, MAX_GRID);

  // initial set of pre-calculated points
  int nr = grid->nr;
  const double R = 5.0;
  const double R3 = R * R * R;

  for (int i = 1; i <= nr; i++) {
    double xi = i / (nr + 1.0);
    double xi3 = xi * xi * xi;
    double ln_xi = log(1 - xi3);
    double wi = 3 * R3 * xi * xi * ln_xi * ln_xi / ((1 - xi3) * (nr + 1.0));
    double ri = -R * ln_xi;

    grid->r[i - 1] = ri;
    grid->w[i - 1] = wi;
    grid->pot_values[i - 1] = params->potential(ri, params->potential_params);
    grid->gto_values[i - 1] = radial_type1_integrand_fun(ri, params);

    for (int lambda = 0; lambda <= lambda_max; lambda++) {
      grid->mod_bessel[lambda][i - 1] =
          libgrpp_modified_bessel_scaled(lambda, ri * params->k);
    }

    for (int n = 0; n <= n_max; n++) {
      grid->r_N[n][i - 1] = pow(ri, n);
    }
  }

  return grid;
}

static void delete_radial_type1_grid(radial_type1_grid_t *grid) {
  free(grid->r);
  free(grid->w);
  free(grid->pot_values);
  free(grid->gto_values);
  free_2d(grid->r_N, grid->n_max + 1);
  free_2d(grid->mod_bessel, grid->lambda_max + 1);
  free(grid);
}

static void expand_radial_type1_grid(radial_type1_grid_t *grid, int nr) {
  const double R = 5.0;
  const double R3 = R * R * R;

  if (nr > MAX_GRID) {
    return;
  }

  if (nr <= grid->nr) { // nothing to do
    return;
  }

  int idx = grid->nr;
  for (int i = 1; i <= nr; i += 2) {
    double xi = i / (nr + 1.0);
    double xi3 = xi * xi * xi;
    double ln_xi = log(1 - xi3);
    double wi = 3 * R3 * xi * xi * ln_xi * ln_xi / ((1 - xi3) * (nr + 1.0));
    double ri = -R * ln_xi;

    grid->r[idx] = ri;
    grid->w[idx] = wi;
    grid->pot_values[idx] =
        grid->params->potential(ri, grid->params->potential_params);
    grid->gto_values[idx] = radial_type1_integrand_fun(ri, grid->params);

    for (int lambda = 0; lambda <= grid->lambda_max; lambda++) {
      double kr = grid->params->k * ri;
      double bessel = libgrpp_modified_bessel_scaled(lambda, kr);
      grid->mod_bessel[lambda][idx] = bessel;
    }

    for (int n = 0; n <= grid->n_max; n++) {
      grid->r_N[n][idx] = pow(ri, n);
    }
    idx++;
  }

  grid->nr = nr;
}

static double calculate_radial_type1_integral(radial_type1_grid_t *grid, int n,
                                              int lambda, double tolerance,
                                              int *converged) {
  int nr = MIN_GRID;

  *converged = 0;
  double prev_sum = 0.0;
  double sum = 0.0;

  double *w = grid->w;
  // double *r = grid->r;
  double *pot_values = grid->pot_values;
  double *gto_values = grid->gto_values;
  double *r_N = grid->r_N[n];
  double *mod_bessel = grid->mod_bessel[lambda];

  /*
   * first step: screening of an integral
   */
  /*double screened = 0.0;
  int screened_success = screening_radial_type1(
          lambda,
          n,
          grid->params->CA_2,
          grid->params->CB_2,
          grid->params->alpha_A,
          grid->params->alpha_B,
          grid->params->k,
          grid->params->prefactor,
          grid->params->potential_params,
          &screened
  );

  if (screened_success == EXIT_SUCCESS && fabs(screened) < tolerance) {
      *converged = 1;
      return screened;
  }*/

  /*
   * second step: calculation on the smallest possible grid
   */
  for (int i = 0; i < nr; i++) {
    sum += w[i] * pot_values[i] * gto_values[i] * r_N[i] * mod_bessel[i];
  }

  /*
   * third step: adaptive integration, refinement of the result
   */
  do {
    int idx = nr;
    nr = 2 * nr + 1;

    if (nr > MAX_GRID) {
      break;
    }

    prev_sum = sum;
    sum = 0.5 * sum;

    expand_radial_type1_grid(grid, nr);

    for (int i = idx; i < nr; i++) {
      sum += w[i] * pot_values[i] * gto_values[i] * r_N[i] * mod_bessel[i];
    }

    /*if (screened_success == EXIT_SUCCESS && (fabs(sum) / fabs(screened) <
    0.001)) { *converged = 0; continue;
    }*/

    *converged = fabs(sum - prev_sum) <= tolerance;

  } while (!(*converged));

  return sum;
}
