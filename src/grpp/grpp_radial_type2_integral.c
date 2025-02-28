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
 * Evaluation of type 2 radial integrals.
 *
 * The procedure in general follows that described in:
 * R. Flores-Moreno et al. Half-numerical evaluation of pseudopotential
 * integrals. J. Comp. Chem. 27, 1009 (2006) (see formulas (12) and (13) for
 * radial integrals invoking contracted Gaussian functions and RPPs)
 *
 * The Log3 integration scheme used here is detailed in:
 * C.-K. Skylaris et al. An efficient method for calculating effective core
 * potential integrals which involve projection operators. Chem. Phys. Lett.
 * 296, 445 (1998)
 */

#include <math.h>
#include <stdlib.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "grpp_radial_type2_integral.h"

#include "grpp_norm_gaussian.h"
#include "grpp_screening.h"
#include "grpp_specfunc.h"
#include "grpp_utils.h"

#define MIN_GRID 31
#define MAX_GRID 10000

typedef struct {
  double CA;
  double CB;
  libgrpp_potential_t *potential;
  libgrpp_shell_t *bra;
  libgrpp_shell_t *ket;
} radial_type2_params_t;

/**
 * RPP and radial contracted Gaussians are pre-calculated on a grid,
 * and then combined into radial integrals
 */
typedef struct {
  int nr;
  int n_max;
  int lambda1_max;
  int lambda2_max;
  double *r;
  double *w;
  double *rpp_values;
  double **r_N;
  double **F1;
  double **F2;
  radial_type2_params_t *params;
} radial_type2_grid_t;

/**
 * pre-definitions of the functions used below
 */

static double calculate_radial_type2_integral(radial_type2_grid_t *grid, int n,
                                              int lambda1, int lambda2,
                                              double tolerance, int *converged);

static radial_type2_grid_t *
create_radial_type2_grid(int lambda1_max, int lambda2_max, int n_max,
                         radial_type2_params_t *params);

static void delete_radial_type2_grid(radial_type2_grid_t *grid);

static double radial_type2_integrand_fun_contracted(double r, int lambda,
                                                    double *k, double CA,
                                                    libgrpp_shell_t *gauss_fun);

static void expand_radial_type2_grid(radial_type2_grid_t *grid, int nr);

static void calc_k_values(int nprim, const double *alpha, double CA, double *k);

double libgrpp_gaussian_integral(int n, double a);

/**
 * Creates table with pre-calculated radial type 2 integrals.
 */
radial_type2_table_t *libgrpp_tabulate_radial_type2_integrals(
    int lambda1_max, int lambda2_max, int n_max, double CA_2, double CB_2,
    libgrpp_potential_t *potential, libgrpp_shell_t *bra,
    libgrpp_shell_t *ket) {
  /*
   * create empty table containing pre-tabulated radial type 2 integrals
   */
  radial_type2_table_t *table;
  table = (radial_type2_table_t *)calloc(1, sizeof(radial_type2_table_t));
  table->lambda1_max = lambda1_max;
  table->lambda2_max = lambda2_max;
  table->n_max = n_max;
  table->radial_integrals = (double *)calloc(
      (lambda1_max + 1) * (lambda2_max + 1) * (n_max + 1), sizeof(double));

  /*
   * the special case of one-center RPP integrals
   */
  if (CA_2 < 1e-14 && CB_2 < 1e-14) {

    for (int i = 0; i < bra->num_primitives; i++) {
      double alpha_A = bra->alpha[i];
      double coef_i =
          bra->coeffs[i] * libgrpp_gaussian_norm_factor(bra->L, 0, 0, alpha_A);

      for (int j = 0; j < ket->num_primitives; j++) {
        double alpha_B = ket->alpha[j];
        double coef_j = ket->coeffs[j] *
                        libgrpp_gaussian_norm_factor(ket->L, 0, 0, alpha_B);

        for (int k = 0; k < potential->num_primitives; k++) {
          double eta = potential->alpha[k];
          int ni = potential->powers[k];
          double coef_k = potential->coeffs[k];

          double p = alpha_A + alpha_B + eta;
          double factor = coef_i * coef_j * coef_k;

          for (int n = 0; n <= n_max; n++) {

            double val_ijk = libgrpp_gaussian_integral(ni + n, p);
            table->radial_integrals[n] += factor * val_ijk;
            ;
          }
        }
      }
    }

    return table;
  }

  /*
   * for numerical integration on the grid
   */
  radial_type2_params_t params;
  params.CA = sqrt(CA_2);
  params.CB = sqrt(CB_2);
  params.potential = libgrpp_shrink_potential(potential);

  params.bra = libgrpp_shell_deep_copy(bra);
  libgrpp_shell_shrink(params.bra);
  libgrpp_shell_mult_normcoef(params.bra);

  params.ket = libgrpp_shell_deep_copy(ket);
  libgrpp_shell_shrink(params.ket);
  libgrpp_shell_mult_normcoef(params.ket);

  /*
   * create radial grid
   */
  radial_type2_grid_t *grid =
      create_radial_type2_grid(lambda1_max, lambda2_max, n_max, &params);

  /*
   * calculate radial integrals and store them into the table
   */
  for (int lambda_1 = 0; lambda_1 <= lambda1_max; lambda_1++) {
    for (int lambda_2 = 0; lambda_2 <= lambda2_max; lambda_2++) {
      for (int n = 0; n <= n_max; n++) {

        int converged;
        double Q = calculate_radial_type2_integral(grid, n, lambda_1, lambda_2,
                                                   1e-16, &converged);

        //        int dim1 = lambda1_max + 1;
        int dim2 = lambda2_max + 1;
        int dimn = n_max + 1;
        table->radial_integrals[dim2 * dimn * lambda_1 + dimn * lambda_2 + n] =
            Q;
      }
    }
  }

  /*
   * clean-up
   */
  libgrpp_delete_potential(params.potential);
  libgrpp_delete_shell(params.bra);
  libgrpp_delete_shell(params.ket);
  delete_radial_type2_grid(grid);

  return table;
}

/**
 * destructor for the table of radial type 2 integrals
 */
void libgrpp_delete_radial_type2_integrals(radial_type2_table_t *table) {
  free(table->radial_integrals);
  free(table);
}

/**
 * Returns radial integral at complex index (lambda1,lambda2,N)
 */
double libgrpp_get_radial_type2_integral(radial_type2_table_t *table,
                                         int lambda1, int lambda2, int n) {
  // int lambda1_max = table->lambda1_max;
  int lambda2_max = table->lambda2_max;
  int n_max = table->n_max;
  // int dim1 = lambda1_max + 1;
  int dim2 = lambda2_max + 1;
  int dimn = n_max + 1;

  double Q =
      table->radial_integrals[dim2 * dimn * lambda1 + dimn * lambda2 + n];

  return Q;
}

/**
 * calculates type 2 radial integral T^N_{lambda1,lambda2}
 * for the two given contracted gaussian functions and the contracted potential
 */
static double calculate_radial_type2_integral(radial_type2_grid_t *grid, int n,
                                              int lambda1, int lambda2,
                                              double tolerance,
                                              int *converged) {
  int nr = MIN_GRID;

  *converged = 0;
  double prev_sum = 0.0;
  double sum = 0.0;

  double *w = grid->w;
  // double *r = grid->r;
  double *pot_values = grid->rpp_values;
  double *F1 = grid->F1[lambda1];
  double *F2 = grid->F2[lambda2];
  double *r_N = grid->r_N[n];

  /*
   * first step: integral screening
   */
  double CA = grid->params->CA;
  double CB = grid->params->CB;

  double screened = 0.0;
  int screen_success = libgrpp_screening_radial_type2(
      lambda1, lambda2, n, CA * CA, CB * CB, grid->params->bra,
      grid->params->ket, grid->params->potential, &screened);

  if (screen_success == EXIT_SUCCESS && fabs(screened) < tolerance) {
    *converged = 1;
    return screened;
  }

  /*
   * second step: calculation on the smallest possible grid
   */
  for (int i = 0; i < nr; i++) {
    sum += w[i] * pot_values[i] * F1[i] * F2[i] * r_N[i];
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

    expand_radial_type2_grid(grid, nr);

    for (int i = idx; i < nr; i++) {
      sum += w[i] * pot_values[i] * F1[i] * F2[i] * r_N[i];
    }

    if (screen_success == EXIT_SUCCESS &&
        (fabs(sum) / fabs(screened) < 0.001)) {
      *converged = 0;
      continue;
    }

    *converged = fabs(sum - prev_sum) <= tolerance;

  } while (!(*converged));

  return sum;
}

/**
 * Numerical integration on the Log3 grid
 */
static radial_type2_grid_t *
create_radial_type2_grid(int lambda1_max, int lambda2_max, int n_max,
                         radial_type2_params_t *params) {
  radial_type2_grid_t *grid =
      (radial_type2_grid_t *)calloc(1, sizeof(radial_type2_grid_t));

  grid->nr = MIN_GRID;
  grid->n_max = n_max;
  grid->lambda1_max = lambda1_max;
  grid->lambda2_max = lambda2_max;
  grid->params = params;

  grid->r = alloc_zeros_1d(MAX_GRID);
  grid->w = alloc_zeros_1d(MAX_GRID);
  grid->rpp_values = alloc_zeros_1d(MAX_GRID);
  grid->F1 = alloc_zeros_2d(lambda1_max + 1, MAX_GRID);
  grid->F2 = alloc_zeros_2d(lambda2_max + 1, MAX_GRID);
  grid->r_N = alloc_zeros_2d(n_max + 1, MAX_GRID);

  // vectors 'k': k = - 2 * alpha * |CA|
  double *bra_k = alloc_zeros_1d(params->bra->num_primitives);
  double *ket_k = alloc_zeros_1d(params->ket->num_primitives);
  calc_k_values(params->bra->num_primitives, params->bra->alpha, params->CA,
                bra_k);
  calc_k_values(params->ket->num_primitives, params->ket->alpha, params->CB,
                ket_k);

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
    grid->rpp_values[i - 1] =
        libgrpp_potential_value(grid->params->potential, ri);
    for (int n = 0; n <= n_max; n++) {
      grid->r_N[n][i - 1] = pow(ri, n);
    }
    for (int lambda1 = 0; lambda1 <= lambda1_max; lambda1++) {
      grid->F1[lambda1][i - 1] = radial_type2_integrand_fun_contracted(
          ri, lambda1, bra_k, params->CA, params->bra);
    }
    for (int lambda2 = 0; lambda2 <= lambda2_max; lambda2++) {
      grid->F2[lambda2][i - 1] = radial_type2_integrand_fun_contracted(
          ri, lambda2, ket_k, params->CB, params->ket);
    }
  }

  free(bra_k);
  free(ket_k);

  return grid;
}

/**
 * constructs new radial grid points
 */
static void expand_radial_type2_grid(radial_type2_grid_t *grid, int nr) {
  const double R = 5.0;
  const double R3 = R * R * R;

  if (nr > MAX_GRID) {
    return;
  }

  if (nr <= grid->nr) { // nothing to do
    return;
  }

  radial_type2_params_t *params = grid->params;

  // vectors 'k': k = - 2 * alpha * |CA|
  double *bra_k = alloc_zeros_1d(params->bra->num_primitives);
  double *ket_k = alloc_zeros_1d(params->ket->num_primitives);
  calc_k_values(params->bra->num_primitives, params->bra->alpha, params->CA,
                bra_k);
  calc_k_values(params->ket->num_primitives, params->ket->alpha, params->CB,
                ket_k);

  // additional set of grid points
  int idx = grid->nr;
  for (int i = 1; i <= nr; i += 2) {
    double xi = i / (nr + 1.0);
    double xi3 = xi * xi * xi;
    double ln_xi = log(1 - xi3);
    double wi = 3 * R3 * xi * xi * ln_xi * ln_xi / ((1 - xi3) * (nr + 1.0));
    double ri = -R * ln_xi;

    grid->r[idx] = ri;
    grid->w[idx] = wi;
    grid->rpp_values[idx] =
        libgrpp_potential_value(grid->params->potential, ri);

    for (int n = 0; n <= grid->n_max; n++) {
      grid->r_N[n][idx] = pow(ri, n);
    }

    for (int lambda1 = 0; lambda1 <= grid->lambda1_max; lambda1++) {
      grid->F1[lambda1][idx] = radial_type2_integrand_fun_contracted(
          ri, lambda1, bra_k, grid->params->CA, params->bra);
    }
    for (int lambda2 = 0; lambda2 <= grid->lambda2_max; lambda2++) {
      grid->F2[lambda2][idx] = radial_type2_integrand_fun_contracted(
          ri, lambda2, ket_k, grid->params->CB, params->ket);
    }

    idx++;
  }

  grid->nr = nr;

  free(bra_k);
  free(ket_k);
}

/**
 * deallocates memory used for the radial grid
 */
static void delete_radial_type2_grid(radial_type2_grid_t *grid) {
  free(grid->r);
  free(grid->w);
  free(grid->rpp_values);
  free_2d(grid->F1, grid->lambda1_max + 1);
  free_2d(grid->F2, grid->lambda2_max + 1);
  free_2d(grid->r_N, grid->n_max + 1);
  free(grid);
}

/**
 * Calculate the value of the integrand function
 */
static double
radial_type2_integrand_fun_contracted(double r, int lambda, double *k,
                                      double CA, libgrpp_shell_t *gauss_fun) {
  double F = 0.0;
  double r_CA_2 = (r - CA) * (r - CA);

  int nprim = gauss_fun->num_primitives;
  double *alpha = gauss_fun->alpha;
  double *coeffs = gauss_fun->coeffs;

  for (int i = 0; i < nprim; i++) {
    double power = -alpha[i] * r_CA_2;
    F += coeffs[i] * exp(power) *
         libgrpp_modified_bessel_scaled(lambda, k[i] * r);
  }

  return F;
}

static void calc_k_values(int nprim, const double *alpha, double CA,
                          double *k) {
  for (int i = 0; i < nprim; i++) {
    k[i] = 2.0 * alpha[i] * CA;
  }
}
