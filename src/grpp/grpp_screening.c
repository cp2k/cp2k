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
 * Screening of radial integrals.
 *
 * The technique of screening is adopted from:
 *   R. A. Shaw, J. G. Hill. Prescreening and efficiency in the evaluation
 *   of integrals over ab initio effective core potentials.
 *   J. Chem. Phys. 147, 074108 (2017). doi: 10.1063/1.4986887
 * (see also Supplementary Material for this article).
 *
 * Note that in this publication the transcendental equation (2) for
 * type 2 integrals is not correct.
 */

#include <math.h>
#include <stdlib.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "grpp_factorial.h"
#include "grpp_screening.h"
#include "grpp_specfunc.h"
#include "libgrpp.h"

/*
 * functions defined below in the file
 */

static int screening_radial_type1_integral_primitive(
    int lambda, int n, double CA_2, double CB_2, double alpha_A, double alpha_B,
    double k, double eta, double *screened_value);

static double screening_type1_equation_for_maximum(double r, int n, int lambda,
                                                   double p, double k);

static int screening_radial_type2_integral_primitive(
    int lambda1, int lambda2, int n, double CA_2, double CB_2, double alpha_A,
    double alpha_B, double k1, double k2, double eta, double *screened_value);

static double screening_type2_equation_for_maximum(double r, int n, int lambda1,
                                                   int lambda2, double p,
                                                   double k1, double k2);

// static double analytic_one_center_rpp_integral_primitive(int L, double
// alpha1,
//                                                          double alpha2, int
//                                                          n, double zeta);

/**
 * screening for the type 1 radial integrals
 * for the pair of contracted gaussian functions.
 */
int libgrpp_screening_radial_type1(int lambda, int n, double CA_2, double CB_2,
                                   double alpha_A, double alpha_B, double k,
                                   double prefactor,
                                   libgrpp_potential_t *potential,
                                   double *screened_value) {
  *screened_value = 0.0;

  if (lambda >= 1 && fabs(k) <= LIBGRPP_ZERO_THRESH) {
    return EXIT_SUCCESS;
  }

  /*
   * loop over RPP primitives
   */
  for (int iprim = 0; iprim < potential->num_primitives; iprim++) {
    double eta = potential->alpha[iprim];
    int ni = n + potential->powers[iprim];
    double coef = potential->coeffs[iprim];

    double val_i = 0.0;
    int err_code = screening_radial_type1_integral_primitive(
        lambda, ni, CA_2, CB_2, alpha_A, alpha_B, k, eta, &val_i);
    if (err_code == EXIT_FAILURE) {
      return EXIT_FAILURE;
    }

    *screened_value += prefactor * coef * val_i;
  }

  return EXIT_SUCCESS;
}

/**
 * screening for the type 1 radial integrals
 * for the pair of primitive gaussian functions.
 */
static int screening_radial_type1_integral_primitive(
    int lambda, int n, double CA_2, double CB_2, double alpha_A, double alpha_B,
    double k, double eta, double *screened_value) {
  double p = alpha_A + alpha_B + eta;
  double CA = sqrt(CA_2);
  double CB = sqrt(CB_2);

  /*
   * find position of the maximum of the integrand
   */
  const double tol = 1e-2;
  double r0 = (alpha_A * CA + alpha_B * CB) / p;
  double r0_prev = 0.0;

  int nsteps = 0;
  do {
    nsteps++;
    if (nsteps == 10) {
      *screened_value = 0.0;
      return EXIT_FAILURE;
    }

    r0_prev = r0;
    r0 = screening_type1_equation_for_maximum(r0, n, lambda, p, k);
  } while (fabs(r0 - r0_prev) > tol);

  /*
   * envelope function for the integrand
   */
  *screened_value =
      sqrt(M_PI / p) * pow(r0, n) *
      libgrpp_modified_bessel_scaled(lambda, k * r0) *
      exp(-p * r0 * r0 - alpha_A * CA_2 - alpha_B * CB_2 + k * r0) * 0.5 *
      (1 + erf(sqrt(p) * r0));

  return EXIT_SUCCESS;
}

/**
 * ratio of two modified scaled Bessel functions guarded against divide by zero
 */
static double modified_bessel_scaled_ratio(int n, double x) {
  double numerator, denominator;
  if (n == 0) {
    numerator = libgrpp_modified_bessel_scaled(1, x);
    denominator = libgrpp_modified_bessel_scaled(0, x);
  } else {
    numerator = libgrpp_modified_bessel_scaled(n - 1, x);
    denominator = libgrpp_modified_bessel_scaled(n, x);
  }
  if (denominator == 0.0) {
    return (2 * n + 1) / x; // asymptote for x->0
  } else {
    return numerator / denominator;
  }
}

/**
 * transcendental equation for finding maximum of the type 1 integrand
 */
static double screening_type1_equation_for_maximum(double r, int n, int lambda,
                                                   double p, double k) {
  double k_r = k * r;
  double K_ratio = modified_bessel_scaled_ratio(lambda, k_r);
  double a = n + K_ratio * k_r;
  if (lambda > 0) {
    a = a - lambda - 1;
  }

  return sqrt(a / (2.0 * p));
}

/**
 * screening for the type 2 radial integrals
 * for the pair of contracted gaussian functions.
 */
int libgrpp_screening_radial_type2(int lambda1, int lambda2, int n, double CA_2,
                                   double CB_2, libgrpp_shell_t *bra,
                                   libgrpp_shell_t *ket,
                                   libgrpp_potential_t *potential,
                                   double *screened_value) {
  *screened_value = 0.0;

  double CA = sqrt(CA_2);
  double CB = sqrt(CB_2);

  /*
   * loop over 'bra' contracted function
   */
  for (int i = 0; i < bra->num_primitives; i++) {
    double alpha_A = bra->alpha[i];
    double coef_i = bra->coeffs[i];
    double k1 = 2 * alpha_A * CA;

    /*
     * loop over 'ket' contracted function
     */
    for (int j = 0; j < ket->num_primitives; j++) {
      double alpha_B = ket->alpha[j];
      double coef_j = ket->coeffs[j];
      double k2 = 2 * alpha_B * CB;

      /*
       * loop over RPP primitives
       */
      for (int k = 0; k < potential->num_primitives; k++) {
        double eta = potential->alpha[k];
        int ni = n + potential->powers[k];
        double coef_k = potential->coeffs[k];

        double val_ijk = 0.0;
        int err_code = screening_radial_type2_integral_primitive(
            lambda1, lambda2, ni, CA_2, CB_2, alpha_A, alpha_B, k1, k2, eta,
            &val_ijk);
        if (err_code == EXIT_FAILURE) {
          return EXIT_FAILURE;
        }

        *screened_value += coef_i * coef_j * coef_k * val_ijk;
      }
    }
  }

  return EXIT_SUCCESS;
}

/**
 * Analytically evaluates Gaussian integral:
 * \int_0^\infty r^n e^(-a r^2) dr
 */
double libgrpp_gaussian_integral(int n, double a) {
  if (n % 2 == 0) {
    int k = n / 2;
    return libgrpp_double_factorial(2 * k - 1) / (pow(2.0, k + 1) * pow(a, k)) *
           sqrt(M_PI / a);
  } else {
    int k = (n - 1) / 2;
    return libgrpp_factorial(k) / (2.0 * pow(a, k + 1));
  }
}

/**
 * screening for the type 2 radial integrals
 * for the pair of primitive gaussian functions.
 */
static int screening_radial_type2_integral_primitive(
    int lambda1, int lambda2, int n, double CA_2, double CB_2, double alpha_A,
    double alpha_B, double k1, double k2, double eta, double *screened_value) {
  *screened_value = 0.0;

  if (lambda1 >= 1 && fabs(k1) <= LIBGRPP_ZERO_THRESH) {
    return EXIT_SUCCESS;
  }
  if (lambda2 >= 1 && fabs(k2) <= LIBGRPP_ZERO_THRESH) {
    return EXIT_SUCCESS;
  }

  double p = alpha_A + alpha_B + eta;
  double CA = sqrt(CA_2);
  double CB = sqrt(CB_2);

  /*
   * special case:
   * lambda1 = lambda2 = 0,
   * k1 = k2 = 0.
   * => M_0(0) = 1
   * => we have one-center integral which can be evaluated analytically
   */
  if (lambda1 == 0 && lambda2 == 0) {
    if (fabs(k1) <= LIBGRPP_ZERO_THRESH && fabs(k2) <= LIBGRPP_ZERO_THRESH) {
      *screened_value = exp(-alpha_A * CA * CA - alpha_B * CB * CB) *
                        libgrpp_gaussian_integral(n, p);
      return EXIT_SUCCESS;
    }
  }

  /*
   * find position of the maximum of the integrand
   */
  const double tol = 1e-2;
  double r0 = (alpha_A * CA + alpha_B * CB) / p;
  double r0_prev = 0.0;
  int nsteps = 0;

  do {
    nsteps++;
    if (nsteps == 5) {
      *screened_value = 0.0;
      return EXIT_FAILURE;
    }
    r0_prev = r0;
    r0 = screening_type2_equation_for_maximum(r0, n, lambda1, lambda2, p, k1,
                                              k2);
  } while (fabs(r0 - r0_prev) > tol);

  /*
   * envelope function for the integrand
   */
  *screened_value = sqrt(M_PI / p) * pow(r0, n) *
                    libgrpp_modified_bessel_scaled(lambda1, k1 * r0) *
                    libgrpp_modified_bessel_scaled(lambda2, k2 * r0) *
                    exp(-eta * r0 * r0 - alpha_A * (r0 - CA) * (r0 - CA) -
                        alpha_B * (r0 - CB) * (r0 - CB)) *
                    0.5 * (1 + erf(sqrt(p) * r0));

  return EXIT_SUCCESS;
}

/**
 * transcendental equation for finding maximum of the type 2 integrand
 */
static double screening_type2_equation_for_maximum(double r, int n, int lambda1,
                                                   int lambda2, double p,
                                                   double k1, double k2) {
  double k1_r = k1 * r;
  double K1_ratio = modified_bessel_scaled_ratio(lambda1, k1_r);

  double k2_r = k2 * r;
  double K2_ratio = modified_bessel_scaled_ratio(lambda2, k2_r);

  double a = K1_ratio * k1_r + K2_ratio * k2_r + n;

  if (lambda1 > 0) {
    a = a - lambda1 - 1;
  }
  if (lambda2 > 0) {
    a = a - lambda2 - 1;
  }

  return sqrt(a / (2.0 * p));
}
