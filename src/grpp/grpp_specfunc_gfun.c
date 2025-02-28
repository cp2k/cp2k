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
 * Implementation of the Gn(x) auxiliary function.
 * This function is required to calculate integrals over the 1/r^2 operator.
 *
 * More on the Gn(x) function evaluation:
 * (1) J. 0. Jensen, A. H. Cameri, C. P. Vlahacos, D. Zeroka, H. F. Hameka, C.
 * N. Merrow, Evaluation of one-electron integrals for arbitrary operators V(r)
 * over cartesian Gaussians: Application to inverse-square distance and Yukawa
 * operators. J. Comput. Chem. 14(8), 986 (1993). doi: 10.1002/jcc.540140814 (2)
 * B. Gao, A. J. Thorvaldsen, K. Ruud, GEN1INT: A unified procedure for the
 * evaluation of one-electron integrals over Gaussian basis functions and their
 * geometric derivatives. Int. J. Quantum Chem. 111(4), 858 (2011).
 *     doi: 10.1002/qua.22886
 */
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433
#endif

#include "grpp_factorial.h"
#include "grpp_specfunc.h"

static double gfun_taylor(int n, double x);

/**
 * Calculates values of the Gn(x) auxiliary function for n = 0, ..., nmax
 * and stores them into the g[] array.
 */
void libgrpp_gfun_values(double x, int nmax, double *g) {
  memset(g, 0, (nmax + 1) * sizeof(double));

  if (x <= 12.0) {
    /*
     * downward recursion
     */
    g[nmax] = gfun_taylor(nmax, x);

    for (int n = nmax; n > 0; n--) {
      g[n - 1] = (1.0 - 2.0 * x * g[n]) / (2.0 * n - 1.0);
    }
  } else {
    /*
     * upward recursion
     */
    double sqrt_x = sqrt(x);
    g[0] = libgrpp_Dawsons_Integral(sqrt_x) / sqrt_x;

    for (int n = 0; n < nmax; n++) {
      g[n + 1] = (1.0 - (2 * n + 1) * g[n]) / (2.0 * x);
    }
  }
}

/**
 * Calculates value of the Gn(x) auxiliary function using the Taylor expansion.
 * The Taylor series converges for x <= 30.
 */
static double gfun_taylor(int n, double x) {
  const double thresh = 1e-15;
  double sum = 0.0;

  for (int k = 0; k < 100; k++) {

    double y_exp = exp(-x);
    double y_pow = pow(x, k);
    double y_fac = libgrpp_factorial(k);
    double y_nk1 = 2.0 * n + 2.0 * k + 1.0;

    double contrib = y_exp * y_pow / y_fac / y_nk1;
    sum += contrib;

    if (fabs(contrib) < thresh) {
      break;
    }
  }

  return sum;
}
