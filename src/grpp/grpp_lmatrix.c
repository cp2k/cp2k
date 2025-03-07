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
 * Functions for construction of matrices of the angular momentum operator L
 * in the bases of either real or complex spherical harmonics.
 */

#include "grpp_lmatrix.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static void get_transformation_coeffs_csh_to_rsh(int m, double complex *a,
                                                 double complex *b);

/**
 * Constructs matrices of the Lx, Ly, Lz operators for the given angular
 * momentum L in the basis of real spherical harmonics (rsh).
 */
void libgrpp_construct_angular_momentum_matrices_rsh(int L, double *lx_matrix,
                                                     double *ly_matrix,
                                                     double *lz_matrix) {
  int dim = 2 * L + 1;

  // set all matrices to zero
  memset(lx_matrix, 0, dim * dim * sizeof(double));
  memset(ly_matrix, 0, dim * dim * sizeof(double));
  memset(lz_matrix, 0, dim * dim * sizeof(double));

  double *lx_matrix_csh = calloc(dim * dim, sizeof(double));
  double *ly_matrix_csh = calloc(dim * dim, sizeof(double));
  double *lz_matrix_csh = calloc(dim * dim, sizeof(double));

  libgrpp_construct_angular_momentum_matrices_csh(L, lx_matrix_csh,
                                                  ly_matrix_csh, lz_matrix_csh);

  for (int m1 = -L; m1 <= L; m1++) {
    for (int m2 = -L; m2 <= L; m2++) {

      // coefficients: S_lm = a * Y_{l,-m} + b * Y_{l,m}
      double complex a1, b1;
      double complex a2, b2;
      get_transformation_coeffs_csh_to_rsh(m1, &a1, &b1); // bra
      get_transformation_coeffs_csh_to_rsh(m2, &a2, &b2); // ket

      int m1m = -abs(m1);
      int m1p = +abs(m1);
      int m2m = -abs(m2);
      int m2p = +abs(m2);

      double complex lx = 0.0 + 0.0 * I;
      lx += conj(a1) * a2 * lx_matrix_csh[dim * (m1m + L) + (m2m + L)];
      lx += conj(a1) * b2 * lx_matrix_csh[dim * (m1m + L) + (m2p + L)];
      lx += conj(b1) * a2 * lx_matrix_csh[dim * (m1p + L) + (m2m + L)];
      lx += conj(b1) * b2 * lx_matrix_csh[dim * (m1p + L) + (m2p + L)];

      double complex ly = 0.0 + 0.0 * I;
      ly += conj(a1) * a2 * ly_matrix_csh[dim * (m1m + L) + (m2m + L)];
      ly += conj(a1) * b2 * ly_matrix_csh[dim * (m1m + L) + (m2p + L)];
      ly += conj(b1) * a2 * ly_matrix_csh[dim * (m1p + L) + (m2m + L)];
      ly += conj(b1) * b2 * ly_matrix_csh[dim * (m1p + L) + (m2p + L)];

      double complex lz = 0.0 + 0.0 * I;
      lz += conj(a1) * a2 * lz_matrix_csh[dim * (m1m + L) + (m2m + L)];
      lz += conj(a1) * b2 * lz_matrix_csh[dim * (m1m + L) + (m2p + L)];
      lz += conj(b1) * a2 * lz_matrix_csh[dim * (m1p + L) + (m2m + L)];
      lz += conj(b1) * b2 * lz_matrix_csh[dim * (m1p + L) + (m2p + L)];

      lx_matrix[(m1 + L) * dim + (m2 + L)] = cimag(lx);
      ly_matrix[(m1 + L) * dim + (m2 + L)] = -creal(ly);
      lz_matrix[(m1 + L) * dim + (m2 + L)] = cimag(lz);
    }
  }

  free(lx_matrix_csh);
  free(ly_matrix_csh);
  free(lz_matrix_csh);
}

/**
 * Constructs matrices of the Lx, Ly, Lz operators in the basis of
 * complex spherical harmonics (csh) |Y_lm> for angular momentum l=L.
 * Matrices of size (2*L+1) x (2*L+1) must be pre-allocated.
 */
void libgrpp_construct_angular_momentum_matrices_csh(int L, double *lx_matrix,
                                                     double *ly_matrix,
                                                     double *lz_matrix) {
  int dim = 2 * L + 1;

  // set all matrices to zero
  memset(lx_matrix, 0, dim * dim * sizeof(double));
  memset(ly_matrix, 0, dim * dim * sizeof(double));
  memset(lz_matrix, 0, dim * dim * sizeof(double));

  for (int m1 = -L; m1 <= L; m1++) {
    for (int m2 = -L; m2 <= L; m2++) {

      double lz = m2 * (m1 == m2);
      double lp = sqrt((L - m2) * (L + m2 + 1)) * (m1 == m2 + 1); // L+
      double lm = sqrt((L + m2) * (L - m2 + 1)) * (m1 == m2 - 1); // L-
      double lx = 0.5 * (lp + lm);
      double ly = 0.5 * (lp - lm);

      lx_matrix[(m1 + L) * dim + (m2 + L)] = lx;
      ly_matrix[(m1 + L) * dim + (m2 + L)] = ly;
      lz_matrix[(m1 + L) * dim + (m2 + L)] = lz;
    }
  }
}

/**
 * Real spherical harmonic S_{l,m} can be represented as a linear combination
 * of two complex spherical harmonics:
 * S_{l,m} = a * Y_{l,-m} + b * Y_{l,m}
 * (except for the case m=0, where S_{l,0} = Y_{l,0})
 *
 * coefficients can be found elsewhere, see, for example,
 * https://en.wikipedia.org/wiki/Table_of_spherical_harmonics
 */
static void get_transformation_coeffs_csh_to_rsh(int m, double complex *a,
                                                 double complex *b) {
  if (m == 0) {
    *a = 0.5;
    *b = 0.5;
  } else if (m < 0) {
    *a = +1.0 * I / sqrt(2);
    *b = -1.0 * I / sqrt(2) * pow(-1, abs(m));
  } else { // m > 0
    *a = +1.0 / sqrt(2);
    *b = +1.0 / sqrt(2) * pow(-1, m);
  }
}
