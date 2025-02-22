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

#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "grpp_norm_gaussian.h"

/**
 * Calculates normalization factor for the cartesian Gaussian x^n y^l z^m
 * exp(-alpha*r^2)
 */
double libgrpp_gaussian_norm_factor(int n, int l, int m, double alpha) {
  return pow(2 * alpha / M_PI, 0.75) *
         pow(4 * alpha, 0.5 * (n + l + m)); /* /
          sqrt((double) double_factorial(2 * n - 1) *
               (double) double_factorial(2 * l - 1) *
               (double) double_factorial(2 * m - 1));*/
}
