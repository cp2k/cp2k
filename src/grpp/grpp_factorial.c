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

#include "grpp_factorial.h"

#include <assert.h>
#include <stdint.h>

static uint64_t pretabulated_factorials[] = {1,
                                             1,
                                             2,
                                             6,
                                             24,
                                             120,
                                             720,
                                             5040,
                                             40320,
                                             362880,
                                             3628800,
                                             39916800,
                                             479001600,
                                             6227020800,
                                             87178291200,
                                             1307674368000,
                                             20922789888000,
                                             355687428096000,
                                             6402373705728000,
                                             121645100408832000,
                                             2432902008176640000};

double libgrpp_factorial(int n) {
  if (n < 0) {
    return 1;
  } else if (n <= 20) {
    return (double)pretabulated_factorials[n];
  } else {
    return n * libgrpp_factorial(n - 1);
  }
}

/*
 * Calculates ratio of two factorials:
 *   n!
 *  ----
 *   m!
 */
double libgrpp_factorial_ratio(int n, int m) {
  if (n == m) {
    return 1.0;
  }
  if (n < m) {
    return 1.0 / libgrpp_factorial_ratio(m, n);
  } else { // n > m
    double prod = 1.0;
    for (int i = m + 1; i <= n; i++) {
      prod *= i;
    }
    return prod;
  }
}

static uint64_t pretabulated_double_factorials[] = {
    1, // 0!!
    1,
    2,
    3,
    8,
    15, // 5!!
    48,
    105,
    384,
    945,
    3840, // 10!!
    10395,
    46080,
    135135,
    645120,
    2027025, // 15!!
    10321920,
    34459425,
    185794560,
    654729075,
    3715891200, // 20!!
    13749310575,
    81749606400,
    316234143225,
    1961990553600,
    7905853580625, // 25!!
    51011754393600,
    213458046676875,
    1428329123020800,
    6190283353629375,
    42849873690624000 // 30!!
};

double libgrpp_double_factorial(int n) {
  assert(n >= -1 && n <= 30);
  if (n == -1) {
    return 1;
  } else if (n <= 30) {
    return (double)pretabulated_double_factorials[n];
  } else {
    return n * libgrpp_double_factorial(n - 2);
  }
}
