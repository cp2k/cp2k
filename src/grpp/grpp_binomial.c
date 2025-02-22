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

#include "grpp_binomial.h"

#include <stdint.h>

/* The code is borrowed from RosettaCode:
 * https://rosettacode.org/wiki/Evaluate_binomial_coefficients#C
 * We go to some effort to handle overflow situations.
 */

static uint64_t gcd_ui(uint64_t x, uint64_t y);

/*
 * returns binomial coefficient:
 * ( n )
 * ( k )
 */
uint64_t libgrpp_binomial(uint64_t n, uint64_t k) {
  uint64_t d, g, r = 1;

  if (k == 0) {
    return 1;
  }
  if (k == 1) {
    return n;
  }
  if (k >= n) {
    return (k == n);
  }
  if (k > n / 2) {
    k = n - k;
  }
  for (d = 1; d <= k; d++) {
    if (r >= UINT64_MAX / n) { /* Possible overflow */
      uint64_t nr, dr;         /* reduced numerator / denominator */
      g = gcd_ui(n, d);
      nr = n / g;
      dr = d / g;
      g = gcd_ui(r, dr);
      r = r / g;
      dr = dr / g;
      if (r >= UINT64_MAX / nr)
        return 0; /* Unavoidable overflow */
      r *= nr;
      r /= dr;
      n--;
    } else {
      r *= n--;
      r /= d;
    }
  }
  return r;
}

static uint64_t gcd_ui(uint64_t x, uint64_t y) {
  uint64_t t;

  if (y < x) {
    t = x;
    x = y;
    y = t;
  }
  while (y > 0) {
    t = y;
    y = x % y;
    x = t; /* y1 <- x0 % y0 ; x1 <- y0 */
  }
  return x;
}
