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
#include <assert.h>

#include "grpp_parameters.h"

static int cartesian_generator_dirac(int L, int *cart_list);

static int cartesian_generator_turbomole(int L, int *cart_list);

libgrpp_parameters_t libgrpp_params = {
    // tolerance of radial integration
    1e-16,
    // tolerance of angular integral screening
    1e-16,
    // tolerance of modified Bessel functions evaluation
    // is it really needed?
    1e-16,
    // subroutine to generate cartesian components with given ang momentum L
    cartesian_generator_dirac};

void libgrpp_set_default_parameters() {
  // #pragma omp critical
  {
    libgrpp_set_radial_tolerance(1e-16);
    libgrpp_set_angular_screening_tolerance(1e-16);
    libgrpp_set_modified_bessel_tolerance(1e-16);
    libgrpp_set_cartesian_generator(cartesian_generator_dirac);
  }
}

void libgrpp_set_radial_tolerance(double tolerance) {
  // #pragma omp critical
  { libgrpp_params.radial_tolerance = tolerance; }
}

void libgrpp_set_angular_screening_tolerance(double tolerance) {
  // #pragma omp critical
  { libgrpp_params.angular_screening_tolerance = tolerance; }
}

void libgrpp_set_modified_bessel_tolerance(double tolerance) {
  // #pragma omp critical
  { libgrpp_params.modified_bessel_tolerance = tolerance; }
}

void libgrpp_set_cartesian_order(int order) {
  // #pragma omp critical
  {
    assert(order == LIBGRPP_CART_ORDER_DIRAC ||
           order == LIBGRPP_CART_ORDER_TURBOMOLE);

    if (order == LIBGRPP_CART_ORDER_DIRAC) {
      libgrpp_set_cartesian_generator(cartesian_generator_dirac);
    } else if (order == LIBGRPP_CART_ORDER_TURBOMOLE) {
      libgrpp_set_cartesian_generator(cartesian_generator_turbomole);
    }
  }
}

void libgrpp_set_cartesian_generator(
    int (*cartesian_generator)(int L, int *cart_list)) {
  // #pragma omp critical
  { libgrpp_params.cartesian_generator = cartesian_generator; }
}

static int cartesian_generator_dirac(int L, int *cart_list) {
  int count = 0;
  int n_cart = (L + 1) * (L + 2) / 2;

  for (int r = L; r >= 0; r--) {
    for (int s = L; s >= 0; s--) {
      for (int t = L; t >= 0; t--) {
        if (r + s + t == L) {
          cart_list[3 * count + 0] = r;
          cart_list[3 * count + 1] = s;
          cart_list[3 * count + 2] = t;
          count++;
        }
      }
    }
  }

  return n_cart;
}

static int cartesian_generator_turbomole(int L, int *cart_list) {
  int count = 0;
  int n_cart = (L + 1) * (L + 2) / 2;

  for (int r = L; r >= 0; r--) {
    for (int s = L - r; s >= 0; s--) {
      int t = L - r - s;

      cart_list[3 * count + 0] = r;
      cart_list[3 * count + 1] = s;
      cart_list[3 * count + 2] = t;

      count++;
    }
  }

  return n_cart;
}
