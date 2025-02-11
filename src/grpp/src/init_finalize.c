/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2023 Alexander Oleynichenko
 */

#include "libgrpp.h"

void create_real_spherical_harmonic_coeffs_tables(int Lmax);
void pretabulate_bessel();

static int libgrpp_initialized = 0;

/**
 * thread-safe initialization
 */
void libgrpp_init() {
#pragma omp critical
  {
    if (libgrpp_initialized == 0) {
      create_real_spherical_harmonic_coeffs_tables(40);
      pretabulate_bessel();

      libgrpp_initialized = 1;
    }
  }
}

int libgrpp_is_initialized() { return libgrpp_initialized; }

/**
 * thread-safe finalization
 */
void libgrpp_finalize() {
#pragma omp critical
  {
    if (libgrpp_initialized == 1) {
      libgrpp_initialized = 0;
    }
  }
}
