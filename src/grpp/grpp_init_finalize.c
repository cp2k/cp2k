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

#include "libgrpp.h"

void libgrpp_create_real_spherical_harmonic_coeffs_tables(int Lmax);
void libgrpp_pretabulate_bessel();

static int libgrpp_initialized = 0;

/**
 * thread-safe initialization
 */
void libgrpp_init() {
#pragma omp critical
  {
    if (libgrpp_initialized == 0) {
      libgrpp_create_real_spherical_harmonic_coeffs_tables(40);
      libgrpp_pretabulate_bessel();

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
