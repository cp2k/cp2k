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

#ifndef LIBGRPP_SCREENING_H
#define LIBGRPP_SCREENING_H

#include "libgrpp.h"

int libgrpp_screening_radial_type1(int lambda, int n, double CA_2, double CB_2,
                                   double alpha_A, double alpha_B, double k,
                                   double prefactor,
                                   libgrpp_potential_t *potential,
                                   double *screened_value);

int libgrpp_screening_radial_type2(int lambda1, int lambda2, int n, double CA_2,
                                   double CB_2, libgrpp_shell_t *bra,
                                   libgrpp_shell_t *ket,
                                   libgrpp_potential_t *potential,
                                   double *screened_value);

#endif // LIBGRPP_SCREENING_H
