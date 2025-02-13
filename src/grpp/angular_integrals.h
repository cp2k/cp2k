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

#ifndef LIBGRPP_ANGULAR_INTEGRALS_H
#define LIBGRPP_ANGULAR_INTEGRALS_H

double libgrpp_angular_type1_integral(int lambda, int II, int JJ, int KK,
                                      double *k);

double libgrpp_angular_type2_integral(int lambda, int L, int m, int a, int b,
                                      int c, const double *rsh_values);

#endif // LIBGRPP_ANGULAR_INTEGRALS_H
