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

#ifndef LIBGRPP_FACTORIAL_H
#define LIBGRPP_FACTORIAL_H

double libgrpp_factorial(int n);

double libgrpp_factorial_ratio(int n, int m);

double libgrpp_double_factorial(int n);

#endif // LIBGRPP_FACTORIAL_H
