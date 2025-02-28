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
 * Interface to the MyMathLib library:
 * http://www.mymathlib.com/
 */

#ifndef LIBGRPP_SPECFUNC_H
#define LIBGRPP_SPECFUNC_H

double libgrpp_modified_bessel_scaled(int n, double x);

void libgrpp_gfun_values(double x, int nmax, double *g);

double libgrpp_boys_function(int n, double x);

void libgrpp_boys_values(double x, int nmax, double *b);

double libgrpp_specfunc_fermi_sk(int k, double x);

double libgrpp_Dawsons_Integral(double x);

#endif // LIBGRPP_SPECFUNC_H
