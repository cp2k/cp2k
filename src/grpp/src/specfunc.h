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

double modified_bessel_scaled(int n, double x);

void gfun_values(double x, int nmax, double *g);

double boys_function(int n, double x);

void boys_values(double x, int nmax, double *b);

double specfunc_fermi_sk(int k, double x);

double Dawsons_Integral(double x);

#endif // LIBGRPP_SPECFUNC_H
