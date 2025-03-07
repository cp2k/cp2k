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

#ifndef LIBGRPP_UTILS_H
#define LIBGRPP_UTILS_H

int int_max2(int x, int y);

int int_max3(int x, int y, int z);

double *alloc_zeros_1d(int n);

double **alloc_zeros_2d(int n, int m);

void free_2d(double **array, int n);

void libgrpp_daxpy(int n, double a, double *x, double *y);

void libgrpp_multiply_matrices(int M, int N, int K, double *A, double *B,
                               double *C);

double distance_squared(double *A, double *B);

double distance(double *A, double *B);

int points_are_equal(double *a, double *b);

#endif // LIBGRPP_UTILS_H
