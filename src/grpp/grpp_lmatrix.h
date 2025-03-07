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

#ifndef LIBGRPP_LMATRIX_H
#define LIBGRPP_LMATRIX_H

void libgrpp_construct_angular_momentum_matrices_rsh(int L, double *lx_matrix,
                                                     double *ly_matrix,
                                                     double *lz_matrix);

void libgrpp_construct_angular_momentum_matrices_csh(int L, double *lx_matrix,
                                                     double *ly_matrix,
                                                     double *lz_matrix);

#endif // LIBGRPP_LMATRIX_H
