/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2023 Alexander Oleynichenko
 */

#ifndef LIBGRPP_LMATRIX_H
#define LIBGRPP_LMATRIX_H

void construct_angular_momentum_matrices_rsh(int L, double *lx_matrix,
                                             double *ly_matrix,
                                             double *lz_matrix);

void construct_angular_momentum_matrices_csh(int L, double *lx_matrix,
                                             double *ly_matrix,
                                             double *lz_matrix);

#endif // LIBGRPP_LMATRIX_H
