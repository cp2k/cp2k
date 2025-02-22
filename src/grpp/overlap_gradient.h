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

#ifndef LIBGRPP_OVERLAP_GRADIENT_H
#define LIBGRPP_OVERLAP_GRADIENT_H

#include "libgrpp.h"

void libgrpp_overlap_integrals_gradient(libgrpp_shell_t *shell_A,
                                        libgrpp_shell_t *shell_B,
                                        double *point_3d, double **grad);

#endif // LIBGRPP_OVERLAP_GRADIENT_H
