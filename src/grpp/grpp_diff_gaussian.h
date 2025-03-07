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
 * Differentiation of contracted Gaussian functions.
 * Derivatives are then used to calculate analytic gradients of 1-el integrals.
 */

#ifndef LIBGRPP_DIFF_GAUSSIAN_H
#define LIBGRPP_DIFF_GAUSSIAN_H

#include "libgrpp_types.h"

void libgrpp_differentiate_shell(libgrpp_shell_t *shell,
                                 libgrpp_shell_t **shell_minus,
                                 libgrpp_shell_t **shell_plus);

#endif // LIBGRPP_DIFF_GAUSSIAN_H
