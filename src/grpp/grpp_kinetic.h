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

#ifndef LIBGRPP_KINETIC_H
#define LIBGRPP_KINETIC_H

#include "libgrpp_types.h"

void libgrpp_kinetic_energy_integrals(libgrpp_shell_t *shell_A,
                                      libgrpp_shell_t *shell_B,
                                      double *kinetic_matrix);

#endif // LIBGRPP_KINETIC_H
