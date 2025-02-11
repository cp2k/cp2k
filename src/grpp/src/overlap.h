/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2023 Alexander Oleynichenko
 */

#ifndef LIBGRPP_OVERLAP_H
#define LIBGRPP_OVERLAP_H

#include "libgrpp.h"

void libgrpp_overlap_integrals(libgrpp_shell_t *shell_A,
                               libgrpp_shell_t *shell_B,
                               double *overlap_matrix);

#endif // LIBGRPP_OVERLAP_H
