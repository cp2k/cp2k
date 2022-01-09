/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#ifndef DBM_MULTIPLY_H
#define DBM_MULTIPLY_H

#include <stdbool.h>

#include "dbm_matrix.h"

/*******************************************************************************
 * \brief Performs a multiplication of two dbm_matrix_t matrices,
          as  C := alpha * op( A ) * op( B ) + beta * C.

          The filter_eps parameter is used to filter the resulting matrix.
          The filtering criterion is whether the block-frobenius norm is less
          than the specified epsilon. One-the-fly filtering is done such that
          individual multiplications are skipped if the product of the frobenius
          norms of the left- and right-matrix blocks are less than the specified
          epsilon divided by the maximum number of possible multiplies in each
          row. In addition a final filtering is done as well with the same
          epsilon value.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply(const bool transa, const bool transb, const double alpha,
                  dbm_matrix_t *matrix_a, dbm_matrix_t *matrix_b,
                  const double beta, dbm_matrix_t *matrix_c,
                  const bool retain_sparsity, const double filter_eps,
                  long *flop);

#endif

// EOF
