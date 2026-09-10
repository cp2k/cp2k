/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola */
#ifndef WIGNERNJ_EXACT_H
#define WIGNERNJ_EXACT_H

#include "bigint.h"

/*
 * Exact symbolic result of a Wigner symbol computation.
 *
 * The floating-point value is:
 *   sign * sum_sign * bigint_to_T(sum) * bigint_to_T(int_num)
 *                   / bigint_to_T(int_den)
 *                   * sqrt(bigint_to_T(sqrt_num) / bigint_to_T(sqrt_den))
 *
 * where T is float, double, or long double.
 *
 * int_den absorbs both the denominator from the outer sqrt factor and
 * the LCM of all Racah sum term denominators.
 *
 * When is_zero is set the symbol vanishes by selection rules; all bigints
 * hold their init-time zero state and carry no useful data.
 */
typedef struct {
  bigint_t sum;      /* magnitude of Racah integer sum                */
  int sum_sign;      /* sign of Racah sum: +1 or -1 (1 when sum == 0) */
  bigint_t int_num;  /* integer-part numerator of outer sqrt factor   */
  bigint_t int_den;  /* integer-part denominator (includes LCM)       */
  bigint_t sqrt_num; /* squarefree sqrt numerator                     */
  bigint_t sqrt_den; /* squarefree sqrt denominator                   */
  int sign;          /* overall phase: +1 or -1                       */
  int is_zero;       /* 1 if symbol is identically zero               */
} wignernj_exact_t;

void wignernj_exact_init(wignernj_exact_t *e);
/* Reset an already-initialised wignernj_exact_t to a fresh-result state
 * without freeing its bigint buffers (used by the scratch cache). */
void wignernj_exact_reset(wignernj_exact_t *e);
void wigner3j_exact(int tj1, int tj2, int tj3, int tm1, int tm2, int tm3,
                    wignernj_exact_t *out);
void wigner6j_exact(int tj1, int tj2, int tj3, int tj4, int tj5, int tj6,
                    wignernj_exact_t *out);
void wigner9j_exact(int tj11, int tj12, int tj13, int tj21, int tj22, int tj23,
                    int tj31, int tj32, int tj33, wignernj_exact_t *out);
void wignernj_exact_free(wignernj_exact_t *e);

/* Multiply the value carried by `e` by sqrt(k), in place, by factoring
 * the positive integer k and distributing its prime powers across the
 * outer-sqrt tuple:
 *   even prime power p^(2a) → p^a into int_num (outside sqrt)
 *   odd  prime power p^(2a+1) → p^a into int_num, one p into sqrt_num
 * No-op when k == 1.  Caller is responsible for k > 0.
 *
 * Used by clebsch_gordan_exact to fold sqrt(2J+1) and by fano_x_exact
 * to fold sqrt[(2j12+1)(2j34+1)(2j13+1)(2j24+1)] into the exact tuple
 * produced by the underlying Wigner symbol. */
void wignernj_exact_mul_sqrt_int(wignernj_exact_t *e, int k);

float wignernj_exact_to_float(const wignernj_exact_t *e);
double wignernj_exact_to_double(const wignernj_exact_t *e);
long double wignernj_exact_to_long_double(const wignernj_exact_t *e);

#ifdef WIGNERNJ_HAVE_QUADMATH
#include <quadmath.h>
__float128 wignernj_exact_to_float128(const wignernj_exact_t *e);
#endif

#ifdef WIGNERNJ_HAVE_MPFR
#include <mpfr.h>
void wignernj_exact_to_mpfr(mpfr_t rop, const wignernj_exact_t *e,
                            mpfr_rnd_t rnd);
#endif

#endif /* WIGNERNJ_EXACT_H */
