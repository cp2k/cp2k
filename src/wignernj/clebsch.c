/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include "cp2k_wignernj_prefix.h"

/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * Clebsch-Gordan coefficients via the Wigner 3j symbol.
 *
 * Definition:
 *   <j1 m1; j2 m2 | J M> = (-1)^(j1-j2+M) * sqrt(2J+1)
 *                           * (j1 j2  J )
 *                             (m1 m2 -M)
 *
 * Arguments are passed as twice their value (tj = 2*j, tm = 2*m).
 * Requires m1+m2 = M (else zero).
 */
#include "scratch.h"
#include "wignernj.h"
#include "wignernj_exact.h"

/*
 * Internal: compute the exact CG coefficient.
 * The 3j symbol is computed exactly; sqrt(2J+1) is folded in by factoring
 * (tJ+1) and distributing prime powers into int_num and sqrt_num.
 */
static void clebsch_gordan_exact(int tj1, int tm1, int tj2, int tm2, int tJ,
                                 int tM, wignernj_exact_t *out) {
  int phase;
  /* m conservation -- short-circuit before wigner3j_exact (which would
   * also reset out, but only after a selection-rule check that doesn't
   * cover the CG-specific m1+m2 = M constraint).  Reset here so that
   * out is in a known-zero state when we set is_zero=1 below. */
  if (tm1 + tm2 != tM) {
    wignernj_exact_reset(out);
    out->is_zero = 1;
    return;
  }

  /* Compute the underlying 3j symbol exactly.  wigner3j_exact resets
   * out itself, so no init is needed here. */
  wigner3j_exact(tj1, tj2, tJ, tm1, tm2, -tM, out);
  if (out->is_zero)
    return;

  /* Phase: (-1)^(j1-j2+M) = (-1)^((tj1-tj2+tM)/2) */
  phase = (((tj1 - tj2 + tM) / 2) & 1) ? -1 : 1;
  out->sign *= phase;

  /* Multiply by sqrt(2J+1). */
  wignernj_exact_mul_sqrt_int(out, tJ + 1);
}

int clebsch_gordan_max_factorial(int tj1, int tm1, int tj2, int tm2, int tJ,
                                 int tM) {
  /* CG = (-1)^(...) * sqrt(2J+1) * 3j(j1, j2, J; m1, m2, -M).  The
   * sqrt(2J+1) factor is folded in via bigint_mul_prime_pow without
   * touching the factorial cache, so the bound is exactly
   * wigner3j_max_factorial of the underlying 3j. */
  return wigner3j_max_factorial(tj1, tj2, tJ, tm1, tm2, -tM);
}

float clebsch_gordan_f(int tj1, int tm1, int tj2, int tm2, int tJ, int tM) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  float r;
  clebsch_gordan_exact(tj1, tm1, tj2, tm2, tJ, tM, &s->exact);
  r = wignernj_exact_to_float(&s->exact);
  wignernj_scratch_relinquish(s);
  return r;
}

double clebsch_gordan(int tj1, int tm1, int tj2, int tm2, int tJ, int tM) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  double r;
  clebsch_gordan_exact(tj1, tm1, tj2, tm2, tJ, tM, &s->exact);
  r = wignernj_exact_to_double(&s->exact);
  wignernj_scratch_relinquish(s);
  return r;
}

long double clebsch_gordan_l(int tj1, int tm1, int tj2, int tm2, int tJ,
                             int tM) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  long double r;
  clebsch_gordan_exact(tj1, tm1, tj2, tm2, tJ, tM, &s->exact);
  r = wignernj_exact_to_long_double(&s->exact);
  wignernj_scratch_relinquish(s);
  return r;
}

#ifdef WIGNERNJ_HAVE_QUADMATH
__float128 clebsch_gordan_q(int tj1, int tm1, int tj2, int tm2, int tJ,
                            int tM) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  __float128 r;
  clebsch_gordan_exact(tj1, tm1, tj2, tm2, tJ, tM, &s->exact);
  r = wignernj_exact_to_float128(&s->exact);
  wignernj_scratch_relinquish(s);
  return r;
}
#endif

#ifdef WIGNERNJ_HAVE_MPFR
#include <mpfr.h>
void clebsch_gordan_mpfr(mpfr_t rop, int tj1, int tm1, int tj2, int tm2, int tJ,
                         int tM, mpfr_rnd_t rnd) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  clebsch_gordan_exact(tj1, tm1, tj2, tm2, tJ, tM, &s->exact);
  wignernj_exact_to_mpfr(rop, &s->exact, rnd);
  wignernj_scratch_relinquish(s);
}
#endif
