/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * Fano X-coefficient via the Wigner 9j symbol.
 *
 * Definition (Fano 1953; Edmonds 1957 §6.4; Tamura 1970):
 *   X(j1 j2 j12; j3 j4 j34; j13 j24 J)
 *     = sqrt[(2j12+1)(2j34+1)(2j13+1)(2j24+1)]
 *       * { j1   j2   j12 }
 *         { j3   j4   j34 }
 *         { j13  j24  J   }
 *
 * The X-coefficient is a normalisation variant of the 9j symbol used in
 * the analysis of polarisation correlations and in the Russell--Saunders
 * coupling literature; it carries no information beyond the underlying
 * 9j symbol but is produced directly because some downstream
 * applications expect that normalisation.
 *
 * Arguments are passed as twice their value (tj = 2*j); the four
 * (2j+1) factors above appear in the exact pipeline as (tj12+1),
 * (tj34+1), (tj13+1), and (tj24+1) and are folded into the existing
 * wignernj_exact_t in the same way the (2J+1) factor is folded into the
 * Clebsch--Gordan coefficient (see src/clebsch.c).
 */
#include "scratch.h"
#include "wignernj.h"
#include "wignernj_exact.h"

/*
 * Internal: compute the exact Fano X coefficient.  The 9j symbol is
 * computed exactly by the existing wigner9j_exact pipeline; each of the
 * four (2j+1) factors above is then folded in via
 * wignernj_exact_mul_sqrt_int, which is the same helper that
 * src/clebsch.c uses for its single sqrt(2J+1) factor.
 */
static void fano_x_exact(int tj1, int tj2, int tj12, int tj3, int tj4, int tj34,
                         int tj13, int tj24, int tJ, wignernj_exact_t *out) {
  /* The underlying 9j */
  wigner9j_exact(tj1, tj2, tj12, tj3, tj4, tj34, tj13, tj24, tJ, out);
  if (out->is_zero)
    return;

  /* Fold sqrt[(2j12+1)(2j34+1)(2j13+1)(2j24+1)] into the exact tuple. */
  wignernj_exact_mul_sqrt_int(out, tj12 + 1);
  wignernj_exact_mul_sqrt_int(out, tj34 + 1);
  wignernj_exact_mul_sqrt_int(out, tj13 + 1);
  wignernj_exact_mul_sqrt_int(out, tj24 + 1);
}

/* ── public API ──────────────────────────────────────────────────────────── */

int fano_x_max_factorial(int tj1, int tj2, int tj12, int tj3, int tj4, int tj34,
                         int tj13, int tj24, int tJ) {
  /* X = sqrt[...] * 9j; the sqrt prefactor is folded in via
   * bigint_mul_prime_pow without touching the factorial cache. */
  return wigner9j_max_factorial(tj1, tj2, tj12, tj3, tj4, tj34, tj13, tj24, tJ);
}

float fano_x_f(int tj1, int tj2, int tj12, int tj3, int tj4, int tj34, int tj13,
               int tj24, int tJ) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  float r;
  fano_x_exact(tj1, tj2, tj12, tj3, tj4, tj34, tj13, tj24, tJ, &s->exact);
  r = wignernj_exact_to_float(&s->exact);
  wignernj_scratch_relinquish(s);
  return r;
}

double fano_x(int tj1, int tj2, int tj12, int tj3, int tj4, int tj34, int tj13,
              int tj24, int tJ) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  double r;
  fano_x_exact(tj1, tj2, tj12, tj3, tj4, tj34, tj13, tj24, tJ, &s->exact);
  r = wignernj_exact_to_double(&s->exact);
  wignernj_scratch_relinquish(s);
  return r;
}

long double fano_x_l(int tj1, int tj2, int tj12, int tj3, int tj4, int tj34,
                     int tj13, int tj24, int tJ) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  long double r;
  fano_x_exact(tj1, tj2, tj12, tj3, tj4, tj34, tj13, tj24, tJ, &s->exact);
  r = wignernj_exact_to_long_double(&s->exact);
  wignernj_scratch_relinquish(s);
  return r;
}

#ifdef WIGNERNJ_HAVE_QUADMATH
__float128 fano_x_q(int tj1, int tj2, int tj12, int tj3, int tj4, int tj34,
                    int tj13, int tj24, int tJ) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  __float128 r;
  fano_x_exact(tj1, tj2, tj12, tj3, tj4, tj34, tj13, tj24, tJ, &s->exact);
  r = wignernj_exact_to_float128(&s->exact);
  wignernj_scratch_relinquish(s);
  return r;
}
#endif

#ifdef WIGNERNJ_HAVE_MPFR
#include "wignernj_mpfr.h"
void fano_x_mpfr(mpfr_t rop, int tj1, int tj2, int tj12, int tj3, int tj4,
                 int tj34, int tj13, int tj24, int tJ, mpfr_rnd_t rnd) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  fano_x_exact(tj1, tj2, tj12, tj3, tj4, tj34, tj13, tj24, tJ, &s->exact);
  wignernj_exact_to_mpfr(rop, &s->exact, rnd);
  wignernj_scratch_relinquish(s);
}
#endif
