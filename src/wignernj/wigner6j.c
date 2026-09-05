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
 * Wigner 6j symbol via the Racah formula with prime-factorized
 * factorials (Dodds & Wiechers, CPC 4, 268 (1972)) and the multiword-
 * integer Racah sum of Johansson & Forssén (doi:10.1137/15M1021908).
 *
 * Racah formula:
 *   {j1 j2 j3}  =
 * Delta(j1,j2,j3)*Delta(j1,j5,j6)*Delta(j4,j2,j6)*Delta(j4,j5,j3) {j4 j5 j6}
 *     * sum_s (-1)^s (s+1)! / [(s-j1-j2-j3)!(s-j1-j5-j6)!(s-j4-j2-j6)!
 *                               (s-j4-j5-j3)!(j1+j2+j4+j5-s)!
 *                               (j2+j3+j5+j6-s)!(j1+j3+j4+j6-s)!]
 */
#include "pfrac.h"
#include "primes.h"
#include "scratch.h"
#include "wignernj.h"
#include "wignernj_exact.h"
#include <stdlib.h>

static int triangle_ok(int ta, int tb, int tc) {
  if ((ta + tb + tc) & 1)
    return 0;
  if (tc < abs(ta - tb))
    return 0;
  if (tc > ta + tb)
    return 0;
  return 1;
}

static int selection_rules_6j(int tj1, int tj2, int tj3, int tj4, int tj5,
                              int tj6) {
  return triangle_ok(tj1, tj2, tj3) && triangle_ok(tj1, tj5, tj6) &&
         triangle_ok(tj4, tj2, tj6) && triangle_ok(tj4, tj5, tj3);
}

/* Summation bounds for the 6j Racah sum.
 * s_min = max(j1+j2+j3, j1+j5+j6, j4+j2+j6, j4+j5+j3)
 * s_max = min(j1+j2+j4+j5, j2+j3+j5+j6, j1+j3+j4+j6)
 * All in half-integer 2x representation.
 */
static void bounds_6j(int tj1, int tj2, int tj3, int tj4, int tj5, int tj6,
                      int *s_min, int *s_max) {
  int a, b, c, d;
  a = (tj1 + tj2 + tj3) / 2;
  b = (tj1 + tj5 + tj6) / 2;
  c = (tj4 + tj2 + tj6) / 2;
  d = (tj4 + tj5 + tj3) / 2;
  *s_min = a;
  if (b > *s_min)
    *s_min = b;
  if (c > *s_min)
    *s_min = c;
  if (d > *s_min)
    *s_min = d;

  a = (tj1 + tj2 + tj4 + tj5) / 2;
  b = (tj2 + tj3 + tj5 + tj6) / 2;
  c = (tj1 + tj3 + tj4 + tj6) / 2;
  *s_max = a;
  if (b < *s_max)
    *s_max = b;
  if (c < *s_max)
    *s_max = c;
}

/*
 * Add all four Delta triangle-coefficient sqrt factors to outer pfrac.
 * Each Delta(a,b,c) contributes:
 *   numerator:   (a+b-c)! (a-b+c)! (-a+b+c)!
 *   denominator: (a+b+c+1)!
 * (all divided by 2 since ta,tb,tc are 2x values)
 */
static void add_delta_sqrt(pfrac_t *outer, int ta, int tb, int tc) {
  pfrac_mul_factorial(outer, (ta + tb - tc) / 2);
  pfrac_mul_factorial(outer, (ta - tb + tc) / 2);
  pfrac_mul_factorial(outer, (-ta + tb + tc) / 2);
  pfrac_div_factorial(outer, (ta + tb + tc) / 2 + 1);
}

void wigner6j_exact(int tj1, int tj2, int tj3, int tj4, int tj5, int tj6,
                    wignernj_exact_t *out) {
  int s, s_min, s_max, pi;
  wignernj_scratch_t *scratch;
  bigint_ws_t *ws;
  pfrac_t *outer, *term;
  int *lcm_exp;
  bigint_t *sum_pos, *sum_neg, *scaled;
  size_t mw;
  int lcm_max_idx;
  int n_terms;

  wignernj_exact_reset(out);

  if (!selection_rules_6j(tj1, tj2, tj3, tj4, tj5, tj6)) {
    out->is_zero = 1;
    return;
  }

  bounds_6j(tj1, tj2, tj3, tj4, tj5, tj6, &s_min, &s_max);
  if (s_min > s_max) {
    out->is_zero = 1;
    return;
  }

  out->sign = 1; /* phase is carried by (-1)^s in the sum */

  /* Largest factorial argument: max β_u + 1 (the (s+1)! numerator).  Use a
   * safe upper bound on N for sizing the workspace and bigints. */
  mw = bigint_words_for_factorial((tj1 + tj2 + tj3 + tj4 + tj5 + tj6) / 2 + 5);

  /* Acquire the thread's cached scratch and bring every buffer up to
   * the size required by the current call. */
  scratch = wignernj_scratch_acquire();
  ws = &scratch->ws;
  outer = &scratch->pfracs[0];
  lcm_exp = scratch->lcm_exp[0];
  sum_pos = &scratch->bigints[0];
  sum_neg = &scratch->bigints[1];
  scaled = &scratch->bigints[2];

  bigint_ws_reserve(ws, mw);
  pfrac_zero(outer);
  wignernj_scratch_lcm_clear(scratch, 0);
  bigint_set_zero(sum_pos);
  bigint_set_zero(sum_neg);
  bigint_set_zero(scaled);
  bigint_reserve(sum_pos, mw);
  bigint_reserve(sum_neg, mw);
  bigint_reserve(scaled, mw);

  /* Outer sqrt = product of 4 Delta triangle coefficients */
  add_delta_sqrt(outer, tj1, tj2, tj3);
  add_delta_sqrt(outer, tj1, tj5, tj6);
  add_delta_sqrt(outer, tj4, tj2, tj6);
  add_delta_sqrt(outer, tj4, tj5, tj3);

  /* ── Pass 1: build each term pfrac once into the cache, find LCM ── */
  n_terms = s_max - s_min + 1;
  wignernj_scratch_terms_reserve(scratch, n_terms);
  lcm_max_idx = 0;
  for (s = s_min; s <= s_max; s++) {
    /* Denominator factorials of term s:
     * (s-j1-j2-j3)! (s-j1-j5-j6)! (s-j4-j2-j6)! (s-j4-j5-j3)!
     * (j1+j2+j4+j5-s)! (j2+j3+j5+j6-s)! (j1+j3+j4+j6-s)!
     * Numerator factorial: (s+1)!
     */
    term = &scratch->terms[s - s_min];
    pfrac_zero(term);
    pfrac_div_factorial(term, s - (tj1 + tj2 + tj3) / 2);
    pfrac_div_factorial(term, s - (tj1 + tj5 + tj6) / 2);
    pfrac_div_factorial(term, s - (tj4 + tj2 + tj6) / 2);
    pfrac_div_factorial(term, s - (tj4 + tj5 + tj3) / 2);
    pfrac_div_factorial(term, (tj1 + tj2 + tj4 + tj5) / 2 - s);
    pfrac_div_factorial(term, (tj2 + tj3 + tj5 + tj6) / 2 - s);
    pfrac_div_factorial(term, (tj1 + tj3 + tj4 + tj6) / 2 - s);
    pfrac_mul_factorial(term, s + 1);

    /* term.exp[pi] is the net exponent; for LCM we want -min(term.exp[pi]) */
    for (pi = 0; pi < term->max_idx; pi++) {
      int neg = -term->exp[pi]; /* positive = denominator contribution */
      if (neg > lcm_exp[pi])
        lcm_exp[pi] = neg;
    }
    if (term->max_idx > lcm_max_idx)
      lcm_max_idx = term->max_idx;
  }
  wignernj_scratch_lcm_dirty(scratch, 0, lcm_max_idx);

  /* ── Pass 2: incremental walk of LCM-scaled term values ──
   *
   * The 6j Racah term is
   *   T_s = (s+1)! / [(s-α₁)!(s-α₂)!(s-α₃)!(s-α₄)!(β₁-s)!(β₂-s)!(β₃-s)!]
   * with α₁..₄ the four triangle sums and β₁..₃ the three quadruple
   * sums.  The successive-term ratio
   *   T_{s+1}/T_s = (s+2)(β₁-s)(β₂-s)(β₃-s)
   *               / [(s+1-α₁)(s+1-α₂)(s+1-α₃)(s+1-α₄)]
   * gives an integer recurrence on scaled_s = LCM·T_s,
   *   scaled_{s+1} (s+1-α₁)(s+1-α₂)(s+1-α₃)(s+1-α₄)
   *     = scaled_s (s+2)(β₁-s)(β₂-s)(β₃-s).
   *
   * Each factor is bounded above by ~2·MAX_FACTORIAL_ARG: triangle
   * sums αᵢ are at most MAX_FACTORIAL_ARG-1 (else build_outer_sqrt_6j
   * would have aborted in Pass 1), and quadruple sums βⱼ pair two
   * triangles that share one j, so βⱼ ≤ 2·MAX_FACTORIAL_ARG.  The
   * batched quadruple product therefore fits in uint64 when
   * (2·MAX_FACTORIAL_ARG)⁴ ≤ 2⁶⁴, i.e. MAX_FACTORIAL_ARG ≤ 2¹⁵ ·
   * something; we conservatively pin the static_assert at 32767.
   * Algorithm D in bigint_div128 handles the resulting 64-bit
   * divisor on every backend.
   */
#if MAX_FACTORIAL_ARG > 32767
#error                                                                         \
    "MAX_FACTORIAL_ARG too large for batched 6j Pass-2 ratio update (quadruple product would overflow uint64); regenerate the prime table with a smaller --limit, or split this loop into pairwise mul/div"
#endif
  {
    uint64_t a1 = (uint64_t)((tj1 + tj2 + tj3) / 2);
    uint64_t a2 = (uint64_t)((tj1 + tj5 + tj6) / 2);
    uint64_t a3 = (uint64_t)((tj4 + tj2 + tj6) / 2);
    uint64_t a4 = (uint64_t)((tj4 + tj5 + tj3) / 2);
    uint64_t b1 = (uint64_t)((tj1 + tj2 + tj4 + tj5) / 2);
    uint64_t b2 = (uint64_t)((tj2 + tj3 + tj5 + tj6) / 2);
    uint64_t b3 = (uint64_t)((tj1 + tj3 + tj4 + tj6) / 2);

    /* Seed: scaled_{s_min} = LCM·T_{s_min} via the existing
     * prime-power helper from the cached pfrac for term 0. */
    pfrac_lcm_scaled_product(scaled, lcm_exp, scratch->terms[0].exp, +1,
                             lcm_max_idx, ws);
    if ((s_min & 1) == 0)
      bigint_add(sum_pos, sum_pos, scaled);
    else
      bigint_add(sum_neg, sum_neg, scaled);

    /* Step s -> s+1.  At s = s_max we stop without applying the
     * ratio: (β-s) would hit 0 there for the i with βᵢ = s_max. */
    for (s = s_min; s < s_max; s++) {
      uint64_t us = (uint64_t)s;
      uint64_t num = (us + 2) * (b1 - us) * (b2 - us) * (b3 - us);
      uint64_t den =
          (us + 1 - a1) * (us + 1 - a2) * (us + 1 - a3) * (us + 1 - a4);
      bigint_mul_u64(scaled, scaled, num);
      bigint_div_u64_exact(scaled, scaled, den);

      if (((s + 1) & 1) == 0)
        bigint_add(sum_pos, sum_pos, scaled);
      else
        bigint_add(sum_neg, sum_neg, scaled);
    }
  }

  out->sum_sign = bigint_sub_signed(&out->sum, sum_pos, sum_neg);

  /* Build output bigints: reserve before set_u64 so the cap-from-0
   * growth and the subsequent grow-to-mw collapse into a single
   * realloc per output bigint. */
  bigint_reserve(&out->int_num, mw);
  bigint_reserve(&out->int_den, mw);
  bigint_reserve(&out->sqrt_num, mw);
  bigint_reserve(&out->sqrt_den, mw);
  bigint_set_u64(&out->int_num, 1);
  bigint_set_u64(&out->int_den, 1);
  bigint_set_u64(&out->sqrt_num, 1);
  bigint_set_u64(&out->sqrt_den, 1);

  pfrac_to_sqrt_rational_ws(outer, &out->int_num, &out->int_den, &out->sqrt_num,
                            &out->sqrt_den, ws);

  pfrac_bigint_mul_prime_pow_array(&out->int_den, lcm_exp, lcm_max_idx, ws);

  wignernj_scratch_relinquish(scratch);
}

/* ── public API ──────────────────────────────────────────────────────────── */

int wigner6j_max_factorial(int tj1, int tj2, int tj3, int tj4, int tj5,
                           int tj6) {
  /* (s+1)! with s_max = max sum-of-4 across the three combinations
   * the Racah sum's denominator factorials reference dominates the
   * factorial-argument upper bound.  Triangle-Δ denominators are
   * sums-of-3 and are bounded by these. */
  int a = (tj1 + tj2 + tj4 + tj5) / 2;
  int b = (tj2 + tj3 + tj5 + tj6) / 2;
  int c = (tj1 + tj3 + tj4 + tj6) / 2;
  int s = a;
  if (b > s)
    s = b;
  if (c > s)
    s = c;
  return s + 1;
}

float wigner6j_f(int tj1, int tj2, int tj3, int tj4, int tj5, int tj6) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  float result;
  wigner6j_exact(tj1, tj2, tj3, tj4, tj5, tj6, &s->exact);
  result = wignernj_exact_to_float(&s->exact);
  wignernj_scratch_relinquish(s);
  return result;
}

double wigner6j(int tj1, int tj2, int tj3, int tj4, int tj5, int tj6) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  double result;
  wigner6j_exact(tj1, tj2, tj3, tj4, tj5, tj6, &s->exact);
  result = wignernj_exact_to_double(&s->exact);
  wignernj_scratch_relinquish(s);
  return result;
}

long double wigner6j_l(int tj1, int tj2, int tj3, int tj4, int tj5, int tj6) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  long double result;
  wigner6j_exact(tj1, tj2, tj3, tj4, tj5, tj6, &s->exact);
  result = wignernj_exact_to_long_double(&s->exact);
  wignernj_scratch_relinquish(s);
  return result;
}

#ifdef WIGNERNJ_HAVE_QUADMATH
__float128 wigner6j_q(int tj1, int tj2, int tj3, int tj4, int tj5, int tj6) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  __float128 result;
  wigner6j_exact(tj1, tj2, tj3, tj4, tj5, tj6, &s->exact);
  result = wignernj_exact_to_float128(&s->exact);
  wignernj_scratch_relinquish(s);
  return result;
}
#endif

#ifdef WIGNERNJ_HAVE_MPFR
void wigner6j_mpfr(mpfr_t rop, int tj1, int tj2, int tj3, int tj4, int tj5,
                   int tj6, mpfr_rnd_t rnd) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  wigner6j_exact(tj1, tj2, tj3, tj4, tj5, tj6, &s->exact);
  wignernj_exact_to_mpfr(rop, &s->exact, rnd);
  wignernj_scratch_relinquish(s);
}
#endif
