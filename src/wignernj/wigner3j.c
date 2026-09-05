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
 * Wigner 3j symbol via the Racah formula with prime-factorized
 * factorials (Dodds & Wiechers, CPC 4, 268 (1972)) and the multiword-
 * integer Racah sum of Johansson & Forssén (doi:10.1137/15M1021908).
 *
 * All arguments are passed as 2*j integers (e.g. j=3/2 → tj=3).
 *
 * Racah formula:
 *   (j1 j2 j3)  = (-1)^(j1-j2-m3) * Delta(j1,j2,j3)
 *   (m1 m2 m3)    * sqrt[(j1+m1)!(j1-m1)!(j2+m2)!(j2-m2)!(j3+m3)!(j3-m3)!]
 *                 * sum_s (-1)^s / [s!(j1+j2-j3-s)!(j1-m1-s)!
 *                                    (j2+m2-s)!(j3-j2+m1+s)!(j3-j1-m2+s)!]
 *
 * Delta(a,b,c) = sqrt[(a+b-c)!(a-b+c)!(-a+b+c)! / (a+b+c+1)!]
 */
#include "pfrac.h"
#include "primes.h"
#include "scratch.h"
#include "wignernj.h"
#include "wignernj_exact.h"
#include <stdlib.h> /* abs */

/* ── selection rules ─────────────────────────────────────────────────────── */

static int triangle_ok(int ta, int tb, int tc) {
  /* ta+tb+tc must be even (so that j-sums are integers) */
  if ((ta + tb + tc) & 1)
    return 0;
  if (tc < abs(ta - tb))
    return 0;
  if (tc > ta + tb)
    return 0;
  return 1;
}

static int selection_rules_3j(int tj1, int tj2, int tj3, int tm1, int tm2,
                              int tm3) {
  /* m conservation */
  if (tm1 + tm2 + tm3 != 0)
    return 0;
  /* m within j */
  if (abs(tm1) > tj1 || abs(tm2) > tj2 || abs(tm3) > tj3)
    return 0;
  /* same half-integer type: tmi must have same parity as tji */
  if (((tj1 - tm1) & 1) || ((tj2 - tm2) & 1) || ((tj3 - tm3) & 1))
    return 0;
  /* triangle */
  if (!triangle_ok(tj1, tj2, tj3))
    return 0;
  return 1;
}

/* ── Racah summation bounds ──────────────────────────────────────────────── */

static void bounds_3j(int tj1, int tj2, int tj3, int tm1, int tm2, int *s_min,
                      int *s_max) {
  /* All factorial arguments in the Racah sum must be non-negative.
   * Denominator factorials: s!, (j1+j2-j3-s)!, (j1-m1-s)!, (j2+m2-s)!,
   *                         (j3-j2+m1+s)!, (j3-j1-m2+s)!
   * Lower bounds (from last two): s >= j2-j3-m1 and s >= j1+m2-j3
   * Upper bounds (from first four): s <= j1+j2-j3, s <= j1-m1, s <= j2+m2
   * All in 2x representation, divided by 2.
   */
  int a = (tj2 - tj3 - tm1) / 2; /* lower from (j3-j2+m1+s)! */
  int b = (tj1 + tm2 - tj3) / 2; /* lower from (j3-j1-m2+s)! */
  *s_min = (a > b) ? a : b;
  if (*s_min < 0)
    *s_min = 0;

  a = (tj1 + tj2 - tj3) / 2;
  b = (tj1 - tm1) / 2;
  int c = (tj2 + tm2) / 2;
  *s_max = a;
  if (b < *s_max)
    *s_max = b;
  if (c < *s_max)
    *s_max = c;
}

/* ── core exact computation ──────────────────────────────────────────────── */

/*
 * Build the outer sqrt factor pfrac_t for the 3j symbol.
 * This represents the argument under the outer sqrt sign:
 *   (j1+j2-j3)!(j1-j2+j3)!(-j1+j2+j3)! *
 * (j1+m1)!(j1-m1)!(j2+m2)!(j2-m2)!(j3+m3)!(j3-m3)! / (j1+j2+j3+1)!
 */
static void build_outer_sqrt_3j(pfrac_t *outer, int tj1, int tj2, int tj3,
                                int tm1, int tm2, int tm3) {
  pfrac_zero(outer);
  /* Delta numerator factorials */
  pfrac_mul_factorial(outer, (tj1 + tj2 - tj3) / 2);
  pfrac_mul_factorial(outer, (tj1 - tj2 + tj3) / 2);
  pfrac_mul_factorial(outer, (-tj1 + tj2 + tj3) / 2);
  /* m-dependent sqrt numerator factorials */
  pfrac_mul_factorial(outer, (tj1 + tm1) / 2);
  pfrac_mul_factorial(outer, (tj1 - tm1) / 2);
  pfrac_mul_factorial(outer, (tj2 + tm2) / 2);
  pfrac_mul_factorial(outer, (tj2 - tm2) / 2);
  pfrac_mul_factorial(outer, (tj3 + tm3) / 2);
  pfrac_mul_factorial(outer, (tj3 - tm3) / 2);
  /* Delta denominator factorial */
  pfrac_div_factorial(outer, (tj1 + tj2 + tj3) / 2 + 1);
}

void wigner3j_exact(int tj1, int tj2, int tj3, int tm1, int tm2, int tm3,
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

  /* Reset (rather than init) `out`: callers always pass a previously-
   * initialised wignernj_exact_t (either freshly allocated or held in
   * the cached scratch).  Reset preserves the bigint capacities so
   * subsequent reserves are no-ops once they have been sized once. */
  wignernj_exact_reset(out);

  if (!selection_rules_3j(tj1, tj2, tj3, tm1, tm2, tm3)) {
    out->is_zero = 1;
    return;
  }

  /* Closed-form fast path for (j1 j2 j3; 0 0 0): the standard Racah
   * single-sum collapses to the well-known closed form
   *
   *   (j1 j2 j3; 0 0 0)
   *     = (-1)^g * sqrt[Delta^2(j1,j2,j3)]
   *               * g! / [(g-j1)! (g-j2)! (g-j3)!]
   *
   * with g = (j1+j2+j3)/2; the symbol vanishes unless g is integer
   * (i.e. (j1+j2+j3) even).  This skips the entire Racah-sum
   * machinery (Pass~1 LCM build + Pass~2 ratio recurrence) in
   * favour of a single multinomial-coefficient bigint
   * materialisation.  Used directly here and inherited by every
   * Gaunt evaluation (which calls (l1 l2 l3; 0 0 0) once on every
   * call). */
  if (tm1 == 0 && tm2 == 0 && tm3 == 0) {
    int Lhalf = (tj1 + tj2 + tj3) / 4; /* g = (j1+j2+j3)/2 if integer */
    if ((tj1 + tj2 + tj3) & 2) {       /* (j1+j2+j3) odd → vanishes */
      out->is_zero = 1;
      return;
    }
    int j1 = tj1 / 2, j2 = tj2 / 2, j3 = tj3 / 2;
    int g = Lhalf;
    out->sign = (g & 1) ? -1 : +1;

    scratch = wignernj_scratch_acquire();
    ws = &scratch->ws;
    outer = &scratch->pfracs[0];
    pfrac_t *closed = &scratch->pfracs[1];

    mw = bigint_words_for_factorial((tj1 + tj2 + tj3) / 2 + 5);
    bigint_ws_reserve(ws, mw);
    bigint_reserve(&out->sum, mw);
    bigint_reserve(&out->int_num, mw);
    bigint_reserve(&out->int_den, mw);
    bigint_reserve(&out->sqrt_num, mw);
    bigint_reserve(&out->sqrt_den, mw);

    /* Outer sqrt factor: Delta^2(j1,j2,j3) */
    pfrac_zero(outer);
    pfrac_mul_factorial(outer, (tj1 + tj2 - tj3) / 2);
    pfrac_mul_factorial(outer, (tj1 - tj2 + tj3) / 2);
    pfrac_mul_factorial(outer, (-tj1 + tj2 + tj3) / 2);
    pfrac_div_factorial(outer, (tj1 + tj2 + tj3) / 2 + 1);

    /* Closed-form Racah sum: g! / [(g-j1)! (g-j2)! (g-j3)!] */
    pfrac_zero(closed);
    pfrac_mul_factorial(closed, g);
    pfrac_div_factorial(closed, g - j1);
    pfrac_div_factorial(closed, g - j2);
    pfrac_div_factorial(closed, g - j3);

    /* Materialise the multinomial coefficient as out->sum.  Every
     * exp[i] of `closed` is non-negative because the multinomial
     * is a positive integer (triangle inequality forces g >= j_i). */
    bigint_set_u64(&out->sum, 1);
    pfrac_bigint_mul_prime_pow_array(&out->sum, closed->exp, closed->max_idx,
                                     ws);
    out->sum_sign = +1;

    /* Outer sqrt split */
    bigint_set_u64(&out->int_num, 1);
    bigint_set_u64(&out->int_den, 1);
    bigint_set_u64(&out->sqrt_num, 1);
    bigint_set_u64(&out->sqrt_den, 1);
    pfrac_to_sqrt_rational_ws(outer, &out->int_num, &out->int_den,
                              &out->sqrt_num, &out->sqrt_den, ws);

    wignernj_scratch_relinquish(scratch);
    return;
  }

  bounds_3j(tj1, tj2, tj3, tm1, tm2, &s_min, &s_max);
  if (s_min > s_max) {
    out->is_zero = 1;
    return;
  }

  /* Overall phase: (-1)^(j1-j2-m3) */
  out->sign = ((((tj1 - tj2 - tm3) / 2) & 1) == 0) ? 1 : -1;

  /* Largest factorial argument is bounded by (tj1+tj2+tj3)/2+1 plus
   * (tji+|tmi|)/2 ≤ tji.  Use a safe upper bound for sizing. */
  mw = bigint_words_for_factorial((tj1 + tj2 + tj3) / 2 + 5);

  /* Acquire the thread's cached scratch and bring every buffer up to
   * the size required by the current call (no-op if it is already
   * big enough from a prior, larger call). */
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

  /* Compute outer sqrt factor */
  build_outer_sqrt_3j(outer, tj1, tj2, tj3, tm1, tm2, tm3);

  /* ── Pass 1: build each term pfrac once into the cache, find LCM ── */
  n_terms = s_max - s_min + 1;
  wignernj_scratch_terms_reserve(scratch, n_terms);
  lcm_max_idx = 0;
  for (s = s_min; s <= s_max; s++) {
    /* Build denominator pfrac for term s:
     * s! * (j1+j2-j3-s)! * (j1-m1-s)! * (j2+m2-s)! *
     * (j3-j2+m1+s)! * (j3-j1-m2+s)!
     */
    term = &scratch->terms[s - s_min];
    pfrac_zero(term);
    pfrac_mul_factorial(term, s);
    pfrac_mul_factorial(term, (tj1 + tj2 - tj3) / 2 - s);
    pfrac_mul_factorial(term, (tj1 - tm1) / 2 - s);
    pfrac_mul_factorial(term, (tj2 + tm2) / 2 - s);
    pfrac_mul_factorial(term, (tj3 - tj2 + tm1) / 2 + s);
    pfrac_mul_factorial(term, (tj3 - tj1 - tm2) / 2 + s);

    for (pi = 0; pi < term->max_idx; pi++) {
      if (term->exp[pi] > lcm_exp[pi])
        lcm_exp[pi] = term->exp[pi];
    }
    if (term->max_idx > lcm_max_idx)
      lcm_max_idx = term->max_idx;
  }
  wignernj_scratch_lcm_dirty(scratch, 0, lcm_max_idx);

  /* ── Pass 2: incremental walk of LCM-scaled term values ──
   *
   * Define D_s = s! · (a-s)! · (b-s)! · (c-s)! · (d+s)! · (e+s)!
   * with
   *   a = (tj1+tj2-tj3)/2,  b = (tj1-tm1)/2,  c = (tj2+tm2)/2,
   *   d = (tj3-tj2+tm1)/2,  e = (tj3-tj1-tm2)/2.
   * Then scaled_s := LCM/D_s satisfies the integer recurrence
   *
   *   scaled_{s+1} · (s+1)(d+s+1)(e+s+1) = scaled_s · (a-s)(b-s)(c-s).
   *
   * All six factors are positive integers <= MAX_FACTORIAL_ARG within
   * the sum range, and Pass 1 has already established the LCM and
   * stored every term's pfrac in scratch->terms[].  Each step
   * therefore reduces to a single multiplication by the triple
   * numerator product and a single exact division by the triple
   * denominator product on the running scaled bigint -- two
   * single-word bigint operations per term, replacing the O(pi(j))
   * prime-power expansion of pfrac_lcm_scaled_product that profiling
   * identified (62.7% in bigint_mul_u64 inside that helper at
   * j=4000) as the dominant cost at large j.
   *
   * Each triple product is bounded by MAX_FACTORIAL_ARG^3.  The
   * static_assert below pins this to a compile-time check against
   * UINT64_MAX so a regenerated prime table that exceeded the
   * batched-product safe range would break the build with a clear
   * message rather than silently miscompute.  At the default
   * MAX_FACTORIAL_ARG = 20020, MAX_FACTORIAL_ARG^3 ~ 8e12 has more
   * than 11 bits of headroom under UINT64_MAX ~ 1.8e19; the
   * crossover at 2^64 is around MAX_FACTORIAL_ARG ~ 2.64e6.
   * Algorithm D in bigint_div128 handles arbitrary 64-bit divisors
   * on every backend (hardware divq on x86-64, __uint128_t-using
   * Algorithm D on aarch64 etc., pure-C99 Algorithm D on toolchains
   * without __uint128_t), so the batched divisor is correct on
   * every supported target.
   */
#if MAX_FACTORIAL_ARG > 2642245
#error                                                                         \
    "MAX_FACTORIAL_ARG too large for batched 3j Pass-2 ratio update (triple product would overflow uint64); regenerate the prime table with a smaller --limit, or replace the batched mul/div in this loop with six single-factor calls"
#endif
  {
    uint64_t a = (uint64_t)((tj1 + tj2 - tj3) / 2);
    uint64_t b = (uint64_t)((tj1 - tm1) / 2);
    uint64_t c = (uint64_t)((tj2 + tm2) / 2);
    uint64_t d = (uint64_t)((tj3 - tj2 + tm1) / 2);
    uint64_t e = (uint64_t)((tj3 - tj1 - tm2) / 2);

    /* Seed: scaled_{s_min} = LCM/D_{s_min} via the existing
     * prime-power helper, using the cached pfrac for the first
     * term.  This is the only call to that helper per wigner3j;
     * every subsequent term is reached by ratio update below. */
    pfrac_lcm_scaled_product(scaled, lcm_exp, scratch->terms[0].exp, -1,
                             lcm_max_idx, ws);
    if ((s_min & 1) == 0)
      bigint_add(sum_pos, sum_pos, scaled);
    else
      bigint_add(sum_neg, sum_neg, scaled);

    /* Step s -> s+1.  At s = s_max we stop without applying the
     * ratio: a-s would be 0 there (one of the bounds enforces
     * s_max <= a, similarly b, c) and the next scaled_s would be
     * 0 by the well-defined limit, but we never need it. */
    for (s = s_min; s < s_max; s++) {
      uint64_t us = (uint64_t)s;
      uint64_t num = (a - us) * (b - us) * (c - us);
      uint64_t den = (us + 1) * (d + us + 1) * (e + us + 1);
      bigint_mul_u64(scaled, scaled, num);
      bigint_div_u64(scaled, scaled, den);

      if (((s + 1) & 1) == 0)
        bigint_add(sum_pos, sum_pos, scaled);
      else
        bigint_add(sum_neg, sum_neg, scaled);
    }
  }

  out->sum_sign = bigint_sub_signed(&out->sum, sum_pos, sum_neg);

  /* ── Build output bigints from outer sqrt factor and LCM ── */
  /* Reserve before set_u64: bigint_set_u64 calls bigint_ensure(.., 1)
   * which would otherwise grow cap from 0 to 4 first (one realloc),
   * and the subsequent bigint_reserve(.., mw) would then grow
   * 4→mw (a second realloc).  Reserving first collapses the two
   * reallocs into one. */
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

int wigner3j_max_factorial(int tj1, int tj2, int tj3, int tm1, int tm2,
                           int tm3) {
  /* The largest factorial referenced in wigner3j_exact is the Δ
   * denominator (j1+j2+j3+1)!.  All m-dependent and Racah-sum
   * factorials have arguments <= (j1+j2+j3)/2. */
  (void)tm1;
  (void)tm2;
  (void)tm3;
  return (tj1 + tj2 + tj3) / 2 + 1;
}

/* The public wrappers compute through the cached wignernj_exact_t held in
 * the thread's scratch.  wigner3j_exact resets it (preserving bigint
 * capacities) before each call, so on every call after the first the
 * output bigints reuse their previously-allocated buffers and no
 * allocation happens.  wignernj_exact_free is therefore not called -- the
 * scratch's destructor takes care of it on thread exit. */

float wigner3j_f(int tj1, int tj2, int tj3, int tm1, int tm2, int tm3) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  float result;
  wigner3j_exact(tj1, tj2, tj3, tm1, tm2, tm3, &s->exact);
  result = wignernj_exact_to_float(&s->exact);
  wignernj_scratch_relinquish(s);
  return result;
}

double wigner3j(int tj1, int tj2, int tj3, int tm1, int tm2, int tm3) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  double result;
  wigner3j_exact(tj1, tj2, tj3, tm1, tm2, tm3, &s->exact);
  result = wignernj_exact_to_double(&s->exact);
  wignernj_scratch_relinquish(s);
  return result;
}

long double wigner3j_l(int tj1, int tj2, int tj3, int tm1, int tm2, int tm3) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  long double result;
  wigner3j_exact(tj1, tj2, tj3, tm1, tm2, tm3, &s->exact);
  result = wignernj_exact_to_long_double(&s->exact);
  wignernj_scratch_relinquish(s);
  return result;
}

#ifdef WIGNERNJ_HAVE_QUADMATH
__float128 wigner3j_q(int tj1, int tj2, int tj3, int tm1, int tm2, int tm3) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  __float128 result;
  wigner3j_exact(tj1, tj2, tj3, tm1, tm2, tm3, &s->exact);
  result = wignernj_exact_to_float128(&s->exact);
  wignernj_scratch_relinquish(s);
  return result;
}
#endif

#ifdef WIGNERNJ_HAVE_MPFR
void wigner3j_mpfr(mpfr_t rop, int tj1, int tj2, int tj3, int tm1, int tm2,
                   int tm3, mpfr_rnd_t rnd) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  wigner3j_exact(tj1, tj2, tj3, tm1, tm2, tm3, &s->exact);
  wignernj_exact_to_mpfr(rop, &s->exact, rnd);
  wignernj_scratch_relinquish(s);
}
#endif
