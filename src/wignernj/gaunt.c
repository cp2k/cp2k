/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * Gaunt coefficients: integral of three complex spherical harmonics.
 *
 * Phase convention:
 *   Y_l^m is taken in the Condon–Shortley phase, i.e.
 *     Y_l^m(θ,φ) = (-1)^m sqrt[(2l+1)/(4π) · (l-m)!/(l+m)!]
 *                  · P_l^m(cos θ) · e^{i m φ}    (m ≥ 0),
 *     Y_l^{-m}    = (-1)^m (Y_l^m)*,
 *   the convention of Edmonds (1957) and Varshalovich et al. (1988).
 *   The closed-form expression below in terms of two 3j symbols holds
 *   only under this phase convention; users of a different phase must
 *   adjust at the call site.
 *
 * Definition:
 *   G(l1,m1,l2,m2,l3,m3)
 *       = integral Y_{l1}^{m1}(Ω) Y_{l2}^{m2}(Ω) Y_{l3}^{m3}(Ω) dΩ
 *       = sqrt[(2l1+1)(2l2+1)(2l3+1) / (4π)]
 *         * (l1 l2 l3; 0 0 0) * (l1 l2 l3; m1 m2 m3)
 *
 * All l arguments must be non-negative integers; m must satisfy
 * m1+m2+m3=0 and |mi|<=li.  Arguments passed as tl=2*l, tm=2*m.
 *
 * Exact arithmetic path:
 *   The product of the two 3j symbols has the form
 *     (-1)^m3 * sqrt[outer_combined] * R0 * Rm
 *   where outer_combined = [Δ(l1,l2,l3)]² * (l1!)²(l2!)²(l3!)²
 *                          * (l1+m1)!(l1-m1)!(l2+m2)!(l2-m2)!(l3+m3)!(l3-m3)!
 *                          * (2l1+1)(2l2+1)(2l3+1) / 4
 *   and R0, Rm are the Racah integer sums for m=(0,0,0) and m=(m1,m2,m3).
 *   All factorial/integer factors are prime-factored exactly; only 1/sqrt(π)
 *   remains as a floating-point factor at the final conversion step.
 *
 * Phase derivation:
 *   phase_0 = (-1)^(l1-l2),   phase_m = (-1)^(l1-l2-m3)
 *   combined = (-1)^(2(l1-l2)-m3) = (-1)^(-m3) = (-1)^m3.
 */
#include "pfrac.h"
#include "primes.h"
#include "scratch.h"
#include "wignernj.h"
#include "wignernj_exact.h"
#include <math.h>   /* sqrtf, sqrt, sqrtl, acosl */
#include <stdlib.h> /* abs */

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

/* ── helpers ─────────────────────────────────────────────────────────────── */

/* Per-precision finish step: convert the exact tuple to T and divide by
 * sqrt(π).  Shared between {complex, real}-Y Gaunt so that the per-precision
 * π constants live in one place per precision rather than being duplicated
 * across both variants. */
static inline float gaunt_finish_f(const wignernj_exact_t *e) {
  return wignernj_exact_to_float(e) / sqrtf((float)M_PI);
}
static inline double gaunt_finish_d(const wignernj_exact_t *e) {
  return wignernj_exact_to_double(e) / sqrt(M_PI);
}
static inline long double gaunt_finish_l(const wignernj_exact_t *e) {
  return wignernj_exact_to_long_double(e) / sqrtl(acosl(-1.0L));
}

#ifdef WIGNERNJ_HAVE_QUADMATH
static inline __float128 gaunt_finish_q(const wignernj_exact_t *e) {
  return wignernj_exact_to_float128(e) / sqrtq(M_PIq);
}
#endif

#ifdef WIGNERNJ_HAVE_MPFR
/* In-place /sqrt(π) at the precision of rop.  No-op when rop is zero
 * (selection rules already fired). */
static void gaunt_finish_mpfr(mpfr_t rop, mpfr_rnd_t rnd) {
  mpfr_t pi;
  if (mpfr_zero_p(rop))
    return;
  mpfr_init2(pi, mpfr_get_prec(rop));
  mpfr_const_pi(pi, rnd);
  mpfr_sqrt(pi, pi, rnd);
  mpfr_div(rop, rop, pi, rnd);
  mpfr_clear(pi);
}
#endif

static int gaunt_selection_rules(int tl1, int tm1, int tl2, int tm2, int tl3,
                                 int tm3) {
  if (tm1 + tm2 + tm3 != 0)
    return 0;
  if (abs(tm1) > tl1 || abs(tm2) > tl2 || abs(tm3) > tl3)
    return 0;
  if ((tl1 + tl2 + tl3) & 1)
    return 0;
  if (tl3 < abs(tl1 - tl2) || tl3 > tl1 + tl2)
    return 0;
  return 1;
}

/*
 * Summation bounds for the 3j Racah sum (same logic as bounds_3j in
 * wigner3j.c).
 */
static void gaunt_racah_bounds(int tj1, int tj2, int tj3, int tm1, int tm2,
                               int *s_min, int *s_max) {
  int a = (tj2 - tj3 - tm1) / 2;
  int b = (tj1 + tm2 - tj3) / 2;
  int c;
  *s_min = (a > b) ? a : b;
  if (*s_min < 0)
    *s_min = 0;
  a = (tj1 + tj2 - tj3) / 2;
  b = (tj1 - tm1) / 2;
  c = (tj2 + tm2) / 2;
  *s_max = a;
  if (b < *s_max)
    *s_max = b;
  if (c < *s_max)
    *s_max = c;
}

/*
 * Racah integer sum for 3j(j1,j2,j3; m1,m2,m3), no outer sqrt factor.
 * Fills sum_out (LCM-scaled magnitude), sum_sign_out, lcm_exp[0..*lcm_max_io)
 * and the new *lcm_max_io upper bound.  lcm_exp must be zeroed by the caller.
 */
static void gaunt_3j_racah_sum(int tj1, int tj2, int tj3, int tm1, int tm2,
                               bigint_t *sum_out, int *sum_sign_out,
                               int *lcm_exp, int *lcm_max_io, bigint_t *sum_pos,
                               bigint_t *sum_neg, bigint_t *scaled,
                               wignernj_scratch_t *scratch, bigint_ws_t *ws) {
  int s, s_min, s_max, pi;
  int lcm_max = 0;
  pfrac_t *term;

  /* Closed-form fast path for the (j1 j2 j3; 0 0 0) sub-sum that
   * gaunt_exact requests once per call.  Substituting m=0 into the
   * Racah single sum collapses it to the closed rational
   *   S_3j(0,0,0) = (-1)^(g - j1 + j2)
   *               * g! / [j1! j2! j3! (g-j1)! (g-j2)! (g-j3)!]
   * with g = (j1+j2+j3)/2 (the symbol vanishes when 2g = j1+j2+j3
   * is odd).  The exponent of each prime in this rational is
   * computed once and split between lcm_exp[] (the negative-
   * exponent / denominator part that the caller folds into D_int)
   * and sum_out (the positive-exponent / numerator part that goes
   * into the running Σ), exactly as gaunt_exact does for the
   * full-Racah-loop output of this routine.  Skips the entire
   * Pass-1/Pass-2 loop below. */
  if (tm1 == 0 && tm2 == 0) {
    int j1 = tj1 / 2, j2 = tj2 / 2, j3 = tj3 / 2;
    int L = j1 + j2 + j3;
    if (L & 1) {
      bigint_set_zero(sum_out);
      *sum_sign_out = 1;
      *lcm_max_io = 0;
      return;
    }
    int g = L / 2;
    pfrac_t *closed = &scratch->pfracs[1];
    pfrac_zero(closed);
    pfrac_mul_factorial(closed, g);
    pfrac_div_factorial(closed, j1);
    pfrac_div_factorial(closed, j2);
    pfrac_div_factorial(closed, j3);
    pfrac_div_factorial(closed, g - j1);
    pfrac_div_factorial(closed, g - j2);
    pfrac_div_factorial(closed, g - j3);

    /* Split the signed prime exponents into denominator
     * (lcm_exp[]) and numerator (the entries that survive in
     * `closed` after the negative ones are zeroed); then update
     * lcm_max to reflect the highest prime touched and
     * materialise the numerator into sum_out. */
    int lcm_max_local = 0;
    for (pi = 0; pi < closed->max_idx; pi++) {
      int e = closed->exp[pi];
      if (e < 0) {
        int neg = -e;
        if (neg > lcm_exp[pi])
          lcm_exp[pi] = neg;
        closed->exp[pi] = 0;
        if (pi + 1 > lcm_max_local)
          lcm_max_local = pi + 1;
      } else if (e > 0) {
        if (pi + 1 > lcm_max_local)
          lcm_max_local = pi + 1;
      }
    }
    bigint_set_u64(sum_out, 1);
    pfrac_bigint_mul_prime_pow_array(sum_out, closed->exp, closed->max_idx, ws);
    *sum_sign_out = ((g - j1 + j2) & 1) ? -1 : 1;
    *lcm_max_io = lcm_max_local;
    return;
  }

  gaunt_racah_bounds(tj1, tj2, tj3, tm1, tm2, &s_min, &s_max);

  if (s_min > s_max) {
    bigint_set_zero(sum_out);
    *sum_sign_out = 1;
    *lcm_max_io = 0;
    return;
  }

  bigint_set_zero(sum_pos);
  bigint_set_zero(sum_neg);

  /* Pass 1: build each term pfrac, fold into LCM.  Pass 2 uses the
   * incremental ratio recurrence and only consults the pfrac cache
   * once (as a seed for s = s_min), so two slots suffice: slot 0
   * holds the s_min term and slot 1 is a single rolling scratch
   * reused for every s > s_min.  The first iteration is peeled
   * explicitly so the inner loop's term pointer is loop-invariant
   * (a conditional `term = &scratch->terms[s == s_min ? 0 : 1]`
   * inside the loop body measurably regressed small-j calls on the
   * x86-64 divq backend, where per-call overhead is otherwise
   * minimal). */
  wignernj_scratch_terms_reserve(scratch, 2);

  /* Peeled iteration s = s_min, written to slot 0 (Pass-2 seed). */
  term = &scratch->terms[0];
  pfrac_zero(term);
  pfrac_mul_factorial(term, s_min);
  pfrac_mul_factorial(term, (tj1 + tj2 - tj3) / 2 - s_min);
  pfrac_mul_factorial(term, (tj1 - tm1) / 2 - s_min);
  pfrac_mul_factorial(term, (tj2 + tm2) / 2 - s_min);
  pfrac_mul_factorial(term, (tj3 - tj2 + tm1) / 2 + s_min);
  pfrac_mul_factorial(term, (tj3 - tj1 - tm2) / 2 + s_min);
  for (pi = 0; pi < term->max_idx; pi++) {
    if (term->exp[pi] > lcm_exp[pi])
      lcm_exp[pi] = term->exp[pi];
  }
  if (term->max_idx > lcm_max)
    lcm_max = term->max_idx;

  /* Remaining iterations s = s_min+1 .. s_max, all overwriting slot 1. */
  term = &scratch->terms[1];
  for (s = s_min + 1; s <= s_max; s++) {
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
    if (term->max_idx > lcm_max)
      lcm_max = term->max_idx;
  }

  /* Pass 2: incremental walk of LCM-scaled term values, mirroring the
   * 3j Pass 2 in wigner3j.c.  The Gaunt 3j Racah denominator
   *   D_s = s! · (a-s)! · (b-s)! · (c-s)! · (d+s)! · (e+s)!
   * with
   *   a = (tj1+tj2-tj3)/2,  b = (tj1-tm1)/2,  c = (tj2+tm2)/2,
   *   d = (tj3-tj2+tm1)/2,  e = (tj3-tj1-tm2)/2,
   * is structurally identical to the 3j case, so the same triple-
   * product ratio recurrence
   *   scaled_{s+1} · (s+1)(d+s+1)(e+s+1) = scaled_s · (a-s)(b-s)(c-s)
   * applies, with the same MAX_FACTORIAL_ARG ≤ 2642245 safe range
   * for the batched uint64 mul/div primitives.  See wigner3j.c for
   * the full derivation. */
#if MAX_FACTORIAL_ARG > 2642245
#error                                                                         \
    "MAX_FACTORIAL_ARG too large for batched Gaunt-Pass-2 ratio update (triple product would overflow uint64); regenerate the prime table with a smaller --limit, or replace the batched mul/div in this loop with six single-factor calls"
#endif
  {
    uint64_t a = (uint64_t)((tj1 + tj2 - tj3) / 2);
    uint64_t b = (uint64_t)((tj1 - tm1) / 2);
    uint64_t c = (uint64_t)((tj2 + tm2) / 2);
    uint64_t d = (uint64_t)((tj3 - tj2 + tm1) / 2);
    uint64_t e = (uint64_t)((tj3 - tj1 - tm2) / 2);

    /* Seed: scaled_{s_min} = LCM/D_{s_min} via the prime-power
     * helper, using the cached pfrac for the first term.  This is
     * the only call to that helper per gaunt_3j_racah_sum; every
     * subsequent term is reached by ratio update below. */
    pfrac_lcm_scaled_product(scaled, lcm_exp, scratch->terms[0].exp, -1,
                             lcm_max, ws);
    if ((s_min & 1) == 0)
      bigint_add(sum_pos, sum_pos, scaled);
    else
      bigint_add(sum_neg, sum_neg, scaled);

    /* Step s -> s+1.  At s = s_max we stop without applying the
     * ratio: a-s would be 0 there (one of the upper-bound factors
     * is zero by the s_max definition), and the next scaled is 0
     * by the well-defined limit, but we never need it. */
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

  *sum_sign_out = bigint_sub_signed(sum_out, sum_pos, sum_neg);
  *lcm_max_io = lcm_max;
}

/* ── exact Gaunt (without the 1/sqrt(π) factor) ─────────────────────────── */

/*
 * Fills *out such that G(l1,m1,l2,m2,l3,m3)
 *   = wignernj_exact_to_T(out) / sqrt(π).
 *
 * When is_zero is set, the Gaunt coefficient is zero and no further work
 * need be done.
 */
static void gaunt_exact(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3,
                        wignernj_exact_t *out) {
  int pi;
  wignernj_scratch_t *scratch;
  bigint_ws_t *ws;
  pfrac_t *outer;
  int *lcm0, *lcmm;
  bigint_t *sum0, *summ, *sum_pos, *sum_neg, *scaled;
  int ss0, ssm;
  size_t mw, mw_prod;
  int lcm0_max, lcmm_max;

  wignernj_exact_reset(out);

  if (!gaunt_selection_rules(tl1, tm1, tl2, tm2, tl3, tm3)) {
    out->is_zero = 1;
    return;
  }

  /* Phase = (-1)^m3 */
  out->sign = (((tm3 / 2) & 1) == 0) ? 1 : -1;

  /* Largest factorial argument bounded by (tl1+tl2+tl3)/2 + 1.  The
   * combined outer pfrac and the cross product of the two Racah sums
   * push the largest bigint to ~2x that, so allocate accordingly. */
  mw = bigint_words_for_factorial((tl1 + tl2 + tl3) / 2 + 5);
  mw_prod = 2 * mw;

  /* Acquire the thread's cached scratch.  Gaunt uses 5 bigints
   * (sum0, summ + sum_pos/sum_neg/scaled inside the inner Racah
   * sum), 1 pfrac (outer; the per-term pfracs live in the cached
   * terms array), and 2 lcm arrays (lcm0, lcmm). */
  scratch = wignernj_scratch_acquire();
  ws = &scratch->ws;
  outer = &scratch->pfracs[0];
  lcm0 = scratch->lcm_exp[0];
  lcmm = scratch->lcm_exp[1];
  sum0 = &scratch->bigints[0];
  summ = &scratch->bigints[1];
  sum_pos = &scratch->bigints[2];
  sum_neg = &scratch->bigints[3];
  scaled = &scratch->bigints[4];

  bigint_ws_reserve(ws, mw_prod);
  pfrac_zero(outer);
  wignernj_scratch_lcm_clear(scratch, 0);
  wignernj_scratch_lcm_clear(scratch, 1);
  bigint_set_zero(sum0);
  bigint_reserve(sum0, mw);
  bigint_set_zero(summ);
  bigint_reserve(summ, mw);
  bigint_set_zero(sum_pos);
  bigint_reserve(sum_pos, mw);
  bigint_set_zero(sum_neg);
  bigint_reserve(sum_neg, mw);
  bigint_set_zero(scaled);
  bigint_reserve(scaled, mw);

  /*
   * Combined outer pfrac (argument under the outer sqrt, without π):
   *
   *   [Δ(l1,l2,l3)]²:
   *     each factor contributes (l1+l2-l3)!(l1-l2+l3)!(-l1+l2+l3)! /
   * (l1+l2+l3+1)! called twice → all exponents even → folds entirely into
   * int_num/int_den
   *
   *   m=0 factor:  (l1!)²(l2!)²(l3!)²
   *
   *   general-m factor: (l1+m1)!(l1-m1)!(l2+m2)!(l2-m2)!(l3+m3)!(l3-m3)!
   *
   *   normalization: (2l1+1)(2l2+1)(2l3+1)/4
   */

  /* Δ(l1,l2,l3) — first copy */
  pfrac_mul_factorial(outer, (tl1 + tl2 - tl3) / 2);
  pfrac_mul_factorial(outer, (tl1 - tl2 + tl3) / 2);
  pfrac_mul_factorial(outer, (-tl1 + tl2 + tl3) / 2);
  pfrac_div_factorial(outer, (tl1 + tl2 + tl3) / 2 + 1);
  /* Δ(l1,l2,l3) — second copy (makes [Δ]² rational) */
  pfrac_mul_factorial(outer, (tl1 + tl2 - tl3) / 2);
  pfrac_mul_factorial(outer, (tl1 - tl2 + tl3) / 2);
  pfrac_mul_factorial(outer, (-tl1 + tl2 + tl3) / 2);
  pfrac_div_factorial(outer, (tl1 + tl2 + tl3) / 2 + 1);

  /* (l1!)²(l2!)²(l3!)² from the m=0 3j */
  pfrac_mul_factorial(outer, tl1 / 2);
  pfrac_mul_factorial(outer, tl1 / 2);
  pfrac_mul_factorial(outer, tl2 / 2);
  pfrac_mul_factorial(outer, tl2 / 2);
  pfrac_mul_factorial(outer, tl3 / 2);
  pfrac_mul_factorial(outer, tl3 / 2);

  /* (l1+m1)!(l1-m1)!(l2+m2)!(l2-m2)!(l3+m3)!(l3-m3)! from the general-m 3j */
  pfrac_mul_factorial(outer, (tl1 + tm1) / 2);
  pfrac_mul_factorial(outer, (tl1 - tm1) / 2);
  pfrac_mul_factorial(outer, (tl2 + tm2) / 2);
  pfrac_mul_factorial(outer, (tl2 - tm2) / 2);
  pfrac_mul_factorial(outer, (tl3 + tm3) / 2);
  pfrac_mul_factorial(outer, (tl3 - tm3) / 2);

  /* (2l1+1)(2l2+1)(2l3+1) / 4 : integer factors folded in via prime
   * factorization. g_primes[0] = 2, so subtracting 2 from exp[0] divides by 4 =
   * 2². Bump max_idx so the manual write to exp[0] is visible to
   * pfrac_to_sqrt_rational's max_idx-bounded sweep. */
  pfrac_mul_int(outer, tl1 + 1);
  pfrac_mul_int(outer, tl2 + 1);
  pfrac_mul_int(outer, tl3 + 1);
  outer->exp[0] -= 2;
  if (outer->max_idx < 1)
    outer->max_idx = 1;

  /* Racah integer sums for both 3j symbols */
  lcm0_max = 0;
  lcmm_max = 0;

  gaunt_3j_racah_sum(tl1, tl2, tl3, 0, 0, sum0, &ss0, lcm0, &lcm0_max, sum_pos,
                     sum_neg, scaled, scratch, ws);
  if (bigint_is_zero(sum0)) {
    out->is_zero = 1;
    goto cleanup;
  }

  gaunt_3j_racah_sum(tl1, tl2, tl3, tm1, tm2, summ, &ssm, lcmm, &lcmm_max,
                     sum_pos, sum_neg, scaled, scratch, ws);
  if (bigint_is_zero(summ)) {
    out->is_zero = 1;
    goto cleanup;
  }

  /* Product of the two Racah sums */
  bigint_reserve(&out->sum, mw_prod);
  bigint_mul_ws(&out->sum, sum0, summ, ws);
  out->sum_sign = ss0 * ssm;

  /* Convert outer pfrac */
  bigint_reserve(&out->int_num, mw_prod);
  bigint_reserve(&out->int_den, mw_prod);
  bigint_reserve(&out->sqrt_num, mw_prod);
  bigint_reserve(&out->sqrt_den, mw_prod);
  bigint_set_u64(&out->int_num, 1);
  bigint_set_u64(&out->int_den, 1);
  bigint_set_u64(&out->sqrt_num, 1);
  bigint_set_u64(&out->sqrt_den, 1);
  pfrac_to_sqrt_rational_ws(outer, &out->int_num, &out->int_den, &out->sqrt_num,
                            &out->sqrt_den, ws);

  /* Combined LCM denominator → int_den.  Per-prime exponent is
   * lcm0[pi] + lcmm[pi]; route through the uint64-batched
   * accumulator helper instead of one bigint_mul_prime_pow_ws call
   * per prime, mirroring the wignerXj_exact end-of-function paths. */
  {
    int lcm_union_max = (lcm0_max > lcmm_max) ? lcm0_max : lcmm_max;
    uint64_t acc = 1;
    for (pi = 0; pi < lcm_union_max; pi++) {
      int e = lcm0[pi] + lcmm[pi];
      if (e > 0)
        pfrac_mul_pow_into_acc(&out->int_den, &acc, (uint64_t)g_primes[pi], e,
                               ws);
    }
    if (acc > 1)
      bigint_mul_u64(&out->int_den, &out->int_den, acc);
  }

cleanup:
  wignernj_scratch_lcm_dirty(scratch, 0, lcm0_max);
  wignernj_scratch_lcm_dirty(scratch, 1, lcmm_max);
  wignernj_scratch_relinquish(scratch);
}

/* ── Real-spherical-harmonic Gaunt coefficient ──────────────────────────────
 *
 *   G^R(l1,m1,l2,m2,l3,m3)  =  ∫ S_{l1,m1} S_{l2,m2} S_{l3,m3} dΩ
 *
 * with the Condon–Shortley / Wikipedia convention for the real harmonics:
 *
 *   S_{l,0}    = Y_l^0
 *   S_{l,m>0}  = (1/√2)( Y_l^{-m} + (-1)^m Y_l^m )
 *   S_{l,m<0}  = (i/√2)( Y_l^{ m}  - (-1)^|m| Y_l^{-m} )
 *
 * Substituting into the integral expresses the real Gaunt as a small (≤ 2)
 * linear combination of complex Gaunts at the same (l1,l2,l3) but at
 * (s_1 |m_1|, s_2 |m_2|, s_3 |m_3|) for sign assignments s_i ∈ {±1} that
 * satisfy s_1|m_1| + s_2|m_2| + s_3|m_3| = 0.
 *
 * Two further reductions make this run at the cost of *one* complex-Gaunt
 * evaluation:
 *
 *   (a) The (≤ 2) valid sign tuples form a sign-flipped pair (s, -s).  The
 *       complex Gaunt is invariant under simultaneous sign flip of all
 *       m-arguments when l1+l2+l3 is even -- which is required for any
 *       non-zero Gaunt at all -- so both tuples evaluate to the same
 *       G^C.  We therefore compute G^C only once and multiply by the
 *       sum of the two tuple coefficients.
 *
 *   (b) The total coefficient is rational (or rational times √2) with
 *       small numerator and denominator (∈ {0, ±1, ±2} divided by
 *       √2^k where k ∈ {0,2,3} is the number of non-zero |m_i|), and
 *       can be absorbed exactly into the wignernj_exact_t pipeline by
 *       multiplying int_num/int_den/sqrt_den by small integers -- so
 *       last-bit accuracy is preserved.
 *
 * If the ≤ 2 tuple coefficients sum to zero, or if the parity of the
 * |m_i| does not admit a valid sign tuple, the real Gaunt is zero.
 */

/* The integer "tau" part of T(σ, s, a):  T = τ * φ / √2^[a > 0],
 * with φ = 1 (real) or φ = i.  Returns τ ∈ {±1}; sets *has_i to 1 iff φ = i.
 * Argument a is the integer m (not 2m); a >= 0 always. */
static int real_gaunt_T_tau(int sigma, int s, int a, int *has_i) {
  if (a == 0) {
    *has_i = 0;
    return 1;
  }
  int parity = a & 1; /* (-1)^a = +1 if a even, -1 if a odd */
  if (sigma > 0) {    /* m_i > 0 */
    *has_i = 0;
    if (s > 0)
      return parity ? -1 : +1; /*  T = (-1)^a/√2  */
    else
      return +1; /*  T =      1/√2  */
  } else {       /* m_i < 0 */
    *has_i = 1;
    if (s > 0)
      return parity ? +1 : -1; /*  T = -i(-1)^a/√2  */
    else
      return +1; /*  T = +i      /√2  */
  }
}

static void gaunt_real_exact(int tl1, int tm1, int tl2, int tm2, int tl3,
                             int tm3, wignernj_exact_t *out) {
  wignernj_exact_reset(out);

  /* ℓ and m must be integer-valued, i.e., tl, tm even. */
  if ((tl1 | tl2 | tl3 | tm1 | tm2 | tm3) & 1) {
    out->is_zero = 1;
    return;
  }
  if (tl1 < 0 || tl2 < 0 || tl3 < 0) {
    out->is_zero = 1;
    return;
  }

  int ta[3] = {tm1 < 0 ? -tm1 : tm1, tm2 < 0 ? -tm2 : tm2,
               tm3 < 0 ? -tm3 : tm3};
  int sigma[3] = {(tm1 > 0) - (tm1 < 0), (tm2 > 0) - (tm2 < 0),
                  (tm3 > 0) - (tm3 < 0)};
  int a[3] = {ta[0] / 2, ta[1] / 2, ta[2] / 2};
  int tl[3] = {tl1, tl2, tl3};

  /* The phase i^{n_-} (n_- = number of m_i < 0) is the same for every
   * sign tuple.  If n_- is odd the contribution is purely imaginary,
   * hence the real Gaunt vanishes by the reality of the integral. */
  int n_minus = (sigma[0] < 0) + (sigma[1] < 0) + (sigma[2] < 0);
  if (n_minus & 1) {
    out->is_zero = 1;
    return;
  }
  int phase_sign = (n_minus == 2) ? -1 : +1;

  /* Sweep the (≤ 8) sign-tuple candidates.  For each that satisfies
   * Σ s_i a_i = 0, accumulate phase_sign * τ_1 τ_2 τ_3 into total_tau,
   * and remember the first such tuple as the representative for the
   * complex-Gaunt evaluation below. */
  int total_tau = 0;
  int has_repr = 0;
  int s_repr[3] = {1, 1, 1};
  for (int s1_idx = 0; s1_idx < 2; s1_idx++) {
    if (s1_idx == 1 && a[0] == 0)
      break;
    int s1 = (s1_idx == 0) ? +1 : -1;
    for (int s2_idx = 0; s2_idx < 2; s2_idx++) {
      if (s2_idx == 1 && a[1] == 0)
        break;
      int s2 = (s2_idx == 0) ? +1 : -1;
      for (int s3_idx = 0; s3_idx < 2; s3_idx++) {
        if (s3_idx == 1 && a[2] == 0)
          break;
        int s3 = (s3_idx == 0) ? +1 : -1;
        if (s1 * a[0] + s2 * a[1] + s3 * a[2] != 0)
          continue;

        int s_tup[3] = {s1, s2, s3};
        int tau_prod = 1;
        int has_i_dummy;
        for (int i = 0; i < 3; i++)
          tau_prod *= real_gaunt_T_tau(sigma[i], s_tup[i], a[i], &has_i_dummy);
        total_tau += phase_sign * tau_prod;
        if (!has_repr) {
          s_repr[0] = s1;
          s_repr[1] = s2;
          s_repr[2] = s3;
          has_repr = 1;
        }
      }
    }
  }
  if (!has_repr || total_tau == 0) {
    out->is_zero = 1;
    return;
  }

  /* Evaluate the complex Gaunt at the representative tuple. */
  int p_repr[3] = {s_repr[0] * ta[0], s_repr[1] * ta[1], s_repr[2] * ta[2]};
  gaunt_exact(tl[0], p_repr[0], tl[1], p_repr[1], tl[2], p_repr[2], out);
  if (out->is_zero)
    return;

  /* Multiply by  total_tau / √2^k  with k = #{i : a_i > 0}.
   * Even k:  pure rational, fold into int_num/int_den.
   * Odd k :  fold the integer part into int_den, the residual 1/√2 into
   *          sqrt_den (representing an extra factor of √(1/2)). */
  int abs_tau = total_tau < 0 ? -total_tau : total_tau;
  if (total_tau < 0)
    out->sign = -out->sign;
  if (abs_tau != 1)
    bigint_mul_u64(&out->int_num, &out->int_num, (uint64_t)abs_tau);
  int k = (a[0] > 0) + (a[1] > 0) + (a[2] > 0);
  int half_k = k / 2;
  for (int i = 0; i < half_k; i++)
    bigint_mul_u64(&out->int_den, &out->int_den, 2);
  if (k & 1)
    bigint_mul_u64(&out->sqrt_den, &out->sqrt_den, 2);
}

/* ── public API ──────────────────────────────────────────────────────────── */

int gaunt_max_factorial(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3) {
  /* The combined outer pfrac and the inner 3j Racah sums all touch
   * factorials whose argument is bounded by (l1+l2+l3+1)! (the
   * triangle-Δ denominator).  m-dependent factorials and the inner
   * Racah sums are <= (l1+l2+l3)/2. */
  (void)tm1;
  (void)tm2;
  (void)tm3;
  return (tl1 + tl2 + tl3) / 2 + 1;
}

int gaunt_real_max_factorial(int tl1, int tm1, int tl2, int tm2, int tl3,
                             int tm3) {
  /* gaunt_real_exact dispatches to one complex-Gaunt evaluation at
   * (sign-adjusted) m's; same factorial bound. */
  return gaunt_max_factorial(tl1, tm1, tl2, tm2, tl3, tm3);
}

float gaunt_f(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  float result;
  gaunt_exact(tl1, tm1, tl2, tm2, tl3, tm3, &s->exact);
  result = gaunt_finish_f(&s->exact);
  wignernj_scratch_relinquish(s);
  return result;
}

double gaunt(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  double result;
  gaunt_exact(tl1, tm1, tl2, tm2, tl3, tm3, &s->exact);
  result = gaunt_finish_d(&s->exact);
  wignernj_scratch_relinquish(s);
  return result;
}

long double gaunt_l(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  long double result;
  gaunt_exact(tl1, tm1, tl2, tm2, tl3, tm3, &s->exact);
  result = gaunt_finish_l(&s->exact);
  wignernj_scratch_relinquish(s);
  return result;
}

#ifdef WIGNERNJ_HAVE_QUADMATH
__float128 gaunt_q(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  __float128 result;
  gaunt_exact(tl1, tm1, tl2, tm2, tl3, tm3, &s->exact);
  result = gaunt_finish_q(&s->exact);
  wignernj_scratch_relinquish(s);
  return result;
}
#endif

#ifdef WIGNERNJ_HAVE_MPFR
#include "wignernj_mpfr.h"
void gaunt_mpfr(mpfr_t rop, int tl1, int tm1, int tl2, int tm2, int tl3,
                int tm3, mpfr_rnd_t rnd) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  gaunt_exact(tl1, tm1, tl2, tm2, tl3, tm3, &s->exact);
  wignernj_exact_to_mpfr(rop, &s->exact, rnd);
  wignernj_scratch_relinquish(s);
  gaunt_finish_mpfr(rop, rnd);
}
#endif

/* ── Real-spherical-harmonic Gaunt: public API ──────────────────────────── */

float gaunt_real_f(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  float result;
  gaunt_real_exact(tl1, tm1, tl2, tm2, tl3, tm3, &s->exact);
  result = gaunt_finish_f(&s->exact);
  wignernj_scratch_relinquish(s);
  return result;
}

double gaunt_real(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  double result;
  gaunt_real_exact(tl1, tm1, tl2, tm2, tl3, tm3, &s->exact);
  result = gaunt_finish_d(&s->exact);
  wignernj_scratch_relinquish(s);
  return result;
}

long double gaunt_real_l(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  long double result;
  gaunt_real_exact(tl1, tm1, tl2, tm2, tl3, tm3, &s->exact);
  result = gaunt_finish_l(&s->exact);
  wignernj_scratch_relinquish(s);
  return result;
}

#ifdef WIGNERNJ_HAVE_QUADMATH
__float128 gaunt_real_q(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  __float128 result;
  gaunt_real_exact(tl1, tm1, tl2, tm2, tl3, tm3, &s->exact);
  result = gaunt_finish_q(&s->exact);
  wignernj_scratch_relinquish(s);
  return result;
}
#endif

#ifdef WIGNERNJ_HAVE_MPFR
void gaunt_real_mpfr(mpfr_t rop, int tl1, int tm1, int tl2, int tm2, int tl3,
                     int tm3, mpfr_rnd_t rnd) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  gaunt_real_exact(tl1, tm1, tl2, tm2, tl3, tm3, &s->exact);
  wignernj_exact_to_mpfr(rop, &s->exact, rnd);
  wignernj_scratch_relinquish(s);
  gaunt_finish_mpfr(rop, rnd);
}
#endif
