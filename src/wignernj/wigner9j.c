/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * Wigner 9j symbol via exact prime-factorization arithmetic.
 *
 * Formula (Racah W-coefficient representation):
 *   {j11 j12 j13}
 *   {j21 j22 j23} = (-1)^(tj13+tj22+tj31) * sum_k (2k+1)
 *   {j31 j32 j33}
 *       * {j11 j21 j31} * {j11 j12 j13} * {j22 j21 j23}
 *         {j32 j33  k }   {j23 j33  k }   { k  j12 j32}
 *
 * (Phase exponent is in 2j units; tj values are always integers so the phase
 * is always ±1.)
 *
 * Exactness argument: the product of the three 6j symbols at each k is
 *   sqrt(C) * rational(k) * R1(k) * R2(k) * R3(k)
 * where sqrt(C) is formed by the six k-independent Delta factors (two per 6j):
 *   Δ(j11,j21,j31), Δ(j32,j33,j31), Δ(j11,j12,j13), Δ(j23,j33,j13),
 *   Δ(j22,j21,j23), Δ(j22,j12,j32)
 * and rational(k) = Δ²(j11,j33,k)·Δ²(j21,j32,k)·Δ²(j12,j23,k) is purely
 * rational (each k-dependent Δ appears exactly twice across the three 6j
 * symbols, making it a perfect square and thus rational).  Ri(k) is the
 * Racah integer sum for the i-th 6j.  The 9j therefore has the form
 * sqrt(C) × exact_integer, enabling the same exact pipeline as 3j and 6j.
 *
 * k range: max(|j11-j33|,|j21-j32|,|j12-j23|) <= k <=
 * min(j11+j33,j21+j32,j12+j23)
 */
#include "pfrac.h"
#include "primes.h"
#include "scratch.h"
#include "wignernj.h"
#include "wignernj_exact.h"
#include <stdlib.h>
#include <string.h>

/* ── helpers ─────────────────────────────────────────────────────────────── */

static int triangle_ok(int ta, int tb, int tc) {
  if ((ta + tb + tc) & 1)
    return 0;
  if (tc < abs(ta - tb))
    return 0;
  if (tc > ta + tb)
    return 0;
  return 1;
}

static int max3(int a, int b, int c) {
  int m = a;
  if (b > m)
    m = b;
  if (c > m)
    m = c;
  return m;
}

static int min3(int a, int b, int c) {
  int m = a;
  if (b < m)
    m = b;
  if (c < m)
    m = c;
  return m;
}

static int selection_rules_9j(int tj11, int tj12, int tj13, int tj21, int tj22,
                              int tj23, int tj31, int tj32, int tj33) {
  return triangle_ok(tj11, tj12, tj13) && triangle_ok(tj21, tj22, tj23) &&
         triangle_ok(tj31, tj32, tj33) && triangle_ok(tj11, tj21, tj31) &&
         triangle_ok(tj12, tj22, tj32) && triangle_ok(tj13, tj23, tj33);
}

/* Δ(a,b,c) = sqrt[(a+b-c)!(a-b+c)!(-a+b+c)!/(a+b+c+1)!]; add to pfrac. */
static void add_delta_sqrt(pfrac_t *f, int ta, int tb, int tc) {
  pfrac_mul_factorial(f, (ta + tb - tc) / 2);
  pfrac_mul_factorial(f, (ta - tb + tc) / 2);
  pfrac_mul_factorial(f, (-ta + tb + tc) / 2);
  pfrac_div_factorial(f, (ta + tb + tc) / 2 + 1);
}

/*
 * Compute the Racah integer sum for 6j {a b c; d e f} (2j-integer args).
 * Fills sum_out with the LCM-scaled magnitude, sum_sign_out with its sign,
 * and lcm_exp[0..g_nprimes-1] with the LCM denominator exponents.
 * lcm_exp must be a caller-allocated array of g_nprimes ints (need not be
 * zeroed; this function initialises it).
 */
/* Inner Racah sum for one of the three 6j symbols.  Caller-owned scratch
 * bigints (sum_pos, sum_neg, scaled) are reused across all calls of one 9j
 * evaluation, eliminating per-call bigint allocations. */
/* lcm_max_idx is in/out: caller passes the previous call's value (or 0 on
 * first use); function zeros lcm_exp[0..*lcm_max_idx_io) and overwrites it
 * with the new bound on return.  This avoids the 2.2 KB memset on every
 * inner-6j call. */
static void racah_6j_sum(int tj1, int tj2, int tj3, int tj4, int tj5, int tj6,
                         bigint_t *sum_out, int *sum_sign_out, int *lcm_exp,
                         int *lcm_max_idx_io, bigint_t *sum_pos,
                         bigint_t *sum_neg, bigint_t *scaled,
                         wignernj_scratch_t *scratch, bigint_ws_t *ws) {
  int s, s_min, s_max, pi;
  int lcm_max_idx = 0;
  int n_terms;
  pfrac_t *term;

  {
    int a = (tj1 + tj2 + tj3) / 2, b = (tj1 + tj5 + tj6) / 2;
    int c = (tj4 + tj2 + tj6) / 2, d = (tj4 + tj5 + tj3) / 2;
    s_min = a;
    if (b > s_min)
      s_min = b;
    if (c > s_min)
      s_min = c;
    if (d > s_min)
      s_min = d;
  }
  {
    int a = (tj1 + tj2 + tj4 + tj5) / 2;
    int b = (tj2 + tj3 + tj5 + tj6) / 2;
    int c = (tj1 + tj3 + tj4 + tj6) / 2;
    s_max = a;
    if (b < s_max)
      s_max = b;
    if (c < s_max)
      s_max = c;
  }

  /* Clear only the [0, prev_max) range that the previous call could have
   * dirtied; the suffix beyond that bound is still zeroed from xcalloc or
   * from an earlier zero-fill. */
  if (*lcm_max_idx_io > 0)
    memset(lcm_exp, 0, (size_t)(*lcm_max_idx_io) * sizeof(int));

  if (s_min > s_max) {
    bigint_set_zero(sum_out);
    *sum_sign_out = 1;
    *lcm_max_idx_io = 0;
    return;
  }

  bigint_set_zero(sum_pos);
  bigint_set_zero(sum_neg);

  /* Pass 1: build each term pfrac once into the cache, find LCM */
  n_terms = s_max - s_min + 1;
  wignernj_scratch_terms_reserve(scratch, n_terms);
  for (s = s_min; s <= s_max; s++) {
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
    for (pi = 0; pi < term->max_idx; pi++) {
      int neg = -term->exp[pi];
      if (neg > lcm_exp[pi])
        lcm_exp[pi] = neg;
    }
    if (term->max_idx > lcm_max_idx)
      lcm_max_idx = term->max_idx;
  }

  /* Pass 2: incremental walk via the 6j Racah-term ratio recurrence
   * (see the analogous block in wigner6j.c for the derivation and
   * the static_assert bound on MAX_FACTORIAL_ARG; the recurrence is
   * identical here since racah_6j_sum is structurally a 6j inner
   * sum).  Replaces the per-term O(pi(j)) prime-power expansion
   * with a single uint64-batched mul/div per term, called O(j)
   * times per outer-k iteration of the 9j accumulator. */
#if MAX_FACTORIAL_ARG > 32767
#error                                                                         \
    "MAX_FACTORIAL_ARG too large for batched 9j inner-Racah ratio update; the wigner6j.c bound applies here too"
#endif
  {
    uint64_t a1 = (uint64_t)((tj1 + tj2 + tj3) / 2);
    uint64_t a2 = (uint64_t)((tj1 + tj5 + tj6) / 2);
    uint64_t a3 = (uint64_t)((tj4 + tj2 + tj6) / 2);
    uint64_t a4 = (uint64_t)((tj4 + tj5 + tj3) / 2);
    uint64_t b1 = (uint64_t)((tj1 + tj2 + tj4 + tj5) / 2);
    uint64_t b2 = (uint64_t)((tj2 + tj3 + tj5 + tj6) / 2);
    uint64_t b3 = (uint64_t)((tj1 + tj3 + tj4 + tj6) / 2);

    pfrac_lcm_scaled_product(scaled, lcm_exp, scratch->terms[0].exp, +1,
                             lcm_max_idx, ws);
    if ((s_min & 1) == 0)
      bigint_add(sum_pos, sum_pos, scaled);
    else
      bigint_add(sum_neg, sum_neg, scaled);

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

  *sum_sign_out = bigint_sub_signed(sum_out, sum_pos, sum_neg);
  *lcm_max_idx_io = lcm_max_idx;
}

/* ── exact 9j ────────────────────────────────────────────────────────────── */

void wigner9j_exact(int tj11, int tj12, int tj13, int tj21, int tj22, int tj23,
                    int tj31, int tj32, int tj33, wignernj_exact_t *out) {
  int tk, tk_min, tk_max, pi;
  wignernj_scratch_t *scratch;
  bigint_ws_t *ws;
  pfrac_t *outer, *kdep;
  int *global_lcm, *lcm1, *lcm2, *lcm3;
  bigint_t *sum_pos, *sum_neg, *scaled, *s1, *s2, *s3, *tmp;
  bigint_t *inner_pos, *inner_neg, *inner_scaled;
  int ss1, ss2, ss3;
  int lcm1_max, lcm2_max, lcm3_max, global_lcm_max;
  size_t mw, mw_prod;

  wignernj_exact_reset(out);

  if (!selection_rules_9j(tj11, tj12, tj13, tj21, tj22, tj23, tj31, tj32,
                          tj33)) {
    out->is_zero = 1;
    return;
  }

  tk_min = max3(abs(tj11 - tj33), abs(tj12 - tj23), abs(tj21 - tj32));
  tk_max = min3(tj11 + tj33, tj12 + tj23, tj21 + tj32);
  if (((tk_min + tj11 + tj33) & 1) != 0)
    tk_min++;
  if (((tk_max + tj11 + tj33) & 1) != 0)
    tk_max--;
  if (tk_min > tk_max) {
    out->is_zero = 1;
    return;
  }

  out->sign = (((tj13 + tj22 + tj31) & 1) == 0) ? 1 : -1;

  /* Largest factorial argument is bounded by max(β_u + 1) over all three
   * inner 6j calls, which in turn is at most (sum of all 9 tj plus 2*tk_max)/2
   * + 1.  We use a uniform safe upper bound for sizing. */
  mw = bigint_words_for_factorial((tj11 + tj12 + tj13 + tj21 + tj22 + tj23 +
                                   tj31 + tj32 + tj33 + 2 * tk_max) /
                                      2 +
                                  5);
  /* The 9j accumulator also holds products of three 6j-sized bigints,
   * so allocate up to ~5x that for the long-lived bigints. */
  mw_prod = 5 * mw;

  /* Acquire the thread's cached scratch and bring every buffer up to
   * the size required by the current call. */
  scratch = wignernj_scratch_acquire();
  ws = &scratch->ws;
  outer = &scratch->pfracs[0];
  kdep = &scratch->pfracs[1];
  global_lcm = scratch->lcm_exp[0];
  lcm1 = scratch->lcm_exp[1];
  lcm2 = scratch->lcm_exp[2];
  lcm3 = scratch->lcm_exp[3];
  s1 = &scratch->bigints[0];
  s2 = &scratch->bigints[1];
  s3 = &scratch->bigints[2];
  tmp = &scratch->bigints[3];
  scaled = &scratch->bigints[4];
  sum_pos = &scratch->bigints[5];
  sum_neg = &scratch->bigints[6];
  inner_pos = &scratch->bigints[7];
  inner_neg = &scratch->bigints[8];
  inner_scaled = &scratch->bigints[9];

  bigint_ws_reserve(ws, mw_prod);
  pfrac_zero(outer);
  pfrac_zero(kdep);
  /* Clear the dirty prefix of every lcm array.  Note: racah_6j_sum
   * tracks its own per-call dirty bound through lcm{1,2,3}_max, so we
   * only need to handle global_lcm here. */
  wignernj_scratch_lcm_clear(scratch, 0);
  /* lcm{1,2,3} are zeroed inside racah_6j_sum via *lcm_max_idx_io;
   * pass in the previous call's value (if any) so it knows how much
   * of the suffix is still clean. */
  bigint_set_zero(s1);
  bigint_reserve(s1, mw);
  bigint_set_zero(s2);
  bigint_reserve(s2, mw);
  bigint_set_zero(s3);
  bigint_reserve(s3, mw);
  bigint_set_zero(tmp);
  bigint_reserve(tmp, mw_prod);
  bigint_set_zero(scaled);
  bigint_reserve(scaled, mw_prod);
  bigint_set_zero(sum_pos);
  bigint_reserve(sum_pos, mw_prod);
  bigint_set_zero(sum_neg);
  bigint_reserve(sum_neg, mw_prod);
  bigint_set_zero(inner_pos);
  bigint_reserve(inner_pos, mw);
  bigint_set_zero(inner_neg);
  bigint_reserve(inner_neg, mw);
  bigint_set_zero(inner_scaled);
  bigint_reserve(inner_scaled, mw);

  /*
   * Outer sqrt: product of the six k-independent Δ factors, two per 6j:
   *   6j_1 = {j11 j21 j31; j32 j33 k}: Δ(j11,j21,j31)  Δ(j32,j33,j31)
   *   6j_2 = {j11 j12 j13; j23 j33 k}: Δ(j11,j12,j13)  Δ(j23,j33,j13)
   *   6j_3 = {j22 j21 j23; k  j12 j32}: Δ(j22,j21,j23)  Δ(j22,j12,j32)
   */
  add_delta_sqrt(outer, tj11, tj21, tj31);
  add_delta_sqrt(outer, tj32, tj33, tj31);
  add_delta_sqrt(outer, tj11, tj12, tj13);
  add_delta_sqrt(outer, tj23, tj33, tj13);
  add_delta_sqrt(outer, tj22, tj21, tj23);
  add_delta_sqrt(outer, tj22, tj12, tj32);

  /* lcm{1,2,3}_max start at the previous call's dirty bound so
   * racah_6j_sum knows how much of those arrays it can leave alone. */
  lcm1_max = scratch->lcm_max_dirty[1];
  lcm2_max = scratch->lcm_max_dirty[2];
  lcm3_max = scratch->lcm_max_dirty[3];
  global_lcm_max = 0;

  /* ── Pass 1: global LCM ───────────────────────────────────────────────── */
  for (tk = tk_min; tk <= tk_max; tk += 2) {
    /*
     * k-dependent rational: Δ²(j11,j33,k)·Δ²(j21,j32,k)·Δ²(j12,j23,k).
     * Two add_delta_sqrt calls per pair → all pfrac exponents are even →
     * purely rational, no sqrt residual.
     */
    pfrac_zero(kdep);
    add_delta_sqrt(kdep, tj11, tj33, tk);
    add_delta_sqrt(kdep, tj11, tj33, tk);
    add_delta_sqrt(kdep, tj21, tj32, tk);
    add_delta_sqrt(kdep, tj21, tj32, tk);
    add_delta_sqrt(kdep, tj12, tj23, tk);
    add_delta_sqrt(kdep, tj12, tj23, tk);

    racah_6j_sum(tj11, tj21, tj31, tj32, tj33, tk, tmp, &ss1, lcm1, &lcm1_max,
                 inner_pos, inner_neg, inner_scaled, scratch, ws);
    racah_6j_sum(tj11, tj12, tj13, tj23, tj33, tk, tmp, &ss2, lcm2, &lcm2_max,
                 inner_pos, inner_neg, inner_scaled, scratch, ws);
    racah_6j_sum(tj22, tj21, tj23, tk, tj12, tj32, tmp, &ss3, lcm3, &lcm3_max,
                 inner_pos, inner_neg, inner_scaled, scratch, ws);

    /* Iterate over the union of contributing prime indices: kdep,
     * lcm1, lcm2, lcm3.  Anything beyond this bound is zero. */
    int union_max = kdep->max_idx;
    if (lcm1_max > union_max)
      union_max = lcm1_max;
    if (lcm2_max > union_max)
      union_max = lcm2_max;
    if (lcm3_max > union_max)
      union_max = lcm3_max;
    for (pi = 0; pi < union_max; pi++) {
      /* kdep.exp[pi] is the sqrt-argument exponent; divide by 2 for value */
      int net = kdep->exp[pi] / 2 - lcm1[pi] - lcm2[pi] - lcm3[pi];
      if (-net > global_lcm[pi])
        global_lcm[pi] = -net;
    }
    if (union_max > global_lcm_max)
      global_lcm_max = union_max;
  }

  /* ── Pass 2: accumulate ───────────────────────────────────────────────── */
  for (tk = tk_min; tk <= tk_max; tk += 2) {
    pfrac_zero(kdep);
    add_delta_sqrt(kdep, tj11, tj33, tk);
    add_delta_sqrt(kdep, tj11, tj33, tk);
    add_delta_sqrt(kdep, tj21, tj32, tk);
    add_delta_sqrt(kdep, tj21, tj32, tk);
    add_delta_sqrt(kdep, tj12, tj23, tk);
    add_delta_sqrt(kdep, tj12, tj23, tk);

    racah_6j_sum(tj11, tj21, tj31, tj32, tj33, tk, s1, &ss1, lcm1, &lcm1_max,
                 inner_pos, inner_neg, inner_scaled, scratch, ws);
    racah_6j_sum(tj11, tj12, tj13, tj23, tj33, tk, s2, &ss2, lcm2, &lcm2_max,
                 inner_pos, inner_neg, inner_scaled, scratch, ws);
    racah_6j_sum(tj22, tj21, tj23, tk, tj12, tj32, s3, &ss3, lcm3, &lcm3_max,
                 inner_pos, inner_neg, inner_scaled, scratch, ws);

    if (bigint_is_zero(s1) || bigint_is_zero(s2) || bigint_is_zero(s3))
      continue;

    /* scale[pi] = global_lcm[pi] + net[pi] >= 0 by construction.
     * Bound iteration by the union max of all contributing arrays. */
    bigint_set_u64(scaled, (uint64_t)(tk + 1)); /* factor (2k+1) */
    int union_max = global_lcm_max;
    if (kdep->max_idx > union_max)
      union_max = kdep->max_idx;
    if (lcm1_max > union_max)
      union_max = lcm1_max;
    if (lcm2_max > union_max)
      union_max = lcm2_max;
    if (lcm3_max > union_max)
      union_max = lcm3_max;
    for (pi = 0; pi < union_max; pi++) {
      int e =
          global_lcm[pi] + kdep->exp[pi] / 2 - lcm1[pi] - lcm2[pi] - lcm3[pi];
      if (e > 0)
        bigint_mul_prime_pow_ws(scaled, (uint64_t)g_primes[pi], e, ws);
    }

    bigint_mul_ws(tmp, scaled, s1, ws);
    bigint_mul_ws(scaled, tmp, s2, ws);
    bigint_mul_ws(tmp, scaled, s3, ws);

    /* outer_phase absorbed into out->sign; accumulate by ss1*ss2*ss3 */
    if (ss1 * ss2 * ss3 > 0)
      bigint_add(sum_pos, sum_pos, tmp);
    else
      bigint_add(sum_neg, sum_neg, tmp);
  }

  out->sum_sign = bigint_sub_signed(&out->sum, sum_pos, sum_neg);

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
  pfrac_bigint_mul_prime_pow_array(&out->int_den, global_lcm, global_lcm_max,
                                   ws);

  /* Update dirty bounds for next call's tail-zeroing. */
  wignernj_scratch_lcm_dirty(scratch, 0, global_lcm_max);
  wignernj_scratch_lcm_dirty(scratch, 1, lcm1_max);
  wignernj_scratch_lcm_dirty(scratch, 2, lcm2_max);
  wignernj_scratch_lcm_dirty(scratch, 3, lcm3_max);

  wignernj_scratch_relinquish(scratch);
}

/* ── public API ──────────────────────────────────────────────────────────── */

int wigner9j_max_factorial(int tj11, int tj12, int tj13, int tj21, int tj22,
                           int tj23, int tj31, int tj32, int tj33) {
  /* The k-loop reaches tk = tk_max = min(tj11+tj33, tj12+tj23,
   * tj21+tj32).  At that tk the inner racah_6j_sum factorials touch
   * (s+1)! with s_max bounded by the sum of all 9 outer quantum
   * numbers plus 2*tk_max divided by 2 (a generous upper bound that
   * covers every triangle-denominator and every sum-of-4 inside any
   * of the three inner 6j's). */
  int tkmax_a = tj11 + tj33;
  int tkmax_b = tj12 + tj23;
  int tkmax_c = tj21 + tj32;
  int tk_max = (tkmax_a < tkmax_b) ? tkmax_a : tkmax_b;
  if (tkmax_c < tk_max)
    tk_max = tkmax_c;
  return (tj11 + tj12 + tj13 + tj21 + tj22 + tj23 + tj31 + tj32 + tj33 +
          2 * tk_max) /
             2 +
         1;
}

float wigner9j_f(int tj11, int tj12, int tj13, int tj21, int tj22, int tj23,
                 int tj31, int tj32, int tj33) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  float result;
  wigner9j_exact(tj11, tj12, tj13, tj21, tj22, tj23, tj31, tj32, tj33,
                 &s->exact);
  result = wignernj_exact_to_float(&s->exact);
  wignernj_scratch_relinquish(s);
  return result;
}

double wigner9j(int tj11, int tj12, int tj13, int tj21, int tj22, int tj23,
                int tj31, int tj32, int tj33) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  double result;
  wigner9j_exact(tj11, tj12, tj13, tj21, tj22, tj23, tj31, tj32, tj33,
                 &s->exact);
  result = wignernj_exact_to_double(&s->exact);
  wignernj_scratch_relinquish(s);
  return result;
}

long double wigner9j_l(int tj11, int tj12, int tj13, int tj21, int tj22,
                       int tj23, int tj31, int tj32, int tj33) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  long double result;
  wigner9j_exact(tj11, tj12, tj13, tj21, tj22, tj23, tj31, tj32, tj33,
                 &s->exact);
  result = wignernj_exact_to_long_double(&s->exact);
  wignernj_scratch_relinquish(s);
  return result;
}

#ifdef WIGNERNJ_HAVE_QUADMATH
__float128 wigner9j_q(int tj11, int tj12, int tj13, int tj21, int tj22,
                      int tj23, int tj31, int tj32, int tj33) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  __float128 result;
  wigner9j_exact(tj11, tj12, tj13, tj21, tj22, tj23, tj31, tj32, tj33,
                 &s->exact);
  result = wignernj_exact_to_float128(&s->exact);
  wignernj_scratch_relinquish(s);
  return result;
}
#endif

#ifdef WIGNERNJ_HAVE_MPFR
void wigner9j_mpfr(mpfr_t rop, int tj11, int tj12, int tj13, int tj21, int tj22,
                   int tj23, int tj31, int tj32, int tj33, mpfr_rnd_t rnd) {
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  wigner9j_exact(tj11, tj12, tj13, tj21, tj22, tj23, tj31, tj32, tj33,
                 &s->exact);
  wignernj_exact_to_mpfr(rop, &s->exact, rnd);
  wignernj_scratch_relinquish(s);
}
#endif
