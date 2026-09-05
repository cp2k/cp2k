/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * 64-bit arithmetic primitives used by bigint.c.
 *
 * Two implementations are selected at compile time:
 *
 *   - When __SIZEOF_INT128__ is defined (GCC, Clang, ICC, and most other
 *     Unix compilers), the primitives use the __uint128_t extension for
 *     64x64 -> 128 products and 128/64 divisions.
 *
 *   - Otherwise, a pure-C99 fallback is used: 64x64 products are formed
 *     from four 32x32 -> 64 partial products combined with explicit carry
 *     tracking, and 128/64 divisions (whose divisor in this code is always
 *     a prime <= PRIME_SIEVE_LIMIT and therefore fits in 32 bits) are
 *     implemented as two 64/32 long-division steps.
 *
 * Both code paths are bit-exact and produce identical results.
 */
#ifndef BIGINT_ARITH_H
#define BIGINT_ARITH_H

#include <stdint.h>

/* Defining BIGINT_FORCE_PORTABLE at compile time forces use of the pure-C99
 * fallback path even on compilers that support __uint128_t.  Useful for
 * exercising the fallback in CI. */
#if defined(__SIZEOF_INT128__) && !defined(BIGINT_FORCE_PORTABLE)
#define BIGINT_HAVE_UINT128 1
typedef unsigned __int128 bigint_u128_t;
#else
#define BIGINT_HAVE_UINT128 0
#endif

/* True when bigint_div128 dispatches to inline-asm hardware divq.
 * Mirrors the conditional inside bigint_div128 below.  Used by
 * bigint_div_u64 to decide whether the Möller-Granlund reciprocal-based
 * limb sweep is worth the per-call setup cost: it isn't on this path,
 * because divq's variable-latency fast path for sub-64-bit divisors
 * already runs at ~2.3 ns/limb and beats GM's mul-add-mul body once
 * the per-call recip_init divq is folded in. */
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__)) &&        \
    !defined(BIGINT_FORCE_PORTABLE) && !defined(BIGINT_NO_DIVQ)
#define BIGINT_HAVE_DIVQ_INLINE 1
#else
#define BIGINT_HAVE_DIVQ_INLINE 0
#endif

/* {*plo, *phi} = a * b */
static inline void bigint_mul64x64(uint64_t a, uint64_t b, uint64_t *plo,
                                   uint64_t *phi) {
#if BIGINT_HAVE_UINT128
  bigint_u128_t p = (bigint_u128_t)a * b;
  *plo = (uint64_t)p;
  *phi = (uint64_t)(p >> 64);
#else
  uint32_t a0 = (uint32_t)a, a1 = (uint32_t)(a >> 32);
  uint32_t b0 = (uint32_t)b, b1 = (uint32_t)(b >> 32);
  uint64_t p00 = (uint64_t)a0 * b0;
  uint64_t p01 = (uint64_t)a0 * b1;
  uint64_t p10 = (uint64_t)a1 * b0;
  uint64_t p11 = (uint64_t)a1 * b1;
  /* Combine partial products: a*b = p11 * 2^64 + (p10+p01) * 2^32 + p00 */
  uint64_t mid = (p00 >> 32) + (uint32_t)p01 + (uint32_t)p10;
  *plo = ((uint32_t)p00) | (mid << 32);
  *phi = p11 + (p01 >> 32) + (p10 >> 32) + (mid >> 32);
#endif
}

/* {*plo, *phi} = a * b + c */
static inline void bigint_mul_add(uint64_t a, uint64_t b, uint64_t c,
                                  uint64_t *plo, uint64_t *phi) {
#if BIGINT_HAVE_UINT128
  bigint_u128_t p = (bigint_u128_t)a * b + c;
  *plo = (uint64_t)p;
  *phi = (uint64_t)(p >> 64);
#else
  uint64_t lo, hi;
  bigint_mul64x64(a, b, &lo, &hi);
  *plo = lo + c;
  *phi = hi + (*plo < lo);
#endif
}

/* {*plo, *phi} = a * b + c + d */
static inline void bigint_mul_add2(uint64_t a, uint64_t b, uint64_t c,
                                   uint64_t d, uint64_t *plo, uint64_t *phi) {
#if BIGINT_HAVE_UINT128
  bigint_u128_t p = (bigint_u128_t)a * b + c + d;
  *plo = (uint64_t)p;
  *phi = (uint64_t)(p >> 64);
#else
  uint64_t lo, hi, t;
  bigint_mul64x64(a, b, &lo, &hi);
  t = lo + c;
  hi += (t < lo);
  *plo = t + d;
  *phi = hi + (*plo < t);
#endif
}

/* a + b + cin.  Returns the sum and writes the carry-out (0 or 1) to *cout. */
static inline uint64_t bigint_addc(uint64_t a, uint64_t b, uint64_t cin,
                                   uint64_t *cout) {
  uint64_t s = a + b;
  uint64_t c1 = (s < a);
  uint64_t r = s + cin;
  uint64_t c2 = (r < s);
  *cout = c1 | c2;
  return r;
}

/* a - b - bin.  Returns the difference and writes the borrow-out (0 or 1)
 * to *bout. */
static inline uint64_t bigint_subb(uint64_t a, uint64_t b, uint64_t bin,
                                   uint64_t *bout) {
  uint64_t d = a - b;
  uint64_t b1 = (a < b);
  uint64_t r = d - bin;
  uint64_t b2 = (d < bin);
  *bout = b1 | b2;
  return r;
}

/*
 * Knuth's Algorithm D for 128/64 unsigned division at base B = 2^32.
 * Caller guarantees hi < d (so the quotient fits in 64 bits) and
 * d >= 2^32 (the simpler 32-bit two-step long division below handles
 * the d < 2^32 case at lower constant cost).
 *
 * Two long-division steps each produce a 32-bit quotient digit,
 * recombined into the 64-bit result.  Uses bigint_mul_add /
 * bigint_subb for the inner 64x64->128 product and 64-bit subtraction-
 * with-borrow, so the same routine compiles correctly on the
 * BIGINT_HAVE_UINT128 path (where the helpers fold to native
 * __uint128_t arithmetic) and on the pure-C99 fallback (where they
 * fall through to four 32x32->64 partial products with explicit carry
 * tracking).
 *
 * Reference: Knuth, TAOCP Vol. 2, Algorithm D (Section 4.3.1).  We
 * specialise to a fixed n=2 (divisor digits) and m=2 (extra
 * dividend digits over the divisor) and inline the trial-quotient
 * refinement (at most two iterations after divisor normalisation)
 * and the post-subtraction add-back correction (which fires with
 * probability ~2^-32 per step).
 */
static inline uint64_t bigint_div128_alg_d(uint64_t hi, uint64_t lo, uint64_t d,
                                           uint64_t *rem) {
  int s, step;
  uint64_t dn, uhi, ulo, d1, d0, mid_digits[2];
  uint64_t cur, qhat, rhat, prod_lo, prod_hi, num_lo, num_hi;
  uint64_t bsub, badd;
  uint32_t qdigit[2];

  /* Step 1: normalise divisor so its top bit is set. */
  s = 0;
  {
    uint64_t t = d;
    while ((t & (1ULL << 63)) == 0) {
      t <<= 1;
      s++;
    }
  }
  dn = d << s;
  /* Apply the same shift to the dividend.  hi < d implies the shifted
   * top half stays within 64 bits. */
  if (s == 0) {
    uhi = hi;
    ulo = lo;
  } else {
    uhi = (hi << s) | (lo >> (64 - s));
    ulo = lo << s;
  }

  /* Decompose into 32-bit digits. */
  d1 = dn >> 32;
  d0 = dn & 0xFFFFFFFFULL;
  mid_digits[0] = ulo >> 32;           /* fed in for step 0 */
  mid_digits[1] = ulo & 0xFFFFFFFFULL; /* fed in for step 1 */

  cur = uhi;
  for (step = 0; step < 2; step++) {
    uint64_t mid = mid_digits[step];

    /* Trial quotient digit q_hat = floor(cur / d1), capped at B-1. */
    if ((cur >> 32) >= d1) {
      qhat = 0xFFFFFFFFULL;
      /* rhat = cur - qhat * d1 = (cur - d1*B) + d1, where
       * cur - d1*B >= 0 since (cur >> 32) >= d1. */
      rhat = (cur - (d1 << 32)) + d1;
    } else {
      qhat = cur / d1;
      rhat = cur - qhat * d1;
    }

    /* Refine: while q_hat * d0 > B*r_hat + mid, decrement q_hat.
     * Skip once r_hat overflows 32 bits (the comparison is then
     * trivially false; computing the RHS would also overflow). */
    while (rhat < (1ULL << 32) && qhat * d0 > (rhat << 32) + mid) {
      qhat--;
      rhat += d1;
    }

    /* Subtract qhat * dn from the 96-bit numerator (cur << 32 | mid).
     * prod = qhat * dn fits in 96 bits; numer96 fits in 96 bits;
     * the difference fits in 64 bits when qhat is correct, with a
     * possible 1-step add-back if qhat overshot by 1. */
    bigint_mul_add(qhat, dn, 0, &prod_lo, &prod_hi);
    num_lo = (cur << 32) | mid;
    num_hi = cur >> 32;
    cur = bigint_subb(num_lo, prod_lo, 0, &bsub);
    num_hi = bigint_subb(num_hi, prod_hi, bsub, &bsub);
    if (bsub) {
      /* qhat was 1 too large (rare: probability ~2^-32 per step
       * after refinement).  Add back one copy of dn. */
      qhat--;
      cur = bigint_addc(cur, dn, 0, &badd);
      num_hi = num_hi + badd; /* clears num_hi back to 0 */
    }

    qdigit[step] = (uint32_t)qhat;
  }

  /* Recombine quotient and denormalise remainder. */
  *rem = cur >> s;
  return ((uint64_t)qdigit[0] << 32) | qdigit[1];
}

/*
 * Divide [hi:lo] (128 bits, with hi < d) by d.  Returns the 64-bit
 * quotient and writes the remainder (< d) to *rem.  Handles arbitrary
 * 64-bit divisors so the library scales correctly when the prime
 * table is regenerated for an enlarged sieve limit (where individual
 * primes themselves can exceed 2^32).
 */
static inline uint64_t bigint_div128(uint64_t hi, uint64_t lo, uint64_t d,
                                     uint64_t *rem) {
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__)) &&        \
    !defined(BIGINT_FORCE_PORTABLE) && !defined(BIGINT_NO_DIVQ)
  /* x86-64 has a hardware divq that does 128/64 -> 64+remainder in
   * ~20 cycles when hi < d, which is always the case in bigint code
   * (bigint_div_u64 carries the remainder of the previous limb as
   * `hi`, and that remainder is < d by definition).  Letting the
   * compiler emit __uint128_t division calls libgcc's generic
   * __udivmodti4 (~50-100 cycles software emulation) instead --
   * profile-confirmed as a top-3 hot spot at j=4000.  Inline asm
   * forces the hardware instruction.  The -DBIGINT_NO_DIVQ build
   * override disables this asm path while keeping __uint128_t
   * available, so a dedicated CI cell exercises the Algorithm-D
   * code below on x86-64 hardware -- that path is the one selected
   * on aarch64, ppc64le, and any other target with __uint128_t but
   * no hardware 128/64 instruction. */
  uint64_t q, r;
  __asm__("divq %4" : "=a"(q), "=d"(r) : "a"(lo), "d"(hi), "rm"(d));
  *rem = r;
  return q;
#else
  /* Software path.  Two-step 64/32 long division for d < 2^32 (the
   * common case in libwignernj: every prime <= PRIME_SIEVE_LIMIT
   * with the default sieve fits comfortably in 32 bits), Knuth's
   * Algorithm D for d >= 2^32 (the case that arises with batched
   * small-integer denominators in the Racah-sum recurrence and with
   * regenerated prime tables whose sieve limit exceeds 2^32).  Both
   * branches use only 64-bit arithmetic primitives, so the routine
   * compiles correctly on toolchains without __uint128_t. */
  if (d < (1ULL << 32)) {
    uint64_t numer, q1, r1, q2;
    numer = (hi << 32) | (lo >> 32);
    q1 = numer / d;
    r1 = numer % d;
    numer = (r1 << 32) | (lo & 0xFFFFFFFFu);
    q2 = numer / d;
    *rem = numer % d;
    return (q1 << 32) | (q2 & 0xFFFFFFFFu);
  } else {
    return bigint_div128_alg_d(hi, lo, d, rem);
  }
#endif
}

/* ── Möller-Granlund "improved division by invariant integers" ──────────────
 *
 * Reference: Niels Möller and Torbjörn Granlund, "Improved division by
 * invariant integers", IEEE Trans. Comput. 60(2), pp. 165-175 (2011),
 * doi:10.1109/TC.2010.143.
 *
 * Used by bigint_div_u64 when the dividend has more than one limb, so that
 * the per-call reciprocal precompute (one bigint_div128 call plus a
 * leading-zero count) amortises across the limb sweep.  Per-limb work is
 * one 64x64->128 mul-high, one 128-bit add, one 64x64->64 mul, and a pair
 * of conditional fixups -- about 6-8 cycles on x86-64 and aarch64, vs
 * ~22 cycles for hardware divq on x86-64 or ~30-50 cycles for the
 * Algorithm-D fallback above.
 *
 * On the pure-C99 path (no __uint128_t), the existing two-step 64/32 long
 * division in bigint_div128 is faster than Möller-Granlund whenever
 * d < 2^32 (the common case for prime-power scalars), so the reciprocal
 * carries a `use_gm` flag and bigint_div128_recip dispatches accordingly.
 */

typedef struct {
  uint64_t v;      /* floor((B^2 - 1) / d_norm) - B,  B = 2^64       */
  uint64_t d_norm; /* d << s  (top bit of d_norm is set)             */
  int s;           /* leading-zero count of d, in [0, 63]            */
  int use_gm;      /* 1 = Möller-Granlund, 0 = portable d<2^32 path  */
} bigint_recip64_t;

/* Count leading zeros of a non-zero 64-bit value. */
static inline int bigint_clz64(uint64_t x) {
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_clzll(x);
#else
  int n = 0;
  if ((x >> 32) == 0) {
    n += 32;
    x <<= 32;
  }
  if ((x >> 48) == 0) {
    n += 16;
    x <<= 16;
  }
  if ((x >> 56) == 0) {
    n += 8;
    x <<= 8;
  }
  if ((x >> 60) == 0) {
    n += 4;
    x <<= 4;
  }
  if ((x >> 62) == 0) {
    n += 2;
    x <<= 2;
  }
  if ((x >> 63) == 0) {
    n += 1;
  }
  return n;
#endif
}

/* Build a reciprocal of d for use by bigint_div128_recip.  d must be > 0.
 * The returned struct is independent of dividend, so one call before the
 * limb sweep is sufficient. */
static inline bigint_recip64_t bigint_recip64_init(uint64_t d) {
  bigint_recip64_t r;
  uint64_t rem;
  r.s = bigint_clz64(d);
  r.d_norm = d << r.s;

  /* Pure-C99 portable path keeps its tight 64/32 fast-path for small
   * divisors -- the Möller-Granlund mul-add-mul body is slightly more
   * arithmetic on that backend than the two-step long division, so
   * for d < 2^32 we tell bigint_div128_recip to fall through to
   * bigint_div128 on every limb.  On the hardware-divq and
   * __uint128_t backends, GM is at-least-as-fast as divq even for
   * small d, so the flag is set unconditionally. */
#if !BIGINT_HAVE_UINT128
  if (d < (1ULL << 32)) {
    r.v = 0;
    r.use_gm = 0;
    return r;
  }
#endif
  /* v = ((B - d_norm) * B + (B - 1)) / d_norm
   *   = bigint_div128(~0 - d_norm, ~0, d_norm)
   * (hi = ~0 - d_norm < d_norm because d_norm has its MSB set, so the
   * bigint_div128 hi<d precondition holds). */
  r.v = bigint_div128(~(uint64_t)0 - r.d_norm, ~(uint64_t)0, r.d_norm, &rem);
  r.use_gm = 1;
  (void)rem;
  return r;
}

/* Divide [hi:lo] by d using a precomputed reciprocal.  hi < d required.
 * Returns the 64-bit quotient and writes the remainder (< d) to *rem. */
static inline uint64_t bigint_div128_recip(uint64_t hi, uint64_t lo, uint64_t d,
                                           bigint_recip64_t r, uint64_t *rem) {
  uint64_t u1, u0, q1, q0, r_n, carry;

  if (!r.use_gm) {
    /* Pure-C99 fallback: defer to the two-step 64/32 long division
     * inside bigint_div128, which beats Möller-Granlund for the
     * d < 2^32 regime on this backend. */
    return bigint_div128(hi, lo, d, rem);
  }

  /* Renormalise dividend by the same shift applied to the divisor.
   * After shifting, the (logical) 128-bit dividend has u1 < d_norm. */
  if (r.s == 0) {
    u1 = hi;
    u0 = lo;
  } else {
    u1 = (hi << r.s) | (lo >> (64 - r.s));
    u0 = lo << r.s;
  }

  /* Möller-Granlund Algorithm 4 step on normalised inputs:
   *   {q1, q0} = v * u1
   *   {q1, q0} += {u1, u0}
   *   q1 += 1
   *   r_n = u0 - q1 * d_norm
   *   if r_n > q0: q1 -= 1; r_n += d_norm
   *   if r_n >= d_norm: q1 += 1; r_n -= d_norm  (rare, prob ~2^-64)
   */
  bigint_mul64x64(r.v, u1, &q0, &q1);
  q0 = bigint_addc(q0, u0, 0, &carry);
  q1 = bigint_addc(q1, u1, carry, &carry);
  (void)carry;
  q1 += 1;
  r_n = u0 - q1 * r.d_norm;
  if (r_n > q0) {
    q1 -= 1;
    r_n += r.d_norm;
  }
  if (r_n >= r.d_norm) {
    q1 += 1;
    r_n -= r.d_norm;
  }

  /* Denormalise the remainder. */
  *rem = r_n >> r.s;
  return q1;
}

#endif /* BIGINT_ARITH_H */
