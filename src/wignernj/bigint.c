/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * 64-bit and 128-bit arithmetic primitives are abstracted into
 * bigint_arith.h, which provides a native __uint128_t-based path on
 * GCC/Clang/ICC and a pure-C99 fallback on other compilers (e.g. MSVC).
 */
#include "bigint.h"
#include "bigint_arith.h"
#include "xalloc.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ── internal helpers ────────────────────────────────────────────────────── */

static void bigint_ensure(bigint_t *a, size_t need) {
  if (need > a->cap) {
    size_t nc = a->cap ? a->cap * 2 : 4;
    while (nc < need)
      nc *= 2;
    a->words = (uint64_t *)xrealloc(a->words, nc * sizeof(uint64_t));
    memset(a->words + a->cap, 0, (nc - a->cap) * sizeof(uint64_t));
    a->cap = nc;
  }
}

void bigint_reserve(bigint_t *a, size_t cap) { bigint_ensure(a, cap); }

static void bigint_trim(bigint_t *a) {
  while (a->size > 0 && a->words[a->size - 1] == 0)
    a->size--;
}

/* ── lifecycle ───────────────────────────────────────────────────────────── */

void bigint_init(bigint_t *a) {
  a->words = NULL;
  a->size = a->cap = 0;
}
void bigint_free(bigint_t *a) {
  free(a->words);
  a->words = NULL;
  a->size = a->cap = 0;
}
void bigint_set_zero(bigint_t *a) { a->size = 0; }
int bigint_is_zero(const bigint_t *a) { return a->size == 0; }

void bigint_copy(bigint_t *dst, const bigint_t *src) {
  if (src->size == 0) {
    dst->size = 0;
    return;
  }
  bigint_ensure(dst, src->size);
  memcpy(dst->words, src->words, src->size * sizeof(uint64_t));
  dst->size = src->size;
}

void bigint_set_u64(bigint_t *a, uint64_t v) {
  bigint_ensure(a, 1);
  a->words[0] = v;
  a->size = (v != 0) ? 1 : 0;
}

/* ── comparison ──────────────────────────────────────────────────────────── */

int bigint_cmp(const bigint_t *a, const bigint_t *b) {
  size_t i;
  if (a->size != b->size)
    return (a->size > b->size) ? 1 : -1;
  for (i = a->size; i-- > 0;)
    if (a->words[i] != b->words[i])
      return (a->words[i] > b->words[i]) ? 1 : -1;
  return 0;
}

/*
 * The four "linear" routines below (add, sub_uu, mul_u64, div_u64) all use
 * a read-then-write pattern at index i: each iteration first reads the
 * operand word(s) at i and only then writes r->words[i].  This makes them
 * correct even when r aliases an operand, so we can write the result
 * directly into r and avoid the per-call temp + copy that the previous
 * implementation always performed.
 *
 * bigint_mul cannot be aliased away as easily: its inner loop reads
 * operand words ahead of where it writes, so when r aliases a or b a
 * separate destination is still required.  We use one in that case only.
 */

/* ── addition ────────────────────────────────────────────────────────────── */

void bigint_add(bigint_t *r, const bigint_t *a, const bigint_t *b) {
  const bigint_t *big = (a->size >= b->size) ? a : b;
  const bigint_t *sml = (a->size >= b->size) ? b : a;
  size_t n = big->size, i;
  uint64_t carry = 0;
  if (n == 0) {
    bigint_set_zero(r);
    return;
  }
  bigint_ensure(r, n + 1);
  for (i = 0; i < sml->size; i++)
    r->words[i] = bigint_addc(big->words[i], sml->words[i], carry, &carry);
  for (; i < n; i++)
    r->words[i] = bigint_addc(big->words[i], 0, carry, &carry);
  r->words[n] = carry;
  r->size = n + 1;
  bigint_trim(r);
}

/* ── subtraction (unsigned, a >= b assumed) ─────────────────────────────── */

static void bigint_sub_uu(bigint_t *r, const bigint_t *a, const bigint_t *b) {
  uint64_t borrow = 0;
  size_t i;
  if (a->size == 0) {
    bigint_set_zero(r);
    return;
  }
  bigint_ensure(r, a->size);
  for (i = 0; i < b->size; i++)
    r->words[i] = bigint_subb(a->words[i], b->words[i], borrow, &borrow);
  for (; i < a->size; i++)
    r->words[i] = bigint_subb(a->words[i], 0, borrow, &borrow);
  r->size = a->size;
  bigint_trim(r);
}

int bigint_sub_signed(bigint_t *r, const bigint_t *a, const bigint_t *b) {
  int c = bigint_cmp(a, b);
  if (c == 0) {
    bigint_set_zero(r);
    return 1;
  }
  if (c > 0) {
    bigint_sub_uu(r, a, b);
    return 1;
  }
  bigint_sub_uu(r, b, a);
  return -1;
}

/* ── multiplication ──────────────────────────────────────────────────────── */

void bigint_mul_u64(bigint_t *r, const bigint_t *a, uint64_t b) {
  size_t i;
  uint64_t carry = 0;
  if (b == 0 || a->size == 0) {
    bigint_set_zero(r);
    return;
  }
  bigint_ensure(r, a->size + 1);
  for (i = 0; i < a->size; i++) {
    uint64_t lo, hi;
    bigint_mul_add(a->words[i], b, carry, &lo, &hi);
    r->words[i] = lo;
    carry = hi;
  }
  r->words[a->size] = carry;
  r->size = a->size + 1;
  bigint_trim(r);
}

/*
 * Crossover below which schoolbook multiplication is faster than
 * Karatsuba.  Tuned via tools/bench_karatsuba.c on x86-64 GCC; the
 * exact value is not critical -- anywhere from ~16 to ~64 limbs gives
 * almost the same performance because the Karatsuba constant factor
 * is dominated by the schoolbook leaf at this size.  Overridable via
 * -DKARATSUBA_THRESHOLD=N at compile time so the benchmark can sweep
 * both code paths.
 */
#ifndef KARATSUBA_THRESHOLD
#define KARATSUBA_THRESHOLD 32
#endif

/* r[0..n-1] += b[0..nb-1]; returns final carry-out (0 or 1). */
static uint64_t add_inplace_words(uint64_t *r, size_t n, const uint64_t *b,
                                  size_t nb) {
  uint64_t carry = 0;
  size_t i;
  for (i = 0; i < nb && i < n; i++)
    r[i] = bigint_addc(r[i], b[i], carry, &carry);
  for (; i < n && carry; i++)
    r[i] = bigint_addc(r[i], 0, carry, &carry);
  return carry;
}

/* r[0..n-1] -= b[0..nb-1]; returns final borrow-out (0 or 1). */
static uint64_t sub_inplace_words(uint64_t *r, size_t n, const uint64_t *b,
                                  size_t nb) {
  uint64_t borrow = 0;
  size_t i;
  for (i = 0; i < nb && i < n; i++)
    r[i] = bigint_subb(r[i], b[i], borrow, &borrow);
  for (; i < n && borrow; i++)
    r[i] = bigint_subb(r[i], 0, borrow, &borrow);
  return borrow;
}

/*
 * Schoolbook leaf:  r[0 .. na+nb-1] += a[0..na-1] * b[0..nb-1].
 * The += is for bit-exact reuse from inside the Karatsuba combine, where
 * the destination region may already hold partial results.  For pure
 * multiplication the caller zeros r first; nothing in this routine reads
 * r before its first write at any given index.  r must be disjoint from
 * a and b.
 */
static void school_mul_words(uint64_t *r, const uint64_t *a, size_t na,
                             const uint64_t *b, size_t nb) {
  size_t i, j;
  for (i = 0; i < na; i++) {
    uint64_t carry = 0;
    for (j = 0; j < nb; j++) {
      uint64_t lo, hi;
      bigint_mul_add2(a[i], b[j], r[i + j], carry, &lo, &hi);
      r[i + j] = lo;
      carry = hi;
    }
    r[i + nb] += carry;
  }
}

/*
 * Karatsuba multiplication on word arrays:  r[0..na+nb-1] = a*b.
 *
 *   - r is zero-initialized and disjoint from a and b.
 *   - scratch points to at least kar_scratch_words(max(na,nb)) words of
 *     working memory, also disjoint from r, a, b.
 *
 * Splits at k = min(na,nb)/2:
 *   a = a1*B^k + a0,  b = b1*B^k + b0,  B = 2^64
 *   z0 = a0*b0,  z2 = a1*b1,  zm = (a0+a1)(b0+b1) - z0 - z2
 *   a*b = z2*B^(2k) + zm*B^k + z0
 *
 * Below KARATSUBA_THRESHOLD the routine falls through to schoolbook.
 */
static void kar_mul_words(uint64_t *r, const uint64_t *a, size_t na,
                          const uint64_t *b, size_t nb, uint64_t *scratch) {
  size_t mn = (na < nb) ? na : nb;
  size_t k, na1, nb1;
  size_t suma_alloc, sumb_alloc, z1_alloc;
  uint64_t *sum_a, *sum_b, *z1, *sub;
  size_t suma_size, sumb_size, z1_size, z0_size, z2_size;

  if (mn < KARATSUBA_THRESHOLD) {
    school_mul_words(r, a, na, b, nb);
    return;
  }

  k = mn / 2;
  na1 = na - k;
  nb1 = nb - k;

  /* Slice the scratch buffer.  na1 >= k and nb1 >= k since na,nb >= 2k. */
  suma_alloc = na1 + 1;
  sumb_alloc = nb1 + 1;
  z1_alloc = na1 + nb1 + 2;
  sum_a = scratch;
  sum_b = sum_a + suma_alloc;
  z1 = sum_b + sumb_alloc;
  sub = z1 + z1_alloc; /* sub-call scratch starts here */

  /* z0 = a0 * b0  ->  r[0 .. 2k-1].  r[2k..] is still zero. */
  kar_mul_words(r, a, k, b, k, sub);
  /* z2 = a1 * b1  ->  r[2k .. 2k + na1+nb1 - 1] */
  kar_mul_words(r + 2 * k, a + k, na1, b + k, nb1, sub);

  /* sum_a = a0 + a1.  Start with the longer half (a1, na1 >= k limbs)
   * and add a0 (k limbs); the result fits in na1+1 limbs. */
  memcpy(sum_a, a + k, na1 * sizeof(uint64_t));
  sum_a[na1] = 0;
  add_inplace_words(sum_a, na1 + 1, a, k);
  suma_size = na1 + 1;
  while (suma_size > 0 && sum_a[suma_size - 1] == 0)
    suma_size--;

  /* sum_b = b0 + b1, same shape. */
  memcpy(sum_b, b + k, nb1 * sizeof(uint64_t));
  sum_b[nb1] = 0;
  add_inplace_words(sum_b, nb1 + 1, b, k);
  sumb_size = nb1 + 1;
  while (sumb_size > 0 && sum_b[sumb_size - 1] == 0)
    sumb_size--;

  /* z1 = sum_a * sum_b (recursive). */
  memset(z1, 0, z1_alloc * sizeof(uint64_t));
  if (suma_size > 0 && sumb_size > 0)
    kar_mul_words(z1, sum_a, suma_size, sum_b, sumb_size, sub);
  z1_size = z1_alloc;
  while (z1_size > 0 && z1[z1_size - 1] == 0)
    z1_size--;

  /* z1 -= z0  (z0 occupies r[0 .. 2k-1]). */
  z0_size = 2 * k;
  while (z0_size > 0 && r[z0_size - 1] == 0)
    z0_size--;
  if (z0_size > 0)
    sub_inplace_words(z1, z1_size, r, z0_size);
  while (z1_size > 0 && z1[z1_size - 1] == 0)
    z1_size--;

  /* z1 -= z2  (z2 occupies r[2k .. 2k + na1+nb1 - 1]). */
  z2_size = na1 + nb1;
  while (z2_size > 0 && r[2 * k + z2_size - 1] == 0)
    z2_size--;
  if (z2_size > 0)
    sub_inplace_words(z1, z1_size, r + 2 * k, z2_size);
  while (z1_size > 0 && z1[z1_size - 1] == 0)
    z1_size--;

  /* r += z1 << (k*64), i.e. r[k .. k+z1_size-1] += z1. */
  if (z1_size > 0)
    add_inplace_words(r + k, na + nb - k, z1, z1_size);
}

/*
 * Upper bound on Karatsuba scratch (in 64-bit words) for a multiplication
 * with operand sizes (na, nb).  Recurrence at level (na, nb):
 *   need(na, nb) = 2*(na1+nb1) + 4  +  max{need over sub-multiplies}
 * The three sub-multiplies are evaluated sequentially and reuse the same
 * sub-scratch slice, so we only need the largest of them; that is the
 * (sum_a, sum_b) cross product on operands of sizes up to (na1+1, nb1+1).
 * Imbalanced operands (na << nb) shrink the smaller dimension by half but
 * leave the larger one nearly intact, so a single-argument estimator is
 * not safe -- both sizes must be tracked.
 */
static size_t kar_scratch_words(size_t na, size_t nb) {
  size_t total = 0;
  while (((na < nb) ? na : nb) >= KARATSUBA_THRESHOLD) {
    size_t mn = (na < nb) ? na : nb;
    size_t k = mn / 2;
    size_t na1 = na - k;
    size_t nb1 = nb - k;
    total += 2 * (na1 + nb1) + 4;
    na = na1 + 1;
    nb = nb1 + 1;
  }
  return total + 16;
}

/*
 * Multiplication with caller-provided scratch:  r[0..na+nb-1] = a*b.
 * r is assumed zero-initialized and disjoint from a, b, ws_scratch.
 * Below the Karatsuba threshold ws_scratch is unused.
 */
static void mul_words_ws(uint64_t *r, const uint64_t *a, size_t na,
                         const uint64_t *b, size_t nb, bigint_t *ws_scratch) {
  size_t mn = (na < nb) ? na : nb;
  if (mn < KARATSUBA_THRESHOLD) {
    school_mul_words(r, a, na, b, nb);
    return;
  }
  bigint_ensure(ws_scratch, kar_scratch_words(na, nb));
  kar_mul_words(r, a, na, b, nb, ws_scratch->words);
}

/*
 * Top-level dispatch without a workspace:  r[0..na+nb-1] = a*b.  Used by
 * the bigint_mul (no-ws) entry point; allocates a one-shot scratch when
 * Karatsuba is selected.
 */
static void mul_words(uint64_t *r, const uint64_t *a, size_t na,
                      const uint64_t *b, size_t nb) {
  size_t mn = (na < nb) ? na : nb;
  uint64_t *scratch;
  if (mn < KARATSUBA_THRESHOLD) {
    school_mul_words(r, a, na, b, nb);
    return;
  }
  scratch = (uint64_t *)xmalloc(kar_scratch_words(na, nb) * sizeof(uint64_t));
  kar_mul_words(r, a, na, b, nb, scratch);
  free(scratch);
}

void bigint_mul(bigint_t *r, const bigint_t *a, const bigint_t *b) {
  bigint_t local;
  bigint_t *dest;
  int aliased;
  if (a->size == 0 || b->size == 0) {
    bigint_set_zero(r);
    return;
  }
  /* Karatsuba reads ahead of where it writes via the recursive
   * combine, so destination aliasing with an operand still requires
   * a separate temp when r == a or r == b. */
  aliased = (r == a || r == b);
  if (aliased) {
    bigint_init(&local);
    bigint_ensure(&local, a->size + b->size);
    memset(local.words, 0, (a->size + b->size) * sizeof(uint64_t));
    local.size = a->size + b->size;
    dest = &local;
  } else {
    bigint_ensure(r, a->size + b->size);
    memset(r->words, 0, (a->size + b->size) * sizeof(uint64_t));
    r->size = a->size + b->size;
    dest = r;
  }
  mul_words(dest->words, a->words, a->size, b->words, b->size);
  bigint_trim(dest);
  if (aliased) {
    bigint_copy(r, &local);
    bigint_free(&local);
  }
}

/* ── division by scalar ──────────────────────────────────────────────────── */

/* Forward declarations of shift helpers used by bigint_div_u64_exact;
 * the definitions live further down in this file. */
static void bigint_shr_inplace(bigint_t *a, int k);

void bigint_div_u64(bigint_t *r, const bigint_t *a, uint64_t b) {
  size_t i;
  uint64_t rem = 0;
  if (a->size == 0) {
    bigint_set_zero(r);
    return;
  }
  bigint_ensure(r, a->size);
#if BIGINT_HAVE_DIVQ_INLINE
  /* Hardware-divq path: keep the tight per-limb loop.  The Möller-
   * Granlund reciprocal-based sweep is a net regression here because
   * divq's variable-latency fast path for sub-64-bit divisors already
   * runs at ~2.3 ns/limb, and the per-call recip_init (one divq +
   * shifts) amortises poorly across the small (m <= ~30) limb counts
   * that arise in libwignernj's Racah sums. */
  for (i = a->size; i-- > 0;)
    r->words[i] = bigint_div128(rem, a->words[i], b, &rem);
#else
  /* No hardware divq (Algorithm-D fallback).  Precompute a Möller-
   * Granlund reciprocal once per call and sweep limbs with mul-high
   * + adjustment -- 1.5-5x faster per limb than the trial-quotient
   * Algorithm D used by bigint_div128 on this backend.  Threshold
   * profile-tuned: at m <= 2 the recip_init cost (one bigint_div128
   * + a leading-zero count) is comparable to the per-limb work it
   * saves; only at m >= 3 does the precompute amortise into a clean
   * win (1.13x on tiny d, 1.61x on full 64-bit d, climbing to 5.8x
   * by m = 250). */
  if (a->size <= 2) {
    for (i = a->size; i-- > 0;)
      r->words[i] = bigint_div128(rem, a->words[i], b, &rem);
  } else {
    bigint_recip64_t recip = bigint_recip64_init(b);
    for (i = a->size; i-- > 0;)
      r->words[i] = bigint_div128_recip(rem, a->words[i], b, recip, &rem);
  }
#endif
  r->size = a->size;
  bigint_trim(r);
}

/*
 * Modular inverse of an odd 64-bit value mod 2^64, computed via
 * Newton's iteration.  Each step doubles the number of correct bits;
 * the seed (3*d) ^ 2 is correct mod 2^5, and four iterations of
 * x <- x*(2 - d*x) take that to >= 2^80, ample headroom for the
 * 64-bit result.
 *
 * d must be odd; the inverse only exists in that case.
 */
static uint64_t modinv_u64_odd(uint64_t d) {
  uint64_t x = (3u * d) ^ 2u;
  int i;
  for (i = 0; i < 5; i++)
    x = x * (2u - d * x);
  return x;
}

/*
 * Exact division of a multi-word non-negative integer by a 64-bit
 * divisor.  Caller guarantees that d divides a evenly; if it does not
 * the result is silently corrupt.  Faster than bigint_div_u64 because
 * each per-limb operation is a multiplication by a precomputed
 * modular inverse plus one carry-propagation multiplication, never a
 * 128/64 hardware division (or its Algorithm-D fallback).
 *
 * Used by the Racah-sum Pass-2 ratio recurrence in wigner3j_exact /
 * wigner6j_exact / racah_6j_sum, where the running scaled bigint is
 * provably divisible by the per-step denominator product.  At
 * wigner6j_exact j=4000 profiling identified bigint_div_u64 as 42% of
 * total runtime; the Hensel-style implementation here is roughly 5x
 * faster per limb, so the per-call wall time drops accordingly.
 *
 * r may alias a.  For d == 1 the function copies a to r unchanged;
 * for a == 0 it sets r to zero.
 *
 * Algorithm:
 *   1. Factor d = 2^s * d_odd (s = trailing zero count).  2^s | a is
 *      implied by the exact-divisibility precondition, so a/d =
 *      (a >> s) / d_odd.
 *   2. Compute m = d_odd^{-1} mod 2^64 via Newton's iteration.
 *   3. Apply Hensel-style exact division with d_odd's inverse:
 *      iterate over limbs low-to-high, computing the quotient digit
 *      q_i = (a_i - borrow) * m mod 2^64, then propagating
 *      borrow_{i+1} = high_64(q_i * d_odd) + (a_i underflow flag).
 */
void bigint_div_u64_exact(bigint_t *r, const bigint_t *a, uint64_t d) {
  int s;
  uint64_t d_odd, m_inv, borrow;
  size_t i;

  if (a->size == 0) {
    bigint_set_zero(r);
    return;
  }
  if (d == 1) {
    bigint_copy(r, a);
    return;
  }
  /* For small dividends, the per-limb divq is faster than the
   * fixed-overhead Hensel scheme; fall back to the generic
   * remainder-tracking division.  Profile-tuned for the 9j inner
   * Racah sum where the running scaled bigint stays in the 1-3 limb
   * range across most s values. */
  if (a->size < 4) {
    bigint_div_u64(r, a, d);
    return;
  }

  /* Factor 2^s out of d. */
  s = 0;
  {
    uint64_t t = d;
    while ((t & 1u) == 0u) {
      t >>= 1;
      s++;
    }
    d_odd = t;
  }

  /* Copy / shift so r holds (a >> s).  bigint_shr_inplace handles
   * s == 0 as a no-op. */
  if (r != a)
    bigint_copy(r, a);
  if (s > 0)
    bigint_shr_inplace(r, s);

  /* If d_odd is 1, we're done after the shift. */
  if (d_odd == 1u)
    return;

  m_inv = modinv_u64_odd(d_odd);

  /* Hensel exact division by d_odd, in-place over r->words. */
  borrow = 0;
  for (i = 0; i < r->size; i++) {
    uint64_t under = (r->words[i] < borrow) ? 1u : 0u;
    uint64_t cur = r->words[i] - borrow;
    uint64_t q = cur * m_inv; /* mod 2^64 */
    uint64_t lo, hi;
    bigint_mul64x64(q, d_odd, &lo, &hi); /* lo == cur (mod 2^64) */
    (void)lo;
    r->words[i] = q;
    borrow = hi + under;
  }
  /* If d divides a exactly, the final borrow is 0 and r->words is
   * the quotient.  We trust the caller's contract here -- a non-zero
   * residual would only surface as a wrong test result, caught
   * downstream by the ctest sweep. */
  bigint_trim(r);
}

/*
 * In-place left-shift by k bits.  a <<= k.  Equivalent to a *= 2^k.
 * Used by bigint_mul_prime_pow{,_ws} as the p == 2 fast path: a single
 * O(words(a)) shift replaces k iterations of O(words(a)) bigint_mul_u64
 * (whose multiplier happens to be 2 but whose code path is the
 * generic per-word multiply-add).
 */
static void bigint_shl_inplace(bigint_t *a, int k) {
  int bw, bb;
  size_t old_size, i;
  if (a->size == 0 || k <= 0)
    return;
  bw = k / 64;
  bb = k % 64;
  old_size = a->size;
  /* Reserve up to bw + 1 extra words at the top.  bigint_trim below
   * removes any unused top word. */
  bigint_ensure(a, old_size + (size_t)bw + 1);
  if (bb == 0) {
    /* Pure word shift. */
    memmove(a->words + bw, a->words, old_size * sizeof(uint64_t));
    memset(a->words, 0, (size_t)bw * sizeof(uint64_t));
    a->size = old_size + (size_t)bw;
    bigint_trim(a);
    return;
  }
  /* Mixed bit/word shift.  Iterate top-down so each write is to an
   * index above every read, avoiding any read-after-write hazard
   * (the same-index aliasing in the bw == 0 case is well-defined
   * because the read happens before the write within a single
   * compound expression). */
  a->words[old_size + (size_t)bw] = a->words[old_size - 1] >> (64 - bb);
  for (i = old_size + (size_t)bw - 1; i > (size_t)bw; i--) {
    uint64_t hi = a->words[i - (size_t)bw];
    uint64_t lo = a->words[i - (size_t)bw - 1];
    a->words[i] = (hi << bb) | (lo >> (64 - bb));
  }
  a->words[bw] = a->words[0] << bb;
  if (bw > 0)
    memset(a->words, 0, (size_t)bw * sizeof(uint64_t));
  a->size = old_size + (size_t)bw + 1;
  bigint_trim(a);
}

/*
 * In-place right-shift by k bits.  a >>= k.  Equivalent to a /= 2^k
 * when a is exactly divisible (the only case that arises in
 * bigint_div_prime_pow's caller -- the floor-division semantics of an
 * unsigned right-shift agree with exact division here).
 */
static void bigint_shr_inplace(bigint_t *a, int k) {
  int bw, bb;
  size_t new_size, i;
  if (a->size == 0 || k <= 0)
    return;
  bw = k / 64;
  bb = k % 64;
  if ((size_t)bw >= a->size) {
    a->size = 0;
    return;
  }
  new_size = a->size - (size_t)bw;
  if (bb == 0) {
    memmove(a->words, a->words + bw, new_size * sizeof(uint64_t));
    a->size = new_size;
    bigint_trim(a);
    return;
  }
  /* Mixed bit/word right-shift; bottom-up so writes go to indices
   * <= the indices we read. */
  for (i = 0; i + 1 < new_size; i++) {
    uint64_t lo = a->words[i + (size_t)bw];
    uint64_t hi = a->words[i + (size_t)bw + 1];
    a->words[i] = (lo >> bb) | (hi << (64 - bb));
  }
  a->words[new_size - 1] = a->words[a->size - 1] >> bb;
  a->size = new_size;
  bigint_trim(a);
}

void bigint_mul_prime_pow(bigint_t *a, uint64_t p, int k) {
  bigint_t pw, base;
  int i;
  if (k <= 0)
    return;
  /* p == 2 fast path: a single in-place shift replaces k generic
   * bigint_mul_u64 multiplications.  In the LCM-scaling loop of the
   * Racah sum the p=2 exponent is the largest by far (Legendre's
   * formula gives v_2(N!) ~ N for large N), so this special case
   * carries most of the weight. */
  if (p == 2) {
    bigint_shl_inplace(a, k);
    return;
  }
  if (k <= 8) {
    for (i = 0; i < k; i++)
      bigint_mul_u64(a, a, p);
    return;
  }
  /* Binary exponentiation: compute p^k, then a *= p^k. */
  bigint_init(&pw);
  bigint_init(&base);
  bigint_set_u64(&pw, 1);
  bigint_set_u64(&base, p);
  while (k > 0) {
    if (k & 1)
      bigint_mul(&pw, &pw, &base);
    k >>= 1;
    if (k > 0)
      bigint_mul(&base, &base, &base);
  }
  bigint_mul(a, a, &pw);
  bigint_free(&pw);
  bigint_free(&base);
}

void bigint_div_prime_pow(bigint_t *a, uint64_t p, int k) {
  int i;
  if (p == 2) {
    bigint_shr_inplace(a, k);
    return;
  }
  for (i = 0; i < k; i++)
    bigint_div_u64(a, a, p);
}

/* ── workspace-aware variants ────────────────────────────────────────────── */

void bigint_ws_init(bigint_ws_t *ws) {
  bigint_init(&ws->mul_temp);
  bigint_init(&ws->pp_pw);
  bigint_init(&ws->pp_base);
  bigint_init(&ws->kar_scratch);
}

void bigint_ws_free(bigint_ws_t *ws) {
  bigint_free(&ws->mul_temp);
  bigint_free(&ws->pp_pw);
  bigint_free(&ws->pp_base);
  bigint_free(&ws->kar_scratch);
}

void bigint_ws_reserve(bigint_ws_t *ws, size_t max_words) {
  bigint_ensure(&ws->mul_temp, max_words);
  bigint_ensure(&ws->pp_pw, max_words);
  bigint_ensure(&ws->pp_base, max_words);
  /* Karatsuba scratch needs to cover the worst-case multiplication of
   * two max_words bigints; that bound dominates every other use. */
  bigint_ensure(&ws->kar_scratch, kar_scratch_words(max_words, max_words));
}

void bigint_mul_ws(bigint_t *r, const bigint_t *a, const bigint_t *b,
                   bigint_ws_t *ws) {
  bigint_t *dest;
  int aliased;
  if (a->size == 0 || b->size == 0) {
    bigint_set_zero(r);
    return;
  }
  aliased = (r == a || r == b);
  if (aliased) {
    bigint_ensure(&ws->mul_temp, a->size + b->size);
    memset(ws->mul_temp.words, 0, (a->size + b->size) * sizeof(uint64_t));
    ws->mul_temp.size = a->size + b->size;
    dest = &ws->mul_temp;
  } else {
    bigint_ensure(r, a->size + b->size);
    memset(r->words, 0, (a->size + b->size) * sizeof(uint64_t));
    r->size = a->size + b->size;
    dest = r;
  }
  mul_words_ws(dest->words, a->words, a->size, b->words, b->size,
               &ws->kar_scratch);
  bigint_trim(dest);
  if (aliased)
    bigint_copy(r, dest);
}

void bigint_mul_prime_pow_ws(bigint_t *a, uint64_t p, int k, bigint_ws_t *ws) {
  int i;
  if (k <= 0)
    return;
  if (p == 2) {
    bigint_shl_inplace(a, k);
    return;
  }
  if (k <= 8) {
    for (i = 0; i < k; i++)
      bigint_mul_u64(a, a, p);
    return;
  }
  /* Binary exponentiation using ws->pp_pw and ws->pp_base as scratch.
   * No allocations occur if ws was reserved to the worst-case size. */
  bigint_set_u64(&ws->pp_pw, 1);
  bigint_set_u64(&ws->pp_base, p);
  while (k > 0) {
    if (k & 1)
      bigint_mul_ws(&ws->pp_pw, &ws->pp_pw, &ws->pp_base, ws);
    k >>= 1;
    if (k > 0)
      bigint_mul_ws(&ws->pp_base, &ws->pp_base, &ws->pp_base, ws);
  }
  bigint_mul_ws(a, a, &ws->pp_pw, ws);
}

size_t bigint_words_for_factorial(int N) {
  /* Stirling: log2(N!) ≈ N*log2(N) - N*log2(e) + 0.5*log2(N) + log2(sqrt(2π))
   *
   * The returned size is generous: we double the Stirling bound so that
   * any LCM-sized bigint, plus the cross-product temporaries it produces
   * in bigint_mul_ws (whose size is operand1.size + operand2.size words),
   * fits without forcing the workspace to realloc inside the hot loop.
   */
  double bits, log2N;
  size_t base;
  if (N <= 1)
    return 8;
  log2N = log((double)N) * 1.4426950408889634;
  bits = (double)N * log2N - (double)N * 1.4426950408889634 + 0.5 * log2N +
         1.33 + 64.0;
  if (bits < 64.0)
    bits = 64.0;
  base = (size_t)(bits / 64.0) + 4;
  return 2 * base; /* headroom for mid-multiplication temporaries */
}

/* ── bit length ──────────────────────────────────────────────────────────── */

int bigint_bit_length(const bigint_t *a) {
  uint64_t top;
  int bits;
  if (a->size == 0)
    return 0;
  top = a->words[a->size - 1];
  bits = (int)(a->size - 1) * 64;
  while (top) {
    bits++;
    top >>= 1;
  }
  return bits;
}

/* ── IEEE 754 float conversion ───────────────────────────────────────────── */

/*
 * Extract the top `mant` bits (mant <= 64) of `a` as a uint64_t, applying
 * round-to-nearest-even.  Sets *out_exp so that result == mantissa * 2^out_exp.
 */
static uint64_t extract_mantissa(const bigint_t *a, int mant, int *out_exp) {
  int n, shift, ws, bs, tsticky_word, tsticky_off, i;
  uint64_t m, lo, hi, mask;
  int round_bit, sticky;

  n = bigint_bit_length(a);
  if (n == 0) {
    *out_exp = 0;
    return 0;
  }

  shift = n - mant; /* bits to discard; negative means number is short */
  *out_exp = shift;

  if (shift <= 0) {
    /* Number has fewer bits than mantissa; left-shift to fill mant slots.
     * -shift < 64 is guaranteed since mant <= 64 and n >= 1. */
    m = a->words[0];
    if (shift < 0)
      m <<= -shift;
    return m;
  }

  /* Extract bits [shift .. shift+mant-1] from a */
  ws = shift / 64;
  bs = shift % 64;
  if (bs == 0) {
    lo = (ws < (int)a->size) ? a->words[ws] : 0;
    m = lo;
  } else {
    lo = (ws < (int)a->size) ? a->words[ws] : 0;
    hi = (ws + 1 < (int)a->size) ? a->words[ws + 1] : 0;
    m = (lo >> bs) | (hi << (64 - bs));
  }
  if (mant < 64)
    m &= ((uint64_t)1 << mant) - 1;

  /* Round bit: bit at position shift-1 */
  {
    int rw = (shift - 1) / 64, rb = (shift - 1) % 64;
    round_bit = (rw < (int)a->size) ? (int)((a->words[rw] >> rb) & 1) : 0;
  }

  /* Sticky bit: OR of bits 0 .. shift-2 */
  sticky = 0;
  if (shift >= 2) {
    tsticky_word = (shift - 2) / 64;
    tsticky_off = (shift - 2) % 64;
    /* Bits 0..tsticky_off in word tsticky_word */
    mask = (tsticky_off == 63) ? UINT64_MAX
                               : (((uint64_t)1 << (tsticky_off + 1)) - 1);
    if ((int)a->size > tsticky_word)
      sticky = (a->words[tsticky_word] & mask) != 0;
    /* All lower words */
    for (i = 0; i < tsticky_word && !sticky; i++)
      sticky = (a->words[i] != 0);
  }

  /* Round-to-nearest-even */
  if (round_bit && (sticky || (m & 1))) {
    m++;
    if (mant < 64 && (m >> mant)) {
      m >>= 1;
      (*out_exp)++;
    } else if (mant == 64 && m == 0) {
      m = (uint64_t)1 << 63;
      (*out_exp)++;
    }
  }
  return m;
}

#if LDBL_MANT_DIG > 64
/*
 * Extract the 64 bits of `a` immediately below the top 64 bits, normalized
 * so that bit 63 of the result corresponds to bit (n-65) of `a`.
 * The positional value of the result is:  return_value * 2^(n-128).
 * Requires n > 64 and no rounding overflow in the caller's extract_mantissa.
 */
static uint64_t lo64_from_words(const bigint_t *a, int n) {
  int lo_top =
      n - 65; /* 0-indexed bit position of the top of the second group */
  if (lo_top < 64) {
    uint64_t w0 = (a->size > 0) ? a->words[0] : 0;
    return (lo_top < 63) ? (w0 << (63 - lo_top)) : w0;
  } else {
    int word_hi = lo_top / 64;
    int off_hi = lo_top % 64;
    uint64_t wh = (word_hi < (int)a->size) ? a->words[word_hi] : 0;
    uint64_t wl = (word_hi - 1 < (int)a->size) ? a->words[word_hi - 1] : 0;
    return (off_hi == 63) ? wh : ((wh << (63 - off_hi)) | (wl >> (off_hi + 1)));
  }
}
#endif

float bigint_to_float(const bigint_t *a) {
  int exp2;
  uint64_t m;
  if (a->size == 0)
    return 0.0f;
  m = extract_mantissa(a, FLT_MANT_DIG, &exp2);
  return ldexpf((float)m, exp2);
}

float bigint_frexpf(const bigint_t *a, int *out_exp) {
  int exp2;
  uint64_t m;
  if (a->size == 0) {
    *out_exp = 0;
    return 0.0f;
  }
  m = extract_mantissa(a, FLT_MANT_DIG, &exp2);
  *out_exp = exp2 + FLT_MANT_DIG;
  return ldexpf((float)m, -FLT_MANT_DIG);
}

double bigint_frexp(const bigint_t *a, int *out_exp) {
  int exp2;
  uint64_t m;
  if (a->size == 0) {
    *out_exp = 0;
    return 0.0;
  }
  m = extract_mantissa(a, DBL_MANT_DIG, &exp2);
  *out_exp = exp2 + DBL_MANT_DIG;
  return ldexp((double)m, -DBL_MANT_DIG);
}

long double bigint_frexpl(const bigint_t *a, int *out_exp) {
  int mant = (LDBL_MANT_DIG <= 64) ? LDBL_MANT_DIG : 64;
  int exp2;
  uint64_t m;
  if (a->size == 0) {
    *out_exp = 0;
    return 0.0L;
  }
  m = extract_mantissa(a, mant, &exp2);
  *out_exp = exp2 + mant;
#if LDBL_MANT_DIG > 64
  {
    /* For quad-precision long double: add a second 64-bit chunk below the
     * top 64. Contribution to the normalized mantissa is always m_lo * 2^-128
     * (since *out_exp = exp2+64 = n, so dividing by 2^n shifts by -n,
     *  and the lower bits sit at 2^(n-128), giving 2^(n-128)/2^n = 2^-128). */
    int n = bigint_bit_length(a);
    long double result = ldexpl((long double)m, -mant);
    if (n > 64 && exp2 == n - 64) {
      result += ldexpl((long double)lo64_from_words(a, n), -128);
    }
    return result;
  }
#else
  return ldexpl((long double)m, -mant);
#endif
}

double bigint_to_double(const bigint_t *a) {
  int exp2;
  uint64_t m;
  if (a->size == 0)
    return 0.0;
  m = extract_mantissa(a, DBL_MANT_DIG, &exp2);
  return ldexp((double)m, exp2);
}

#ifdef WIGNERNJ_HAVE_QUADMATH
/*
 * Build the float128 value by Horner-style evaluation of the top up to
 * three 64-bit words: r = w[size-1]*2^128 + w[size-2]*2^64 + w[size-3], scaled
 * by 2^((size-3)*64).  192 bits of input feed a 113-bit mantissa, so the
 * top-bit retention is unconditional and the bottom ~79 bits are truncated
 * by the float128 add itself (round-to-nearest-even).  Sticky bits below
 * the top three words are not preserved, so a tied half-ulp result may
 * round to even either way; this is within the 2-ulp tolerance the
 * verification harness uses for every precision.
 */
static __float128 bigint_top_float128(const bigint_t *a, int *out_exp) {
  if (a->size == 0) {
    *out_exp = 0;
    return 0.0Q;
  }
  size_t k = (a->size >= 3) ? (a->size - 3) : 0;
  __float128 r = 0.0Q;
  for (size_t i = a->size; i-- > k;) {
    r = r * 0x1.0p64Q + (__float128)a->words[i];
  }
  *out_exp = (int)k * 64;
  return r;
}

__float128 bigint_to_float128(const bigint_t *a) {
  int e;
  __float128 r = bigint_top_float128(a, &e);
  return ldexpq(r, e);
}

__float128 bigint_frexp_q(const bigint_t *a, int *out_exp) {
  int e;
  __float128 r = bigint_top_float128(a, &e);
  if (a->size == 0) {
    *out_exp = 0;
    return 0.0Q;
  }
  int e2;
  __float128 mant = frexpq(r, &e2);
  *out_exp = e + e2;
  return mant;
}
#endif /* WIGNERNJ_HAVE_QUADMATH */

long double bigint_to_long_double(const bigint_t *a) {
  int mant = (LDBL_MANT_DIG <= 64) ? LDBL_MANT_DIG : 64;
  int exp2;
  uint64_t m;
  if (a->size == 0)
    return 0.0L;
  m = extract_mantissa(a, mant, &exp2);
#if LDBL_MANT_DIG > 64
  {
    /* For quad-precision long double: add a second 64-bit chunk.
     * m's positional value is m * 2^exp2; the lower chunk adds m_lo *
     * 2^(n-128). */
    int n = bigint_bit_length(a);
    long double hi = ldexpl((long double)m, exp2);
    if (n > 64 && exp2 == n - 64) {
      hi += ldexpl((long double)lo64_from_words(a, n), n - 128);
    }
    return hi;
  }
#else
  return ldexpl((long double)m, exp2);
#endif
}

/* ── MPFR conversion ─────────────────────────────────────────────────────── */

#ifdef WIGNERNJ_HAVE_MPFR
void bigint_to_mpfr(mpfr_t rop, const bigint_t *a, mpfr_rnd_t rnd) {
  size_t i;
  if (a->size == 0) {
    mpfr_set_zero(rop, +1);
    return;
  }
  /* Accumulate most-significant word first: rop = rop*2^32 + half */
  mpfr_set_zero(rop, +1);
  for (i = a->size; i-- > 0;) {
    mpfr_mul_2ui(rop, rop, 32, rnd);
    mpfr_add_ui(rop, rop, (unsigned long)(a->words[i] >> 32), rnd);
    mpfr_mul_2ui(rop, rop, 32, rnd);
    mpfr_add_ui(rop, rop, (unsigned long)(a->words[i] & 0xFFFFFFFFUL), rnd);
  }
}
#endif
