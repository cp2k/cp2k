/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola */
#ifndef PFRAC_H
#define PFRAC_H

#include "bigint.h"
#include "primes.h"

/*
 * Prime-factored exact representation.
 *
 * A pfrac_t stores, for each prime p_i in g_primes[], the signed exponent
 * of p_i in the represented value.  Positive exponents are in the numerator;
 * negative exponents are in the denominator.
 *
 * Two distinct uses:
 *   - "rational" pfrac: represents an exact rational number q = prod
 * p_i^exp[i].
 *   - "sqrt" pfrac:     represents sqrt(q) = prod p_i^(exp[i]/2).
 *     In both cases the same struct is used; the interpretation is the
 * caller's.
 *
 * exp[] is sized to g_nprimes and zero-initialised at allocation.
 *
 * max_idx is an upper bound on the largest prime index that has been
 * touched by any pfrac_mul_factorial / pfrac_div_factorial / pfrac_mul_int
 * / pfrac_copy operation since the last pfrac_zero.  All entries with
 * i >= max_idx are guaranteed to be zero, so loops over the pfrac may
 * iterate to max_idx instead of g_nprimes.  For angular momenta well below
 * MAX_FACTORIAL_ARG this is a substantial reduction (max_idx ~ pi(j) for
 * moderate j, vs g_nprimes = MAX_PRIME_COUNT ~ 2.3k for the default
 * sieve limit).
 */
typedef struct {
  int *exp;    /* exp[i] = exponent of g_primes[i]; len = g_nprimes */
  int max_idx; /* exp[i] == 0 for all i >= max_idx (invariant) */
} pfrac_t;

void pfrac_init(pfrac_t *f);
void pfrac_free(pfrac_t *f);
void pfrac_zero(pfrac_t *f); /* set all exponents to 0 */
void pfrac_copy(pfrac_t *dst, const pfrac_t *src);

/* Multiply the rational by n! (adds Legendre valuations to exp[]).  */
void pfrac_mul_factorial(pfrac_t *f, int n);

/* Divide the rational by n! (subtracts Legendre valuations from exp[]). */
void pfrac_div_factorial(pfrac_t *f, int n);

/* Multiply the rational by the positive integer k.
 * k must have no prime factor exceeding PRIME_SIEVE_LIMIT. */
void pfrac_mul_int(pfrac_t *f, int k);

/*
 * Convert a pfrac_t representing sqrt(Q) into four bigint components such
 * that   sqrt(Q) = (int_num / int_den) * sqrt(sqrt_num / sqrt_den).
 *
 * For each prime p_i with exponent e = f->exp[i]:
 *   If e > 0:  p_i^floor(e/2) → int_num;  if e is odd, p_i → sqrt_num.
 *   If e < 0:  p_i^floor(|e|/2) → int_den; if |e| is odd, p_i → sqrt_den.
 *   If e == 0: skipped.
 *
 * All four bigints are initialised by the caller; this function multiplies
 * into them (does not reset to 1).  The caller must initialise them to 1
 * (bigint_set_u64(x, 1)) before calling.
 */
void pfrac_to_sqrt_rational(const pfrac_t *f, bigint_t *int_num,
                            bigint_t *int_den, bigint_t *sqrt_num,
                            bigint_t *sqrt_den);

/* Workspace-aware variant: identical output to pfrac_to_sqrt_rational, but
 * uses ws->pp_pw / ws->pp_base / ws->mul_temp as scratch for prime-power
 * binary exponentiation, avoiding the per-call allocations of the original. */
void pfrac_to_sqrt_rational_ws(const pfrac_t *f, bigint_t *int_num,
                               bigint_t *int_den, bigint_t *sqrt_num,
                               bigint_t *sqrt_den, bigint_ws_t *ws);

/*
 * Multiply prime power g_primes[pi]^exp into a bigint accumulator
 * pattern.  A uint64_t *acc collects small contributions and is
 * flushed to *dst via bigint_mul_u64 only when overflow is imminent;
 * larger contributions fall back to bigint_mul_prime_pow_ws.  Exposed
 * for callers (gaunt.c) whose per-prime exponent is computed from a
 * combination of arrays that doesn't fit pfrac_bigint_mul_prime_pow_array's
 * single-array signature; the simpler-array case should call that
 * helper instead.
 *
 * Caller must initialise *acc = 1 before the loop and flush a final
 * `if (*acc > 1) bigint_mul_u64(dst, dst, *acc)` after.
 */
void pfrac_mul_pow_into_acc(bigint_t *dst, uint64_t *acc, uint64_t p, int exp,
                            bigint_ws_t *ws);

/*
 * Multiply a (in place) by prod_{i: exp[i] > 0} g_primes[i]^exp[i].
 * Uses the same uint64-batched accumulator pattern as
 * pfrac_lcm_scaled_product: small contributions are folded into a
 * uint64_t and only flushed to the bigint when the accumulator would
 * overflow.  Falls back to bigint_mul_prime_pow_ws (with shift fast
 * path for p == 2 and binary-exponentiation for larger exponents)
 * when an individual p^exp doesn't fit in uint64_t.
 *
 * Replaces the simpler `for (pi) if (exp[pi] > 0) bigint_mul_prime_pow_ws(...)`
 * pattern at the end-of-function int_den construction in every Racah
 * exact() routine -- profile identified the per-prime call chain
 * there as the second-largest residual cost at large j.
 */
void pfrac_bigint_mul_prime_pow_array(bigint_t *a, const int *exp, int max_idx,
                                      bigint_ws_t *ws);

/*
 * Build the per-term scaling bigint scaled = LCM / term_denominator that
 * the Racah-sum pass-2 inner loop multiplies into the running sum.  For
 * each prime g_primes[i] with effective exponent
 *
 *   diff = lcm[i] + sign * term_exp[i]   (sign = -1 for 3j-style sums where
 *                                          term_exp tracks denominator
 *                                          factorials directly; sign = +1
 *                                          for 6j/9j-style sums where
 *                                          term_exp is the signed net
 *                                          exponent of the rational term)
 *
 * the contribution g_primes[i]^diff is folded into a uint64_t accumulator
 * when it fits, flushing to `scaled` via a single bigint_mul_u64 only when
 * the accumulator would overflow.  Replaces a chain of up to max_idx
 * separate bigint_mul_prime_pow_ws calls with a small handful of
 * bigint_mul_u64s, eliminating the per-call binary-exponentiation setup
 * and bigint sweep that dominated the inner loop at small to moderate j.
 *
 * `scaled` is set to 1 first.  Output is identical bit-for-bit to the
 * unbatched loop.
 */
void pfrac_lcm_scaled_product(bigint_t *scaled, const int *lcm,
                              const int *term_exp, int sign, int max_idx,
                              bigint_ws_t *ws);

#endif /* PFRAC_H */
