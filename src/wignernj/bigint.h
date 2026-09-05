/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola */
#ifndef BIGINT_H
#define BIGINT_H

#include <stddef.h>
#include <stdint.h>

/*
 * Arbitrary-precision non-negative integer.  Sign is tracked externally
 * by callers (wignernj_exact_t).  Two interchangeable backends are
 * provided:
 *
 * - The default schoolbook backend (src/bigint.c) stores the value as
 *   a little-endian array of uint64_t base-2^64 words and uses
 *   schoolbook multiplication.  ``Little-endian'' here refers to the
 *   logical word ordering---words[0] holds the least-significant
 *   digit, words[size-1] the most-significant---and is independent of
 *   the host CPU's byte order, so the library behaves identically on
 *   little- and big-endian machines.  Invariant: words[size-1] != 0
 *   unless size == 0 (representing zero).
 *
 * - The optional FLINT backend (src/bigint_flint.c, selected by the
 *   CMake option WIGNERNJ_BUILD_FLINT) stores the value as a single fmpz_t and
 *   delegates every arithmetic operation to FLINT/GMP.  This trades a
 *   runtime dependency for sub-quadratic multiplication
 *   (Karatsuba/Toom-Cook/Schönhage--Strassen via FLINT) at very large
 *   angular momenta.  The two backends produce bit-identical
 *   floating-point output, verified in CI.
 *
 * The struct definition switches at compile time on
 * WIGNERNJ_USE_FLINT, which is set by the build system.  No public C
 * API consumer should depend on the layout of bigint_t.
 */
#ifdef WIGNERNJ_USE_FLINT
/* fmpz.h conditionally exposes fmpz_get_mpz (and other mpz/mpfr
 * conversion helpers) only when <gmp.h> has already been included,
 * so we pull GMP in first. */
#include <flint/fmpz.h>
#include <gmp.h>
typedef struct {
  fmpz_t v;
} bigint_t;
#else
typedef struct {
  uint64_t *words; /* logical little-endian base-2^64 digits */
  size_t size;     /* number of significant words            */
  size_t cap;      /* allocated capacity                     */
} bigint_t;
#endif

/* Lifecycle */
void bigint_init(bigint_t *a);
void bigint_free(bigint_t *a);
void bigint_copy(bigint_t *dst, const bigint_t *src);
void bigint_set_u64(bigint_t *a, uint64_t v);
void bigint_set_zero(bigint_t *a);
int bigint_is_zero(const bigint_t *a);
/* Pre-allocate at least `cap` words of capacity, without changing the
 * stored value or its size.  Use to avoid realloc churn when the maximum
 * size is known in advance. */
void bigint_reserve(bigint_t *a, size_t cap);

/* Comparison (unsigned magnitudes) */
int bigint_cmp(const bigint_t *a, const bigint_t *b); /* -1, 0, +1 */

/* Arithmetic (all operands are non-negative; sign is caller's problem) */
void bigint_add(bigint_t *r, const bigint_t *a, const bigint_t *b);
/* r = |a| - |b|, returns sign (+1 if a>=b, -1 if a<b). r may alias a or b. */
int bigint_sub_signed(bigint_t *r, const bigint_t *a, const bigint_t *b);
void bigint_mul(bigint_t *r, const bigint_t *a, const bigint_t *b);
void bigint_mul_u64(bigint_t *r, const bigint_t *a, uint64_t b);
/* Exact (no remainder) division by a small divisor. */
void bigint_div_u64(bigint_t *r, const bigint_t *a, uint64_t b);
/* Exact division by a 64-bit divisor: the dividend a must be exactly
 * divisible by d.  Uses Hensel-style modular-inverse arithmetic
 * (multiply by precomputed d^{-1} mod 2^64), avoiding the 128/64
 * hardware division per limb that bigint_div_u64 does -- typically
 * ~5x faster per limb when both inputs are multi-word.  r may alias a. */
void bigint_div_u64_exact(bigint_t *r, const bigint_t *a, uint64_t d);
/* Multiply a by p^k in-place (p is a prime that fits in uint64_t). */
void bigint_mul_prime_pow(bigint_t *a, uint64_t p, int k);
/* Exact division by p^k in-place. */
void bigint_div_prime_pow(bigint_t *a, uint64_t p, int k);

/*
 * Workspace for the workspace-aware multiplication helpers below.
 * Pre-allocating once at the start of a wignerXj_exact call eliminates the
 * malloc/realloc traffic that otherwise dominates the per-call cost.
 *
 * Lifecycle:  bigint_ws_init  ->  bigint_ws_reserve  ->  ...use...  ->
 * bigint_ws_free
 */
typedef struct {
  bigint_t mul_temp;    /* temp for bigint_mul_ws when destination is aliased */
  bigint_t pp_pw;       /* p^k accumulator in the binary-exponentiation loop */
  bigint_t pp_base;     /* current base in the binary-exponentiation loop  */
  bigint_t kar_scratch; /* Karatsuba scratch (schoolbook backend); unused
                           by the FLINT backend, where multiplication is
                           delegated to FLINT's own scratch management.   */
} bigint_ws_t;

void bigint_ws_init(bigint_ws_t *ws);
void bigint_ws_free(bigint_ws_t *ws);
/* Pre-allocate every workspace bigint to at least `max_words` words. */
void bigint_ws_reserve(bigint_ws_t *ws, size_t max_words);

/* In-place product: r = a * b, using ws->mul_temp as scratch when r aliases
 * either operand.  Bit-identical output to bigint_mul. */
void bigint_mul_ws(bigint_t *r, const bigint_t *a, const bigint_t *b,
                   bigint_ws_t *ws);
/* In-place a *= p^k, using ws->pp_pw and ws->pp_base as scratch.
 * Bit-identical output to bigint_mul_prime_pow. */
void bigint_mul_prime_pow_ws(bigint_t *a, uint64_t p, int k, bigint_ws_t *ws);

/*
 * Upper bound on the number of 64-bit words needed to represent N!, with
 * enough slack to also cover the LCM of (N+1) sum terms and one cross
 * product of two such bigints.  Use this to pre-size every long-lived
 * bigint inside a wignerXj_exact call so that no realloc occurs in the
 * inner Racah-sum loop.
 */
size_t bigint_words_for_factorial(int N);

/* Bit operations */
int bigint_bit_length(const bigint_t *a); /* floor(log2(a))+1, 0 for zero */

/* Conversion to floating-point.  IEEE 754 round-to-nearest-even. */
float bigint_to_float(const bigint_t *a);
double bigint_to_double(const bigint_t *a);
long double bigint_to_long_double(const bigint_t *a);

/*
 * Normalised frexp-style decomposition: returns m ∈ [0.5, 1) and sets
 * *out_exp such that  a ≈ m * 2^(*out_exp).
 * Returns 0.0 and *out_exp=0 for zero.  Rounds the mantissa to the
 * precision of the respective type (FLT_MANT_DIG / DBL_MANT_DIG /
 * LDBL_MANT_DIG). Used by wignernj_exact_to_* to avoid intermediate
 * overflow/underflow.
 */
float bigint_frexpf(const bigint_t *a, int *out_exp);
double bigint_frexp(const bigint_t *a, int *out_exp);
long double bigint_frexpl(const bigint_t *a, int *out_exp);

/* IEEE 754 binary128 ("quadruple precision", __float128).  Built only when
 * the compiler provides the type -- detected at configure time and exposed
 * through WIGNERNJ_HAVE_QUADMATH.  The conversion combines the top three
 * 64-bit words of `a` in __float128 arithmetic, which gives 192 bits of
 * raw precision for a 113-bit mantissa: the bottom ~79 bits are dropped,
 * which is enough headroom for the 2-ulp test tolerance used everywhere
 * in the verification harness. */
#ifdef WIGNERNJ_HAVE_QUADMATH
#include <quadmath.h>
__float128 bigint_to_float128(const bigint_t *a);
__float128 bigint_frexp_q(const bigint_t *a, int *out_exp);
#endif

/* MPFR conversion (only available when compiled with WIGNERNJ_HAVE_MPFR). */
#ifdef WIGNERNJ_HAVE_MPFR
#include <mpfr.h>
void bigint_to_mpfr(mpfr_t rop, const bigint_t *a, mpfr_rnd_t rnd);
#endif

#endif /* BIGINT_H */
