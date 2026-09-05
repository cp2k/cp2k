/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola */
#include "pfrac.h"
#include "wignernj.h"
#include "wignernj_tls.h"
#include "xalloc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Public accessor: absolute compile-time ceiling on factorial arguments.
 * Determined by the prime sieve in src/prime_table.inc. */
int wignernj_max_factorial_arg(void) { return MAX_FACTORIAL_ARG; }

/* ── factorial-decomposition cache ────────────────────────────────────────────
 *
 * pfrac_mul_factorial / pfrac_div_factorial used to call
 * legendre_valuation(N, p_i) for every prime p_i <= N on every entry,
 * paying a per-call cost proportional to the number of primes below N
 * times the per-Legendre-valuation divide loop.  Because the prime
 * decomposition of N! depends only on N, this is wasteful: at the
 * scale of a single 3j evaluation the same N is decomposed dozens of
 * times across the Racah sum.
 *
 * The cache is a per-thread table indexed by N: g_fact_cache[N], when
 * non-NULL, points to a malloc'd array of fact_width(N) ints holding
 * the exponent of each g_primes[i] (i = 0 .. fact_width(N)-1) in N!.
 * pfrac_mul_factorial / pfrac_div_factorial are now a vector add /
 * subtract over that row instead of a per-prime divide loop.
 *
 * Lifecycle:  lazy-allocated on first use, persists for the lifetime
 * of the calling thread, never freed (the OS reclaims at thread exit).
 * Memory grows with the largest N that thread sees; at the absolute
 * default-build ceiling (MAX_FACTORIAL_ARG = 20020) the per-thread
 * cost is bounded by sum_{N=1}^{20020} fact_width(N) ints ~ 80 MB.
 * For typical workloads (j <~ a few thousand) it is a small fraction
 * of that.
 *
 * Toolchains without thread-local storage skip the cache entirely and
 * fall back to the per-call legendre_valuation loop -- still
 * thread-safe (no shared state), just slower.
 */

#if WIGNERNJ_HAVE_TLS

/* Number of primes <= n in g_primes.  O(1) lookup into the
 * compile-time table generated alongside g_primes.  Only used by the
 * TLS-cached factorial-decomposition cache; the no-TLS fallback walks
 * the prime list directly via g_primes[i] <= n. */
static int fact_width(int n) { return (int)g_pi_table[n]; }

/* Pointer table: indexed by N, lazy-allocated rows.  The table itself
 * is allocated lazily on first cache access. */
static WIGNERNJ_TLS int **g_fact_cache = NULL;

static void fact_cache_ensure_table(void) {
  if (!g_fact_cache)
    g_fact_cache =
        (int **)xcalloc((size_t)(MAX_FACTORIAL_ARG + 1), sizeof(int *));
}

static const int *fact_cache_get(int n) {
  fact_cache_ensure_table();
  if (g_fact_cache[n])
    return g_fact_cache[n];

  int width = fact_width(n);
  int *row = (int *)xmalloc((size_t)width * sizeof(int));
  int i;
  for (i = 0; i < width; i++)
    row[i] = legendre_valuation(n, i);
  g_fact_cache[n] = row;
  return row;
}

/* Internal: invoked by wignernj_warmup_to() in src/scratch.c.
 * N_max <= 0 means "use the absolute prime-table ceiling". */
void wignernj_factorial_cache_fill(int N_max) {
  int n;
  if (N_max <= 0 || N_max > MAX_FACTORIAL_ARG)
    N_max = MAX_FACTORIAL_ARG;
  fact_cache_ensure_table();
  for (n = 2; n <= N_max; n++) {
    if (!g_fact_cache[n])
      (void)fact_cache_get(n);
  }
}

/* Release the calling thread's factorial-decomposition cache.  Frees
 * every cached row plus the pointer table itself.  After return, the
 * next pfrac_mul_factorial / pfrac_div_factorial call in this thread
 * starts fresh (re-allocates the table on first use, rebuilds rows
 * on demand). */
void wignernj_factorial_cache_release(void) {
  int n;
  if (!g_fact_cache)
    return;
  for (n = 0; n <= MAX_FACTORIAL_ARG; n++) {
    free(g_fact_cache[n]);
  }
  free(g_fact_cache);
  g_fact_cache = NULL;
}

#else /* !WIGNERNJ_HAVE_TLS */

void wignernj_factorial_cache_fill(int N_max) {
  /* No persistent cache to populate -- the no-TLS path already
   * computes legendre_valuation per call.  No-op for source
   * compatibility with the cached path. */
  (void)N_max;
}

void wignernj_factorial_cache_release(void) {
  /* No persistent cache; nothing to free. */
}

#endif

/* ───────────────────────────────────────────────────────────────────────────
 */

static void pfrac_fatal(const char *msg) {
  fprintf(stderr, "libwignernj: %s\n", msg);
  abort();
}

void pfrac_init(pfrac_t *f) {
  f->exp = (int *)xcalloc((size_t)g_nprimes, sizeof(int));
  f->max_idx = 0;
}

void pfrac_free(pfrac_t *f) {
  free(f->exp);
  f->exp = NULL;
  f->max_idx = 0;
}

void pfrac_zero(pfrac_t *f) {
  /* Zero only [0, max_idx); the invariant guarantees the rest is zero. */
  if (f->max_idx > 0)
    memset(f->exp, 0, (size_t)f->max_idx * sizeof(int));
  f->max_idx = 0;
}

void pfrac_copy(pfrac_t *dst, const pfrac_t *src) {
  /* If dst was wider than src, zero the excess so the invariant holds. */
  if (dst->max_idx > src->max_idx) {
    memset(dst->exp + src->max_idx, 0,
           (size_t)(dst->max_idx - src->max_idx) * sizeof(int));
  }
  if (src->max_idx > 0)
    memcpy(dst->exp, src->exp, (size_t)src->max_idx * sizeof(int));
  dst->max_idx = src->max_idx;
}

void pfrac_mul_factorial(pfrac_t *f, int n) {
  int i, width;
  if (n <= 1)
    return;
  if (n > MAX_FACTORIAL_ARG) {
    fprintf(stderr,
            "libwignernj: factorial argument %d exceeds MAX_FACTORIAL_ARG=%d "
            "(angular momenta too large for prime table)\n",
            n, MAX_FACTORIAL_ARG);
    pfrac_fatal("aborting");
  }
#if WIGNERNJ_HAVE_TLS
  {
    /* `restrict`: row points into the per-thread factorial-
     * decomposition cache, ex points to the term's own buffer;
     * they are disjoint by construction.  Without restrict the
     * compiler must assume aliasing through `f->exp` (a struct
     * field reachable by pointer), which blocks the auto-
     * vectoriser at -O3.  Lifting f->exp to a restrict-qualified
     * local lets GCC emit 16-byte SSE2 vpaddd / vpsubd. */
    const int *restrict row = fact_cache_get(n);
    int *restrict ex = f->exp;
    width = fact_width(n);
    for (i = 0; i < width; i++)
      ex[i] += row[i];
  }
#else
  for (i = 0; i < g_nprimes && g_primes[i] <= n; i++)
    f->exp[i] += legendre_valuation(n, i);
  width = i;
#endif
  if (width > f->max_idx)
    f->max_idx = width;
}

void pfrac_div_factorial(pfrac_t *f, int n) {
  int i, width;
  if (n <= 1)
    return;
  if (n > MAX_FACTORIAL_ARG) {
    fprintf(stderr,
            "libwignernj: factorial argument %d exceeds MAX_FACTORIAL_ARG=%d "
            "(angular momenta too large for prime table)\n",
            n, MAX_FACTORIAL_ARG);
    pfrac_fatal("aborting");
  }
#if WIGNERNJ_HAVE_TLS
  {
    const int *restrict row = fact_cache_get(n);
    int *restrict ex = f->exp;
    width = fact_width(n);
    for (i = 0; i < width; i++)
      ex[i] -= row[i];
  }
#else
  for (i = 0; i < g_nprimes && g_primes[i] <= n; i++)
    f->exp[i] -= legendre_valuation(n, i);
  width = i;
#endif
  if (width > f->max_idx)
    f->max_idx = width;
}

void pfrac_mul_int(pfrac_t *f, int k) {
  int i;
  int last = -1;
  if (k <= 1)
    return;
  if (k > PRIME_SIEVE_LIMIT) {
    fprintf(
        stderr,
        "libwignernj: pfrac_mul_int argument %d exceeds PRIME_SIEVE_LIMIT=%d "
        "(angular momenta too large for prime table)\n",
        k, PRIME_SIEVE_LIMIT);
    pfrac_fatal("aborting");
  }
  /*
   * Trial-divide k by primes up to floor(sqrt(k)).  When the loop
   * exits with k > 1, the residual factor must itself be prime (it
   * has no divisor below sqrt(k_original)), so its index can be
   * looked up in O(1) via prime_index_of() rather than continuing
   * the trial-division loop.  For prime or near-prime k this turns
   * an O(pi(k)) sweep into an O(pi(sqrt(k))) one -- e.g. for
   * k = 1031 (prime, just above 32^2), 173 prime tests become 11.
   *
   * The bound check uses the *current* k (which decreases as we
   * divide it), so smooth k -- which collapse to 1 quickly -- exit
   * the loop as soon as no remaining factor is small enough to
   * matter, matching the speed of the original tight loop.
   */
  for (i = 0; i < g_nprimes; i++) {
    long long p = g_primes[i];
    if (p * p > (long long)k)
      break;
    while (k % g_primes[i] == 0) {
      f->exp[i]++;
      last = i;
      k /= g_primes[i];
    }
    if (k == 1)
      break;
  }
  if (k > 1) {
    int idx = prime_index_of(k);
    if (idx < 0)
      pfrac_fatal("pfrac_mul_int: residual factor not in prime table "
                  "(internal error)");
    f->exp[idx]++;
    if (idx > last)
      last = idx;
  }
  if (last + 1 > f->max_idx)
    f->max_idx = last + 1;
}

void pfrac_to_sqrt_rational(const pfrac_t *f, bigint_t *int_num,
                            bigint_t *int_den, bigint_t *sqrt_num,
                            bigint_t *sqrt_den) {
  /* Forward to the workspace-aware version with a one-shot stack
   * workspace.  The split logic and uint64-batched accumulators live
   * exclusively in pfrac_to_sqrt_rational_ws to avoid duplicating the
   * four-accumulator dance. */
  bigint_ws_t ws;
  bigint_ws_init(&ws);
  pfrac_to_sqrt_rational_ws(f, int_num, int_den, sqrt_num, sqrt_den, &ws);
  bigint_ws_free(&ws);
}

/*
 * Multiply prime power g_primes[pi]^exp into a bigint accumulator
 * pattern: a uint64_t acc collects small contributions and is flushed
 * to *dst via bigint_mul_u64 only when adding the next contribution
 * would overflow.  Falls back to bigint_mul_prime_pow_ws (which has a
 * shift fast path for p == 2 and a binary-exponentiation path for
 * larger exponents) when p^exp itself doesn't fit in uint64_t.
 *
 * Used by every routine in this file that needs to multiply a bigint
 * by a product of prime powers indexed by g_primes -- the per-call
 * function-call overhead of bigint_mul_u64 / bigint_mul_prime_pow_ws
 * was the dominant Pass-2 cost in wigner3j_exact before this batched
 * scheme replaced it (62.7% in profile at j=4000, profile-confirmed).
 */
void pfrac_mul_pow_into_acc(bigint_t *dst, uint64_t *acc, uint64_t p, int exp,
                            bigint_ws_t *ws) {
  uint64_t pp = 1;
  int e, fits = 1;
  if (exp <= 0)
    return;
  for (e = 0; e < exp; e++) {
    if (pp > UINT64_MAX / p) {
      fits = 0;
      break;
    }
    pp *= p;
  }
  if (!fits) {
    if (*acc > 1) {
      bigint_mul_u64(dst, dst, *acc);
      *acc = 1;
    }
    bigint_mul_prime_pow_ws(dst, p, exp, ws);
    return;
  }
  if (*acc > UINT64_MAX / pp) {
    bigint_mul_u64(dst, dst, *acc);
    *acc = pp;
  } else {
    *acc *= pp;
  }
}

void pfrac_lcm_scaled_product(bigint_t *scaled, const int *lcm,
                              const int *term_exp, int sign, int max_idx,
                              bigint_ws_t *ws) {
  int pi;
  uint64_t acc = 1;

  bigint_set_u64(scaled, 1);
  for (pi = 0; pi < max_idx; pi++) {
    int diff = lcm[pi] + sign * term_exp[pi];
    if (diff > 0)
      pfrac_mul_pow_into_acc(scaled, &acc, (uint64_t)g_primes[pi], diff, ws);
  }
  if (acc > 1)
    bigint_mul_u64(scaled, scaled, acc);
}

void pfrac_bigint_mul_prime_pow_array(bigint_t *a, const int *exp, int max_idx,
                                      bigint_ws_t *ws) {
  int pi;
  uint64_t acc = 1;
  for (pi = 0; pi < max_idx; pi++) {
    if (exp[pi] > 0)
      pfrac_mul_pow_into_acc(a, &acc, (uint64_t)g_primes[pi], exp[pi], ws);
  }
  if (acc > 1)
    bigint_mul_u64(a, a, acc);
}

void pfrac_to_sqrt_rational_ws(const pfrac_t *f, bigint_t *int_num,
                               bigint_t *int_den, bigint_t *sqrt_num,
                               bigint_t *sqrt_den, bigint_ws_t *ws) {
  /* Four parallel uint64 accumulators -- one per output bigint --
   * each batched and flushed via the same scheme as
   * pfrac_lcm_scaled_product.  Replaces a chain of up to 4*pi(N)
   * separate bigint_mul_prime_pow_ws calls (each with its own
   * binary-exponentiation setup, plus big*big multiplies whenever
   * p^k overflows uint64) with a small handful of bigint_mul_u64s
   * per output -- profile-confirmed as ~25% of total time at large
   * j before this batching landed. */
  int i, e, ae, half;
  uint64_t acc_in = 1, acc_id = 1, acc_sn = 1, acc_sd = 1;
  for (i = 0; i < f->max_idx; i++) {
    e = f->exp[i];
    if (e == 0)
      continue;
    ae = (e > 0) ? e : -e;
    half = ae / 2;
    uint64_t p = (uint64_t)g_primes[i];
    if (e > 0) {
      if (half > 0)
        pfrac_mul_pow_into_acc(int_num, &acc_in, p, half, ws);
      if (ae & 1)
        pfrac_mul_pow_into_acc(sqrt_num, &acc_sn, p, 1, ws);
    } else {
      if (half > 0)
        pfrac_mul_pow_into_acc(int_den, &acc_id, p, half, ws);
      if (ae & 1)
        pfrac_mul_pow_into_acc(sqrt_den, &acc_sd, p, 1, ws);
    }
  }
  if (acc_in > 1)
    bigint_mul_u64(int_num, int_num, acc_in);
  if (acc_id > 1)
    bigint_mul_u64(int_den, int_den, acc_id);
  if (acc_sn > 1)
    bigint_mul_u64(sqrt_num, sqrt_num, acc_sn);
  if (acc_sd > 1)
    bigint_mul_u64(sqrt_den, sqrt_den, acc_sd);
}
