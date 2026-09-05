/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include "cp2k_wignernj_prefix.h"

/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola */
#include "primes.h"
#include "scratch.h"
#include "wignernj.h"
#include "wignernj_tls.h"
#include "xalloc.h"
#include <stdlib.h>
#include <string.h>

/* ── shared helpers ──────────────────────────────────────────────────────── */

static void scratch_init(wignernj_scratch_t *s) {
  int i;
  bigint_ws_init(&s->ws);
  for (i = 0; i < WIGNERNJ_SCRATCH_BIGINTS; i++)
    bigint_init(&s->bigints[i]);
  for (i = 0; i < WIGNERNJ_SCRATCH_PFRACS; i++)
    pfrac_init(&s->pfracs[i]);
  for (i = 0; i < WIGNERNJ_SCRATCH_LCMEXP; i++) {
    s->lcm_exp[i] = (int *)xcalloc((size_t)g_nprimes, sizeof(int));
    s->lcm_max_dirty[i] = 0;
  }
  s->terms = NULL;
  s->terms_cap = 0;
  wignernj_exact_init(&s->exact);
}

static void scratch_destroy(wignernj_scratch_t *s) {
  int i;
  bigint_ws_free(&s->ws);
  for (i = 0; i < WIGNERNJ_SCRATCH_BIGINTS; i++)
    bigint_free(&s->bigints[i]);
  for (i = 0; i < WIGNERNJ_SCRATCH_PFRACS; i++)
    pfrac_free(&s->pfracs[i]);
  for (i = 0; i < WIGNERNJ_SCRATCH_LCMEXP; i++)
    free(s->lcm_exp[i]);
  for (i = 0; i < s->terms_cap; i++)
    pfrac_free(&s->terms[i]);
  free(s->terms);
  s->terms = NULL;
  s->terms_cap = 0;
  /* Mark the cached exact as never-zero so wignernj_exact_free actually
   * frees its bigints (the is_zero short-circuit would otherwise leak). */
  s->exact.is_zero = 0;
  wignernj_exact_free(&s->exact);
}

void wignernj_scratch_lcm_clear(wignernj_scratch_t *s, int idx) {
  if (s->lcm_max_dirty[idx] > 0)
    memset(s->lcm_exp[idx], 0, (size_t)s->lcm_max_dirty[idx] * sizeof(int));
  s->lcm_max_dirty[idx] = 0;
}

void wignernj_scratch_lcm_dirty(wignernj_scratch_t *s, int idx, int new_max) {
  if (new_max > s->lcm_max_dirty[idx])
    s->lcm_max_dirty[idx] = new_max;
}

void wignernj_scratch_terms_reserve(wignernj_scratch_t *s, int n_terms) {
  int i;
  if (n_terms <= s->terms_cap)
    return;
  s->terms = (pfrac_t *)xrealloc(s->terms, (size_t)n_terms * sizeof(pfrac_t));
  for (i = s->terms_cap; i < n_terms; i++)
    pfrac_init(&s->terms[i]);
  s->terms_cap = n_terms;
}

/* ── TLS-cached path (the common case) ──────────────────────────────────────
 */

#if WIGNERNJ_HAVE_TLS

static WIGNERNJ_TLS wignernj_scratch_t *g_scratch = NULL;

wignernj_scratch_t *wignernj_scratch_acquire(void) {
  wignernj_scratch_t *s = g_scratch;
  if (s)
    return s;
  s = (wignernj_scratch_t *)xcalloc(1, sizeof(*s));
  scratch_init(s);
  g_scratch = s;
  return s;
}

void wignernj_scratch_relinquish(wignernj_scratch_t *s) {
  /* Cached path: nothing to do, the scratch persists for the
   * lifetime of the calling thread.  Argument unused. */
  (void)s;
}

void wignernj_scratch_release(void) {
  wignernj_scratch_t *s = g_scratch;
  if (!s)
    return;
  scratch_destroy(s);
  free(s);
  g_scratch = NULL;
}

/* ── per-call fallback (no TLS available) ───────────────────────────────────
 */

#else /* !WIGNERNJ_HAVE_TLS */

wignernj_scratch_t *wignernj_scratch_acquire(void) {
  wignernj_scratch_t *s = (wignernj_scratch_t *)xcalloc(1, sizeof(*s));
  scratch_init(s);
  return s;
}

void wignernj_scratch_relinquish(wignernj_scratch_t *s) {
  if (!s)
    return;
  scratch_destroy(s);
  free(s);
}

void wignernj_scratch_release(void) {
  /* No persistent state to drop.  Each acquire/relinquish pair is
   * already independent and thread-safe by construction. */
}

#endif /* WIGNERNJ_HAVE_TLS */

/* ── public warmup ──────────────────────────────────────────────────────── */

/* Internal: pre-fill the factorial-decomposition cache for arguments
 * 2..N_max.  Implemented in src/pfrac.c.  Not part of the public API. */
extern void wignernj_factorial_cache_fill(int N_max);

void wignernj_warmup_to(int N_max) {
#if WIGNERNJ_HAVE_TLS
  /* Pre-grow every cached buffer to fit factorial arguments up to
   * N_max.  Pass 0 to size to the absolute prime-table ceiling.
   *
   * The dominant factor is mw = bigint_words_for_factorial(N), which
   * is the size of the largest factorial the request will index.
   * 9j and Gaunt size their long-lived bigints to mw_prod = 5*mw to
   * cover triple-product cross terms; we use that as the universal
   * upper bound.
   *
   * After this call, every subsequent symbol evaluation in the
   * calling thread whose worst-case factorial argument is <= N is
   * allocation-free in this thread, regardless of which symbol
   * family is invoked. */
  int N = (N_max <= 0 || N_max > MAX_FACTORIAL_ARG) ? MAX_FACTORIAL_ARG : N_max;
  wignernj_scratch_t *s = wignernj_scratch_acquire();
  size_t mw = bigint_words_for_factorial(N);
  size_t mw_prod = 5 * mw;
  int i;

  bigint_ws_reserve(&s->ws, mw_prod);
  for (i = 0; i < WIGNERNJ_SCRATCH_BIGINTS; i++)
    bigint_reserve(&s->bigints[i], mw_prod);
  bigint_reserve(&s->exact.sum, mw_prod);
  bigint_reserve(&s->exact.int_num, mw_prod);
  bigint_reserve(&s->exact.int_den, mw_prod);
  bigint_reserve(&s->exact.sqrt_num, mw_prod);
  bigint_reserve(&s->exact.sqrt_den, mw_prod);

  /* Term-pfrac cache: the longest Racah sum any symbol with
   * factorials bounded by N can produce has at most N+1 terms (the
   * 6j sum range, bounded by half the sum of four j's).  Pre-
   * allocating that many pfrac_t slots costs ~g_nprimes ints per
   * slot. */
  wignernj_scratch_terms_reserve(s, N + 1);

  /* lcm_exp arrays are already sized to g_nprimes ints in
   * scratch_init, which is the worst-case for the prime table.
   * Pfracs are similarly sized once on the first
   * pfrac_mul_factorial call. */

  wignernj_scratch_relinquish(s);

  /* Pre-fill the factorial-decomposition cache, the second
   * per-thread cache that the public-API call path consults. */
  wignernj_factorial_cache_fill(N_max);
#else
  /* No-TLS fallback: every public-API call already allocates a fresh
   * scratch and frees it on return, so there is no persistent state
   * for the warmup to populate.  The function is therefore a no-op on
   * toolchains without thread-local storage; correctness is unaffected,
   * but the caller will continue to pay the per-call allocation cost. */
  (void)N_max;
#endif
}

int wignernj_thread_local_scratch_available(void) {
#if WIGNERNJ_HAVE_TLS
  return 1;
#else
  return 0;
#endif
}

/* Forward declaration -- implemented in src/pfrac.c.  Not exposed in
 * any internal header because it is a pure cleanup hook. */
extern void wignernj_factorial_cache_release(void);

void wignernj_thread_cleanup(void) {
  /* Release every per-thread cache held by the calling thread:
   * the per-symbol scratch (bigint workspace, exact-output bigints,
   * pfrac scratch, lcm-exponent arrays) and the factorial-decomposition
   * cache (per-N prime-exponent rows + the pointer table).  No-op
   * on no-TLS toolchains. */
  wignernj_scratch_release();
  wignernj_factorial_cache_release();
}
