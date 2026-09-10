/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * Thread-safe scratch reused across every public-API call.
 *
 * Background: each public entry point (wigner3j, wigner6j, ...) used to
 * call bigint_init/reserve/free, pfrac_init/free, and a couple of
 * xcalloc-backed LCM-exponent arrays per call.  Across all symbol
 * families this added up to ~18-19 allocations per call, none of which
 * scaled with the angular-momentum input.  At small j (j <~ 10) this
 * allocation overhead dominated the wall-clock cost.
 *
 * Design: every public entry point fetches a wignernj_scratch_t via
 * wignernj_scratch_acquire(), uses it, and returns it via
 * wignernj_scratch_relinquish().  Two implementations of this acquire /
 * relinquish pair coexist:
 *
 *   - When the compiler provides thread-local storage (essentially
 *     every modern toolchain: GCC, Clang, MSVC, Intel, plus any
 *     C11-conforming compiler), the scratch is cached in a TLS slot.
 *     acquire returns the calling thread's cached scratch (lazy-
 *     allocated on first call, kept resident afterwards) and
 *     relinquish is a no-op.  Allocations after the first call are
 *     zero.  This is the fast path.
 *
 *   - When TLS is unavailable, acquire xcalloc's a fresh scratch and
 *     relinquish frees it.  This regresses to the historical per-call
 *     allocation cost, but is thread-safe by construction (each call
 *     owns its own scratch with no shared mutable state).
 *
 * Both code paths are selected entirely at compile time -- the call
 * sites in wigner3j.c et al. are identical.
 *
 * No public API change: wigner3j(), wigner6j(), etc. signatures and
 * behavior are unchanged.
 *
 * The xcalloc'd lcm_exp[] arrays are sized once to g_nprimes ints
 * each; lcm_max_dirty[] tracks how much of each array the previous
 * caller wrote, so subsequent callers only memset the dirty prefix.
 */
#ifndef WIGNERNJ_SCRATCH_H
#define WIGNERNJ_SCRATCH_H

#include "bigint.h"
#include "pfrac.h"
#include "wignernj_exact.h"

/* Slot counts sized for the most demanding caller in each category
 * (wigner9j for bigints, lcm_exp, and pfracs). */
#define WIGNERNJ_SCRATCH_BIGINTS 10
#define WIGNERNJ_SCRATCH_PFRACS 2
#define WIGNERNJ_SCRATCH_LCMEXP 4

typedef struct {
  bigint_ws_t ws;
  bigint_t bigints[WIGNERNJ_SCRATCH_BIGINTS];
  pfrac_t pfracs[WIGNERNJ_SCRATCH_PFRACS];
  int *lcm_exp[WIGNERNJ_SCRATCH_LCMEXP];
  int lcm_max_dirty[WIGNERNJ_SCRATCH_LCMEXP];
  /* Variable-sized cache of per-Racah-sum-term pfracs.  Lazy-grown to
   * the largest sum a given thread has seen.  Each Racah-sum loop
   * builds its term pfracs into [0, n_terms) once during pass 1, then
   * pass 2 walks them without rebuilding.  Capacity persists across
   * calls; pfracs at indices < terms_cap are pfrac_init'd and ready
   * to be pfrac_zero'd by the caller. */
  pfrac_t *terms;
  int terms_cap;
  /* Cached wignernj_exact_t reused across calls.  Lifecycle: init'd
   * once at scratch creation, the public wrappers reset it before
   * each call, and it is destroyed when the scratch is released. */
  wignernj_exact_t exact;
} wignernj_scratch_t;

/* Acquire a scratch for the duration of one public-API call.  Always
 * returns a valid pointer.  Pair every acquire with one relinquish. */
wignernj_scratch_t *wignernj_scratch_acquire(void);
void wignernj_scratch_relinquish(wignernj_scratch_t *s);

/* Helpers: zero only the dirty prefix of lcm_exp[idx], and update the
 * dirty bound when a new caller writes further than the last one. */
void wignernj_scratch_lcm_clear(wignernj_scratch_t *s, int idx);
void wignernj_scratch_lcm_dirty(wignernj_scratch_t *s, int idx, int new_max);

/* Ensure the scratch's term-pfrac cache holds at least n_terms slots,
 * growing in place if necessary.  Newly created slots are pfrac_init'd
 * (caller is responsible for pfrac_zero before use).  Capacity persists
 * across calls; this is how Racah-sum loops avoid rebuilding the
 * per-term pfrac in pass 2. */
void wignernj_scratch_terms_reserve(wignernj_scratch_t *s, int n_terms);

/* Test-only.  In the TLS-cached build, drops the calling thread's
 * cached scratch so the next acquire goes through the lazy-init path
 * again -- needed by test_oom for malloc-failure injection.  In the
 * no-TLS fallback every acquire already reallocates, so this is a
 * no-op. */
void wignernj_scratch_release(void);

/* Internal.  Returns 1 when the calling thread has access to the
 * per-thread scratch cache (TLS path) and 0 when each public-API
 * call allocates fresh (no-TLS fallback).  The build-time macro
 * WIGNERNJ_HAVE_TLS in src/wignernj_tls.h is the same answer at
 * compile time; this runtime accessor exists for tests that want to
 * branch without a preprocessor guard. */
int wignernj_thread_local_scratch_available(void);

#endif /* WIGNERNJ_SCRATCH_H */
