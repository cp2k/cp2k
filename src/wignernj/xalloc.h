/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola */
#ifndef XALLOC_H
#define XALLOC_H

#include <stddef.h>

/*
 * Allocation wrappers that abort with a diagnostic on failure.
 *
 * The library has no recovery path for OOM: every allocation feeds an
 * exact-arithmetic computation whose intermediate state is invalid if any
 * single allocation is missing.  Centralising the failure in xmalloc/xcalloc/
 * xrealloc lets the rest of the code dereference the result unconditionally
 * and gives the user a stable, testable failure mode (write to stderr,
 * abort()) instead of a NULL deref segfault.
 *
 * The diagnostic is a single line of the form
 *     libwignernj: out of memory (xmalloc, NN bytes)
 * written to stderr before abort().
 *
 * For testing, a malloc-failure injector can be installed via
 * xalloc_set_test_failure_countdown(N): the next N allocations succeed and
 * the (N+1)-th returns NULL, triggering the abort path.  N=-1 disables the
 * injector (the default).  This entry point is intended only for the test
 * suite; production code must not call it.
 */

void *xmalloc(size_t n);
void *xcalloc(size_t nmemb, size_t size);
void *xrealloc(void *p, size_t n);

/* Test-only allocation-failure injector.  See header comment. */
void xalloc_set_test_failure_countdown(long n);

#endif /* XALLOC_H */
