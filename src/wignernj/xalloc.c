/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola */
#include "xalloc.h"

#include <stdio.h>
#include <stdlib.h>

/* -1 means "injector disabled".  Otherwise the value is the number of
 * allocations that will still succeed before the next one is forced to fail. */
static long g_fail_countdown = -1;

void xalloc_set_test_failure_countdown(long n) { g_fail_countdown = n; }

static int xalloc_should_fail(void) {
  if (g_fail_countdown < 0)
    return 0;
  if (g_fail_countdown == 0) {
    /* Latch failure so every subsequent call also fails until the test
     * resets the countdown.  This makes behaviour deterministic when one
     * library call performs many allocations after the trigger. */
    return 1;
  }
  g_fail_countdown--;
  return 0;
}

static void xalloc_die(const char *what, size_t n) {
  fprintf(stderr, "libwignernj: out of memory (%s, %zu bytes)\n", what, n);
  abort();
}

void *xmalloc(size_t n) {
  void *p;
  if (xalloc_should_fail())
    xalloc_die("xmalloc", n);
  p = malloc(n);
  if (!p)
    xalloc_die("xmalloc", n);
  return p;
}

void *xcalloc(size_t nmemb, size_t size) {
  void *p;
  if (xalloc_should_fail())
    xalloc_die("xcalloc", nmemb * size);
  p = calloc(nmemb, size);
  if (!p)
    xalloc_die("xcalloc", nmemb * size);
  return p;
}

void *xrealloc(void *p, size_t n) {
  void *q;
  if (xalloc_should_fail())
    xalloc_die("xrealloc", n);
  q = realloc(p, n);
  if (!q)
    xalloc_die("xrealloc", n);
  return q;
}
