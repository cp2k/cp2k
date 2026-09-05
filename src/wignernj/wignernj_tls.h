/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * Compiler-portable thread-local-storage keyword detection.
 *
 * Used by both src/scratch.c (the per-call cached scratch) and
 * src/pfrac.c (the per-thread factorial-decomposition cache).
 *
 *   WIGNERNJ_TLS         -- the storage-class keyword to use, or empty
 *                           if no TLS is available.
 *   WIGNERNJ_HAVE_TLS    -- 1 when the calling toolchain supports TLS
 *                           via one of the detected keywords, 0 when
 *                           the empty fallback was selected.
 *
 * The detection prefers compiler-specific keywords (which work pre-C11)
 * over the C11 _Thread_local for compatibility with older toolchains.
 * Every modern compiler on every libwignernj-supported target falls
 * through one of the first three branches; the empty fallback is dead
 * code in practice but is retained so the library still compiles in
 * thread-unaware mode.
 */
#ifndef WIGNERNJ_TLS_H
#define WIGNERNJ_TLS_H

#if defined(WIGNERNJ_FORCE_NO_TLS)
/* Build override: pretend the toolchain has no TLS keyword.  Used by
 * the dedicated CI cell that exercises the per-call-allocation
 * fallback path on a compiler that would otherwise pick TLS. */
#define WIGNERNJ_TLS
#define WIGNERNJ_HAVE_TLS 0
#elif defined(_MSC_VER)
#define WIGNERNJ_TLS __declspec(thread)
#define WIGNERNJ_HAVE_TLS 1
#elif defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER) ||  \
    defined(__INTEL_LLVM_COMPILER)
#define WIGNERNJ_TLS __thread
#define WIGNERNJ_HAVE_TLS 1
#elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L &&              \
    !defined(__STDC_NO_THREADS__)
#define WIGNERNJ_TLS _Thread_local
#define WIGNERNJ_HAVE_TLS 1
#else
#define WIGNERNJ_TLS
#define WIGNERNJ_HAVE_TLS 0
#endif

#endif /* WIGNERNJ_TLS_H */
