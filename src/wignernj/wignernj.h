/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * Public C API for libwignernj.
 *
 * Citation: if libwignernj contributes to published work, please cite
 *   S. Lehtola, "libwignernj: a reusable C/C++/Fortran/Python library
 *   for exact Wigner symbols and related coefficients", arXiv:2605.06634
 *   (2026), doi:10.48550/arXiv.2605.06634.
 *
 * Every coupling-coefficient routine (3j, 6j, 9j, Clebsch-Gordan,
 * Racah W, Fano X, Gaunt, real-Gaunt) takes its angular-momentum
 * arguments as TWICE their value so half-integer values stay
 * representable as plain integers:
 *   j = 3/2  →  tj = 3
 *   m = -1/2 →  tm = -1
 * The single exception is wignernj_real_ylm_in_complex_ylm, whose
 * argument is always an integer orbital angular momentum and is
 * therefore taken as plain `l` rather than `2*l`.
 *
 * Functions returning 0 / 0.0f / 0.0L for symbols that vanish by selection
 * rules are not errors; only programmer-visible violations (wrong parity,
 * |m| > j) are treated as returning zero without any diagnostics.
 *
 * Phase conventions
 * -----------------
 * The Wigner 3j, 6j, and 9j symbols, the Clebsch–Gordan coefficient, the
 * Racah W coefficient, and the Fano X-coefficient are pure algebraic
 * SU(2) objects.  Their values
 * are fixed entirely by the Racah/Wigner combinatorial formulas; no
 * spherical-harmonic phase convention enters anywhere.  The Clebsch–Gordan
 * coefficient is defined with the Condon–Shortley sign convention used by
 * Edmonds (1957) and Varshalovich et al. (1988):
 *   <j1 m1; j2 m2 | J M> = (-1)^(j1-j2+M) sqrt(2J+1) * 3j(j1,j2,J;m1,m2,-M),
 * which makes every Clebsch–Gordan coefficient real.
 *
 * The phase convention for Y_l^m enters the Gaunt routines below
 * (gaunt and gaunt_real), since both are defined as integrals of
 * three spherical harmonics, and the real ↔ complex Y_lm
 * basis-overlap matrix wignernj_real_ylm_in_complex_ylm, whose
 * entries are the inner products <Y_l^{m_c} | S_{l, m_r}>.  We
 * adopt the Condon–Shortley phase for Y_l^m,
 *   Y_l^m(θ,φ) = (-1)^m sqrt[(2l+1)/(4π) · (l-m)!/(l+m)!] P_l^m(cos θ) e^{i m
 * φ} (m ≥ 0), Y_l^{-m}   = (-1)^m (Y_l^m)*, which is the convention of Edmonds,
 * Varshalovich, and the bulk of the atomic-, molecular-, and nuclear-physics
 * literature.  The real spherical harmonics used by gaunt_real and
 * wignernj_real_ylm_in_complex_ylm follow the Wikipedia/Condon–Shortley
 * construction; see the comment block above gaunt_real() below.
 */
#ifndef WIGNERNJ_H
#define WIGNERNJ_H

/* ── Complex-output element type ────────────────────────────────────────────
 *
 * `wignernj_real_ylm_in_complex_ylm*` below fills a complex-valued matrix.  The
 * C99 `_Complex` keyword is recognised by gcc/clang/Apple-Clang/Intel-icx
 * but not by Microsoft Visual C (which provides `_Dcomplex`/`_Fcomplex`/
 * `_Lcomplex` typedef'd structs in <complex.h> instead) and not by C++
 * (which uses `std::complex<T>` from <complex>).  To keep the C function
 * signatures type-safe across every supported toolchain, we typedef the
 * three element types here so each compiler sees its native form:
 *
 *   gcc/clang/Intel C99:  wignernj_cfloat_t == float       _Complex
 *                         wignernj_cdouble_t == double      _Complex
 *                         wignernj_cldouble_t == long double _Complex
 *   MSVC C:               wignernj_cfloat_t == _Fcomplex (struct)
 *                         wignernj_cdouble_t == _Dcomplex (struct)
 *                         wignernj_cldouble_t == _Lcomplex (struct)
 *   any C++:              wignernj_c{float,double,ldouble}_t == layout-compat
 *                         struct{T[2]}.  C++ users should go through the
 *                         std::complex<T> overloads in <wignernj.hpp>; the
 *                         shim type exists only so the C-side declarations
 *                         parse in C++ translation units.
 *
 * All three forms have identical memory layout (real then imaginary,
 * interleaved) per C99 §6.2.5/13 and C++11 [complex.numbers]/4, so the
 * library implementation reads/writes the buffer as a plain `T *` of
 * length 2*(2l+1)^2 internally. */
#ifndef __cplusplus
#include <complex.h>
#if defined(_MSC_VER)
typedef _Fcomplex wignernj_cfloat_t;
typedef _Dcomplex wignernj_cdouble_t;
typedef _Lcomplex wignernj_cldouble_t;
#else
typedef float _Complex wignernj_cfloat_t;
typedef double _Complex wignernj_cdouble_t;
typedef long double _Complex wignernj_cldouble_t;
#endif
#else
typedef struct {
  float _pair[2];
} wignernj_cfloat_t;
typedef struct {
  double _pair[2];
} wignernj_cdouble_t;
typedef struct {
  long double _pair[2];
} wignernj_cldouble_t;
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* ── Wigner 3j symbol ────────────────────────────────────────────────────── */
/* ( j1  j2  j3 )                                                            */
/* ( m1  m2  m3 )                                                            */

float wigner3j_f(int tj1, int tj2, int tj3, int tm1, int tm2, int tm3);
double wigner3j(int tj1, int tj2, int tj3, int tm1, int tm2, int tm3);
long double wigner3j_l(int tj1, int tj2, int tj3, int tm1, int tm2, int tm3);

/* ── Wigner 6j symbol ────────────────────────────────────────────────────── */
/* { j1 j2 j3 }                                                              */
/* { j4 j5 j6 }                                                              */

float wigner6j_f(int tj1, int tj2, int tj3, int tj4, int tj5, int tj6);
double wigner6j(int tj1, int tj2, int tj3, int tj4, int tj5, int tj6);
long double wigner6j_l(int tj1, int tj2, int tj3, int tj4, int tj5, int tj6);

/* ── Wigner 9j symbol ────────────────────────────────────────────────────── */
/* { j11 j12 j13 }  (row-major order)                                        */
/* { j21 j22 j23 }                                                           */
/* { j31 j32 j33 }                                                           */

float wigner9j_f(int tj11, int tj12, int tj13, int tj21, int tj22, int tj23,
                 int tj31, int tj32, int tj33);
double wigner9j(int tj11, int tj12, int tj13, int tj21, int tj22, int tj23,
                int tj31, int tj32, int tj33);
long double wigner9j_l(int tj11, int tj12, int tj13, int tj21, int tj22,
                       int tj23, int tj31, int tj32, int tj33);

/* ── Clebsch-Gordan coefficient ──────────────────────────────────────────── */
/* <j1 m1; j2 m2 | J M>  =  (-1)^(j1-j2+M) sqrt(2J+1) * (j1 j2  J )       */
/*                                                          (m1 m2 -M)       */

float clebsch_gordan_f(int tj1, int tm1, int tj2, int tm2, int tJ, int tM);
double clebsch_gordan(int tj1, int tm1, int tj2, int tm2, int tJ, int tM);
long double clebsch_gordan_l(int tj1, int tm1, int tj2, int tm2, int tJ,
                             int tM);

/* ── Racah W-coefficient ─────────────────────────────────────────────────── */
/* W(j1 j2 J j3; j12 j23)  =  (-1)^(j1+j2+J+j3) * {j1  j2  j12}           */
/*                                                    {j3  J   j23}           */

float racah_w_f(int tj1, int tj2, int tJ, int tj3, int tj12, int tj23);
double racah_w(int tj1, int tj2, int tJ, int tj3, int tj12, int tj23);
long double racah_w_l(int tj1, int tj2, int tJ, int tj3, int tj12, int tj23);

/* ── Fano X-coefficient ──────────────────────────────────────────────────── */
/* X(j1 j2 j12; j3 j4 j34; j13 j24 J)                                        */
/*   = sqrt[(2j12+1)(2j34+1)(2j13+1)(2j24+1)]                                */
/*     * { j1   j2   j12 }                                                   */
/*       { j3   j4   j34 }                                                   */
/*       { j13  j24  J   }                                                   */
/* This is a normalisation variant of the 9j symbol (Fano 1953;             */
/* Edmonds 1957 §6.4) used in the analysis of polarisation correlations.   */

float fano_x_f(int tj1, int tj2, int tj12, int tj3, int tj4, int tj34, int tj13,
               int tj24, int tJ);
double fano_x(int tj1, int tj2, int tj12, int tj3, int tj4, int tj34, int tj13,
              int tj24, int tJ);
long double fano_x_l(int tj1, int tj2, int tj12, int tj3, int tj4, int tj34,
                     int tj13, int tj24, int tJ);

/* ── Gaunt coefficient ───────────────────────────────────────────────────── */
/* G(l1,m1,l2,m2,l3,m3)                                                      */
/*   = integral Y_{l1}^{m1} Y_{l2}^{m2} Y_{l3}^{m3} dΩ                     */
/*   = sqrt[(2l1+1)(2l2+1)(2l3+1)/(4π)]                                     */
/*     * (l1 l2 l3; 0 0 0) * (l1 l2 l3; m1 m2 m3)                          */
/* The closed-form expression in terms of two 3j symbols assumes the         */
/* Condon–Shortley phase for Y_l^m (see file header above).                  */
/* l arguments must be non-negative integers; tl = 2*l (always even).       */

float gaunt_f(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3);
double gaunt(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3);
long double gaunt_l(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3);

/* ── Real-spherical-harmonic Gaunt coefficient ──────────────────────────── */
/* G^R(l1,m1,l2,m2,l3,m3) = integral S_{l1,m1} S_{l2,m2} S_{l3,m3} dΩ        */
/*                                                                           */
/* with the Condon–Shortley / Wikipedia convention for the real spherical    */
/* harmonics S_{l,m}:                                                        */
/*   S_{l, 0}    = Y_l^0                                                     */
/*   S_{l,  m>0} = (1/√2)(Y_l^{-m} + (-1)^m Y_l^m)                           */
/*   S_{l,  m<0} = (i/√2)(Y_l^{ m} - (-1)^|m| Y_l^{-m})                      */
/*                                                                           */
/* ℓ and |m| arguments must be non-negative integers; tl = 2*l, tm = 2*m,    */
/* both always even.  The implementation runs at the cost of a single        */
/* complex-Gaunt evaluation per call: see src/gaunt.c for the algorithm.     */

float gaunt_real_f(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3);
double gaunt_real(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3);
long double gaunt_real_l(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3);

/* ── Real <-> complex spherical-harmonic basis transformation ───────────────
 *
 * Fill the (2l+1) x (2l+1) unitary matrix C such that
 *
 *     S_{l,m_r}  =  sum_{m_c}  C[m_r, m_c]  Y_l^{m_c}
 *
 * with the same real-spherical-harmonic convention used by gaunt_real
 * (S_{l,0} = Y_l^0, S_{l,m>0} = (1/sqrt2)(Y_l^{-m} + (-1)^m Y_l^m),
 *  S_{l,m<0} = (i/sqrt2)(Y_l^{m} - (-1)^|m| Y_l^{-m})).  Both indices
 * use the physics convention m = -l, ..., -1, 0, 1, ..., l, mapped to
 * memory index m + l.
 *
 * Storage convention: COLUMN-MAJOR (Fortran / BLAS / LAPACK), with
 * leading dimension (2l+1).  Element C[m_r, m_c] lives at flat index
 *     col*(2l+1) + row,    col = m_c + l,    row = m_r + l.
 * Each complex slot stores real and imaginary parts in the native
 * compiler representation (`T _Complex` on gcc/clang/Intel, `_Tcomplex`
 * on MSVC, layout-compat struct in C++); all three are interchangeable
 * at the ABI level.
 *
 * The matrix is unitary (C C^H = I).  Each row has at most two
 * non-zero entries (1 at m_r = 0; two of magnitude 1/sqrt(2)
 * otherwise).  Returns silently for l < 0 or NULL output (no
 * diagnostics).
 *
 * Use case: downstream consumers building rotation-equivariant
 * operator matrix elements in the real-Y basis (orbital angular
 * momentum, dipole, momentum, ...) from a complex-basis form.
 * Because C encodes the basis-vector relation S = C Y, the
 * corresponding operator similarity transform is
 *
 *     O_real  =  conj(C)  @  O_complex  @  transpose(C),
 *
 * i.e. two ZGEMM calls; equivalently U^H @ O_complex @ U with
 * U = transpose(C). */

void wignernj_real_ylm_in_complex_ylm_f(int l, wignernj_cfloat_t *C_out);
void wignernj_real_ylm_in_complex_ylm(int l, wignernj_cdouble_t *C_out);
void wignernj_real_ylm_in_complex_ylm_l(int l, wignernj_cldouble_t *C_out);

/* ── Optional warmup and scratch introspection ──────────────────────────────
 *
 * libwignernj caches a per-thread scratch buffer that is lazy-allocated on
 * the first call from a given thread and reused for every subsequent call.
 * After the first call, no further allocations happen unless the angular
 * momenta exceed what the cached buffer was last sized for; in that case
 * the buffer grows in place.  The lazy growth is invisible to the caller.
 *
 * The cache requires thread-local storage and is therefore only active on
 * toolchains that provide one of `__thread` (GCC/Clang/Intel),
 * `__declspec(thread)` (MSVC), or C11 `_Thread_local`.  On a toolchain
 * that exposes none of these, libwignernj falls back to allocating a
 * fresh scratch on every public-API call and freeing it on return.  That
 * fallback is thread-safe by construction (each call owns its own
 * scratch) but pays the historical per-call allocation cost.
 *
 * `wignernj_warmup_to(N_max)` is an optional convenience that pre-
 * grows both per-thread caches---the Racah-pipeline scratch buffers
 * and the factorial-decomposition cache---to fit any symbol evaluation
 * whose worst-case factorial argument is bounded by `N_max`.  After it
 * returns, every subsequent symbol evaluation in this thread within
 * that bound is guaranteed allocation-free.  Pass `N_max = 0` to size
 * the caches to the absolute prime-table ceiling
 * (`wignernj_max_factorial_arg() = 20020` in the default build, with
 * a memory cost of ~6 MB scratch + ~80 MB factorial cache per thread).
 * Useful in benchmarks where the first-call lazy-init cost would
 * otherwise be observed, in hot loops where allocation overhead
 * matters, or in applications that prefer predictable up-front memory
 * usage to lazy growth.
 *
 * The wigner*_max_factorial helpers below return the largest factorial
 * argument that the corresponding symbol would reference for the given
 * inputs.  Use them to size the warmup precisely:
 *
 *   int N = wigner3j_max_factorial(tj1, tj2, tj3, tm1, tm2, tm3);
 *   wignernj_warmup_to(N);
 *
 * To cover an entire workload, take the maximum across calls or pass
 * 0 to populate up to the absolute prime-table ceiling.
 *
 * The function is thread-safe and only touches the calling thread's
 * caches.  Calling it is purely an optimisation; correctness is
 * unchanged whether or not it is called, and it is harmless to call
 * more than once (subsequent calls grow only what is missing).
 *
 * On a toolchain without thread-local storage, `wignernj_warmup_to` is
 * a no-op: there is no persistent scratch to pre-grow, so the call
 * returns immediately and subsequent symbol evaluations continue to
 * allocate per-call. */
void wignernj_warmup_to(int N_max);
int wignernj_max_factorial_arg(void);

/* Release every per-thread cache held by the calling thread (the
 * scratch buffers used by the per-symbol Racah pipelines and the
 * factorial-decomposition cache used by the prime-factorisation
 * helpers).  After this call, the next symbol evaluation in the
 * calling thread re-enters the lazy-init path: scratch buffers are
 * re-allocated on first use, and the factorial cache is rebuilt as
 * each new N is encountered.
 *
 * Useful for applications that want to reclaim memory after a long-
 * running thread is done computing symbols, or before forking a worker.
 * Calling this is purely a memory-management hint; correctness is
 * unaffected.  Each thread can call it independently; threads other
 * than the caller are not affected.  No-op on toolchains without
 * thread-local storage (there is no persistent state to release). */
void wignernj_thread_cleanup(void);

int wigner3j_max_factorial(int tj1, int tj2, int tj3, int tm1, int tm2,
                           int tm3);
int wigner6j_max_factorial(int tj1, int tj2, int tj3, int tj4, int tj5,
                           int tj6);
int wigner9j_max_factorial(int tj11, int tj12, int tj13, int tj21, int tj22,
                           int tj23, int tj31, int tj32, int tj33);
int clebsch_gordan_max_factorial(int tj1, int tm1, int tj2, int tm2, int tJ,
                                 int tM);
int racah_w_max_factorial(int tj1, int tj2, int tJ, int tj3, int tj12,
                          int tj23);
int fano_x_max_factorial(int tj1, int tj2, int tj12, int tj3, int tj4, int tj34,
                         int tj13, int tj24, int tJ);
int gaunt_max_factorial(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3);
int gaunt_real_max_factorial(int tl1, int tm1, int tl2, int tm2, int tl3,
                             int tm3);

#ifdef __cplusplus
}
#endif

#endif /* WIGNERNJ_H */
