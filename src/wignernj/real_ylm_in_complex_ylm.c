/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * Real <-> complex spherical-harmonic basis transformation.
 *
 * Builds the (2l+1) x (2l+1) unitary matrix C such that
 *
 *     S_{l,m_r}  =  sum_{m_c}  C[m_r, m_c]  Y_l^{m_c}
 *
 * with the same real-spherical-harmonic convention used by gaunt_real
 * (see the long comment block above gaunt_real_exact() in src/gaunt.c):
 *
 *     S_{l, 0}     = Y_l^0
 *     S_{l,  m>0}  = (1/sqrt2) ( Y_l^{-m} + (-1)^m Y_l^m )
 *     S_{l,  m<0}  = (i/sqrt2) ( Y_l^{ m} - (-1)^|m| Y_l^{-m} )
 *
 * Each row has at most two non-zero entries; |entry| in { 1, 1/sqrt2 }.
 * Index ordering for both rows (m_real) and columns (m_cplx) is the
 * physics convention  m = -l, -l+1, ..., l-1, l;  row/column index in
 * memory is  m + l.
 *
 * Output is column-major (Fortran / BLAS / LAPACK convention) with
 * leading dimension (2l+1) and interleaved (re, im) complex layout
 * inside each slot.  Element C[m_r, m_c] lives at flat index
 *     col * (2l+1) + row,   col = m_c + l,   row = m_r + l.
 * The C entry points take a typed pointer
 * `wignernj_c{float,double,ldouble,float128}_t *` -- a typedef shim
 * defined in wignernj.h / wignernj_quadmath.h that maps to the
 * native complex element type for each toolchain (`T _Complex` on
 * gcc/clang/Intel, `_Tcomplex` on MSVC C, a layout-compat struct in
 * C++).  All three representations have identical (re, im)
 * interleaved memory layout per C99 §6.2.5/13 and C++11
 * [complex.numbers]/4.  Column-major matches the natural Fortran
 * array layout and the BLAS/LAPACK call convention.
 *
 * The matrix is unitary: C C^H = I.  Because the basis-vector
 * relation is  S = C Y  (and not the coefficient-vector relation),
 * the corresponding operator-matrix similarity transform is
 *
 *     O_real[m_r, m_r']  =  sum_{m_c, m_c'}  bar(C[m_r, m_c])
 *                             * O_complex[m_c, m_c'] * C[m_r', m_c']
 *
 * i.e.  O_real = conj(C) @ O_complex @ transpose(C).  In BLAS terms:
 *
 *     zgemm 'N','T':  T   = O_complex @ transpose(C)
 *     zgemm 'N','N':  O_r = conj(C)   @ T            (apply conj
 *                                                       elementwise)
 *
 * One frequently sees the alternate "U^H @ O @ U" form in textbooks;
 * that requires the matrix U = transpose(C), which is C with its row
 * and column indices swapped (U[m_c, m_r] = C[m_r, m_c]).  The two
 * conventions are equivalent; libwignernj returns C, not U, to match
 * the basis-vector formula stated above (and used by gaunt_real).
 */
#include "wignernj.h"
#include <math.h>
#include <stddef.h>
#include <string.h>

#ifdef WIGNERNJ_HAVE_QUADMATH
#include "wignernj_quadmath.h"
#endif

#ifdef WIGNERNJ_HAVE_MPFR
#include "wignernj_mpfr.h"
#endif

/* Column-major index into the interleaved (re, im) flat array.  N = 2l+1. */
#define IDX_RE(N, r, c) (2u * ((size_t)(c) * (size_t)(N) + (size_t)(r)))
#define IDX_IM(N, r, c) (IDX_RE(N, r, c) + 1u)

/* Per-precision body.  Each instantiation rounds 1/sqrt(2) at native
 * precision; entries with magnitude 1 are exact in every precision. */
#define FILL_RTC_BODY(T, INV_SQRT2)                                            \
  do {                                                                         \
    const int dim = 2 * l + 1;                                                 \
    memset(out, 0, sizeof(T) * 2u * (size_t)dim * (size_t)dim);                \
    for (int m = -l; m <= l; m++) {                                            \
      const int row = m + l;                                                   \
      if (m == 0) {                                                            \
        out[IDX_RE(dim, row, l)] = (T)1;                                       \
      } else if (m > 0) {                                                      \
        const int k = m;                                                       \
        const T sign_k = (k & 1) ? (T)(-1) : (T)(+1);                          \
        out[IDX_RE(dim, row, -k + l)] = (T)(INV_SQRT2);                        \
        out[IDX_RE(dim, row, +k + l)] = sign_k * (T)(INV_SQRT2);               \
      } else {                                                                 \
        const int k = -m;                                                      \
        const T sign_k = (k & 1) ? (T)(-1) : (T)(+1);                          \
        out[IDX_IM(dim, row, -k + l)] = (T)(INV_SQRT2);                        \
        out[IDX_IM(dim, row, +k + l)] = -sign_k * (T)(INV_SQRT2);              \
      }                                                                        \
    }                                                                          \
  } while (0)

void wignernj_real_ylm_in_complex_ylm_f(int l, wignernj_cfloat_t *out_c) {
  /* The wignernj_cfloat_t typedef is layout-compatible with float[2]
   * on every supported toolchain (C99 §6.2.5/13 for `_Complex`,
   * MSVC's `_Fcomplex` is `struct { float _Val[2]; }`, and the C++
   * shim is `struct { float _pair[2]; }`), so the cast to `float *`
   * is well-defined and lets us fill the interleaved (re, im)
   * scalars uniformly. */
  float *out = (float *)out_c;
  if (l < 0 || out_c == NULL)
    return;
  const float inv_sqrt2 = 0.70710678118654752440f;
  FILL_RTC_BODY(float, inv_sqrt2);
}

void wignernj_real_ylm_in_complex_ylm(int l, wignernj_cdouble_t *out_c) {
  double *out = (double *)out_c;
  if (l < 0 || out_c == NULL)
    return;
  const double inv_sqrt2 = 0.70710678118654752440;
  FILL_RTC_BODY(double, inv_sqrt2);
}

void wignernj_real_ylm_in_complex_ylm_l(int l, wignernj_cldouble_t *out_c) {
  long double *out = (long double *)out_c;
  if (l < 0 || out_c == NULL)
    return;
  /* Cover binary64 long double (MSVC/ARM), binary80 (x86-64), and
   * binary128 (aarch64/POWER) with a single correctly-rounded
   * computation rather than a literal. */
  const long double inv_sqrt2 = 1.0L / sqrtl(2.0L);
  FILL_RTC_BODY(long double, inv_sqrt2);
}

#ifdef WIGNERNJ_HAVE_QUADMATH
void wignernj_real_ylm_in_complex_ylm_q(int l, wignernj_cfloat128_t *out_c) {
  __float128 *out = (__float128 *)out_c;
  if (l < 0 || out_c == NULL)
    return;
  const __float128 inv_sqrt2 = 1.0Q / sqrtq(2.0Q);
  FILL_RTC_BODY(__float128, inv_sqrt2);
}
#endif

#ifdef WIGNERNJ_HAVE_MPFR
/* MPFR variant: write the real and imaginary parts into two parallel
 * column-major mpfr_t arrays, each of length (2l+1)^2.  Element
 * C[m_r, m_c] lives at flat index (m_c+l)*(2l+1) + (m_r+l) in both
 * C_re and C_im.  The caller is responsible for mpfr_init2'ing every
 * entry before the call and mpfr_clear'ing them afterwards.  Entries
 * set to zero or +/- 1 are exact at every precision; the 1/sqrt(2)
 * entries are correctly rounded to the per-entry precision under the
 * chosen rounding mode. */
void wignernj_real_ylm_in_complex_ylm_mpfr(int l, mpfr_t *C_re, mpfr_t *C_im,
                                           mpfr_rnd_t rnd) {
  if (l < 0 || C_re == NULL || C_im == NULL)
    return;
  const int dim = 2 * l + 1;
  const size_t n = (size_t)dim * (size_t)dim;

  /* Zero-fill both parts. */
  for (size_t i = 0; i < n; i++) {
    mpfr_set_zero(C_re[i], +1);
    mpfr_set_zero(C_im[i], +1);
  }

  for (int m = -l; m <= l; m++) {
    const int row = m + l;
    if (m == 0) {
      /* C[0, 0]: row=l, col=l */
      const size_t idx = (size_t)l * (size_t)dim + (size_t)l;
      mpfr_set_si(C_re[idx], 1, rnd);
    } else if (m > 0) {
      const int k = m;
      const int sign_k = (k & 1) ? -1 : +1;
      const size_t cn = (size_t)(-k + l) * (size_t)dim + (size_t)row;
      const size_t cp = (size_t)(+k + l) * (size_t)dim + (size_t)row;
      /* re[cn] = 1/sqrt(2);  re[cp] = sign_k / sqrt(2) */
      mpfr_set_ui(C_re[cn], 2, rnd);
      mpfr_rec_sqrt(C_re[cn], C_re[cn], rnd);
      mpfr_set(C_re[cp], C_re[cn], rnd);
      if (sign_k < 0)
        mpfr_neg(C_re[cp], C_re[cp], rnd);
    } else {
      const int k = -m;
      const int sign_k = (k & 1) ? -1 : +1;
      const size_t cn = (size_t)(-k + l) * (size_t)dim + (size_t)row;
      const size_t cp = (size_t)(+k + l) * (size_t)dim + (size_t)row;
      /* im[cn] = 1/sqrt(2);  im[cp] = -sign_k / sqrt(2) */
      mpfr_set_ui(C_im[cn], 2, rnd);
      mpfr_rec_sqrt(C_im[cn], C_im[cn], rnd);
      mpfr_set(C_im[cp], C_im[cn], rnd);
      if (sign_k > 0)
        mpfr_neg(C_im[cp], C_im[cp], rnd);
    }
  }
}
#endif

#undef FILL_RTC_BODY
#undef IDX_RE
#undef IDX_IM
