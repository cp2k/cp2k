/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include "cp2k_wignernj_prefix.h"

/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * Racah W-coefficient via the Wigner 6j symbol.
 *
 * Definition:
 *   W(j1 j2 J j3; j12 j23) = (-1)^(j1+j2+J+j3) * {j1  j2  j12}
 *                                                    {j3  J   j23}
 *
 * Arguments are passed as twice their value (tj = 2*j).
 */
#include "wignernj.h"

/* Racah W phase: (-1)^(j1+j2+J+j3) = (-1)^((tj1+tj2+tJ+tj3)/2). */
static inline int racah_w_phase(int tj1, int tj2, int tJ, int tj3) {
  return ((((tj1 + tj2 + tJ + tj3) / 2) & 1) == 0) ? 1 : -1;
}

int racah_w_max_factorial(int tj1, int tj2, int tJ, int tj3, int tj12,
                          int tj23) {
  /* W = phase * 6j{j1, j2, j12; j3, J, j23}. */
  return wigner6j_max_factorial(tj1, tj2, tj12, tj3, tJ, tj23);
}

float racah_w_f(int tj1, int tj2, int tJ, int tj3, int tj12, int tj23) {
  return (float)racah_w_phase(tj1, tj2, tJ, tj3) *
         wigner6j_f(tj1, tj2, tj12, tj3, tJ, tj23);
}

double racah_w(int tj1, int tj2, int tJ, int tj3, int tj12, int tj23) {
  return (double)racah_w_phase(tj1, tj2, tJ, tj3) *
         wigner6j(tj1, tj2, tj12, tj3, tJ, tj23);
}

long double racah_w_l(int tj1, int tj2, int tJ, int tj3, int tj12, int tj23) {
  return (long double)racah_w_phase(tj1, tj2, tJ, tj3) *
         wigner6j_l(tj1, tj2, tj12, tj3, tJ, tj23);
}

#ifdef WIGNERNJ_HAVE_QUADMATH
#include "wignernj_quadmath.h"
__float128 racah_w_q(int tj1, int tj2, int tJ, int tj3, int tj12, int tj23) {
  return (__float128)racah_w_phase(tj1, tj2, tJ, tj3) *
         wigner6j_q(tj1, tj2, tj12, tj3, tJ, tj23);
}
#endif

#ifdef WIGNERNJ_HAVE_MPFR
#include "wignernj_mpfr.h"
void racah_w_mpfr(mpfr_t rop, int tj1, int tj2, int tJ, int tj3, int tj12,
                  int tj23, mpfr_rnd_t rnd) {
  wigner6j_mpfr(rop, tj1, tj2, tj12, tj3, tJ, tj23, rnd);
  if (racah_w_phase(tj1, tj2, tJ, tj3) < 0)
    mpfr_neg(rop, rop, MPFR_RNDN);
}
#endif
