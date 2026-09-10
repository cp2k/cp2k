/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

/* Thin C wrapper over libwignernj.
 *
 * CP2K uses libwignernj either from an installation found by CMake, which
 * exports the plain upstream names, or from the copy bundled in src/wignernj,
 * whose symbols are namespaced to cp2k_wignernj_* so that they cannot collide
 * with the other libraries linked into libcp2k (see
 * src/wignernj/cp2k_wignernj_prefix.h). The two cases would otherwise need
 * different BIND(C, name=) strings in src/common/wignernj_interface.F, which
 * cannot be selected by the preprocessor without splitting the interface body.
 * These wrappers give the Fortran side one set of names that works either way.
 */

#if defined(__LIBWIGNERNJ_BUNDLED)
#include "wignernj/cp2k_wignernj_prefix.h"
#include "wignernj/wignernj.h"
#else
#include <wignernj.h>
#endif

/* Gaunt coefficient of three complex spherical harmonics. All angular momenta
 * are passed as twice their value, as libwignernj expects. */
double cp2k_gaunt(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3) {
  return gaunt(tl1, tm1, tl2, tm2, tl3, tm3);
}

/* Gaunt coefficient of three real spherical harmonics. */
double cp2k_gaunt_real(int tl1, int tm1, int tl2, int tm2, int tl3, int tm3) {
  return gaunt_real(tl1, tm1, tl2, tm2, tl3, tm3);
}
