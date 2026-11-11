/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#if defined(__LIBINT)
#include <libint2.h>

/* Alternative to C_F_PROCPOINTER followed by a Fortran call. */
void cp2k_libint_build(void (*build)(const Libint_t *),
                       const Libint_t *libint) {
  build(libint);
}
#endif
