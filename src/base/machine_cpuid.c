/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

/* shared between C and Fortran */
#include "machine_cpuid.h"

#if defined(__cplusplus)
extern "C" {
#endif

/*******************************************************************************
 * \brief This routine determines the CPUID according to the given compiler
 *        flags (expected to be similar to Fortran). Similar to other Fortran
 *        compilers, "gfortran -E -dM -mavx - < /dev/null | grep AVX" defines a
 *        variety of predefined macros (also similar to C). However, with a
 *        Fortran translation unit only a subset of these definitions disappears
 *        ("gfortran -E -dM -mavx my.F | grep AVX")
 *        hence an implementation in C is used.
 ******************************************************************************/
int m_cpuid_static(void); /* avoid pedantic warning about missing prototype */
int m_cpuid_static(void) {
#if (__AVX512F__ && __AVX512CD__ && __AVX2__ && __FMA__ && __AVX__ &&          \
     __SSE4_2__ && __SSE4_1__ && __SSE3__)
  return CP_MACHINE_X86_AVX512;
#elif (__AVX2__ && __FMA__ && __AVX__ && __SSE4_2__ && __SSE4_1__ && __SSE3__)
  return CP_MACHINE_X86_AVX2;
#elif (__AVX__ && __SSE4_2__ && __SSE4_1__ && __SSE3__)
  return CP_MACHINE_X86_AVX;
#elif (__SSE4_2__ && __SSE4_1__ && __SSE3__)
  return CP_MACHINE_X86_SSE4;
#else
  return CP_MACHINE_CPU_GENERIC;
#endif
}

#if defined(__cplusplus)
}
#endif
