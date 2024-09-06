/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2024 CP2K developers group <https://cp2k.org>              */
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
#if (defined(__x86_64__) && 0 != (__x86_64__)) ||                              \
    (defined(__amd64__) && 0 != (__amd64__)) ||                                \
    (defined(_M_X64) || defined(_M_AMD64)) ||                                  \
    (defined(__i386__) && 0 != (__i386__)) || (defined(_M_IX86))
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
#elif defined(__aarch64__) || defined(__arm64__) || defined(_M_ARM64)
  /* only care about SVE, embrace "thermometer approach" like for x86 */
#if (512 <= __ARM_FEATURE_SVE_BITS)
  return CP_MACHINE_ARM_SVE512;
#elif (256 <= __ARM_FEATURE_SVE_BITS)
  return CP_MACHINE_ARM_SVE256;
#elif (128 <= __ARM_FEATURE_SVE_BITS)
  return CP_MACHINE_ARM_SVE128;
#else
  return CP_MACHINE_ARM_ARCH64;
#endif
#elif defined(__ARM_ARCH)
  return CP_MACHINE_CPU_GENERIC;
#else
  return CP_MACHINE_UNKNOWN;
#endif
}

#if defined(__cplusplus)
}
#endif
