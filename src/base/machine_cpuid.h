/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2024 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#define CP_MACHINE_CPU_GENERIC 0

#define CP_MACHINE_X86_SSE4 1000
#define CP_MACHINE_X86_AVX 1001
#define CP_MACHINE_X86_AVX2 1002
#define CP_MACHINE_X86_AVX512 1003

#define CP_MACHINE_ARM_ARCH64 2000
#define CP_MACHINE_ARM_SVE128 2100
#define CP_MACHINE_ARM_SVE256 2200
#define CP_MACHINE_ARM_SVE512 2300

#define CP_MACHINE_UNKNOWN 3000
