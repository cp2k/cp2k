/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2021 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#ifndef FFT_CUDA_UTILS_H
#define FFT_CUDA_UTILS_H
#include <cufft.h>

extern void cufft_error_check(cufftResult_t cufftError, int line);

#endif
