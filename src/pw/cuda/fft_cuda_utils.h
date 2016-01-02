/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2016  CP2K developers group                         *
 *****************************************************************************/

#ifndef FFT_CUDA_UTILS_H
#define FFT_CUDA_UTILS_H
#include <cufft.h>

extern void cufft_error_check (cufftResult_t cufftError, int line);

#endif
