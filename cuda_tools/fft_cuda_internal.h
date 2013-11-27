#ifndef FFT_CUDA_INTERNAL_H
#define FFT_CUDA_INTERNAL_H
/******************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations
 *  Copyright (C) 2000 - 2013  CP2K developers group
 *
 *  Authors: Benjamin G Levine, Andreas Gloess
 *
 *  2012/05/18                 Refacturing - original files:
 *                              - cuda_tools/cufft.h
 *                              - cuda_tools/cufft.cu
 *
 *****************************************************************************/
#if defined ( __PW_CUDA )

/******************************************************************************
 * \brief   Sets up static data for FFT plan storage and reuse.
 * \author  Andreas Gloess
 * \date    2012-05-18
 * \version 0.01
 *****************************************************************************/
static const int   max_3d_plans = 30;
static const int   max_2d_plans = 0;
static const int   max_1d_plans = 30;


static const int   sum_plans = max_3d_plans + max_2d_plans + max_1d_plans; 
static const int   max_plans = sum_plans > 1? sum_plans : 1; 

// configuration(s)
#define FFT_ALIGNMENT CUFFT_COMPATIBILITY_NATIVE // potentially faster
//#define FFT_ALIGNMENT CUFFT_COMPATIBILITY_FFTW_PADDING // the default

#endif
#endif
