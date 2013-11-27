#ifndef FFT_CUDA_H
#define FFT_CUDA_H
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

// Double precision complex procedures
extern "C" void fftcu_run_3d_z_  (const int                 fsign,
                                  const int                *n,
                                  const double              scale,
                                        cufftDoubleComplex *data,
                                  const cudaStream_t        cuda_stream);


extern "C" void fftcu_run_2dm_z_ (const int                 fsign,
                                  const int                *n,
                                  const double              scale,
                                        cufftDoubleComplex *data_in,
                                        cufftDoubleComplex *data_out,
                                  const cudaStream_t        cuda_stream);


extern "C" void fftcu_run_1dm_z_ (const int                 fsign,
                                  const int                 n,
                                  const int                 m,
                                  const double              scale,
                                        cufftDoubleComplex *data_in,
                                        cufftDoubleComplex *data_out,
                                  const cudaStream_t        cuda_stream);


extern "C" void fftcu_release_   ();
#endif
#endif
