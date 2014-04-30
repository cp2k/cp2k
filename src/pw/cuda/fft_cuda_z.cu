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

// global dependencies
#include <cuda_runtime.h>
#include <cufft.h>
#include <cublas.h>
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>

// local dependencies
#include "pw_cuda_utils.h"
#include "fft_cuda_internal.h"

// debug flag
#define CHECK 1
#define VERBOSE 0

// dimensions
static int         n_plans = 0;
static cufftHandle saved_plans[max_plans];
static int         iplandims[max_plans][5];

// --- CODE -------------------------------------------------------------------


extern void cufft_error_check (cufftResult_t cufftError, int line) {
  int         pid;
  size_t      free, total;
  cudaError_t cErr2;

  cErr2 = cudaGetLastError();
  if (cufftError != CUFFT_SUCCESS || cErr2 != cudaSuccess) {
    pid = getpid();
    printf("%d CUDA FFT Error line: %d \n", pid, line);
    switch (cufftError) {
      case CUFFT_INVALID_PLAN:   printf("%d CUDA FFT1 Error (CUFFT_INVALID_PLAN)\n", pid); break;
      case CUFFT_ALLOC_FAILED:   printf("%d CUDA FFT1 Error (CUFFT_ALLOC_FAILED)\n", pid); break;
      case CUFFT_INVALID_VALUE:  printf("%d CUDA FFT1 Error (CUFFT_INVALID_VALUE)\n", pid); break;
      case CUFFT_INTERNAL_ERROR: printf("%d CUDA FFT1 Error (CUFFT_INTERNAL_ERROR)\n", pid); break;
      case CUFFT_EXEC_FAILED:    printf("%d CUDA FFT1 Error (CUFFT_EXEC_FAILED)\n", pid); break;
      case CUFFT_INVALID_SIZE:   printf("%d CUDA FFT1 Error (CUFFT_INVALID_SIZE)\n", pid); break;
      default: printf("%d CUDA FFT1 Error (--unimplemented--) %d %d\n", pid, cufftError, cErr2); break;
    }
    printf("%d CUDA FFT2 Error %s \n", pid, cudaGetErrorString(cErr2));
    cudaMemGetInfo(&free,&total);
    printf("%d Free: %zu , Total: %zu\n", pid, free, total);
    fflush(stdout);
    exit(-1);
  }
}

/******************************************************************************
 * \brief   Sets up and save a double precision complex 3D-FFT plan on the GPU.
 *          Saved plans are reused if they fit the requirements.
 * \author  Andreas Gloess
 * \date    2012-05-18
 * \version 0.01
 *****************************************************************************/
void fftcu_plan3d_z(      cufftHandle  &plan,
                          int          &ioverflow,
                    const int          *n,
                    const cudaStream_t  cuda_stream) {

  int i;
  cufftResult_t cErr;

  ioverflow = 0;
  for (i = 0; i < n_plans; i++) {
    if ( iplandims[i][0] == 3 && 
         iplandims[i][1] == n[0] && 
         iplandims[i][2] == n[1] && 
         iplandims[i][3] == n[2] &&
         iplandims[i][4] == 0 ) {
      plan = saved_plans[i];
      return;
    }
  }

  if (VERBOSE) printf("FFT 3D (%d-%d-%d)\n", n[0], n[1], n[2]);
  cErr = cufftPlan3d(&plan, n[2], n[1], n[0], CUFFT_Z2Z);
  if (CHECK) cufft_error_check(cErr, __LINE__);
  cErr = cufftSetStream(plan, cuda_stream);
  if (CHECK) cufft_error_check(cErr, __LINE__);
  cErr = cufftSetCompatibilityMode(plan, FFT_ALIGNMENT);
  if (CHECK) cufft_error_check(cErr, __LINE__);

  if ( n_plans < max_3d_plans ) {
    saved_plans[n_plans] = plan;
    iplandims[n_plans][0] = 3;
    iplandims[n_plans][1] = n[0];
    iplandims[n_plans][2] = n[1];
    iplandims[n_plans][3] = n[2];
    iplandims[n_plans][4] = 0;
    n_plans++;
    return;
  }
  ioverflow=1;
}

/******************************************************************************
 * \brief   Sets up and save a double precision complex 2D-FFT plan on the GPU.
 *          Saved plans are reused if they fit the requirements.
 * \author  Andreas Gloess
 * \date    2012-07-16
 * \version 0.01
 *****************************************************************************/
void fftcu_plan2dm_z(      cufftHandle  &plan,
                           int          &ioverflow,
                     const int          *n,
                     const int           fsign,
                     const cudaStream_t  cuda_stream) {

  int i, istride, idist, ostride, odist, batch;
  int nsize[2], inembed[2], onembed[2];
  cufftResult_t cErr;

  ioverflow = 0;
  for (i = 0; i < n_plans; i++) {
    if ( iplandims[i][0] == 2 &&
         iplandims[i][1] == n[0] &&
         iplandims[i][2] == n[1] &&
         iplandims[i][3] == n[2] &&
         iplandims[i][4] == fsign ) {
      plan = saved_plans[i];
      return;
    }
  }

  nsize[0] = n[2];
  nsize[1] = n[1];
  inembed[0] = n[2];
  inembed[1] = n[1];
  onembed[0] = n[2];
  onembed[1] = n[1];
  batch = n[0];
  if ( fsign == +1 ) {
    istride = n[0];
    idist = 1;
    ostride = 1;
    odist = n[1]*n[2];
  } else {
    istride = 1;
    idist = n[1]*n[2];
    ostride = n[0];
    odist = 1;
  }

  if (VERBOSE) printf("FFT 2D (%d) (%d-%d-%d) %d %d %d %d\n", fsign, n[0], n[1], n[2], istride, idist, ostride, odist);
  cErr = cufftPlanMany(&plan, 2, nsize, inembed, istride, idist, onembed, ostride, odist, CUFFT_Z2Z, batch);
  if (CHECK) cufft_error_check(cErr, __LINE__);
  cErr = cufftSetStream(plan, cuda_stream);
  if (CHECK) cufft_error_check(cErr, __LINE__);
  cErr = cufftSetCompatibilityMode(plan, FFT_ALIGNMENT);
  if (CHECK) cufft_error_check(cErr, __LINE__);

  if ( n_plans < max_2d_plans ) {
    saved_plans[n_plans] = plan;
    iplandims[n_plans][0] = 2;
    iplandims[n_plans][1] = n[0];
    iplandims[n_plans][2] = n[1];
    iplandims[n_plans][3] = n[2];
    iplandims[n_plans][4] = fsign;
    n_plans++;
    return;
  }

  ioverflow=1;
}

/******************************************************************************
 * \brief   Sets up and save a double precision complex 1D-FFT plan on the GPU.
 *          Saved plans are reused if they fit the requirements.
 * \author  Andreas Gloess
 * \date    2012-07-04
 * \version 0.01
 *****************************************************************************/
void fftcu_plan1dm_z(      cufftHandle  &plan,
                           int          &ioverflow,
                     const int           n,
                     const int           m,
                     const int           fsign,
                     const cudaStream_t  cuda_stream) {

  int i, istride, idist, ostride, odist, batch;
  int nsize[1], inembed[1], onembed[1];
  cufftResult_t cErr;

  ioverflow = 0;
  for (i = 0; i<n_plans; i++) {
    if ( iplandims[i][0] == 1 &&
         iplandims[i][1] == n &&
         iplandims[i][2] == m &&
         iplandims[i][3] == 0 &&
         iplandims[i][4] == fsign ) {
      plan = saved_plans[i];
      return;
    }
  }

  nsize[0] = n;
  inembed[0] = 0; // is ignored, but is not allowed to be NULL pointer (for adv. strided I/O)
  onembed[0] = 0; // is ignored, but is not allowed to be NULL pointer (for adv. strided I/O)
  batch = m;
  if ( fsign == +1 ) {
    istride = m;
    idist = 1;
    ostride = 1;
    odist = n;
  } else {
    istride = 1;
    idist = n;
    ostride = m;
    odist = 1;
  }

  if (VERBOSE) printf("FFT 1D (%d) (%d-%d) %d %d %d %d\n", fsign, n, m, istride, idist, ostride, odist);
  cErr = cufftPlanMany(&plan, 1, nsize, inembed, istride, idist, onembed, ostride, odist, CUFFT_Z2Z, batch);
  if (CHECK) cufft_error_check(cErr, __LINE__);
  cErr = cufftSetStream(plan, cuda_stream);
  if (CHECK) cufft_error_check(cErr, __LINE__);
  cErr = cufftSetCompatibilityMode(plan, FFT_ALIGNMENT);
  if (CHECK) cufft_error_check(cErr, __LINE__);

  if ( n_plans < max_1d_plans ) {
    saved_plans[n_plans] = plan;
    iplandims[n_plans][0] = 1;
    iplandims[n_plans][1] = n;
    iplandims[n_plans][2] = m;
    iplandims[n_plans][3] = 0;
    iplandims[n_plans][4] = fsign;
    n_plans++;
    return;
  }

  ioverflow=1;
}

/******************************************************************************
 * \brief   Performs a scaled double precision complex 3D-FFT on the GPU.
 *          Input/output is a DEVICE pointer (data).
 * \author  Andreas Gloess
 * \date    2012-05-18
 * \version 0.01
 *****************************************************************************/
extern "C" void fftcu_run_3d_z_(const int                 fsign,
                                const int                *n,
                                const double              scale,
                                      cufftDoubleComplex *data,
                                const cudaStream_t        cuda_stream) {

  int ioverflow, lmem;
  cufftHandle   plan;
  cufftResult_t cErr;
  cudaError_t  cuErr;

  lmem = n[0] * n[1] * n[2];

  fftcu_plan3d_z(plan, ioverflow, n, cuda_stream);
  if ( fsign < 0  ) {
    cErr = cufftExecZ2Z(plan, data, data, CUFFT_INVERSE);
    if (CHECK) cufft_error_check(cErr, __LINE__);
  }
  else {
    cErr = cufftExecZ2Z(plan, data, data, CUFFT_FORWARD);
    if (CHECK) cufft_error_check(cErr, __LINE__);
  }

  if (scale != 1.0e0) {
    cuErr = cudaStreamSynchronize(cuda_stream);
    if (CHECK) pw_cuda_error_check(cuErr, __LINE__);
    cublasDscal(2*lmem, scale, (double *) data, 1);
  }

  if (ioverflow) {
    cuErr = cudaStreamSynchronize(cuda_stream);
    if (CHECK) pw_cuda_error_check(cuErr, __LINE__);
    cErr = cufftDestroy(plan);
    if (CHECK) cufft_error_check(cErr, __LINE__);
  }
}

/******************************************************************************
 * \brief   Performs a scaled double precision complex 2D-FFT many times on
 *          the GPU.
 *          Input/output are DEVICE pointers (data_in, date_out).
 * \author  Andreas Gloess
 * \date    2012-07-16
 * \version 0.01
 *****************************************************************************/
extern "C" void fftcu_run_2dm_z_(const int                 fsign,
                                 const int                *n,
                                 const double              scale,
                                       cufftDoubleComplex *data_in,
                                       cufftDoubleComplex *data_out,
                                 const cudaStream_t        cuda_stream) {

  int ioverflow, lmem;
  cufftHandle   plan;
  cufftResult_t cErr;
  cudaError_t  cuErr;
  
  lmem = n[0] * n[1] * n[2];

  fftcu_plan2dm_z(plan, ioverflow, n, fsign, cuda_stream);
  if ( fsign < 0 ) {
    cErr = cufftExecZ2Z(plan, data_in, data_out, CUFFT_INVERSE);
    if (CHECK) cufft_error_check(cErr, __LINE__);
  }
  else {
    cErr = cufftExecZ2Z(plan, data_in, data_out, CUFFT_FORWARD);
    if (CHECK) cufft_error_check(cErr, __LINE__);
  }

  if (scale != 1.0e0) {
    cuErr = cudaStreamSynchronize(cuda_stream);
    if (CHECK) pw_cuda_error_check(cuErr, __LINE__);
    cublasDscal(2 * lmem, scale, (double *) data_out, 1);
  }

  if (ioverflow) {
    cuErr = cudaStreamSynchronize(cuda_stream);
    if (CHECK) pw_cuda_error_check(cuErr, __LINE__);
    cErr = cufftDestroy(plan);
    if (CHECK) cufft_error_check(cErr, __LINE__);
  }
}

/******************************************************************************
 * \brief   Performs a scaled double precision complex 1D-FFT many times on
 *          the GPU.
 *          Input/output are DEVICE pointers (data_in, date_out).
 * \author  Andreas Gloess
 * \date    2012-05-18
 * \version 0.01
 *****************************************************************************/
extern "C" void fftcu_run_1dm_z_(const int                 fsign,
                                 const int                 n,
                                 const int                 m,
                                 const double              scale,
                                       cufftDoubleComplex *data_in,
                                       cufftDoubleComplex *data_out,
                                 const cudaStream_t        cuda_stream) {

  int ioverflow, lmem;
  cufftHandle   plan;
  cufftResult_t cErr;
  cudaError_t  cuErr;
  
  lmem = n * m;

  fftcu_plan1dm_z(plan, ioverflow, n, m, fsign, cuda_stream);
  if ( fsign < 0 ) {
    cErr = cufftExecZ2Z(plan, data_in, data_out, CUFFT_INVERSE);
    if (CHECK) cufft_error_check(cErr, __LINE__);
  }
  else {
    cErr = cufftExecZ2Z(plan, data_in, data_out, CUFFT_FORWARD);
    if (CHECK) cufft_error_check(cErr, __LINE__);
  }

  if (scale != 1.0e0) {
    cuErr = cudaStreamSynchronize(cuda_stream);
    if (CHECK) pw_cuda_error_check(cuErr, __LINE__);
    cublasDscal(2 * lmem, scale, (double *) data_out, 1);
  }

  if (ioverflow) {
    cuErr = cudaStreamSynchronize(cuda_stream);
    if (CHECK) pw_cuda_error_check(cuErr, __LINE__);
    cErr = cufftDestroy(plan);
    if (CHECK) cufft_error_check(cErr, __LINE__);
  }
}

/******************************************************************************
 * \brief   Release all stored plans.
 *
 * \author  Andreas Gloess
 * \date    2013-04-23
 * \version 0.01
 *****************************************************************************/
extern "C" void fftcu_release_() {
  int           i;
  cufftHandle   plan;
  cufftResult_t cErr;
  cudaError_t   cuErr;

  for (i = 0; i < n_plans; i++) {
    plan = saved_plans[i];
    cuErr = cudaDeviceSynchronize();
    if (CHECK) pw_cuda_error_check(cuErr, __LINE__);
    cErr = cufftDestroy(plan);
    if (CHECK) cufft_error_check(cErr, __LINE__);
  }
  n_plans = 0;
}

#endif
