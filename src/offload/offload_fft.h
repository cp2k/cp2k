/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#ifndef OFFLOAD_FFT_H
#define OFFLOAD_FFT_H

#include <stdio.h>
#include <stdlib.h>

#include "offload_library.h"
#include "offload_operations.h"

#if (defined(__OFFLOAD_CUDA) || defined(__OFFLOAD_HIP))

#if defined(__OFFLOAD_CUDA)
#include <cufft.h>
#elif defined(__OFFLOAD_HIP)
#include <hipfft.h>
#endif

#if defined(__OFFLOAD_CUDA)
typedef cufftHandle offload_fftHandle;
typedef cufftType offload_fftType;
typedef cufftResult offload_fftResult;
#elif defined(__OFFLOAD_HIP)
typedef hipfftHandle offload_fftHandle;
typedef hipfftType offload_fftType;
typedef hipfftResult offload_fftResult;
#endif

#if defined(__OFFLOAD_CUDA)
#define OFFLOAD_FFT_SUCCESS CUFFT_SUCCESS
#define OFFLOAD_FFT_FORWARD CUFFT_FORWARD
#define OFFLOAD_FFT_INVERSE CUFFT_INVERSE
#define OFFLOAD_FFT_Z2Z CUFFT_Z2Z
#elif defined(__OFFLOAD_HIP)
#define OFFLOAD_FFT_SUCCESS HIPFFT_SUCCESS
#define OFFLOAD_FFT_FORWARD HIPFFT_FORWARD
#define OFFLOAD_FFT_INVERSE HIPFFT_BACKWARD // inconsistent
#define OFFLOAD_FFT_Z2Z HIPFFT_Z2Z
#endif

/*******************************************************************************
 * \brief Check given cufft status and upon failure abort with a nice message.
 * \author Ole Schuett
 ******************************************************************************/
#define OFFLOAD_FFT_CHECK(status)                                              \
  if (status != OFFLOAD_FFT_SUCCESS) {                                         \
    fprintf(stderr, "ERROR: %s %s %d\n", offload_fftGetErrorString(status),    \
            __FILE__, __LINE__);                                               \
    abort();                                                                   \
  }

/*******************************************************************************
 * \brief Equivalent to cudaGetErrorString for cufft.
 ******************************************************************************/
static inline const char *offload_fftGetErrorString(offload_fftResult error) {
#if defined(__OFFLOAD_CUDA)
  switch (error) {
  case CUFFT_SUCCESS:
    return "CUFFT_SUCCESS";

  case CUFFT_INVALID_PLAN:
    return "CUFFT_INVALID_PLAN";

  case CUFFT_ALLOC_FAILED:
    return "CUFFT_ALLOC_FAILED";

  case CUFFT_INVALID_TYPE:
    return "CUFFT_INVALID_TYPE";

  case CUFFT_INVALID_VALUE:
    return "CUFFT_INVALID_VALUE";

  case CUFFT_INTERNAL_ERROR:
    return "CUFFT_INTERNAL_ERROR";

  case CUFFT_EXEC_FAILED:
    return "CUFFT_EXEC_FAILED";

  case CUFFT_SETUP_FAILED:
    return "CUFFT_SETUP_FAILED";

  case CUFFT_INVALID_SIZE:
    return "CUFFT_INVALID_SIZE";

  case CUFFT_INCOMPLETE_PARAMETER_LIST:
    return "CUFFT_INCOMPLETE_PARAMETER_LIST";

  case CUFFT_INVALID_DEVICE:
    return "CUFFT_INVALID_DEVICE";

  case CUFFT_PARSE_ERROR:
    return "CUFFT_PARSE_ERROR";

  case CUFFT_NO_WORKSPACE:
    return "CUFFT_NO_WORKSPACE";

  case CUFFT_NOT_IMPLEMENTED:
    return "CUFFT_NOT_IMPLEMENTED";

  case CUFFT_NOT_SUPPORTED:
    return "CUFFT_NOT_SUPPORTED";

  case CUFFT_UNALIGNED_DATA:
    return "CUFFT_UNALIGNED_DATA";

  case CUFFT_LICENSE_ERROR:
    return "CUFFT_LICENSE_ERROR";
  }

#elif defined(__OFFLOAD_HIP)

  switch (error) {
  case HIPFFT_SUCCESS:
    return "HIPFFT_SUCCESS";

  case HIPFFT_INVALID_PLAN:
    return "HIPFFT_INVALID_PLAN";

  case HIPFFT_ALLOC_FAILED:
    return "HIPFFT_ALLOC_FAILED";

  case HIPFFT_INVALID_TYPE:
    return "HIPFFT_INVALID_TYPE";

  case HIPFFT_INVALID_VALUE:
    return "HIPFFT_INVALID_VALUE";

  case HIPFFT_INTERNAL_ERROR:
    return "HIPFFT_INTERNAL_ERROR";

  case HIPFFT_EXEC_FAILED:
    return "HIPFFT_EXEC_FAILED";

  case HIPFFT_SETUP_FAILED:
    return "HIPFFT_SETUP_FAILED";

  case HIPFFT_INVALID_SIZE:
    return "HIPFFT_INVALID_SIZE";

  case HIPFFT_INCOMPLETE_PARAMETER_LIST:
    return "HIPFFT_INCOMPLETE_PARAMETER_LIST";

  case HIPFFT_INVALID_DEVICE:
    return "HIPFFT_INVALID_DEVICE";

  case HIPFFT_PARSE_ERROR:
    return "HIPFFT_PARSE_ERROR";

  case HIPFFT_NO_WORKSPACE:
    return "HIPFFT_NO_WORKSPACE";

  case HIPFFT_NOT_IMPLEMENTED:
    return "HIPFFT_NOT_IMPLEMENTED";

  case HIPFFT_NOT_SUPPORTED:
    return "HIPFFT_NOT_SUPPORTED";

  case HIPFFT_UNALIGNED_DATA:
    return "HIPFFT_UNALIGNED_DATA";

    // case HIPFFT_LICENSE_ERROR:
    //  return "HIPFFT_LICENSE_ERROR";
  }
#endif

  return "<unknown>";
}

/*******************************************************************************
 * \brief Wrapper around cufftPlan3d.
 ******************************************************************************/
static inline void offload_fftPlan3d(offload_fftHandle *plan, int nx, int ny,
                                     int nz, offload_fftType type) {
#if defined(__OFFLOAD_CUDA)
  OFFLOAD_FFT_CHECK(cufftPlan3d(plan, nx, ny, nz, type));
#elif defined(__OFFLOAD_HIP)
  OFFLOAD_FFT_CHECK(hipfftPlan3d(plan, nx, ny, nz, type));
#endif
}

/*******************************************************************************
 * \brief Wrapper around cufftPlanMany.
 ******************************************************************************/
static inline void offload_fftPlanMany(offload_fftHandle *plan, int rank,
                                       int *n, int *inembed, int istride,
                                       int idist, int *onembed, int ostride,
                                       int odist, offload_fftType type,
                                       int batch) {
#if defined(__OFFLOAD_CUDA)
  OFFLOAD_FFT_CHECK(cufftPlanMany(plan, rank, n, inembed, istride, idist,
                                  onembed, ostride, odist, type, batch));
#elif defined(__OFFLOAD_HIP)
  OFFLOAD_FFT_CHECK(hipfftPlanMany(plan, rank, n, inembed, istride, idist,
                                   onembed, ostride, odist, type, batch));
#endif
}

/*******************************************************************************
 * \brief Wrapper around cufftSetStream.
 ******************************************************************************/
static inline void offload_fftSetStream(offload_fftHandle plan,
                                        offloadStream_t stream) {
#if defined(__OFFLOAD_CUDA)
  OFFLOAD_FFT_CHECK(cufftSetStream(plan, stream))
#elif defined(__OFFLOAD_HIP)
  OFFLOAD_FFT_CHECK(hipfftSetStream(plan, stream))
#endif
}

/*******************************************************************************
 * \brief Wrapper around cufftDestroy.
 ******************************************************************************/
static inline void offload_fftDestroy(offload_fftHandle plan) {
#if defined(__OFFLOAD_CUDA)
  OFFLOAD_FFT_CHECK(cufftDestroy(plan));
#elif defined(__OFFLOAD_HIP)
  OFFLOAD_FFT_CHECK(hipfftDestroy(plan));
#endif
}

/*******************************************************************************
 * \brief Wrapper around cufftExecZ2Z.
 ******************************************************************************/
static inline void offload_fftExecZ2Z(offload_fftHandle plan,
                                      const double *idata, double *odata,
                                      offload_fftType type) {
#if defined(__OFFLOAD_CUDA)
  OFFLOAD_FFT_CHECK(cufftExecZ2Z(plan, (cufftDoubleComplex *)idata,
                                 (cufftDoubleComplex *)odata, type));
#elif defined(__OFFLOAD_HIP)
  OFFLOAD_FFT_CHECK(hipfftExecZ2Z(plan, (hipfftDoubleComplex *)idata,
                                  (hipfftDoubleComplex *)odata, type));
#endif
}

#endif // #if (defined(__OFFLOAD_CUDA) || defined(__OFFLOAD_HIP))

#endif

// EOF
