/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#ifndef CUDA_FFT_PRIVATE_HEADER_HPP
#define CUDA_FFT_PRIVATE_HEADER_HPP

#include <array>
#include <cstdio>
#include <cublas_v2.h>
#include <cuda.h>
#include <cufft.h>
#include <map>
#include <string>
#include <unistd.h>
#include <vector>

#include "../../../offload/offload_library.h"
#include "../../../offload/offload_operations.h"

typedef cufftDoubleComplex pw_complex_type;
typedef cublasHandle_t blasHandle_t;

enum fft_direction {
  FFT_FORWARD = CUFFT_FORWARD,
  FFT_BACKWARD = CUFFT_INVERSE,
  FFT_UNKOWN
};
class fft_plan;

static void fft_error_check(cufftResult cufftError, int line,
                            const std::string filename) {
  int pid;

  std::map<cufftResult, std::string> error_messages_ = {
      {CUFFT_SUCCESS, "CUFFT_SUCCESS"},
      {CUFFT_INVALID_PLAN, "CUFFT_INVALID_PLAN"},
      {CUFFT_ALLOC_FAILED, "CUFFT_ALLOC_FAILED"},
      {CUFFT_INVALID_TYPE, "CUFFT_INVALID_TYPE"},
      {CUFFT_INVALID_VALUE, "CUFFT_INVALID_VALUE"},
      {CUFFT_INTERNAL_ERROR, "CUFFT_INTERNAL_ERROR"},
      {CUFFT_EXEC_FAILED, "CUFFT_EXEC_FAILED"},
      {CUFFT_SETUP_FAILED, "CUFFT_SETUP_FAILED"},
      {CUFFT_INVALID_SIZE, "CUFFT_INVALID_SIZE"},
      {CUFFT_INCOMPLETE_PARAMETER_LIST, "CUFFT_INCOMPLETE_PARAMETER_LIST"},
      {CUFFT_INVALID_DEVICE, "CUFFT_INVALID_DEVICE"},
      {CUFFT_PARSE_ERROR, "CUFFT_PARSE_ERROR"},
      {CUFFT_NO_WORKSPACE, "CUFFT_NO_WORKSPACE"},
      {CUFFT_NOT_IMPLEMENTED, "CUFFT_NOT_IMPLEMENTED"},
      {CUFFT_NOT_SUPPORTED, "CUFFT_NOT_SUPPORTED"}};

  if (cufftError != CUFFT_SUCCESS) {
    pid = getpid();
    printf("%d CUDA FFT Error line: %s %d \n", pid, filename.c_str(), line);
    printf("%d CUDA FFT Error (%s)\n", pid,
           error_messages_.at(cufftError).c_str());
    fflush(stdout);
    abort();
  }
}

static void blas_error_check(cublasStatus_t cublasError, int line,
                             const std::string filename) {
  int pid;

  std::map<cublasStatus_t, std::string> error_messages_ = {
      {CUBLAS_STATUS_SUCCESS, "CUBLAS_STATUS_SUCCESS"},
      {CUBLAS_STATUS_NOT_INITIALIZED, "CUBLAS_STATUS_NOT_INITIALIZED"},
      {CUBLAS_STATUS_ALLOC_FAILED, "CUBLAS_STATUS_ALLOC_FAILED"},
      {CUBLAS_STATUS_INVALID_VALUE, "CUBLAS_STATUS_INVALID_VALUE"},
      {CUBLAS_STATUS_MAPPING_ERROR, "CUBLAS_STATUS_MAPPING_ERROR"},
      {CUBLAS_STATUS_EXECUTION_FAILED, "CUBLAS_STATUS_EXECUTION_FAILED"},
      {CUBLAS_STATUS_INTERNAL_ERROR, "CUBLAS_STATUS_INTERNAL_ERROR"},
      {CUBLAS_STATUS_NOT_SUPPORTED, "CUBLAS_STATUS_NOT_SUPPORTED"},
      {CUBLAS_STATUS_ARCH_MISMATCH, "CUBLAS_STATUS_ARCH_MISMATCH"}};

  if (cublasError != CUBLAS_STATUS_SUCCESS) {
    pid = getpid();
    printf("%d CUDA BLAS Error file: %s line: %d \n", pid, filename.c_str(),
           line);
    printf("%d CUDA BLAS Error (%s)\n", pid,
           error_messages_.at(cublasError).c_str());
    fflush(stdout);
    abort();
  }
}

static void error_check(cudaError_t cudaError, int line,
                        const std::string filename) {
  int pid;
  if (!(cudaError == cudaSuccess || cudaError == cudaErrorNotReady)) {
    pid = getpid();
    printf("%d CUDA RT Error line %s %d\n", pid, filename.c_str(), line);
    printf("%d CUDA RT1 Error: %s\n", pid, cudaGetErrorString(cudaError));
    fflush(stdout);
    abort();
  }
}

static void blasCreate(blasHandle_t *handle__) {
  blas_error_check(cublasCreate(handle__), __LINE__, __FILE__);
}

static void blasDestroy(blasHandle_t &handle__) {
  blas_error_check(cublasDestroy(handle__), __LINE__, __FILE__);
}

static void blasSetStream(blasHandle_t &handle__, offloadStream_t &stream__) {
  blas_error_check(cublasSetStream(handle__, stream__), __LINE__, __FILE__);
}

inline void gpuDcopy(blasHandle_t &handle__, int n__, const double *x__,
                     int incx__, double *y__, int incy__) {
  blas_error_check(cublasDcopy(handle__, n__, x__, incx__, y__, incy__),
                   __LINE__, __FILE__);
}

void gpu_gather(offloadStream_t &stream__, const double scale__,
                int num_points__, const int *map_index__,
                const double *dataIn__, double *dataOut__);

void gpu_scatter(offloadStream_t &stream__, const double scale__,
                 int num_points__, const int nmaps__, const int *map_index__,
                 const double *dataIn__, double *dataOut__);
#endif
