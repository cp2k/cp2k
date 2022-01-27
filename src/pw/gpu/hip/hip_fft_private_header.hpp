/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#ifndef HIP_FFT_PRIVATE_HEADER_HPP
#define HIP_FFT_PRIVATE_HEADER_HPP
#include <array>
#include <cstdio>
#include <hipblas.h>
#include <hipfft.h>
#include <map>
#include <string>
#include <unistd.h>
#include <vector>

#include "../../../offload/offload_library.h"
#include "../../../offload/offload_operations.h"

typedef hipfftDoubleComplex pw_complex_type;
typedef hipblasHandle_t blasHandle_t;

class fft_plan;

enum fft_direction {
  FFT_FORWARD = HIPFFT_FORWARD,
  FFT_BACKWARD = HIPFFT_BACKWARD,
  FFT_UNKOWN
};

static void fft_error_check(hipfftResult hipfftError, int line,
                            const std::string filename) {
  int pid;

  std::map<hipfftResult, std::string> error_messages_ = {
      {HIPFFT_SUCCESS, "HIPFFT_SUCCESS"},
      {HIPFFT_INVALID_PLAN, "HIPFFT_INVALID_PLAN"},
      {HIPFFT_ALLOC_FAILED, "HIPFFT_ALLOC_FAILED"},
      {HIPFFT_INVALID_TYPE, "HIPFFT_INVALID_TYPE"},
      {HIPFFT_INVALID_VALUE, "HIPFFT_INVALID_VALUE"},
      {HIPFFT_INTERNAL_ERROR, "HIPFFT_INTERNAL_ERROR"},
      {HIPFFT_EXEC_FAILED, "HIPFFT_EXEC_FAILED"},
      {HIPFFT_SETUP_FAILED, "HIPFFT_SETUP_FAILED"},
      {HIPFFT_INVALID_SIZE, "HIPFFT_INVALID_SIZE"},
      {HIPFFT_INCOMPLETE_PARAMETER_LIST, "HIPFFT_INCOMPLETE_PARAMETER_LIST"},
      {HIPFFT_INVALID_DEVICE, "HIPFFT_INVALID_DEVICE"},
      {HIPFFT_PARSE_ERROR, "HIPFFT_PARSE_ERROR"},
      {HIPFFT_NO_WORKSPACE, "HIPFFT_NO_WORKSPACE"},
      {HIPFFT_NOT_IMPLEMENTED, "HIPFFT_NOT_IMPLEMENTED"},
      {HIPFFT_NOT_SUPPORTED, "HIPFFT_NOT_SUPPORTED"}};

  if (hipfftError != HIPFFT_SUCCESS) {
    pid = getpid();
    printf("%d HIP FFT Error line: %s %d \n", pid, filename.c_str(), line);
    printf("%d HIP FFT Error (%s)\n", pid,
           error_messages_.at(hipfftError).c_str());
    fflush(stdout);
    exit(-1);
  }
}

static void blas_error_check(hipblasStatus_t hipblasError, int line,
                             const std::string filename) {
  int pid;

  std::map<hipblasStatus_t, std::string> error_messages_ = {
      {HIPBLAS_STATUS_SUCCESS, "HIPBLAS_STATUS_SUCCESS"},
      {HIPBLAS_STATUS_NOT_INITIALIZED, "HIPBLAS_STATUS_NOT_INITIALIZED"},
      {HIPBLAS_STATUS_ALLOC_FAILED, "HIPBLAS_STATUS_ALLOC_FAILED"},
      {HIPBLAS_STATUS_INVALID_VALUE, "HIPBLAS_STATUS_INVALID_VALUE"},
      {HIPBLAS_STATUS_MAPPING_ERROR, "HIPBLAS_STATUS_MAPPING_ERROR"},
      {HIPBLAS_STATUS_EXECUTION_FAILED, "HIPBLAS_STATUS_EXECUTION_FAILED"},
      {HIPBLAS_STATUS_INTERNAL_ERROR, "HIPBLAS_STATUS_INTERNAL_ERROR"},
      {HIPBLAS_STATUS_NOT_SUPPORTED, "HIPBLAS_STATUS_NOT_SUPPORTED"},
      {HIPBLAS_STATUS_ARCH_MISMATCH, "HIPBLAS_STATUS_ARCH_MISMATCH"},
      {HIPBLAS_STATUS_HANDLE_IS_NULLPTR, "HIPBLAS_STATUS_HANDLE_IS_NULLPTR"},
      {HIPBLAS_STATUS_INVALID_ENUM, "HIPBLAS_STATUS_INVALID_ENUM"},
      {HIPBLAS_STATUS_UNKNOWN, "HIPBLAS_STATUS_UNKNOWN"}};

  if (hipblasError != HIPBLAS_STATUS_SUCCESS) {
    pid = getpid();
    printf("%d HIP BLAS Error file: %s line: %d \n", pid, filename.c_str(),
           line);
    printf("%d HIP BLAS Error (%s)\n", pid,
           error_messages_.at(hipblasError).c_str());
    fflush(stdout);
    exit(-1);
  }
}

static void blasCreate(blasHandle_t *handle__) {
  blas_error_check(hipblasCreate(handle__), __LINE__, __FILE__);
}

static void blasDestroy(blasHandle_t &handle__) {
  blas_error_check(hipblasDestroy(handle__), __LINE__, __FILE__);
}

static void blasSetStream(blasHandle_t &handle__, offloadStream_t &stream__) {
  blas_error_check(hipblasSetStream(handle__, stream__), __LINE__, __FILE__);
}

static void gpuDcopy(blasHandle_t &handle__, int n__, const double *x__,
                     int incx__, double *y__, int incy__) {
  blas_error_check(hipblasDcopy(handle__, n__, x__, incx__, y__, incy__),
                   __LINE__, __FILE__);
}

void gpu_gather(hipStream_t &stream__, const double scale__, int num_points__,
                const int *map_index__, const double *dataIn__,
                double *dataOut__);

void gpu_scatter(hipStream_t &stream__, const double scale__, int num_points__,
                 const int nmaps__, const int *map_index__,
                 const double *dataIn__, double *dataOut__);

#endif
