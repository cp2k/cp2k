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

class fft_plan {
private:
  cufftHandle plan_;
  int size_[3] = {0, 0, 0};
  int dim_{3};
  int batch_size_ = -1;
  // forward or backward / reverse
  enum fft_direction direction_ { FFT_UNKOWN };
  cudaStream_t stream_;
  bool is_initialized_{false};

public:
  fft_plan() {}

  fft_plan(const std::vector<int> &fft_size__, const int dim__,
           const int batch_size__, const enum fft_direction direction__)
      : dim_(dim__), batch_size_(batch_size__) {
    int inembed[2] = {0, 0};
    int onembed[2] = {0, 0};
    int istride, idist, ostride, odist;

    direction_ = direction__;

    switch (dim__) {
    case 3: {
      size_[0] = fft_size__[0];
      size_[1] = fft_size__[1];
      size_[2] = fft_size__[2];
    } break;
    case 2: {
      size_[0] = fft_size__[1];
      size_[1] = fft_size__[0];
      size_[2] = batch_size__;
      inembed[0] = fft_size__[1];
      inembed[1] = fft_size__[0];
      onembed[0] = fft_size__[1];
      onembed[1] = fft_size__[0];
      batch_size_ = batch_size__;

      if (direction_ == CUFFT_FORWARD) {
        istride = batch_size__;
        idist = 1;
        ostride = 1;
        odist = size_[0] * size_[1];
      } else {
        istride = 1;
        idist = size_[0] * size_[1];
        ostride = batch_size__;
        odist = 1;
      }
    } break;
    case 1: {
      size_[0] = fft_size__[0];
      size_[1] = 1;
      size_[2] = batch_size__;
      batch_size_ = batch_size__;
      if (direction_ == CUFFT_FORWARD) {
        istride = batch_size__;
        idist = 1;
        ostride = 1;
        odist = fft_size__[0];
      } else {
        istride = 1;
        idist = fft_size__[0];
        ostride = batch_size__;
        odist = 1;
      }
      break;
    }
    default:
      abort();
      break;
    }
    if (dim_ == 3) {
      fft_error_check(
          cufftPlan3d(&plan_, size_[2], size_[1], size_[0], CUFFT_Z2Z),
          __LINE__, __FILE__);
    } else {
      fft_error_check(cufftPlanMany(&plan_, dim_, &size_[0], inembed, istride,
                                    idist, onembed, ostride, odist, CUFFT_Z2Z,
                                    batch_size_),
                      __LINE__, __FILE__);
    }

    is_initialized_ = true;
  }

  void set_stream(const cudaStream_t &cuda_stream) {
    stream_ = cuda_stream;
    fft_error_check(cufftSetStream(plan_, stream_), __LINE__, __FILE__);
  }

  ~fft_plan() {
    if (is_initialized_)
      destroy();
  }

  void destroy() {
    error_check(cudaStreamSynchronize(stream_), __LINE__, __FILE__);
    fft_error_check(cufftDestroy(plan_), __LINE__, __FILE__);
    is_initialized_ = false;
  }

  /// run the fft on the data inplace
  void execute_fft(const enum fft_direction direction__,
                   cufftDoubleComplex *data__) {
    fft_error_check(cufftExecZ2Z(plan_, data__, data__, direction__), __LINE__,
                    __FILE__);
  }

  /// run the fft on the data out of place
  void execute_fft(const enum fft_direction direction__,
                   cufftDoubleComplex *dataIn__,
                   cufftDoubleComplex *dataOut__) {
    cufftResult_t cErr;
    // set the stream

    fft_error_check(cufftExecZ2Z(plan_, dataIn__, dataOut__, direction__),
                    __LINE__, __FILE__);
  }

  /// check if this plane can be used to execute the fft
  bool is_it_valid(std::vector<int> size__, const int dim__, const int batch__,
                   const enum fft_direction direction__) const {
    if (dim_ != dim__)
      return false;
    if (batch__ != batch_size_)
      return false;
    switch (dim__) {
    case 3:
      return ((size_[0] != size__[0]) || (size_[1] != size__[1]) ||
              (size_[2] == size__[2]));
      break;
    case 2:
      return ((size_[0] != size__[0]) || (size_[1] != size__[1]));
      break;
    case 1:
      return (size_[0] != size__[0]);
      break;
    default:
      return false;
      break;
    }

    // check for the direction

    if ((direction_ != direction__) && (dim_ != 3)) {
      return false;
    } else {
      return true;
    }
  }

  bool is_initialized() const { return is_initialized_; }
};

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
