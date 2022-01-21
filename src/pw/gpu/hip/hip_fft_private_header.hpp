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

class fft_plan {
private:
  hipfftHandle plan_;
  int size_[3] = {0, 0, 0};
  int sign_{0};
  int dim_{3};
  enum fft_direction direction_ { FFT_UNKOWN };
  int batch_size_ = -1;
  offloadStream_t stream_;
  bool is_initialized_{false};

public:
  fft_plan() {}

  fft_plan(const std::vector<int> &fft_size__, const int dim__,
           const int batch_size__, const enum fft_direction direction__)
      : dim_(dim__), batch_size_(batch_size__) {
    int inembed[2] = {-1, -1};
    int onembed[2] = {-1, -1};
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

      if (direction_ == FFT_FORWARD) {
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
      if (direction_ == FFT_FORWARD) {
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
          hipfftPlan3d(&plan_, size_[2], size_[1], size_[0], HIPFFT_Z2Z),
          __LINE__, __FILE__);
    } else {
      fft_error_check(hipfftPlanMany(&plan_, dim_, &size_[0], inembed, istride,
                                     idist, onembed, ostride, odist, HIPFFT_Z2Z,
                                     batch_size_),
                      __LINE__, __FILE__);
    }

    fft_error_check(hipfftSetAutoAllocation(plan_, 1), __LINE__, __FILE__);
    is_initialized_ = true;
  }

  void set_stream(const offloadStream_t &hip_stream) {
    stream_ = hip_stream;
    fft_error_check(hipfftSetStream(plan_, stream_), __LINE__, __FILE__);
  }

  ~fft_plan() {
    if (is_initialized_)
      destroy();
  }

  void destroy() {
    offloadStreamSynchronize(stream_);
    fft_error_check(hipfftDestroy(plan_), __LINE__, __FILE__);
    is_initialized_ = false;
  }

  /// run the fft on the data inplace
  void execute_fft(const enum fft_direction direction__,
                   hipfftDoubleComplex *data__) {
    fft_error_check(hipfftExecZ2Z(plan_, data__, data__, direction__), __LINE__,
                    __FILE__);
  }

  /// run the fft on the data out of place
  void execute_fft(const enum fft_direction direction__,
                   hipfftDoubleComplex *dataIn__,
                   hipfftDoubleComplex *dataOut__) {
    fft_error_check(hipfftExecZ2Z(plan_, dataIn__, dataOut__, direction__),
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

    if ((direction__ != direction_) && (dim_ != 3)) {
      return false;
    } else {
      return true;
    }
  }

  bool is_initialized() const { return is_initialized_; }
};

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
