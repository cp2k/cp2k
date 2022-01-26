/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#ifndef PW_GPU_INTERNAL_HPP
#define PW_GPU_INTERNAL_HPP

#include "../../offload/offload_library.h"

class fft_plan {
private:
#if defined(__OFFLOAD_HIP)
  hipfftHandle plan_;
#elif defined(__OFFLOAD_CUDA)
  cufftHandle plan_;
#endif
  int size_[3] = {0, 0, 0};
  int sign_{0};
  int dim_{3};
  enum fft_direction direction_ { FFT_UNKOWN };
  int batch_size_ = -1;
  offloadStream_t stream_{nullptr};
  bool is_initialized_{false};
  int num_points_{0};
  int *gmap_{nullptr};
  size_t gmap_size_{0};
  bool should_destroy_{false};
  double *ptr_1_{nullptr};
  double *ptr_2_{nullptr};
  size_t gmap_allocated_size_{0};

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
      num_points_ = size_[0] * size_[1] * size_[2];
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

      num_points_ = batch_size_ * size_[0] * size_[1];
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
      num_points_ = batch_size_ * size_[0];

      break;
    }
    default:
      abort();
      break;
    }

    if (dim_ == 3) {
#if defined(__OFFLOAD_HIP)
      fft_error_check(
          hipfftPlan3d(&plan_, size_[2], size_[1], size_[0], HIPFFT_Z2Z),
          __LINE__, __FILE__);
#elif defined(__OFFLOAD_CUDA)
      fft_error_check(
          cufftPlan3d(&plan_, size_[2], size_[1], size_[0], CUFFT_Z2Z),
          __LINE__, __FILE__);
#endif
    } else {
#if defined(__OFFLOAD_HIP)
      fft_error_check(hipfftPlanMany(&plan_, dim_, &size_[0], inembed, istride,
                                     idist, onembed, ostride, odist, HIPFFT_Z2Z,
                                     batch_size_),
                      __LINE__, __FILE__);
#elif defined(__OFFLOAD_CUDA)
      fft_error_check(cufftPlanMany(&plan_, dim_, &size_[0], inembed, istride,
                                    idist, onembed, ostride, odist, CUFFT_Z2Z,
                                    batch_size_),
                      __LINE__, __FILE__);
#endif
    }

#ifdef __OFFLOAD_HIP
    fft_error_check(hipfftSetAutoAllocation(plan_, 1), __LINE__, __FILE__);
#endif
    is_initialized_ = true;

    // allocate memory for the fft buffers
    offloadMalloc((void **)&ptr_1_, sizeof(double) * 2 * num_points_);
    offloadMalloc((void **)&ptr_2_, sizeof(double) * 2 * num_points_);
  }

  // fft_plan(fft_plan &&src__) {
  //   plan_ = src__.plan_;
  //   size_[0] = src__.size_[0];
  //   size_[1] = src__.size_[1];
  //   size_[2] = src__.size_[2];
  //   sign_ = src__.sign_;
  //   dim_ = src__.dim_;
  //   direction_ = src__.direction_;
  //   batch_size_ = src__.batch_size_;
  //   stream_ = src__.stream_;
  //   is_initialized_ = src__.is_initialized_;
  //   num_points_ = src__.num_points_;
  //   gmap_ = src__.gmap_;
  //   gmap_size_ = src__.gmap_size_;
  //   should_destroy_ = src__.should_destroy_;
  //   ptr_1_ = src__.ptr_1_;
  //   ptr_2_  = src__.ptr_2_;
  //   gmap_allocated_size_ = src__.gmap_allocated_size_;

  //   // prevents from releasing memory for src__
  //   src__.is_initialized_ = false;
  // };

  // fft_plan(const fft_plan &src__) {
  //   plan_ = src__.plan_;
  //   size_[0] = src__.size_[0];
  //   size_[1] = src__.size_[1];
  //   size_[2] = src__.size_[2];
  //   sign_ = src__.sign_;
  //   dim_ = src__.dim_;
  //   direction_ = src__.direction_;
  //   batch_size_ = src__.batch_size_;
  //   stream_ = src__.stream_;
  //   is_initialized_ = src__.is_initialized_;
  //   num_points_ = src__.num_points_;
  //   gmap_ = src__.gmap_;
  //   gmap_size_ = src__.gmap_size_;
  //   should_destroy_ = src__.should_destroy_;
  //   ptr_1_ = src__.ptr_1_;
  //   ptr_2_  = src__.ptr_2_;
  //   gmap_allocated_size_ = src__.gmap_allocated_size_;

  //   // prevents from releasing memory for src__
  // };

  void allocate_gmap(int num_elems__) {
    if ((gmap_allocated_size_ < num_elems__) && (gmap_ != nullptr)) {
      offloadFree(gmap_);
    }

    if (gmap_ == nullptr) {
      offloadMalloc((void **)&gmap_, sizeof(int) * num_elems__);
      gmap_allocated_size_ = num_elems__;
    }

    gmap_size_ = num_elems__;
  }

  void set_stream(const offloadStream_t &hip_stream) {
    stream_ = hip_stream;
#if defined(__OFFLOAD_HIP)
    fft_error_check(hipfftSetStream(plan_, stream_), __LINE__, __FILE__);
#elif defined(__OFFLOAD_CUDA)
    fft_error_check(cufftSetStream(plan_, stream_), __LINE__, __FILE__);
#endif
  }

  ~fft_plan() {
    if (should_destroy_)
      destroy();
  }

  void should_destroy(const bool val__) { should_destroy_ = val__; }

  void should_destroy() {
    if (should_destroy_) {
      destroy();
    }
  }

  void destroy() {
#if defined(__OFFLOAD_HIP)
    fft_error_check(hipfftDestroy(plan_), __LINE__, __FILE__);
#else
    fft_error_check(cufftDestroy(plan_), __LINE__, __FILE__);
#endif
    offloadFree(ptr_1_);
    offloadFree(ptr_2_);
    if (gmap_ != nullptr)
      offloadFree(gmap_);
    is_initialized_ = false;
  }

  double *ptr_1() { return ptr_1_; }

  double *ptr_2() { return ptr_2_; }

  int *ghatmap() { return gmap_; }

  /// run the fft on the data inplace
  void execute_fft(const enum fft_direction direction__,
                   pw_complex_type *data__) {
#if defined(__OFFLOAD_HIP)
    fft_error_check(hipfftExecZ2Z(plan_, data__, data__, direction__), __LINE__,
                    __FILE__);
#elif defined(__OFFLOAD_CUDA)
    fft_error_check(cufftExecZ2Z(plan_, data__, data__, direction__), __LINE__,
                    __FILE__);
#endif
  }

  /// run the fft on the data out of place
  void execute_fft(const enum fft_direction direction__,
                   pw_complex_type *dataIn__, pw_complex_type *dataOut__) {
#if defined(__OFFLOAD_HIP)
    fft_error_check(hipfftExecZ2Z(plan_, dataIn__, dataOut__, direction__),
                    __LINE__, __FILE__);
#elif defined(__OFFLOAD_CUDA)
    fft_error_check(cufftExecZ2Z(plan_, dataIn__, dataOut__, direction__),
                    __LINE__, __FILE__);
#endif
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
      return ((size_[0] != size__[0]) || (size_[1] != size__[1]) ||
              (batch__ != batch_size_));
      break;
    case 1:
      return ((size_[0] != size__[0]) || (batch__ != batch_size_));
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

void real_to_complex(offloadStream_t &stream__, const int length__,
                     const double *src__, double *const dst__);
void complex_to_real(offloadStream_t &stream__, const int length__,
                     const double *src__, double *const dst__);
#endif
