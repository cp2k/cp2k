/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#if defined(__PW_GPU)

#include <array>
#include <cstdio>
#include <map>
#include <string>
#include <unistd.h>
#include <vector>

#include "../../offload/offload_library.h"
#include "../../offload/offload_operations.h"

#ifdef __OFFLOAD_CUDA
#include "cuda/cuda_fft_private_header.hpp"
#endif

#ifdef __OFFLOAD_HIP
#include "hip/hip_fft_private_header.hpp"
#endif

#include "pw_gpu_internal.hpp"

#define MAX_NUM_PLANS 8

static bool is_configured{false};

// static blasHandle_t handle_;
static std::vector<offloadStream_t> streams_;
static std::vector<fft_plan> fft_plans_;
static blasHandle_t handle_;

// INIT/RELEASE
extern "C" int pw_gpu_init() {
  if (!is_configured) {
    // it does not really release memory, just initialize the vector length to
    // zero and call the class destructor which does nothing by default.
    //
    // so doing something like fft_plan.push_back(...) does not cost allocation
    // realloc and memory copy unless we have more than 16 fft plans
    fft_plans_.clear();
    offload_set_device();
    streams_.clear();
    streams_.resize(4);
    /* stream allocation  */

#ifdef __PW_CUDA
    for (auto &stream_ : streams_)
      OFFLOAD_CHECK(cudaStreamCreateWithFlags(&stream_, cudaStreamNonBlocking));
#else
    for (auto &stream_ : streams_)
      OFFLOAD_CHECK(hipStreamCreateWithFlags(&stream_, hipStreamNonBlocking));
#endif
    blasCreate(&handle_);
    blasSetStream(handle_, streams_[0]);
    is_configured = true;
  }
  return 0;
}

extern "C" void pw_gpu_finalize() {
  if (is_configured) {
    offload_set_device();
    for (auto &stream_ : streams_)
      offloadStreamDestroy(stream_);
    blasDestroy(handle_);
    is_configured = false;

    for (auto &plan : fft_plans_)
      plan.should_destroy(true);

    fft_plans_.clear();
  }
}

void search_for_plan(const std::vector<int> &fft_size__,
                     const int dim__,
                     const int batch_size__,
                     const enum fft_direction direction__,
                     fft_plan &plan_)
{
  for (auto &plan__ : fft_plans_) {
    if (plan__.is_it_valid(fft_size__, dim__, batch_size__, direction__))
      plan_ = std::move(plan__);
  }

  // we must create the plan
  plan_ = std::move(fft_plan(fft_size__, dim__, batch_size__, direction__));
  if (fft_plans_.size() == MAX_NUM_PLANS) {
    plan_.should_destroy(true);
  } else {
    plan_.should_destroy(false);
    fft_plans_.push_back(std::move(plan_));
  }
}

template <typename T> void offloadMalloc(T **ptr__, size_t size__) {
  offloadMalloc((void **)ptr__, size__ * sizeof(T));
}


inline void retrieve_3d_plan(const offloadStream_t &stream,
                             const enum fft_direction fsign,
                             const int *n, fft_plan &plan)
{
  std::vector<int> size(3);
  size[0] = n[0];
  size[1] = n[1];
  size[2] = n[2];
  search_for_plan(size, 3, 0, fsign, plan);
  plan.set_stream(stream);
}

inline void retrieve_2d_plan(const offloadStream_t &stream,
                             const enum fft_direction fsign,
                             const int *n, fft_plan &plan) {
  std::vector<int> size(2);
  size[0] = n[1];
  size[1] = n[2];
  const int batch = n[0];
  search_for_plan(size, 2, batch, fsign, plan);
  plan.set_stream(stream);
}

inline void retrieve_1d_plan(const offloadStream_t &stream,
                             const enum fft_direction fsign,
                             const int n, const int m,
                             fft_plan &plan) {
  search_for_plan({n}, 1, m, fsign, plan);
  plan.set_stream(stream);
}


// Double precision complex procedures
inline void fft_gpu_run_3d_z_(const offloadStream_t &stream,
                              const enum fft_direction fsign, const int *n,
                              pw_complex_type *data_in__,
                              pw_complex_type *data_out__) {
  std::vector<int> size(3);
  size[0] = n[0];
  size[1] = n[1];
  size[2] = n[2];
  fft_plan plan;
  search_for_plan(size, 3, 0, fsign, plan);

  plan.set_stream(stream);
  plan.execute_fft(fsign, data_in__, data_out__);
}

inline void fft_gpu_run_2dm_z_(const offloadStream_t &stream,
                               const enum fft_direction fsign, const int *n,
                               pw_complex_type *data_in,
                               pw_complex_type *data_out) {
  std::vector<int> size(2);
  size[0] = n[1];
  size[1] = n[2];
  const int batch = n[0];
  fft_plan plan;
  search_for_plan(size, 2, batch, fsign, plan);

  plan.set_stream(stream);

  plan.execute_fft(fsign, data_in, data_out);
}

inline void fft_gpu_run_1dm_z_(const offloadStream_t &stream,
                               const enum fft_direction fsign, const int n,
                               const int m, pw_complex_type *data_in,
                               pw_complex_type *data_out) {
  fft_plan plan;
  search_for_plan({n}, 1, m, fsign, plan);
  plan.set_stream(stream);

  plan.execute_fft(fsign, data_in, data_out);
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) FFT, followed by a (double
 *          precision complex) gather, on the GPU.
 ******************************************************************************/

extern "C" void pw_gpu_cfffg_z_(const double *din, double *zout,
                                const int *ghatmap, const int *npts,
                                const int ngpts, const double scale) {

  if (!is_configured) {
    fprintf(stderr, "call to pw_gpu_init is missing\n");
    abort();
  }

  double *ptr_1, *ptr_2;
  int *ghatmap_dev;
  // dimensions of double and complex arrays
  const int nrpts = npts[0] * npts[1] * npts[2];
  if (nrpts == 0 || ngpts == 0)
    return;

  // get streams and events
  offload_set_device();

  fft_plan plan;
  retrieve_3d_plan(streams_[0], fft_direction::FFT_FORWARD, npts, plan);

  ptr_2 = plan.ptr_2();
  ptr_1 = plan.ptr_1();

  plan.allocate_gmap(ngpts);
  ghatmap_dev = plan.ghatmap();

  // get device memory pointers
  // offloadMemsetAsync(ptr_2, 0, nrpts * sizeof(pw_complex_type), streams_[0]);
  offloadMemcpyAsyncHtoD(ptr_1, din, sizeof(double) * nrpts, streams_[0]);

  // real to complex blow-up
  real_to_complex(streams_[0], nrpts, ptr_1, ptr_2);

  // copy gather map array from host to the device
  offloadMemcpyAsyncHtoD(ghatmap_dev, ghatmap, sizeof(int) * ngpts, streams_[1]);

  // offloadStreamSynchronize(streams_[2]);
  // gpuDcopy(handle_, nrpts, ptr_1, 1, ptr_2, 2);

  // fft on the GPU

  plan.execute_fft(fft_direction::FFT_FORWARD,
                   (pw_complex_type *)ptr_2,
                   (pw_complex_type *)ptr_1);

  offloadStreamSynchronize(streams_[1]);
  // gather on the GPU
  gpu_gather(streams_[0], scale, ngpts, ghatmap_dev, ptr_1, ptr_2);
  offloadMemcpyAsyncDtoH(zout, ptr_2, sizeof(pw_complex_type) * ngpts,
                         streams_[0]);

  // synchronize with respect to host
  offloadStreamSynchronize(streams_[0]);
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) scatter, followed by a inverse
 *          (double precision complex) FFT, on the GPU.
 ******************************************************************************/

extern "C" void pw_gpu_sfffc_z_(const double *zin, double *dout,
                                const int *ghatmap, const int *npts,
                                const int ngpts, const int nmaps,
                                const double scale) {
  if (!is_configured) {
    fprintf(stderr, "call to pw_gpu_init is missing\n");
    abort();
  }

  double *ptr_1, *ptr_2;
  int *ghatmap_dev;

  // dimensions of double and complex arrays
  const int nrpts = npts[0] * npts[1] * npts[2];
  if (nrpts == 0 || ngpts == 0)
    return;

  // get streams and events
  offload_set_device();

  fft_plan plan;
  retrieve_3d_plan(streams_[0], fft_direction::FFT_BACKWARD, npts, plan);
  ptr_2 = plan.ptr_2();
  ptr_1 = plan.ptr_1();

  plan.allocate_gmap(nmaps * ngpts);
  ghatmap_dev = plan.ghatmap();
  // copy all arrays from host to the device
  offloadMemcpyAsyncHtoD(ghatmap_dev, ghatmap, sizeof(int) * nmaps * ngpts,
                         streams_[2]);

  // get device memory pointers
  offloadMemsetAsync(ptr_2, 0, sizeof(pw_complex_type) * nrpts,
                     streams_[0]); // we need to do this only if spherical cut-off is used!


  offloadMemcpyAsyncHtoD(ptr_1, zin, sizeof(pw_complex_type) * ngpts, streams_[0]);


  offloadStreamSynchronize(streams_[2]);
  gpu_scatter(streams_[0], scale, ngpts, nmaps, ghatmap_dev, ptr_1, ptr_2);

  // fft on the GPU (streams_[0][1])

  plan.execute_fft(fft_direction::FFT_BACKWARD,
                   (pw_complex_type *)ptr_2,
                   (pw_complex_type *)ptr_1);

  // convert the complex (device) pointer 'ptr_2' into a real (host)
  // pointer 'dout' (NOTE: Only first half of ptr_1 is written!)
  // pw_copy_cr_cu_z<<<blocksPerGrid, threadsPerBlock, 0, streams_[0]>>>(ptr_2,
  // ptr_1, nrpts);
  //gpuDcopy(handle_, nrpts, ptr_1, 2, ptr_2, 1);
  complex_to_real(streams_[0], nrpts, ptr_1, ptr_2);

  offloadMemcpyAsyncDtoH(dout, ptr_2, sizeof(double) * nrpts, streams_[0]);

  // synchronize with respect to host
  offloadStreamSynchronize(streams_[0]);
}

/*******************************************************************************
 * \brief   Performs a (double to complex double) blow-up and a (double
 *          precision complex) 2D-FFT on the GPU.
 ******************************************************************************/
extern "C" void pw_gpu_cff_z_(const double *din, double *zout,
                              const int *npts) {
  double *ptr_1, *ptr_2;
  dim3 blocksPerGrid, threadsPerBlock;

  if (!is_configured) {
    fprintf(stderr, "call to pw_gpu_init is missing\n");
    abort();
  }

  // dimensions of double and complex arrays
  const int nrpts = npts[0] * npts[1] * npts[2];
  if (nrpts == 0)
    return;

  fft_plan plan1, plan2;
  retrieve_1d_plan(streams_[0], fft_direction::FFT_FORWARD, npts[2],
                   npts[0] * npts[1], plan1);
  retrieve_1d_plan(streams_[0], fft_direction::FFT_FORWARD, npts[1],
                   npts[0] * npts[2], plan2);

  // get streams and events
  offload_set_device();

  // get device memory pointers for:
  // offloadMalloc<pw_complex_type>((pw_complex_type **)&ptr_2, nrpts);

  ptr_1 = plan1.ptr_1();
  offloadMemcpyAsyncHtoD(ptr_1, din, sizeof(double) * nrpts, streams_[1]);

  ptr_2 = plan1.ptr_2();
  offloadMemsetAsync(ptr_2, 0, sizeof(pw_complex_type) * nrpts, streams_[0]);

  // offloadMalloc<pw_complex_type>((pw_complex_type **)&ptr_1, nrpts);
  offloadStreamSynchronize(streams_[1]);

  gpuDcopy(handle_, nrpts, ptr_1, 1, ptr_2, 2);

  // pw_copy_rc_cu_z<<<blocksPerGrid, threadsPerBlock, 0, streams_[0]>>>(
  //     ptr_1, ptr_2, nrpts);

  // fft on the GPU (streams_[0][1])
  // NOTE: the following works, but CUDA does 2D-FFT in C-shaped (not optimal)
  // order fftcu_run_2dm_z_(1, npts, 1.0e0, (cufftDoubleComplex *) ptr_2,
  // (cufftDoubleComplex *) ptr_1, streams_[0][1]);

  plan1.execute_fft(fft_direction::FFT_FORWARD, (pw_complex_type *)ptr_2,
                    (pw_complex_type *)ptr_1);
  plan2.execute_fft(fft_direction::FFT_FORWARD, (pw_complex_type *)ptr_1,
                    (pw_complex_type *)ptr_2);

  // fft_gpu_run_1dm_z_(streams_[0], fft_direction::FFT_FORWARD, npts[2],
  //                    npts[0] * npts[1], (pw_complex_type *)ptr_2,
  //                    (pw_complex_type *)ptr_1);
  // fft_gpu_run_1dm_z_(streams_[0], fft_direction::FFT_FORWARD, npts[1],
  //                    npts[0] * npts[2], (pw_complex_type *)ptr_1,
  //                    (pw_complex_type *)ptr_2);

  offloadMemcpyAsyncDtoH(zout, ptr_2, sizeof(pw_complex_type) * nrpts,
                         streams_[0]);

  // synchronize with respect to host
  offloadStreamSynchronize(streams_[0]);
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) 2D-FFT and a (double complex
 *          to double) shrink-down on the GPU.
 ******************************************************************************/
extern "C" void pw_gpu_ffc_z_(const double *zin, double *dout,
                              const int *npts) {
  if (!is_configured) {
    fprintf(stderr, "call to pw_gpu_init is missing\n");
    abort();
  }

  // dimensions of double and complex arrays
  const int nrpts = npts[0] * npts[1] * npts[2];
  if (nrpts == 0)
    return;

  fft_plan plan1, plan2;
  retrieve_1d_plan(streams_[0], fft_direction::FFT_BACKWARD, npts[1],
                   npts[0] * npts[2], plan1);
  retrieve_1d_plan(streams_[0], fft_direction::FFT_BACKWARD, npts[2],
                   npts[0] * npts[1], plan2);
  // get streams and events
  offload_set_device();

  // get device memory pointers for:
  double *ptr_1 = plan1.ptr_1();
  double *ptr_2 = plan1.ptr_2();

  // copy input data from host to the device
  offloadMemcpyAsyncHtoD(ptr_1, zin, sizeof(pw_complex_type) * nrpts, streams_[0]);

  // fftcu_run_2dm_z_(-1, npts, 1.0e0, (cufftDoubleComplex *) ptr_1,
  // (cufftDoubleComplex *) ptr_2, streams_[0][1]);
  plan1.execute_fft(fft_direction::FFT_BACKWARD, (pw_complex_type *)ptr_1,
                    (pw_complex_type *)ptr_2);
  plan2.execute_fft(fft_direction::FFT_BACKWARD, (pw_complex_type *)ptr_2,
                    (pw_complex_type *)ptr_1);

  // fft_gpu_run_1dm_z_(streams_[0], fft_direction::FFT_BACKWARD, npts[1],
  //                    npts[0] * npts[2], (pw_complex_type *)ptr_1,
  //                    (pw_complex_type *)ptr_2);
  // fft_gpu_run_1dm_z_(streams_[0], fft_direction::FFT_BACKWARD, npts[2],
  //                    npts[0] * npts[1], (pw_complex_type *)ptr_2,
  //                    (pw_complex_type *)ptr_1);

  // CUDA blocking for pw_copy_cr (currently only 2-D grid)
  // get_grid_params(nrpts, MAXTHREADS, threadsPerBlock, blocksPerGrid);

  // take the real part of the FFT and store the result in dout
  gpuDcopy(handle_, nrpts, ptr_1, 2, ptr_2, 1);

  offloadMemcpyAsyncDtoH(dout, ptr_2, sizeof(double) * nrpts, streams_[0]);
  offloadStreamSynchronize(streams_[0]);
}

/*******************************************************************************
 * \brief   Performs a (double to complex double) blow-up and a (double
 *          precision complex) 1D-FFT on the GPU.
 ******************************************************************************/
extern "C" void pw_gpu_cf_z_(const double *din, double *zout, const int *npts) {

  if (!is_configured) {
    fprintf(stderr, "call to pw_gpu_init is missing\n");
    abort();
  }
  // dimensions of double and complex arrays
  const int nrpts = npts[0] * npts[1] * npts[2];
  if (nrpts == 0)
    return;

  // get streams and events
  offload_set_device();

  fft_plan plan;
  retrieve_1d_plan(streams_[0], fft_direction::FFT_FORWARD, npts[2],
                   npts[0] * npts[1], plan);

  // get device memory pointers for:
  double *ptr_1 = plan.ptr_1();
  double *ptr_2 = plan.ptr_2();

  // convert the real (host) pointer 'din' into a complex (device)
  // pointer 'ptr_2' (NOTE: Only first half of ptr_1 is written!)
  // copy all arrays from host to the device
  offloadMemcpyAsyncHtoD(ptr_1, din, sizeof(double) * nrpts, streams_[1]);

  // get device memory pointers for:
  offloadMemsetAsync(ptr_2, 0, sizeof(pw_complex_type) * nrpts, streams_[0]);

  offloadStreamSynchronize(streams_[1]);

  gpuDcopy(handle_, nrpts, ptr_1, 1, ptr_2, 2);

  plan.execute_fft(fft_direction::FFT_FORWARD, (pw_complex_type *)ptr_2,
                   (pw_complex_type *)ptr_1);

  // fft on the GPU (streams_[0][1])

  // fft_gpu_run_1dm_z_(streams_[0], fft_direction::FFT_FORWARD, npts[2],
  //                    npts[0] * npts[1], (pw_complex_type *)ptr_2,
  //                    (pw_complex_type *)ptr_1);

  offloadMemcpyAsyncDtoH(zout, ptr_1, sizeof(pw_complex_type) * nrpts,
                         streams_[0]);

  // synchronize with respect to host
  offloadStreamSynchronize(streams_[0]);
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) 1D-FFT and a (double complex
 *          to double) shrink-down on the GPU.
 ******************************************************************************/
extern "C" void pw_gpu_fc_z_(const double *zin, double *dout, const int *npts) {

  if (!is_configured) {
    fprintf(stderr, "call to pw_gpu_init is missing\n");
    abort();
  }

  // dimensions of double and complex arrays
  const int nrpts = npts[0] * npts[1] * npts[2];
  if (nrpts == 0)
    return;

  // get streams and events
  offload_set_device();

  fft_plan plan;
  retrieve_1d_plan(streams_[0], fft_direction::FFT_BACKWARD, npts[2],
                   npts[0] * npts[1], plan);

  // get device memory pointers for:
  double *ptr_1 = plan.ptr_1();
  double *ptr_2 = plan.ptr_2();

  // get device memory pointers for:
  // copy input data from host to the device
  offloadMemcpyAsyncHtoD(ptr_1, zin, sizeof(pw_complex_type) * nrpts, streams_[0]);

  plan.execute_fft(fft_direction::FFT_BACKWARD, (pw_complex_type *)ptr_1,
                   (pw_complex_type *)ptr_2);

  // fft_gpu_run_1dm_z_(streams_[0], fft_direction::FFT_BACKWARD, npts[2],
  //                    npts[0] * npts[1], (pw_complex_type *)ptr_1,
  //                    (pw_complex_type *)ptr_2);

  // convert the complex (device) pointer 'ptr_2' into a real (host)
  // pointer 'dout' (NOTE: Only first half of ptr_1 is written!)
  gpuDcopy(handle_, nrpts, ptr_2, 2, ptr_1, 1);

  offloadMemcpyAsyncDtoH(dout, ptr_1, sizeof(double) * nrpts, streams_[0]);

  // synchronize with respect to host
  offloadStreamSynchronize(streams_[0]);
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) 1D-FFT on the GPU.
 ******************************************************************************/
extern "C" void pw_gpu_f_z_(const double *zin, double *zout, const int dir,
                            const int n, const int m) {
  if (!is_configured) {
    fprintf(stderr, "call to pw_gpu_init is missing\n");
    abort();
  }

  const enum fft_direction direction =
      (dir > 0) ? (fft_direction::FFT_FORWARD) : (fft_direction::FFT_BACKWARD);
  // dimensions of complex arrays
  const int nrpts = n * m;
  if (nrpts == 0)
    return;

  // get streams and events
  offload_set_device();

  fft_plan plan;
  retrieve_1d_plan(streams_[0], direction, n, m, plan);

  // get device memory pointers for:
  double *ptr_1 = plan.ptr_1();
  double *ptr_2 = plan.ptr_2();

  // get device memory pointers for:
  // copy input data from host to the device
  offloadMemcpyAsyncHtoD(ptr_1, zin, sizeof(pw_complex_type) * nrpts, streams_[0]);

  plan.execute_fft(direction, (pw_complex_type *)ptr_1,
                   (pw_complex_type *)ptr_2);
  // fft_gpu_run_1dm_z_(streams_[0], direction, n, m, (pw_complex_type *)ptr_1,
  //                    (pw_complex_type *)ptr_2);

  offloadMemcpyAsyncDtoH(zout, ptr_2, sizeof(pw_complex_type) * nrpts,
                         streams_[0]);

  // synchronize with respect to host
  offloadStreamSynchronize(streams_[0]);
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) 1D-FFT, followed by a (double
 *          precision complex) gather, on the GPU.
 ******************************************************************************/
extern "C" void pw_gpu_fg_z_(const double *zin, double *zout,
                             const int *ghatmap, const int *npts,
                             const int mmax, const int ngpts,
                             const double scale) {
  if (!is_configured) {
    fprintf(stderr, "call to pw_gpu_init is missing\n");
    abort();
  }


  // dimensions of double and complex arrays
  const int nrpts = npts[0] * mmax;
  if (nrpts == 0 || ngpts == 0)
    return;

  // get streams and events
  offload_set_device();
  fft_plan plan;
  retrieve_1d_plan(streams_[0], fft_direction::FFT_FORWARD, npts[0], mmax, plan);
  // get device memory pointers for:
  double *ptr_1 = plan.ptr_1();
  double *ptr_2 = plan.ptr_2();
  plan.allocate_gmap(ngpts);
  int *ghatmap_dev = plan.ghatmap();

  // get device memory pointers
  // transfer input data from host to device
  offloadMemcpyAsyncHtoD(ptr_1, zin, sizeof(pw_complex_type) * nrpts, streams_[0]);
  // transfer gather data from host to device
  offloadMemcpyAsyncHtoD(ghatmap_dev, ghatmap, sizeof(int) * ngpts, streams_[1]);

  plan.execute_fft(fft_direction::FFT_FORWARD, (pw_complex_type *)ptr_1, (pw_complex_type *)ptr_2);
  // fft_gpu_run_1dm_z_(streams_[0], fft_direction::FFT_FORWARD, npts[0], mmax,
  //                    (pw_complex_type *)ptr_1, (pw_complex_type *)ptr_2);

  // gather on the GPU
  offloadStreamSynchronize(streams_[1]);
  offloadStreamSynchronize(streams_[0]);

  gpu_gather(streams_[0], scale, ngpts, ghatmap_dev, ptr_2, ptr_1);
  offloadMemcpyAsyncDtoH(zout, ptr_1, sizeof(pw_complex_type) * ngpts,
                         streams_[0]);

  // synchronize with respect to host
  offloadStreamSynchronize(streams_[0]);
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) scatter, followed by a
 *          (double precision complex) 1D-FFT, on the GPU.
 ******************************************************************************/
extern "C" void pw_gpu_sf_z_(const double *zin, double *zout,
                             const int *ghatmap, const int *npts,
                             const int mmax, const int ngpts, const int nmaps,
                             const double scale) {

  if (!is_configured) {
    fprintf(stderr, "call to pw_gpu_init is missing\n");
    abort();
  }

  // dimensions of double and complex arrays
  const int nrpts = npts[0] * mmax;
  if (nrpts == 0 || ngpts == 0)
    return;

  // get streams and events
  offload_set_device();
  fft_plan plan;
  retrieve_1d_plan(streams_[0], fft_direction::FFT_BACKWARD, npts[0], mmax, plan);
  // get device memory pointers for:
  double *ptr_1 = plan.ptr_1();
  double *ptr_2 = plan.ptr_2();
  plan.allocate_gmap(ngpts * nmaps);
  int *ghatmap_dev = plan.ghatmap();
  offloadMemsetAsync(ptr_2, 0, sizeof(pw_complex_type) * nrpts,
                     streams_[2]); // we need to do this only if spherical cut-off is used!

  // transfer input data from host to the device
  offloadMemcpyAsyncHtoD(ptr_1, zin, sizeof(pw_complex_type) * ngpts, streams_[0]);
  // transfer scatter data from host to device
  offloadMemcpyAsyncHtoD(ghatmap_dev, ghatmap, sizeof(int) * nmaps * ngpts,
                         streams_[1]);

  offloadStreamSynchronize(streams_[1]);
  offloadStreamSynchronize(streams_[2]);
  // scatter on the GPU
  gpu_scatter(streams_[0], scale, ngpts, nmaps, ghatmap_dev, ptr_1, ptr_2);
  // fft on the GPU

  plan.execute_fft(fft_direction::FFT_BACKWARD, (pw_complex_type *)ptr_2, (pw_complex_type *)ptr_1);
  // fft_gpu_run_1dm_z_(streams_[0], fft_direction::FFT_BACKWARD, npts[0], mmax,
  //                    (pw_complex_type *)ptr_2, (pw_complex_type *)ptr_1);
  offloadMemcpyAsyncDtoH(zout, ptr_1, sizeof(pw_complex_type) * nrpts,
                         streams_[0]);
  // synchronize with respect to host
  offloadStreamSynchronize(streams_[0]);
  // release memory stack
}

#endif // __PW_GPU
