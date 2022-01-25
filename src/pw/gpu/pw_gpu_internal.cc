/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#if defined (__PW_GPU)

#include <array>
#include <cstdio>
#include <map>
#include <string>
#include <unistd.h>
#include <vector>

#include "../../offload/offload_library.h"
#include "../../offload/offload_operations.h"

#ifdef __PW_CUDA
#include "cuda/cuda_fft_private_header.hpp"
#endif

#ifdef __PW_HIP
#include "hip/hip_fft_private_header.hpp"
#endif

static bool is_configured{false};

static blasHandle_t handle_;
static offloadStream_t streams_;

// INIT/RELEASE
extern "C" int pw_gpu_init() {
  if (!is_configured) {
    offload_set_device();
    /* stream allocation  */
    offloadStreamCreate(&streams_);
    blasCreate(&handle_);
    blasSetStream(handle_, streams_);
    is_configured = true;
  }
  return 0;
}

extern "C" void pw_gpu_finalize() {
  if (is_configured) {
    offload_set_device();

    offloadStreamDestroy(streams_);
    blasDestroy(handle_);
    is_configured = false;
  }
}

template <typename T> void offloadMalloc(T **ptr__, size_t size__) {
  offloadMalloc((void **)ptr__, size__ * sizeof(T));
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
  fft_plan plan = fft_plan(size, 3, 0, fsign);

  plan.set_stream(stream);

  plan.execute_fft(fsign, data_in__, data_out__);

  offloadStreamSynchronize(stream);
}

inline void fft_gpu_run_2dm_z_(const offloadStream_t &stream,
                               const enum fft_direction fsign, const int *n,
                               pw_complex_type *data_in,
                               pw_complex_type *data_out) {
  std::vector<int> size(2);
  size[0] = n[1];
  size[1] = n[2];
  const int batch = n[0];
  fft_plan plan = fft_plan(size, 2, batch, fsign);

  plan.set_stream(stream);

  plan.execute_fft(fsign, data_in, data_out);

  offloadStreamSynchronize(stream);
}

inline void fft_gpu_run_1dm_z_(const offloadStream_t &stream,
                               const enum fft_direction fsign, const int n,
                               const int m, pw_complex_type *data_in,
                               pw_complex_type *data_out) {
  fft_plan plan = fft_plan({n}, 1, m, fsign);
  plan.set_stream(stream);

  plan.execute_fft(fsign, data_in, data_out);

  offloadStreamSynchronize(stream);
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

  // get device memory pointers
  offloadMalloc<pw_complex_type>((pw_complex_type **)&ptr_2, nrpts); //? (ngpts)
  offloadMemsetAsync(ptr_2, 0, nrpts * sizeof(pw_complex_type), streams_);

  offloadMalloc<pw_complex_type>((pw_complex_type **)&ptr_1, nrpts);
  offloadMalloc<int>(&ghatmap_dev, ngpts);
  // copy gather map array from host to the device
  offloadMemcpyAsyncHtoD(ghatmap_dev, ghatmap, sizeof(int) * ngpts, streams_);

  // convert the real (host) pointer 'din' into a complex (device)
  // pointer 'ptr_1'
  // copy to device (NOTE: only first half of ptr_1 is written!)
  offloadMemcpyAsyncHtoD(ptr_1, din, sizeof(double) * nrpts, streams_);

  // real to complex blow-up
  gpuDcopy(handle_, nrpts, ptr_1, 1, ptr_2, 2);

  // fft on the GPU
  fft_gpu_run_3d_z_(streams_, fft_direction::FFT_FORWARD, npts,
                    (pw_complex_type *)ptr_2, (pw_complex_type *)ptr_1);

  // gather on the GPU
  gpu_gather(streams_, scale, ngpts, ghatmap_dev, ptr_1, ptr_2);
  offloadMemcpyAsyncDtoH(zout, ptr_2, sizeof(pw_complex_type) * ngpts,
                         streams_);

  // synchronize with respect to host
  offloadStreamSynchronize(streams_);

  // release memory stack
  offloadFree(ptr_1);
  offloadFree(ptr_2);
  offloadFree(ghatmap_dev);
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

  // get device memory pointers
  offloadMalloc<pw_complex_type>((pw_complex_type **)&ptr_2, nrpts);
  offloadMemsetAsync(
      ptr_2, 0, sizeof(pw_complex_type) * nrpts,
      streams_); // we need to do this only if spherical cut-off is used!
  offloadMalloc<pw_complex_type>((pw_complex_type **)&ptr_1, nrpts);
  offloadMalloc<int>(&ghatmap_dev, nmaps * ngpts);

  // copy all arrays from host to the device
  offloadMemcpyAsyncHtoD(ptr_1, zin, sizeof(pw_complex_type) * ngpts, streams_);
  offloadMemcpyAsyncHtoD(ghatmap_dev, ghatmap, sizeof(int) * nmaps * ngpts,
                         streams_);
  gpu_scatter(streams_, scale, ngpts, nmaps, ghatmap_dev, ptr_1, ptr_2);

  // fft on the GPU (streams_[1])
  fft_gpu_run_3d_z_(streams_, fft_direction::FFT_BACKWARD, npts,
                    (pw_complex_type *)ptr_2, (pw_complex_type *)ptr_1);

  // convert the complex (device) pointer 'ptr_2' into a real (host)
  // pointer 'dout' (NOTE: Only first half of ptr_1 is written!)
  // pw_copy_cr_cu_z<<<blocksPerGrid, threadsPerBlock, 0, streams_>>>(ptr_2,
  // ptr_1, nrpts);
  gpuDcopy(handle_, nrpts, ptr_1, 2, ptr_2, 1);

  offloadMemcpyAsyncDtoH(dout, ptr_2, sizeof(double) * nrpts, streams_);

  // synchronize with respect to host
  offloadStreamSynchronize(streams_);

  // release memory stack
  offloadFree(ptr_1);
  offloadFree(ptr_2);
  offloadFree(ghatmap_dev);
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

  // get streams and events
  offload_set_device();

  // get device memory pointers for:
  offloadMalloc<pw_complex_type>((pw_complex_type **)&ptr_2, nrpts);
  offloadMemsetAsync(ptr_2, 0, sizeof(pw_complex_type) * nrpts, streams_);
  offloadMalloc<pw_complex_type>((pw_complex_type **)&ptr_1, nrpts);

  // convert the real (host) pointer 'din' into a complex (device)
  // pointer 'ptr_in'
  // copy all arrays from host to the device (NOTE: Only first half of ptr_1 is
  // written!)
  offloadMemcpyAsyncHtoD(ptr_1, din, sizeof(double) * nrpts, streams_);
  gpuDcopy(handle_, nrpts, ptr_1, 1, ptr_2, 2);

  // pw_copy_rc_cu_z<<<blocksPerGrid, threadsPerBlock, 0, streams_>>>(
  //     ptr_1, ptr_2, nrpts);

  // fft on the GPU (streams_[1])
  // NOTE: the following works, but CUDA does 2D-FFT in C-shaped (not optimal)
  // order fftcu_run_2dm_z_(1, npts, 1.0e0, (cufftDoubleComplex *) ptr_2,
  // (cufftDoubleComplex *) ptr_1, streams_[1]);

  fft_gpu_run_1dm_z_(streams_, fft_direction::FFT_FORWARD, npts[2],
                     npts[0] * npts[1], (pw_complex_type *)ptr_2,
                     (pw_complex_type *)ptr_1);
  fft_gpu_run_1dm_z_(streams_, fft_direction::FFT_FORWARD, npts[1],
                     npts[0] * npts[2], (pw_complex_type *)ptr_1,
                     (pw_complex_type *)ptr_2);

  offloadMemcpyAsyncDtoH(zout, ptr_2, sizeof(pw_complex_type) * nrpts,
                         streams_);

  // synchronize with respect to host
  offloadStreamSynchronize(streams_);

  // release memory stack
  offloadFree(ptr_1);
  offloadFree(ptr_2);
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

  double *ptr_1, *ptr_2;
  // dimensions of double and complex arrays
  const int nrpts = npts[0] * npts[1] * npts[2];
  if (nrpts == 0)
    return;

  // get streams and events
  offload_set_device();

  // get device memory pointers for:
  offloadMalloc<pw_complex_type>((pw_complex_type **)&ptr_1, nrpts);
  // copy input data from host to the device
  offloadMemcpyAsyncHtoD(ptr_1, zin, sizeof(pw_complex_type) * nrpts, streams_);
  offloadMalloc<pw_complex_type>((pw_complex_type **)&ptr_2, nrpts);

  // fftcu_run_2dm_z_(-1, npts, 1.0e0, (cufftDoubleComplex *) ptr_1,
  // (cufftDoubleComplex *) ptr_2, streams_[1]);
  fft_gpu_run_1dm_z_(streams_, fft_direction::FFT_BACKWARD, npts[1],
                     npts[0] * npts[2], (pw_complex_type *)ptr_1,
                     (pw_complex_type *)ptr_2);
  fft_gpu_run_1dm_z_(streams_, fft_direction::FFT_BACKWARD, npts[2],
                     npts[0] * npts[1], (pw_complex_type *)ptr_2,
                     (pw_complex_type *)ptr_1);

  // CUDA blocking for pw_copy_cr (currently only 2-D grid)
  // get_grid_params(nrpts, MAXTHREADS, threadsPerBlock, blocksPerGrid);

  // take the real part of the FFT and store the result in dout
  gpuDcopy(handle_, nrpts, ptr_1, 2, ptr_2, 1);

  offloadMemcpyAsyncDtoH(dout, ptr_2, sizeof(double) * nrpts, streams_);

  offloadStreamSynchronize(streams_);

  // release memory stack
  offloadFree(ptr_1);
  offloadFree(ptr_2);
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

  double *ptr_1, *ptr_2;
  dim3 blocksPerGrid, threadsPerBlock;

  // dimensions of double and complex arrays
  const int nrpts = npts[0] * npts[1] * npts[2];
  if (nrpts == 0)
    return;

  // get streams and events
  offload_set_device();

  offloadMalloc<pw_complex_type>((pw_complex_type **)&ptr_1, nrpts);
  // convert the real (host) pointer 'din' into a complex (device)
  // pointer 'ptr_2' (NOTE: Only first half of ptr_1 is written!)
  // copy all arrays from host to the device
  offloadMemcpyAsyncHtoD(ptr_1, din, sizeof(double) * nrpts, streams_);

  // get device memory pointers for:
  offloadMalloc<pw_complex_type>((pw_complex_type **)&ptr_2, nrpts);
  offloadMemsetAsync(ptr_2, 0, sizeof(pw_complex_type) * nrpts, streams_);

  gpuDcopy(handle_, nrpts, ptr_1, 1, ptr_2, 2);

  // fft on the GPU (streams_[1])
  fft_gpu_run_1dm_z_(streams_, fft_direction::FFT_FORWARD, npts[2],
                     npts[0] * npts[1], (pw_complex_type *)ptr_2,
                     (pw_complex_type *)ptr_1);

  offloadMemcpyAsyncDtoH(zout, ptr_1, sizeof(pw_complex_type) * nrpts,
                         streams_);

  // synchronize with respect to host
  offloadStreamSynchronize(streams_);

  // release memory stack
  offloadFree(ptr_1);
  offloadFree(ptr_2);
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

  double *ptr_1, *ptr_2;
  dim3 blocksPerGrid, threadsPerBlock;

  // dimensions of double and complex arrays
  const int nrpts = npts[0] * npts[1] * npts[2];
  if (nrpts == 0)
    return;

  // get streams and events
  offload_set_device();

  // get device memory pointers for:
  offloadMalloc<pw_complex_type>((pw_complex_type **)&ptr_1, nrpts);
  // copy input data from host to the device
  offloadMemcpyAsyncHtoD(ptr_1, zin, sizeof(pw_complex_type) * nrpts, streams_);
  offloadMalloc<pw_complex_type>((pw_complex_type **)&ptr_2, nrpts);

  fft_gpu_run_1dm_z_(streams_, fft_direction::FFT_BACKWARD, npts[2],
                     npts[0] * npts[1], (pw_complex_type *)ptr_1,
                     (pw_complex_type *)ptr_2);

  // CUDA blocking for pw_copy_cr (currently only 2-D grid)
  // get_grid_params(nrpts, MAXTHREADS, threadsPerBlock, blocksPerGrid);

  // convert the complex (device) pointer 'ptr_2' into a real (host)
  // pointer 'dout' (NOTE: Only first half of ptr_1 is written!)
  gpuDcopy(handle_, nrpts, ptr_2, 2, ptr_1, 1);

  offloadMemcpyAsyncDtoH(dout, ptr_1, sizeof(double) * nrpts, streams_);

  // synchronize with respect to host
  offloadStreamSynchronize(streams_);

  // release memory stack
  offloadFree(ptr_1);
  offloadFree(ptr_2);
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

  double *ptr_1, *ptr_2;
  const enum fft_direction direction =
      (dir > 0) ? (fft_direction::FFT_FORWARD) : (fft_direction::FFT_BACKWARD);
  // dimensions of complex arrays
  const int nrpts = n * m;
  if (nrpts == 0)
    return;

  // get streams and events
  offload_set_device();

  // get device memory pointers for:
  offloadMalloc<pw_complex_type>((pw_complex_type **)&ptr_1, nrpts);
  // copy input data from host to the device
  offloadMemcpyAsyncHtoD(ptr_1, zin, sizeof(pw_complex_type) * nrpts, streams_);

  offloadMalloc<pw_complex_type>((pw_complex_type **)&ptr_2, nrpts);

  fft_gpu_run_1dm_z_(streams_, direction, n, m, (pw_complex_type *)ptr_1,
                     (pw_complex_type *)ptr_2);

  offloadMemcpyAsyncDtoH(zout, ptr_2, sizeof(pw_complex_type) * nrpts,
                         streams_);

  // synchronize with respect to host
  offloadStreamSynchronize(streams_);

  // release memory stack
  offloadFree(ptr_1);
  offloadFree(ptr_2);
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

  double *ptr_1, *ptr_2;
  int *ghatmap_dev;
  //  int nrpts;
  dim3 blocksPerGrid, threadsPerBlock;

  // dimensions of double and complex arrays
  const int nrpts = npts[0] * mmax;
  if (nrpts == 0 || ngpts == 0)
    return;

  // get streams and events
  offload_set_device();

  // get device memory pointers
  offloadMalloc<pw_complex_type>((pw_complex_type **)&ptr_1, nrpts);
  // transfer input data from host to device
  offloadMemcpyAsyncHtoD(ptr_1, zin, sizeof(pw_complex_type) * nrpts, streams_);
  offloadMalloc<int>(&ghatmap_dev, ngpts);
  // transfer gather data from host to device
  offloadMemcpyAsyncHtoD(ghatmap_dev, ghatmap, sizeof(int) * ngpts, streams_);
  offloadMalloc<pw_complex_type>((pw_complex_type **)&ptr_2, nrpts);

  fft_gpu_run_1dm_z_(streams_, fft_direction::FFT_FORWARD, npts[0], mmax,
                     (pw_complex_type *)ptr_1, (pw_complex_type *)ptr_2);

  // gather on the GPU
  gpu_gather(streams_, scale, ngpts, ghatmap_dev, ptr_2, ptr_1);
  offloadMemcpyAsyncDtoH(zout, ptr_1, sizeof(pw_complex_type) * ngpts,
                         streams_);

  // synchronize with respect to host
  offloadStreamSynchronize(streams_);

  // release memory stack
  offloadFree(ptr_1);
  offloadFree(ptr_2);
  offloadFree(ghatmap_dev);
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

  double *ptr_1, *ptr_2;
  int *ghatmap_dev;

  // dimensions of double and complex arrays
  const int nrpts = npts[0] * mmax;
  if (nrpts == 0 || ngpts == 0)
    return;

  // get streams and events
  offload_set_device();

  // get device memory pointers
  offloadMalloc<pw_complex_type>((pw_complex_type **)&ptr_1, nrpts);
  // transfer input data from host to the device
  offloadMemcpyAsyncHtoD(ptr_1, zin, sizeof(pw_complex_type) * ngpts, streams_);

  offloadMalloc<int>(&ghatmap_dev, nmaps * ngpts);
  // transfer scatter data from host to device
  offloadMemcpyAsyncHtoD(ghatmap_dev, ghatmap, sizeof(int) * nmaps * ngpts,
                         streams_);

  offloadMalloc<pw_complex_type>((pw_complex_type **)&ptr_2, nrpts);
  offloadMemsetAsync(
      ptr_2, 0, sizeof(pw_complex_type) * nrpts,
      streams_); // we need to do this only if spherical cut-off is used!

  // scatter on the GPU
  gpu_scatter(streams_, scale, ngpts, nmaps, ghatmap_dev, ptr_1, ptr_2);

  // fft on the GPU
  fft_gpu_run_1dm_z_(streams_, fft_direction::FFT_BACKWARD, npts[0], mmax,
                     (pw_complex_type *)ptr_2, (pw_complex_type *)ptr_1);

  offloadMemcpyAsyncDtoH(zout, ptr_1, sizeof(pw_complex_type) * nrpts,
                         streams_);

  // synchronize with respect to host
  offloadStreamSynchronize(streams_);

  // release memory stack
  offloadFree(ptr_1);
  offloadFree(ptr_2);
  offloadFree(ghatmap_dev);
}

#endif // __PW_GPU
