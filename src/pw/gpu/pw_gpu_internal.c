/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)

#include <assert.h>
#include <omp.h>
#include <stdbool.h>
#include <stddef.h>
#include <string.h>

#include "../../offload/offload_fft.h"
#include "../../offload/offload_library.h"
#include "../../offload/offload_runtime.h"

#include "pw_gpu_kernels.h"

/*******************************************************************************
 * \brief Static variables for retaining objects that are expensive to create.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  int key[4];
  offload_fftHandle *plan;
} cache_entry;

#define PW_GPU_CACHE_SIZE 32
static cache_entry cache[PW_GPU_CACHE_SIZE];
static int cache_oldest_entry = 0; // used for LRU eviction

static double *buffer_dev_1, *buffer_dev_2;
static int *ghatmap_dev;
static size_t allocated_buffer_size, allocated_map_size;

static offloadStream_t stream;
static bool is_initialized = false;

/*******************************************************************************
 * \brief Initializes the pw_gpu library.
 * \author Ole Schuett
 ******************************************************************************/
void pw_gpu_init(void) {
  assert(omp_get_num_threads() == 1);
  if (is_initialized) {
    // fprintf(stderr, "Error: pw_gpu was already initialized.\n");
    // TODO abort();
    return;
  }
  memset(cache, 0, sizeof(cache_entry) * PW_GPU_CACHE_SIZE);
  cache_oldest_entry = 0;

  allocated_buffer_size = 1; // start small
  allocated_map_size = 1;
  offload_activate_chosen_device();
  offloadMalloc((void **)&buffer_dev_1, allocated_buffer_size);
  offloadMalloc((void **)&buffer_dev_2, allocated_buffer_size);
  offloadMalloc((void **)&ghatmap_dev, allocated_map_size);

  offloadStreamCreate(&stream);
  is_initialized = true;
}

/*******************************************************************************
 * \brief Releases resources held by the pw_gpu library.
 * \author Ole Schuett
 ******************************************************************************/
void pw_gpu_finalize(void) {
  assert(omp_get_num_threads() == 1);
  if (!is_initialized) {
    // fprintf(stderr, "Error: pw_gpu is not initialized.\n");
    // TODO abort();
    return;
  }
  for (int i = 0; i < PW_GPU_CACHE_SIZE; i++) {
    if (cache[i].plan != NULL) {
      offload_fftDestroy(*cache[i].plan);
      free(cache[i].plan);
    }
  }
  offloadFree(buffer_dev_1);
  offloadFree(buffer_dev_2);
  offloadFree(ghatmap_dev);
  offloadStreamDestroy(stream);
  is_initialized = false;
}

/*******************************************************************************
 * \brief Checks size of device buffers and re-allocates them if necessary.
 * \author Ole Schuett
 ******************************************************************************/
static void ensure_memory_sizes(const size_t requested_buffer_size,
                                const size_t requested_map_size) {
  assert(is_initialized);
  if (requested_buffer_size > allocated_buffer_size) {
    offloadFree(buffer_dev_1);
    offloadFree(buffer_dev_2);
    offloadMalloc((void **)&buffer_dev_1, requested_buffer_size);
    offloadMalloc((void **)&buffer_dev_2, requested_buffer_size);
    allocated_buffer_size = requested_buffer_size;
  }
  if (requested_map_size > allocated_map_size) {
    offloadFree(ghatmap_dev);
    offloadMalloc((void **)&ghatmap_dev, requested_map_size);
    allocated_map_size = requested_map_size;
  }
}

/*******************************************************************************
 * \brief Fetches an fft plan from the cache. Returns NULL if not found.
 * \author Ole Schuett
 ******************************************************************************/
static offload_fftHandle *lookup_plan_from_cache(const int key[4]) {
  assert(is_initialized);
  for (int i = 0; i < PW_GPU_CACHE_SIZE; i++) {
    const int *x = cache[i].key;
    if (x[0] == key[0] && x[1] == key[1] && x[2] == key[2] && x[3] == key[3]) {
      return cache[i].plan;
    }
  }
  return NULL;
}

/*******************************************************************************
 * \brief Adds an fft plan to the cache. Assumes ownership of plan's memory.
 * \author Ole Schuett
 ******************************************************************************/
static void add_plan_to_cache(const int key[4], offload_fftHandle *plan) {
  const int i = cache_oldest_entry;
  cache_oldest_entry = (cache_oldest_entry + 1) % PW_GPU_CACHE_SIZE;
  if (cache[i].plan != NULL) {
    offload_fftDestroy(*cache[i].plan);
    free(cache[i].plan);
  }
  cache[i].key[0] = key[0];
  cache[i].key[1] = key[1];
  cache[i].key[2] = key[2];
  cache[i].key[3] = key[3];
  cache[i].plan = plan;
}

/*******************************************************************************
 * \brief   Performs a scaled double precision complex 1D-FFT many times on
 *          the GPU.
 *          Input/output are DEVICE pointers (data_in, date_out).
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
static void fft_1d(const offload_fftType fft_type, const int n, const int m,
                   const double *data_in, double *data_out) {
  const int key[4] = {1, fft_type, n, m}; // first key entry is dimensions
  offload_fftHandle *plan = lookup_plan_from_cache(key);

  if (plan == NULL) {
    int nsize[1] = {n};
    int inembed[1] = {0}; // Is ignored, but is not allowed to be NULL.
    int onembed[1] = {0}; // Is ignored, but is not allowed to be NULL.
    int batch = m;
    int istride, idist, ostride, odist;
    if ((int)fft_type == (int)OFFLOAD_FFT_FORWARD) {
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
    plan = malloc(sizeof(cache_entry));
    offload_fftPlanMany(plan, 1, nsize, inembed, istride, idist, onembed,
                        ostride, odist, OFFLOAD_FFT_Z2Z, batch);
    offload_fftSetStream(*plan, stream);
    add_plan_to_cache(key, plan);
  }

  offload_fftExecZ2Z(*plan, data_in, data_out, fft_type);
}

/*******************************************************************************
 * \brief   Performs a scaled double precision complex 3D-FFT on the GPU.
 *          Input/output is a DEVICE pointer (data).
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
static void fft_3d(const offload_fftType fft_type, const int nx, const int ny,
                   const int nz, double *data) {
  const int key[4] = {3, nx, ny, nz}; // first key entry is dimensions
  offload_fftHandle *plan = lookup_plan_from_cache(key);

  if (plan == NULL) {
    plan = malloc(sizeof(cache_entry));
    offload_fftPlan3d(plan, nx, ny, nz, OFFLOAD_FFT_Z2Z);
    offload_fftSetStream(*plan, stream);
    add_plan_to_cache(key, plan);
  }

  offload_fftExecZ2Z(*plan, data, data, fft_type);
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) FFT, followed by a (double
 *          precision complex) gather, on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void pw_gpu_cfffg(const double *din, double *zout, const int *ghatmap,
                  const int *npts, const int ngpts, const double scale) {
  // Check inputs.
  assert(omp_get_num_threads() == 1);
  const int nrpts = npts[0] * npts[1] * npts[2];
  assert(ngpts <= nrpts);
  if (nrpts == 0 || ngpts == 0) {
    return; // Nothing to do.
  }

  // Allocate device memory.
  offload_activate_chosen_device();
  const size_t buffer_size = 2 * sizeof(double) * nrpts;
  const size_t map_size = sizeof(int) * ngpts;
  ensure_memory_sizes(buffer_size, map_size);

  // Upload REAL input and convert to COMPLEX on device.
  offloadMemcpyAsyncHtoD(buffer_dev_1, din, buffer_size / 2, stream);
  pw_gpu_launch_real_to_complex(buffer_dev_1, buffer_dev_2, nrpts, stream);

  // Run FFT on the device.
  fft_3d(OFFLOAD_FFT_FORWARD, npts[2], npts[1], npts[0], buffer_dev_2);

  // Upload map and run gather on the device.
  offloadMemcpyAsyncHtoD(ghatmap_dev, ghatmap, map_size, stream);
  pw_gpu_launch_gather(buffer_dev_1, buffer_dev_2, scale, ngpts, ghatmap_dev,
                       stream);

  // Download COMPLEX results to host.
  offloadMemcpyAsyncDtoH(zout, buffer_dev_1, 2 * sizeof(double) * ngpts,
                         stream);
  offloadStreamSynchronize(stream);
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) scatter, followed by a
 *          (double precision complex) FFT, on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void pw_gpu_sfffc(const double *zin, double *dout, const int *ghatmap,
                  const int *npts, const int ngpts, const int nmaps,
                  const double scale) {
  // Check inputs.
  assert(omp_get_num_threads() == 1);
  const int nrpts = npts[0] * npts[1] * npts[2];
  assert(ngpts <= nrpts);
  if (nrpts == 0 || ngpts == 0) {
    return; // Nothing to do.
  }

  // Allocate device memory.
  offload_activate_chosen_device();
  const size_t buffer_size = 2 * sizeof(double) * nrpts;
  const size_t map_size = sizeof(int) * nmaps * ngpts;
  ensure_memory_sizes(buffer_size, map_size);

  // Upload COMPLEX inputs to device.
  offloadMemcpyAsyncHtoD(buffer_dev_1, zin, 2 * sizeof(double) * ngpts, stream);

  // Upload map and run scatter on the device.
  offloadMemcpyAsyncHtoD(ghatmap_dev, ghatmap, map_size, stream);
  offloadMemsetAsync(buffer_dev_2, 0, buffer_size, stream);
  pw_gpu_launch_scatter(buffer_dev_2, buffer_dev_1, scale, ngpts, nmaps,
                        ghatmap_dev, stream);

  // Run FFT on the device.
  fft_3d(OFFLOAD_FFT_INVERSE, npts[2], npts[1], npts[0], buffer_dev_2);

  // Convert COMPLEX results to REAL and download to host.
  pw_gpu_launch_complex_to_real(buffer_dev_2, buffer_dev_1, nrpts, stream);
  offloadMemcpyAsyncDtoH(dout, buffer_dev_1, buffer_size / 2, stream);
  offloadStreamSynchronize(stream);
}

/*******************************************************************************
 * \brief   Performs a (double to complex double) blow-up and a (double
 *          precision complex) 2D-FFT on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void pw_gpu_cff(const double *din, double *zout, const int *npts) {
  // Check inputs.
  assert(omp_get_num_threads() == 1);
  const int nrpts = npts[0] * npts[1] * npts[2];
  if (nrpts == 0) {
    return; // Nothing to do.
  }

  // Allocate device memory.
  offload_activate_chosen_device();
  const size_t buffer_size = 2 * sizeof(double) * nrpts;
  ensure_memory_sizes(buffer_size, 0);

  // Upload REAL input and convert to COMPLEX on device.
  offloadMemcpyAsyncHtoD(buffer_dev_1, din, buffer_size / 2, stream);
  pw_gpu_launch_real_to_complex(buffer_dev_1, buffer_dev_2, nrpts, stream);

  // Run FFT on the device.
  // NOTE: Could use 2D-FFT, but CUDA does them C-shaped which is not optimal.
  fft_1d(OFFLOAD_FFT_FORWARD, npts[2], npts[0] * npts[1], buffer_dev_2,
         buffer_dev_1);
  fft_1d(OFFLOAD_FFT_FORWARD, npts[1], npts[0] * npts[2], buffer_dev_1,
         buffer_dev_2);

  // Download COMPLEX results to host.
  offloadMemcpyAsyncDtoH(zout, buffer_dev_2, buffer_size, stream);
  offloadStreamSynchronize(stream);
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) 2D-FFT and a (double complex
 *          to double) shrink-down on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void pw_gpu_ffc(const double *zin, double *dout, const int *npts) {
  // Check inputs.
  assert(omp_get_num_threads() == 1);
  const int nrpts = npts[0] * npts[1] * npts[2];
  if (nrpts == 0) {
    return; // Nothing to do.
  }

  // Allocate device memory.
  offload_activate_chosen_device();
  const size_t buffer_size = 2 * sizeof(double) * nrpts;
  ensure_memory_sizes(buffer_size, 0);

  // Upload COMPLEX input to device.
  offloadMemcpyAsyncHtoD(buffer_dev_1, zin, buffer_size, stream);

  // Run FFT on the device.
  // NOTE: Could use 2D-FFT, but CUDA does them C-shaped which is not optimal.
  fft_1d(OFFLOAD_FFT_INVERSE, npts[1], npts[0] * npts[2], buffer_dev_1,
         buffer_dev_2);
  fft_1d(OFFLOAD_FFT_INVERSE, npts[2], npts[0] * npts[1], buffer_dev_2,
         buffer_dev_1);
  pw_gpu_launch_complex_to_real(buffer_dev_1, buffer_dev_2, nrpts, stream);

  // Download REAL results to host.
  offloadMemcpyAsyncDtoH(dout, buffer_dev_2, buffer_size / 2, stream);
  offloadStreamSynchronize(stream);
}

/*******************************************************************************
 * \brief   Performs a (double to complex double) blow-up and a (double
 *          precision complex) 1D-FFT on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void pw_gpu_cf(const double *din, double *zout, const int *npts) {
  // Check inputs.
  assert(omp_get_num_threads() == 1);
  const int nrpts = npts[0] * npts[1] * npts[2];
  if (nrpts == 0) {
    return; // Nothing to do.
  }

  // Allocate device memory.
  offload_activate_chosen_device();
  const size_t buffer_size = 2 * sizeof(double) * nrpts;
  ensure_memory_sizes(buffer_size, 0);

  // Upload REAL input and convert to COMPLEX on device.
  offloadMemcpyAsyncHtoD(buffer_dev_1, din, buffer_size / 2, stream);
  pw_gpu_launch_real_to_complex(buffer_dev_1, buffer_dev_2, nrpts, stream);

  // Run FFT on the device.
  fft_1d(OFFLOAD_FFT_FORWARD, npts[2], npts[0] * npts[1], buffer_dev_2,
         buffer_dev_1);

  // Download COMPLEX results from device.
  offloadMemcpyAsyncDtoH(zout, buffer_dev_1, buffer_size, stream);
  offloadStreamSynchronize(stream);
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) 1D-FFT and a (double complex
 *          to double) shrink-down on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void pw_gpu_fc(const double *zin, double *dout, const int *npts) {
  // Check inputs.
  assert(omp_get_num_threads() == 1);
  const int nrpts = npts[0] * npts[1] * npts[2];
  if (nrpts == 0) {
    return; // Nothing to do.
  }

  // Allocate device memory.
  offload_activate_chosen_device();
  const size_t buffer_size = 2 * sizeof(double) * nrpts;
  ensure_memory_sizes(buffer_size, 0);

  // Upload COMPLEX input to device.
  offloadMemcpyAsyncHtoD(buffer_dev_1, zin, buffer_size, stream);

  // Run FFT on the device.
  fft_1d(OFFLOAD_FFT_INVERSE, npts[2], npts[0] * npts[1], buffer_dev_1,
         buffer_dev_2);

  // Convert COMPLEX results to REAL and download to host.
  pw_gpu_launch_complex_to_real(buffer_dev_2, buffer_dev_1, nrpts, stream);
  offloadMemcpyAsyncDtoH(dout, buffer_dev_1, buffer_size / 2, stream);
  offloadStreamSynchronize(stream);
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) 1D-FFT on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void pw_gpu_f(const double *zin, double *zout, const int dir, const int n,
              const int m) {
  // Check inputs.
  assert(omp_get_num_threads() == 1);
  const int nrpts = n * m;
  if (nrpts == 0) {
    return; // Nothing to do.
  }

  // Allocate device memory.
  offload_activate_chosen_device();
  const size_t buffer_size = 2 * sizeof(double) * nrpts;
  ensure_memory_sizes(buffer_size, 0);

  // Upload COMPLEX input to device.
  offloadMemcpyAsyncHtoD(buffer_dev_1, zin, buffer_size, stream);

  // Run FFT on the device.
  if (dir > 0) {
    fft_1d(OFFLOAD_FFT_FORWARD, n, m, buffer_dev_1, buffer_dev_2);
  } else {
    fft_1d(OFFLOAD_FFT_INVERSE, n, m, buffer_dev_1, buffer_dev_2);
  }

  // Download COMPLEX results from device.
  offloadMemcpyAsyncDtoH(zout, buffer_dev_2, buffer_size, stream);
  offloadStreamSynchronize(stream);
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) 1D-FFT, followed by a (double
 *          precision complex) gather, on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void pw_gpu_fg(const double *zin, double *zout, const int *ghatmap,
               const int *npts, const int mmax, const int ngpts,
               const double scale) {
  // Check inputs.
  assert(omp_get_num_threads() == 1);
  const int nrpts = npts[0] * mmax;
  assert(ngpts <= nrpts);
  if (nrpts == 0 || ngpts == 0) {
    return; // Nothing to do.
  }

  // Allocate device memory.
  offload_activate_chosen_device();
  const size_t buffer_size = 2 * sizeof(double) * nrpts;
  const size_t map_size = sizeof(int) * ngpts;
  ensure_memory_sizes(buffer_size, map_size);

  // Upload COMPLEX inputs to device.
  offloadMemcpyAsyncHtoD(buffer_dev_1, zin, buffer_size, stream);

  // Run FFT on the device.
  fft_1d(OFFLOAD_FFT_FORWARD, npts[0], mmax, buffer_dev_1, buffer_dev_2);

  // Upload map and run gather on the device.
  offloadMemcpyAsyncHtoD(ghatmap_dev, ghatmap, map_size, stream);
  pw_gpu_launch_gather(buffer_dev_1, buffer_dev_2, scale, ngpts, ghatmap_dev,
                       stream);

  // Download COMPLEX results from device.
  offloadMemcpyAsyncDtoH(zout, buffer_dev_1, 2 * sizeof(double) * ngpts,
                         stream);
  offloadStreamSynchronize(stream);
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) scatter, followed by a
 *          (double precision complex) 1D-FFT, on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void pw_gpu_sf(const double *zin, double *zout, const int *ghatmap,
               const int *npts, const int mmax, const int ngpts,
               const int nmaps, const double scale) {
  // Check inputs.
  assert(omp_get_num_threads() == 1);
  const int nrpts = npts[0] * mmax;
  assert(ngpts <= nrpts);
  if (nrpts == 0 || ngpts == 0) {
    return; // Nothing to do.
  }

  // Allocate device memory.
  offload_activate_chosen_device();
  const size_t buffer_size = 2 * sizeof(double) * nrpts;
  const size_t map_size = sizeof(int) * nmaps * ngpts;
  ensure_memory_sizes(buffer_size, map_size);

  // Upload COMPLEX inputs to device.
  offloadMemcpyAsyncHtoD(buffer_dev_1, zin, 2 * sizeof(double) * ngpts, stream);

  // Upload map and run scatter on the device.
  offloadMemcpyAsyncHtoD(ghatmap_dev, ghatmap, map_size, stream);
  offloadMemsetAsync(buffer_dev_2, 0, buffer_size, stream);
  pw_gpu_launch_scatter(buffer_dev_2, buffer_dev_1, scale, ngpts, nmaps,
                        ghatmap_dev, stream);

  // Run FFT on the device.
  fft_1d(OFFLOAD_FFT_INVERSE, npts[0], mmax, buffer_dev_2, buffer_dev_1);

  // Download COMPLEX results from device.
  offloadMemcpyAsyncDtoH(zout, buffer_dev_1, buffer_size, stream);
  offloadStreamSynchronize(stream);
}

#endif // defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)

// EOF
