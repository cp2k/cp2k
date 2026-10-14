/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#include "fft_gpu.h"
#include "../../offload/offload_runtime.h"
#include "../fft_key.h"
#include "../fft_utils.h"
#include "../mpiwrap/cp_mpi.h"

#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
#include "../../offload/offload_fft.h"
#include "../../offload/offload_library.h"
#include "fft_gpu_kernels.h"
#endif

#include <assert.h>
#include <omp.h>
#include <stdbool.h>
#include <stddef.h>
#include <string.h>

#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
/*******************************************************************************
 * \brief Static variables for retaining objects that are expensive to create.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  fft_key_t key;
  offload_fftHandle *plan;
} cache_entry;

#define FFT_GPU_CACHE_SIZE 64
static cache_entry cache[FFT_GPU_CACHE_SIZE];
static int cache_oldest_entry = 0; // used for LRU eviction

static double *buffer_dev_1, *buffer_dev_2;
static int *ghatmap_dev;
static size_t allocated_buffer_size, allocated_map_size;

static offloadStream_t stream;
static bool is_initialized = false;
#endif

/*******************************************************************************
 * \brief Initializes the fft_gpu library.
 * \author Ole Schuett
 ******************************************************************************/
void fft_gpu_init(void) {
  assert(omp_get_num_threads() == 1);
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
  if (is_initialized) {
    // fprintf(stderr, "Error: fft_gpu was already initialized.\n");
    // TODO abort();
    return;
  }
  memset(cache, 0, sizeof(cache_entry) * FFT_GPU_CACHE_SIZE);
  cache_oldest_entry = 0;

  allocated_buffer_size = 1; // start small
  allocated_map_size = 1;
  offload_activate_chosen_device();
  offloadMalloc((void **)&buffer_dev_1, allocated_buffer_size);
  offloadMalloc((void **)&buffer_dev_2, allocated_buffer_size);
  offloadMalloc((void **)&ghatmap_dev, allocated_map_size);

  offloadStreamCreate(&stream);
  is_initialized = true;
#else
// Nothing to initialize
#endif
}

/*******************************************************************************
 * \brief Releases resources held by the fft_gpu library.
 * \author Ole Schuett
 ******************************************************************************/
void fft_gpu_finalize(void) {
  assert(omp_get_num_threads() == 1);
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
  if (!is_initialized) {
    // fprintf(stderr, "Error: fft_gpu is not initialized.\n");
    // TODO abort();
    return;
  }
  for (int i = 0; i < FFT_GPU_CACHE_SIZE; i++) {
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
#else
// Nothing to finalize
#endif
}

/*******************************************************************************
 * \brief Checks size of device buffers and re-allocates them if necessary.
 * \author Ole Schuett
 ******************************************************************************/
void ensure_memory_sizes(const size_t requested_buffer_size,
                         const size_t requested_map_size) {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
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
#else
  (void)requested_buffer_size;
  (void)requested_map_size;
#endif
}

/*******************************************************************************
 * \brief Allocate buffer of type double.
 * \author Frederick Stein
 ******************************************************************************/
void fft_gpu_allocate_double(const int length, double **buffer) {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
  assert(is_initialized);
  offload_host_malloc((void **)buffer, length * sizeof(double));
#else
  assert(0 && "FFT backend was compiled without GPU support!");
  (void)length;
  (void)buffer;
#endif
}

/*******************************************************************************
 * \brief Allocate buffer of type double complex.
 * \author Frederick Stein
 ******************************************************************************/
void fft_gpu_allocate_complex(const int length, double complex **buffer) {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
  assert(is_initialized);
  offload_host_malloc((void **)buffer, length * sizeof(double complex));
#else
  assert(0 && "FFT backend was compiled without GPU support!");
  (void)length;
  (void)buffer;
#endif
}

/*******************************************************************************
 * \brief Allocate buffer of type double.
 * \author Frederick Stein
 ******************************************************************************/
void fft_gpu_free_double(double *buffer) {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
  assert(is_initialized);
  offload_host_free(buffer);
#else
  assert(0 && "FFT backend was compiled without GPU support!");
  (void)buffer;
#endif
}

/*******************************************************************************
 * \brief Allocate buffer of type double complex.
 * \author Frederick Stein
 ******************************************************************************/
void fft_gpu_free_complex(double complex *buffer) {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
  assert(is_initialized);
  offload_host_free(buffer);
#else
  assert(0 && "FFT backend was compiled without GPU support!");
  (void)buffer;
#endif
}

#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
/*******************************************************************************
 * \brief Fetches an fft plan from the cache. Returns NULL if not found.
 * \author Ole Schuett
 ******************************************************************************/
static offload_fftHandle *lookup_plan_from_cache(const fft_key_t key) {
  assert(is_initialized);
  for (int i = 0; i < FFT_GPU_CACHE_SIZE; i++) {
    const int *x = cache[i].key;
    if (memcmp(key, x, sizeof(fft_key_t)) == 0) {
      return cache[i].plan;
    }
  }
  return NULL;
}

/*******************************************************************************
 * \brief Adds an fft plan to the cache. Assumes ownership of plan's memory.
 * \author Ole Schuett
 ******************************************************************************/
static void add_plan_to_cache(const fft_key_t key, offload_fftHandle *plan) {
  const int i = cache_oldest_entry;
  cache_oldest_entry = (cache_oldest_entry + 1) % FFT_GPU_CACHE_SIZE;
  if (cache[i].plan != NULL) {
    offload_fftDestroy(*cache[i].plan);
    free(cache[i].plan);
  }
  memcpy(cache[i].key, key, sizeof(fft_key_t));
  cache[i].plan = plan;
}

/*******************************************************************************
 * \brief   Performs a scaled double precision complex 1D-FFT many times on
 *          the GPU.
 *          Input/output are DEVICE pointers (data_in, date_out).
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
static void fft_register_1d_gpu(const int direction, const int n, const int m,
                                const bool transpose_in,
                                const bool transpose_out,
                                const int leading_dimension_in,
                                const int leading_dimension_out,
                                const bool inplace) {
  fft_key_t key;
  get_key_1d(direction == OFFLOAD_FFT_FORWARD, n, m, transpose_in,
             transpose_out, leading_dimension_in, leading_dimension_out,
             omp_get_max_threads(), inplace, key);

  if (lookup_plan_from_cache(key) == NULL) {
    int rank, number_of_threads, number_of_ffts, istride, ostride, idist, odist;
    bool dir, inpl;
    int fft_size, inembed, onembed;
    fetch_data_from_key_nd(key, &rank, &fft_size, &number_of_ffts,
                           &number_of_threads, &dir, &inpl, &inembed, &onembed,
                           &idist, &odist, &istride, &ostride);
    offload_fftHandle *plan = malloc(sizeof(cache_entry));
    offload_fftPlanMany(plan, 1, &fft_size, &inembed, istride, idist, &onembed,
                        ostride, odist, OFFLOAD_FFT_Z2Z, number_of_ffts);
    offload_fftSetStream(*plan, stream);
    add_plan_to_cache(key, plan);
  }
}

/*******************************************************************************
 * \brief   Performs a scaled double precision complex 1D-FFT many times on
 *          the GPU.
 *          Input/output are DEVICE pointers (data_in, date_out).
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
static void fft_register_r2c_1d_gpu(const int direction, const int n,
                                    const int m, const bool transpose_in,
                                    const bool transpose_out,
                                    const int leading_dimension_in,
                                    const int leading_dimension_out,
                                    const bool inplace) {
  fft_key_t key;
  get_key_1d_r2c(direction == OFFLOAD_FFT_FORWARD, n, m, transpose_in,
                 transpose_out, leading_dimension_in, leading_dimension_out,
                 omp_get_max_threads(), inplace, key);

  if (lookup_plan_from_cache(key) == NULL) {
    int rank, number_of_threads, number_of_ffts, istride, ostride, idist, odist;
    bool dir, inpl;
    int fft_size, inembed, onembed;
    fetch_data_from_key_nd(key, &rank, &fft_size, &number_of_ffts,
                           &number_of_threads, &dir, &inpl, &inembed, &onembed,
                           &idist, &odist, &istride, &ostride);
    offload_fftHandle *plan = malloc(sizeof(cache_entry));
    if (direction == OFFLOAD_FFT_FORWARD) {
      offload_fftPlanMany(plan, 1, &fft_size, &inembed, istride, idist,
                          &onembed, ostride, odist, OFFLOAD_FFT_D2Z,
                          number_of_ffts);
    } else {
      offload_fftPlanMany(plan, 1, &fft_size, &inembed, istride, idist,
                          &onembed, ostride, odist, OFFLOAD_FFT_Z2D,
                          number_of_ffts);
    }
    offload_fftSetStream(*plan, stream);
    add_plan_to_cache(key, plan);
  }
}

/*******************************************************************************
 * \brief   Performs a scaled double precision complex 1D-FFT many times on
 *          the GPU.
 *          Input/output are DEVICE pointers (data_in, date_out).
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
static void fft_register_2d_gpu(const int direction, const int n[2],
                                const int m, const bool transpose_in,
                                const bool transpose_out, const bool inplace) {
  fft_key_t key;
  get_key_2d(direction == OFFLOAD_FFT_FORWARD, n, m, transpose_in,
             transpose_out, transpose_in ? m : n[0] * n[1],
             transpose_out ? m : n[0] * n[1], omp_get_max_threads(), inplace,
             key);

  if (lookup_plan_from_cache(key) == NULL) {
    int rank, number_of_threads, number_of_ffts, istride, ostride, idist, odist;
    bool dir, inpl;
    int fft_size[2], inembed[2], onembed[2];
    fetch_data_from_key_nd(key, &rank, fft_size, &number_of_ffts,
                           &number_of_threads, &dir, &inpl, inembed, onembed,
                           &idist, &odist, &istride, &ostride);
    offload_fftHandle *plan = malloc(sizeof(cache_entry));
    offload_fftPlanMany(plan, 2, fft_size, inembed, istride, idist, onembed,
                        ostride, odist, OFFLOAD_FFT_Z2Z, number_of_ffts);
    offload_fftSetStream(*plan, stream);
    add_plan_to_cache(key, plan);
  }
}

/*******************************************************************************
 * \brief   Performs a scaled double precision complex R2C/C2R-2D-FFT many times
 * on the GPU. Input/output are DEVICE pointers (data_in, date_out).
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
static void fft_register_r2c_2d_gpu(const int direction, const int n[2],
                                    const int m, const bool transpose_in,
                                    const bool transpose_out,
                                    const bool inplace) {
  fft_key_t key;
  get_key_2d_r2c(
      direction == OFFLOAD_FFT_FORWARD, n, m, transpose_in, transpose_out,
      transpose_in
          ? m
          : n[0] * (direction == OFFLOAD_FFT_FORWARD
                        ? (inplace ? 2 * (n[1] / 2 + 1) : n[1])
                        : n[1] / 2 + 1),
      transpose_in
          ? m
          : n[0] * (direction == OFFLOAD_FFT_FORWARD
                        ? n[1] / 2 + 1
                        : (inplace ? 2 * (n[1] / 2 + 1) : n[1])),
      omp_get_max_threads(), inplace, key);

  if (lookup_plan_from_cache(key) == NULL) {
    int rank, number_of_threads, number_of_ffts, istride, ostride, idist, odist;
    bool dir, inpl;
    int fft_size[2], inembed[2], onembed[2];
    fetch_data_from_key_nd(key, &rank, fft_size, &number_of_ffts,
                           &number_of_threads, &dir, &inpl, inembed, onembed,
                           &idist, &odist, &istride, &ostride);
    offload_fftHandle *plan = malloc(sizeof(cache_entry));
    offload_fftPlanMany(plan, 2, fft_size, inembed, istride, idist, onembed,
                        ostride, odist, dir ? OFFLOAD_FFT_D2Z : OFFLOAD_FFT_Z2D,
                        number_of_ffts);
    offload_fftSetStream(*plan, stream);
    add_plan_to_cache(key, plan);
  }
}

/*******************************************************************************
 * \brief   Performs a scaled double precision complex 3D-FFT on the GPU.
 *          Input/output is a DEVICE pointer (data).
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
static void fft_register_3d_gpu(const int direction, const int nx, const int ny,
                                const int nz, const bool inplace) {
  fft_key_t key;
  get_key_3d(direction == OFFLOAD_FFT_FORWARD, (const int[]){nx, ny, nz},
             omp_get_max_threads(), inplace, key);

  if (lookup_plan_from_cache(key) == NULL) {
    offload_fftHandle *plan = malloc(sizeof(cache_entry));
    offload_fftPlan3d(plan, nx, ny, nz, OFFLOAD_FFT_Z2Z);
    offload_fftSetStream(*plan, stream);
    add_plan_to_cache(key, plan);
  }
}

/*******************************************************************************
 * \brief   Performs a scaled double precision complex 3D-FFT on the GPU.
 *          Input/output is a DEVICE pointer (data).
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
static void fft_register_r2c_3d_gpu(const int direction, const int nx,
                                    const int ny, const int nz,
                                    const bool inplace) {
  fft_key_t key;
  get_key_3d_r2c(direction == OFFLOAD_FFT_FORWARD, (const int[]){nx, ny, nz},
                 omp_get_max_threads(), inplace, key);

  if (lookup_plan_from_cache(key) == NULL) {
    offload_fftHandle *plan = malloc(sizeof(cache_entry));
    offload_fftPlan3d(plan, nx, ny, nz,
                      direction == OFFLOAD_FFT_FORWARD ? OFFLOAD_FFT_D2Z
                                                       : OFFLOAD_FFT_Z2D);
    offload_fftSetStream(*plan, stream);
    add_plan_to_cache(key, plan);
  }
}

/*******************************************************************************
 * \brief   Performs a scaled double precision complex 1D-FFT many times on
 *          the GPU.
 *          Input/output are DEVICE pointers (data_in, date_out).
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
static void fft_1d_gpu(const int direction, const int n, const int m,
                       const bool transpose_in, const bool transpose_out,
                       const int leading_dimension_in,
                       const int leading_dimension_out, const double *data_in,
                       double *data_out) {
  fft_key_t key;
  get_key_1d(direction == OFFLOAD_FFT_FORWARD, n, m, transpose_in,
             transpose_out, leading_dimension_in, leading_dimension_out,
             omp_get_max_threads(), data_in == data_out, key);
  offload_fftHandle *plan = lookup_plan_from_cache(key);

  if (plan == NULL) {
    assert(false);
    int rank, number_of_threads, number_of_ffts, istride, ostride, idist, odist;
    bool dir, inplace;
    int fft_size, inembed, onembed;
    fetch_data_from_key_nd(key, &rank, &fft_size, &number_of_ffts,
                           &number_of_threads, &dir, &inplace, &inembed,
                           &onembed, &idist, &odist, &istride, &ostride);
    plan = malloc(sizeof(cache_entry));
    offload_fftPlanMany(plan, 1, &fft_size, &inembed, istride, idist, &onembed,
                        ostride, odist, OFFLOAD_FFT_Z2Z, number_of_ffts);
    offload_fftSetStream(*plan, stream);
    add_plan_to_cache(key, plan);
  }

  offload_fftExecZ2Z(*plan, data_in, data_out, direction);
}

/*******************************************************************************
 * \brief   Performs a scaled double precision complex 1D-FFT many times on
 *          the GPU.
 *          Input/output are DEVICE pointers (data_in, date_out).
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
static void fft_r2c_1d_gpu(const int direction, const int n, const int m,
                           const bool transpose_in, const bool transpose_out,
                           const int leading_dimension_in,
                           const int leading_dimension_out,
                           const double *data_in, double *data_out) {
  fft_key_t key;
  get_key_1d_r2c(direction == OFFLOAD_FFT_FORWARD, n, m, transpose_in,
                 transpose_out, leading_dimension_in, leading_dimension_out,
                 omp_get_max_threads(), data_in == data_out, key);
  offload_fftHandle *plan = lookup_plan_from_cache(key);

  if (plan == NULL) {
    assert(false);
    int rank, number_of_threads, number_of_ffts, istride, ostride, idist, odist;
    bool dir, inplace;
    int fft_size, inembed, onembed;
    fetch_data_from_key_nd(key, &rank, &fft_size, &number_of_ffts,
                           &number_of_threads, &dir, &inplace, &inembed,
                           &onembed, &idist, &odist, &istride, &ostride);
    plan = malloc(sizeof(cache_entry));
    offload_fftPlanMany(plan, 1, &fft_size, &inembed, istride, idist, &onembed,
                        ostride, odist, dir ? OFFLOAD_FFT_D2Z : OFFLOAD_FFT_Z2D,
                        number_of_ffts);
    offload_fftSetStream(*plan, stream);
    add_plan_to_cache(key, plan);
  }

  if (direction == OFFLOAD_FFT_FORWARD) {
    offload_fftExecD2Z(*plan, data_in, data_out);
  } else {
    offload_fftExecZ2D(*plan, data_in, data_out);
  }
}

/*******************************************************************************
 * \brief   Performs a scaled double precision complex 1D-FFT many times on
 *          the GPU.
 *          Input/output are DEVICE pointers (data_in, date_out).
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
static void fft_2d_gpu(const int direction, const int n[2], const int m,
                       const bool transpose_in, const bool transpose_out,
                       const double *data_in, double *data_out) {
  fft_key_t key;
  get_key_2d(direction == OFFLOAD_FFT_FORWARD, n, m, transpose_in,
             transpose_out, transpose_in ? m : n[0] * n[1],
             transpose_out ? m : n[0] * n[1], omp_get_max_threads(),
             data_in == data_out, key);
  offload_fftHandle *plan = lookup_plan_from_cache(key);

  if (plan == NULL) {
    assert(false);
    int rank, number_of_threads, number_of_ffts, istride, ostride, idist, odist;
    bool dir, inplace;
    int fft_size[2], inembed[2], onembed[2];
    fetch_data_from_key_nd(key, &rank, fft_size, &number_of_ffts,
                           &number_of_threads, &dir, &inplace, inembed, onembed,
                           &idist, &odist, &istride, &ostride);
    plan = malloc(sizeof(cache_entry));
    offload_fftPlanMany(plan, 2, fft_size, inembed, istride, idist, onembed,
                        ostride, odist, OFFLOAD_FFT_Z2Z, number_of_ffts);
    offload_fftSetStream(*plan, stream);
    add_plan_to_cache(key, plan);
  }

  offload_fftExecZ2Z(*plan, data_in, data_out, direction);
}

/*******************************************************************************
 * \brief   Performs a scaled double precision complex R2C/C2R-2D-FFT many times
 * on the GPU. Input/output are DEVICE pointers (data_in, date_out).
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
static void fft_r2c_2d_gpu(const int direction, const int n[2], const int m,
                           const bool transpose_in, const bool transpose_out,
                           const double *data_in, double *data_out) {
  fft_key_t key;
  get_key_2d_r2c(
      direction == OFFLOAD_FFT_FORWARD, n, m, transpose_in, transpose_out,
      transpose_in
          ? m
          : n[0] * (direction == OFFLOAD_FFT_FORWARD
                        ? (data_in == data_out ? 2 * (n[1] / 2 + 1) : n[1])
                        : n[1] / 2 + 1),
      transpose_in
          ? m
          : n[0] * (direction == OFFLOAD_FFT_FORWARD
                        ? n[1] / 2 + 1
                        : (data_in == data_out ? 2 * (n[1] / 2 + 1) : n[1])),
      omp_get_max_threads(), data_in == data_out, key);
  offload_fftHandle *plan = lookup_plan_from_cache(key);

  if (plan == NULL) {
    assert(false);
    int rank, number_of_threads, number_of_ffts, istride, ostride, idist, odist;
    bool dir, inplace;
    int fft_size[2], inembed[2], onembed[2];
    fetch_data_from_key_nd(key, &rank, fft_size, &number_of_ffts,
                           &number_of_threads, &dir, &inplace, inembed, onembed,
                           &idist, &odist, &istride, &ostride);
    plan = malloc(sizeof(cache_entry));
    offload_fftPlanMany(plan, 2, fft_size, inembed, istride, idist, onembed,
                        ostride, odist, dir ? OFFLOAD_FFT_D2Z : OFFLOAD_FFT_Z2D,
                        number_of_ffts);
    offload_fftSetStream(*plan, stream);
    add_plan_to_cache(key, plan);
  }

  if (direction == OFFLOAD_FFT_FORWARD) {
    offload_fftExecD2Z(*plan, data_in, data_out);
  } else {
    offload_fftExecZ2D(*plan, data_in, data_out);
  }
}

/*******************************************************************************
 * \brief   Performs a scaled double precision complex 3D-FFT on the GPU.
 *          Input/output is a DEVICE pointer (data).
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
static void fft_3d_gpu(const int direction, const int nx, const int ny,
                       const int nz, double *data_in, double *data_out) {
  fft_key_t key;
  get_key_3d(direction == OFFLOAD_FFT_FORWARD, (const int[]){nx, ny, nz},
             omp_get_max_threads(), data_in == data_out, key);
  offload_fftHandle *plan = lookup_plan_from_cache(key);

  if (plan == NULL) {
    assert(false);
    plan = malloc(sizeof(cache_entry));
    offload_fftPlan3d(plan, nx, ny, nz, OFFLOAD_FFT_Z2Z);
    offload_fftSetStream(*plan, stream);
    add_plan_to_cache(key, plan);
  }

  offload_fftExecZ2Z(*plan, data_in, data_out, direction);
}

/*******************************************************************************
 * \brief   Performs a scaled double precision complex 3D-FFT on the GPU.
 *          Input/output is a DEVICE pointer (data).
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
static void fft_r2c_3d_gpu(const int direction, const int nx, const int ny,
                           const int nz, const double *data_in,
                           double *data_out) {
  fft_key_t key;
  get_key_3d_r2c(direction == OFFLOAD_FFT_FORWARD, (const int[]){nx, ny, nz},
                 omp_get_max_threads(), data_in == data_out, key);
  offload_fftHandle *plan = lookup_plan_from_cache(key);

  if (plan == NULL) {
    assert(false);
    plan = malloc(sizeof(cache_entry));
    offload_fftPlan3d(plan, nx, ny, nz,
                      direction == OFFLOAD_FFT_FORWARD ? OFFLOAD_FFT_D2Z
                                                       : OFFLOAD_FFT_Z2D);
    offload_fftSetStream(*plan, stream);
    add_plan_to_cache(key, plan);
  }

  if (direction == OFFLOAD_FFT_FORWARD) {
    offload_fftExecD2Z(*plan, data_in, data_out);
  } else {
    offload_fftExecZ2D(*plan, data_in, data_out);
  }
}
#endif

/*******************************************************************************
 * \brief   Performs a (double precision complex) 3D-FFT on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void fft_register_gpu_fff(const bool dir, const int *npts, const bool inplace) {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
  // Check inputs.
  assert(omp_get_num_threads() == 1);
  if (npts[0] == 0 || npts[1] == 0 || npts[2] == 0) {
    return; // Nothing to do.
  }
  // Run FFT on the device.
  fft_register_3d_gpu(dir ? OFFLOAD_FFT_FORWARD : OFFLOAD_FFT_INVERSE, npts[0],
                      npts[1], npts[2], inplace);
#else
  (void)dir;
  (void)npts;
  (void)inplace;
#endif
}

/*******************************************************************************
 * \brief   Performs a 3D-R2C/C2R-FFT, on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void fft_register_r2c_gpu_fff(const bool dir, const int *npts,
                              const bool inplace) {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
  // Check inputs.
  assert(omp_get_num_threads() == 1);
  if (npts[0] == 0 || npts[1] == 0 || npts[2] == 0) {
    return; // Nothing to do.
  }

  // Run FFT on the device.
  fft_register_r2c_3d_gpu(dir ? OFFLOAD_FFT_FORWARD : OFFLOAD_FFT_INVERSE,
                          npts[0], npts[1], npts[2], inplace);
#else
  (void)dir;
  (void)npts;
  (void)inplace;
#endif
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) 1D-FFT on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void fft_register_gpu_f(const bool dir, const int n, const int m,
                        const bool transpose_in, const bool transpose_out,
                        const int leading_dimension_in,
                        const int leading_dimension_out, const bool inplace) {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
  // Check inputs.
  assert(omp_get_num_threads() == 1);
  if (n == 0 || m == 0) {
    return; // Nothing to do.
  }

  // Run FFT on the device.
  fft_register_1d_gpu(dir ? OFFLOAD_FFT_FORWARD : OFFLOAD_FFT_INVERSE, n, m,
                      transpose_in, transpose_out, leading_dimension_in,
                      leading_dimension_out, inplace);
#else
  (void)dir;
  (void)n;
  (void)m;
  (void)transpose_in;
  (void)transpose_out;
  (void)leading_dimension_in;
  (void)leading_dimension_out;
  (void)inplace;
#endif
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) 1D-FFT on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void fft_register_r2c_gpu_f(const bool dir, const int n, const int m,
                            const bool transpose_in, const bool transpose_out,
                            const int leading_dimension_in,
                            const int leading_dimension_out,
                            const bool inplace) {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
  // Check inputs.
  assert(omp_get_num_threads() == 1);
  const int nrpts = n * m;
  if (nrpts == 0) {
    return;
  }

  // Run FFT on the device.
  fft_register_r2c_1d_gpu(dir ? OFFLOAD_FFT_FORWARD : OFFLOAD_FFT_INVERSE, n, m,
                          transpose_in, transpose_out, leading_dimension_in,
                          leading_dimension_out, inplace);
#else
  (void)dir;
  (void)n;
  (void)m;
  (void)transpose_in;
  (void)transpose_out;
  (void)leading_dimension_in;
  (void)leading_dimension_out;
  (void)inplace;
#endif
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) 2D-FFT on the GPU.
 * \author  Frederick Stein
 ******************************************************************************/
void fft_register_gpu_ff(const bool dir, const int n[2], const int m,
                         const bool transpose_in, const bool transpose_out,
                         const bool inplace) {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
  // Check inputs.
  assert(omp_get_num_threads() == 1);
  const int nrpts = n[0] * n[1] * m;
  if (nrpts == 0) {
    return;
  }

  // Run FFT on the device.
  fft_register_2d_gpu(dir ? OFFLOAD_FFT_FORWARD : OFFLOAD_FFT_INVERSE, n, m,
                      transpose_in, transpose_out, inplace);
#else
  (void)dir;
  (void)n;
  (void)m;
  (void)transpose_in;
  (void)transpose_out;
  (void)inplace;
#endif
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) R2C/C2R2D-FFT on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void fft_register_r2c_gpu_ff(const bool dir, const int n[2], const int m,
                             const bool transpose_in, const bool transpose_out,
                             const bool inplace) {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
  // Check inputs.
  assert(omp_get_num_threads() == 1);
  if (n[0] == 0 || n[1] == 0 || m == 0) {
    return;
  }

  // Run FFT on the device.
  fft_register_r2c_2d_gpu(dir ? OFFLOAD_FFT_FORWARD : OFFLOAD_FFT_INVERSE, n, m,
                          transpose_in, transpose_out, inplace);
#else
  (void)dir;
  (void)n;
  (void)m;
  (void)transpose_in;
  (void)transpose_out;
  (void)inplace;
#endif
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) FFT, followed by a (double
 *          precision complex) gather, on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void fft_gpu_cfffg(const double *din, double *zout, const int *ghatmap,
                   const int *npts, const int ngpts, const double scale) {
  // Check inputs.
  assert(omp_get_num_threads() == 1);
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
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
  fft_gpu_launch_real_to_complex(buffer_dev_1, buffer_dev_2, nrpts, stream);

  // Run FFT on the device.
  fft_3d_gpu(OFFLOAD_FFT_FORWARD, npts[2], npts[1], npts[0], buffer_dev_2,
             buffer_dev_2);

  // Upload map and run gather on the device.
  offloadMemcpyAsyncHtoD(ghatmap_dev, ghatmap, map_size, stream);
  fft_gpu_launch_gather(buffer_dev_1, buffer_dev_2, scale, ngpts, ghatmap_dev,
                        stream);

  // Download COMPLEX results to host.
  offloadMemcpyAsyncDtoH(zout, buffer_dev_1, 2 * sizeof(double) * ngpts,
                         stream);
  offloadStreamSynchronize(stream);
#else
  (void)din;
  (void)zout;
  (void)ghatmap;
  (void)npts;
  (void)ngpts;
  (void)scale;
#endif
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) scatter, followed by a
 *          (double precision complex) FFT, on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void fft_gpu_sfffc(const double *zin, double *dout, const int *ghatmap,
                   const int *npts, const int ngpts, const int nmaps,
                   const double scale) {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
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
  fft_gpu_launch_scatter(buffer_dev_2, buffer_dev_1, scale, ngpts, nmaps,
                         ghatmap_dev, stream);

  // Run FFT on the device.
  fft_3d_gpu(OFFLOAD_FFT_INVERSE, npts[2], npts[1], npts[0], buffer_dev_2,
             buffer_dev_2);

  // Convert COMPLEX results to REAL and download to host.
  fft_gpu_launch_complex_to_real(buffer_dev_2, buffer_dev_1, nrpts, stream);
  offloadMemcpyAsyncDtoH(dout, buffer_dev_1, buffer_size / 2, stream);
  offloadStreamSynchronize(stream);
#else
  (void)zin;
  (void)dout;
  (void)ghatmap;
  (void)npts;
  (void)ngpts;
  (void)nmaps;
  (void)scale;
#endif
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) 3D-FFT on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void fft_gpu_fff(const double *zin, double *zout, const bool dir,
                 const int *npts) {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
  // Check inputs.
  assert(omp_get_num_threads() == 1);
  if (npts[0] == 0 || npts[1] == 0 || npts[2] == 0) {
    return; // Nothing to do.
  }

  // Allocate device memory.
  offload_activate_chosen_device();
  const size_t buffer_size = 2 * sizeof(double) * npts[0] * npts[1] * npts[2];
  ensure_memory_sizes(buffer_size, 0);

  // Upload COMPLEX inputs to device.
  offloadMemcpyAsyncHtoD(buffer_dev_1, zin, buffer_size, stream);

  // Run FFT on the device.
  fft_3d_gpu(dir ? OFFLOAD_FFT_FORWARD : OFFLOAD_FFT_INVERSE, npts[0], npts[1],
             npts[2], buffer_dev_1, zin != zout ? buffer_dev_2 : buffer_dev_1);

  // Download to host
  offloadMemcpyAsyncDtoH(zout, zin != zout ? buffer_dev_2 : buffer_dev_1,
                         buffer_size, stream);
  offloadStreamSynchronize(stream);
#else
  (void)zin;
  (void)zout;
  (void)dir;
  (void)npts;
#endif
}

/*******************************************************************************
 * \brief   Performs a 3D-R2C/C2R-FFT, on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void fft_r2c_gpu_fff(const double *zin, double *zout, const bool dir,
                     const int *npts) {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
  // Check inputs.
  assert(omp_get_num_threads() == 1);
  if (npts[0] == 0 || npts[1] == 0 || npts[2] == 0) {
    return; // Nothing to do.
  }

  // Allocate device memory. cuFFT halves the last dimension, and the
  // transform is done out-of-place so that neither side needs padding.
  offload_activate_chosen_device();
  const size_t real_size = sizeof(double) * npts[0] * npts[1] *
                           (zin != zout ? npts[2] : 2 * (npts[2] / 2 + 1));
  const size_t complex_size =
      sizeof(double) * 2 * npts[0] * npts[1] * (npts[2] / 2 + 1);
  ensure_memory_sizes(complex_size > real_size ? complex_size : real_size, 0);

  // Upload inputs to device.
  offloadMemcpyAsyncHtoD(buffer_dev_1, zin, dir ? real_size : complex_size,
                         stream);

  // Run FFT on the device.
  fft_r2c_3d_gpu(dir ? OFFLOAD_FFT_FORWARD : OFFLOAD_FFT_INVERSE, npts[0],
                 npts[1], npts[2], buffer_dev_1,
                 zin != zout ? buffer_dev_2 : buffer_dev_1);

  // Download to host
  offloadMemcpyAsyncDtoH(zout, zin != zout ? buffer_dev_2 : buffer_dev_1,
                         dir ? complex_size : real_size, stream);
  offloadStreamSynchronize(stream);
#else
  (void)zin;
  (void)zout;
  (void)dir;
  (void)npts;
#endif
}

/*******************************************************************************
 * \brief   Performs a (double to complex double) blow-up and a (double
 *          precision complex) 2D-FFT on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void fft_gpu_cff(const double *din, double *zout, const int *npts) {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
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
  fft_gpu_launch_real_to_complex(buffer_dev_1, buffer_dev_2, nrpts, stream);

  // Run FFT on the device.
  // NOTE: Could use 2D-FFT, but CUDA does them C-shaped which is not optimal.
  fft_1d_gpu(OFFLOAD_FFT_FORWARD, npts[2], npts[0] * npts[1], false, false,
             npts[2], npts[2], buffer_dev_2, buffer_dev_1);
  fft_1d_gpu(OFFLOAD_FFT_FORWARD, npts[1], npts[0] * npts[2], false, false,
             npts[1], npts[1], buffer_dev_1, buffer_dev_2);

  // Download COMPLEX results to host.
  offloadMemcpyAsyncDtoH(zout, buffer_dev_2, buffer_size, stream);
  offloadStreamSynchronize(stream);
#else
  (void)din;
  (void)zout;
  (void)npts;
#endif
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) 2D-FFT and a (double complex
 *          to double) shrink-down on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void fft_gpu_ffc(const double *zin, double *dout, const int *npts) {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
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
  fft_1d_gpu(OFFLOAD_FFT_INVERSE, npts[1], npts[0] * npts[2], false, false,
             npts[1], npts[1], buffer_dev_1, buffer_dev_2);
  fft_1d_gpu(OFFLOAD_FFT_INVERSE, npts[2], npts[0] * npts[1], false, false,
             npts[2], npts[2], buffer_dev_2, buffer_dev_1);
  fft_gpu_launch_complex_to_real(buffer_dev_1, buffer_dev_2, nrpts, stream);

  // Download REAL results to host.
  offloadMemcpyAsyncDtoH(dout, buffer_dev_2, buffer_size / 2, stream);
  offloadStreamSynchronize(stream);
#else
  (void)zin;
  (void)dout;
  (void)npts;
#endif
}

/*******************************************************************************
 * \brief   Performs a (double to complex double) blow-up and a (double
 *          precision complex) 1D-FFT on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void fft_gpu_cf(const double *din, double *zout, const int *npts) {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
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
  fft_gpu_launch_real_to_complex(buffer_dev_1, buffer_dev_2, nrpts, stream);

  // Run FFT on the device.
  fft_1d_gpu(OFFLOAD_FFT_FORWARD, npts[2], npts[0] * npts[1], false, false,
             npts[2], npts[2], buffer_dev_2, buffer_dev_1);

  // Download COMPLEX results from device.
  offloadMemcpyAsyncDtoH(zout, buffer_dev_1, buffer_size, stream);
  offloadStreamSynchronize(stream);
#else
  (void)din;
  (void)zout;
  (void)npts;
#endif
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) 1D-FFT and a (double complex
 *          to double) shrink-down on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void fft_gpu_fc(const double *zin, double *dout, const int *npts) {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
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
  fft_1d_gpu(OFFLOAD_FFT_INVERSE, npts[2], npts[0] * npts[1], false, false,
             npts[2], npts[2], buffer_dev_1, buffer_dev_2);

  // Convert COMPLEX results to REAL and download to host.
  fft_gpu_launch_complex_to_real(buffer_dev_2, buffer_dev_1, nrpts, stream);
  offloadMemcpyAsyncDtoH(dout, buffer_dev_1, buffer_size / 2, stream);
  offloadStreamSynchronize(stream);
#else
  (void)zin;
  (void)dout;
  (void)npts;
#endif
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) 1D-FFT on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void fft_gpu_f(const double *zin, double *zout, const bool dir, const int n,
               const int m, const bool transpose_in, const bool transpose_out,
               const int leading_dimension_in,
               const int leading_dimension_out) {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
  // Check inputs.
  assert(omp_get_num_threads() == 1);
  if (n == 0 || m == 0) {
    return; // Nothing to do.
  }

  // Allocate device memory.
  offload_activate_chosen_device();
  const size_t input_size =
      2 * sizeof(double) * leading_dimension_in * (transpose_in ? n : m);
  const size_t output_size =
      2 * sizeof(double) * leading_dimension_out * (transpose_out ? n : m);
  ensure_memory_sizes(imax(input_size, output_size), 0);

  // Upload COMPLEX input to device.
  offloadMemcpyAsyncHtoD(buffer_dev_1, zin, input_size, stream);

  // Run FFT on the device.
  fft_1d_gpu(dir ? OFFLOAD_FFT_FORWARD : OFFLOAD_FFT_INVERSE, n, m,
             transpose_in, transpose_out, leading_dimension_in,
             leading_dimension_out, buffer_dev_1,
             zin != zout ? buffer_dev_2 : buffer_dev_1);

  // Download COMPLEX results from device.
  offloadMemcpyAsyncDtoH(zout, zin != zout ? buffer_dev_2 : buffer_dev_1,
                         output_size, stream);
  offloadStreamSynchronize(stream);
#else
  (void)zin;
  (void)zout;
  (void)dir;
  (void)n;
  (void)m;
  (void)transpose_in;
  (void)transpose_out;
  (void)leading_dimension_in;
  (void)leading_dimension_out;
#endif
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) 1D-FFT on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void fft_r2c_gpu_f(const double *zin, double *zout, const bool dir, const int n,
                   const int m, const bool transpose_in,
                   const bool transpose_out, const int leading_dimension_in,
                   const int leading_dimension_out) {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
  // Check inputs.
  assert(omp_get_num_threads() == 1);
  if (n == 0 || m == 0) {
    return;
  }

  // Allocate device memory.
  offload_activate_chosen_device();
  const size_t input_size =
      sizeof(double) * leading_dimension_in *
      (transpose_in ? (zin != zout && dir ? n : 2 * (n / 2 + 1))
                    : (dir ? m : 2 * m));
  const size_t output_size =
      sizeof(double) * leading_dimension_out *
      (transpose_out ? (zin != zout && !dir ? n : 2 * (n / 2 + 1))
                     : (!dir ? m : 2 * m));
  ensure_memory_sizes(imax(input_size, output_size), 0);

  // Upload COMPLEX input to device.
  offloadMemcpyAsyncHtoD(buffer_dev_1, zin, input_size, stream);

  // Run FFT on the device.
  fft_r2c_1d_gpu(dir ? OFFLOAD_FFT_FORWARD : OFFLOAD_FFT_INVERSE, n, m,
                 transpose_in, transpose_out, leading_dimension_in,
                 leading_dimension_out, buffer_dev_1,
                 zin != zout ? buffer_dev_2 : buffer_dev_1);

  // Download COMPLEX results from device.
  offloadMemcpyAsyncDtoH(zout, zin != zout ? buffer_dev_2 : buffer_dev_1,
                         output_size, stream);
  offloadStreamSynchronize(stream);
#else
  (void)zin;
  (void)zout;
  (void)dir;
  (void)n;
  (void)m;
  (void)transpose_in;
  (void)transpose_out;
  (void)leading_dimension_in;
  (void)leading_dimension_out;
#endif
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) 2D-FFT on the GPU.
 * \author  Frederick Stein
 ******************************************************************************/
void fft_gpu_ff(const double *zin, double *zout, const bool dir, const int n[2],
                const int m, const bool transpose_in,
                const bool transpose_out) {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
  // Check inputs.
  assert(omp_get_num_threads() == 1);
  if (n[0] == 0 || n[1] == 0 || m == 0) {
    return;
  }

  // Allocate device memory.
  offload_activate_chosen_device();
  const size_t buffer_size = 2 * sizeof(double) * n[0] * n[1] * m;
  ensure_memory_sizes(buffer_size, 0);

  // Upload COMPLEX input to device.
  offloadMemcpyAsyncHtoD(buffer_dev_1, zin, buffer_size, stream);

  // Run FFT on the device.
  fft_2d_gpu(dir ? OFFLOAD_FFT_FORWARD : OFFLOAD_FFT_INVERSE, n, m,
             transpose_in, transpose_out, buffer_dev_1,
             zin != zout ? buffer_dev_2 : buffer_dev_1);

  // Download COMPLEX results from device.
  offloadMemcpyAsyncDtoH(zout, zin != zout ? buffer_dev_2 : buffer_dev_1,
                         buffer_size, stream);
  offloadStreamSynchronize(stream);
#else
  (void)zin;
  (void)zout;
  (void)dir;
  (void)n;
  (void)m;
  (void)transpose_in;
  (void)transpose_out;
#endif
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) R2C/C2R2D-FFT on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void fft_r2c_gpu_ff(const double *zin, double *zout, const bool dir,
                    const int n[2], const int m, const bool transpose_in,
                    const bool transpose_out) {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
  // Check inputs.
  assert(omp_get_num_threads() == 1);
  if (n[0] == 0 || n[1] == 0 || m == 0) {
    return;
  }

  // Allocate device memory.
  offload_activate_chosen_device();
  const int input_size = sizeof(double) * n[0] *
                         (zin != zout && dir ? n[1] : 2 * (n[1] / 2 + 1)) * m;
  const int output_size = sizeof(double) * n[0] *
                          (zin != zout && !dir ? n[1] : 2 * (n[1] / 2 + 1)) * m;
  ensure_memory_sizes(imax(input_size, output_size), 0);

  // Upload COMPLEX input to device.
  offloadMemcpyAsyncHtoD(buffer_dev_1, zin, input_size, stream);

  // Run FFT on the device.
  fft_r2c_2d_gpu(dir ? OFFLOAD_FFT_FORWARD : OFFLOAD_FFT_INVERSE, n, m,
                 transpose_in, transpose_out, buffer_dev_1,
                 zin != zout ? buffer_dev_2 : buffer_dev_1);

  // Download COMPLEX results from device.
  offloadMemcpyAsyncDtoH(zout, zin != zout ? buffer_dev_2 : buffer_dev_1,
                         output_size, stream);
  offloadStreamSynchronize(stream);
#else
  (void)zin;
  (void)zout;
  (void)dir;
  (void)n;
  (void)m;
  (void)transpose_in;
  (void)transpose_out;
#endif
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) 1D-FFT, followed by a (double
 *          precision complex) gather, on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void fft_gpu_fg(const double *zin, double *zout, const int *ghatmap,
                const int *npts, const int mmax, const int ngpts,
                const double scale) {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
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
  fft_1d_gpu(OFFLOAD_FFT_FORWARD, npts[0], mmax, false, false, npts[0], npts[0],
             buffer_dev_1, buffer_dev_2);

  // Upload map and run gather on the device.
  offloadMemcpyAsyncHtoD(ghatmap_dev, ghatmap, map_size, stream);
  fft_gpu_launch_gather(buffer_dev_1, buffer_dev_2, scale, ngpts, ghatmap_dev,
                        stream);

  // Download COMPLEX results from device.
  offloadMemcpyAsyncDtoH(zout, buffer_dev_1, 2 * sizeof(double) * ngpts,
                         stream);
  offloadStreamSynchronize(stream);
#else
  (void)zin;
  (void)zout;
  (void)ghatmap;
  (void)npts;
  (void)mmax;
  (void)ngpts;
  (void)scale;
#endif
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) scatter, followed by a
 *          (double precision complex) 1D-FFT, on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void fft_gpu_sf(const double *zin, double *zout, const int *ghatmap,
                const int *npts, const int mmax, const int ngpts,
                const int nmaps, const double scale) {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
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
  fft_gpu_launch_scatter(buffer_dev_2, buffer_dev_1, scale, ngpts, nmaps,
                         ghatmap_dev, stream);

  // Run FFT on the device.
  fft_1d_gpu(OFFLOAD_FFT_INVERSE, npts[0], mmax, false, false, npts[0], npts[0],
             buffer_dev_2, buffer_dev_1);

  // Download COMPLEX results from device.
  offloadMemcpyAsyncDtoH(zout, buffer_dev_1, buffer_size, stream);
  offloadStreamSynchronize(stream);
#else
  (void)zin;
  (void)zout;
  (void)ghatmap;
  (void)npts;
  (void)mmax;
  (void)ngpts;
  (void)nmaps;
  (void)scale;
#endif
}

// EOF
