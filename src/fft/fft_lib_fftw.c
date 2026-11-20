/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include "fft_lib_fftw.h"
#include "fft_key.h"
#include "fft_timer.h"
#include "fft_utils.h"

#include "../base/cp_assert.h"
#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#if defined(__FFTW3)
#include <fftw3.h>

/*******************************************************************************
 * \brief Static variables for retaining objects that are expensive to create.
 * \author Ole Schuett, Frederick Stein
 ******************************************************************************/
typedef struct {
  // The key contains
  // 0: rank (see below), R2C, direction
  // 1: associated Fortran communicator handle (to store it as an integer)
  // 2: Number of associated OpenMP threads
  // 3, 4, 5: FFT sizes (or FFT sizes and number of FFTs)
  // 6, 7, 8: input stride
  // 0, 10, 11: output strides
  fft_key_t key;
  fftw_plan *plan;
} cache_entry;

// We need to reserve more space because of the different combinations
// (local/distributed, C2C/R2C) and potential overlap from OpenMP. This works to
// run all tests.
#define FFTW_CACHE_SIZE 256

// A single FFTW plan cannot be executed concurrently by several OpenMP threads.
// Keep one independent plan cache per thread so that every worker reuses only
// its own plans. The caches are allocated once at initialization.
static cache_entry *caches = NULL;
static int *cache_oldest_entries = NULL;
static int number_of_caches = 0;

static bool is_initialized = false;

static omp_nest_lock_t cache_lock;

static int fftw_planning_mode = -1;

/*******************************************************************************
 * \brief Index of the cache that belongs to the calling thread.
 * \return cache index, clamped to a valid range
 * \author Frederick Stein
 ******************************************************************************/
static inline int current_cache_id(void) {
  int thread_id = omp_get_thread_num();
  if (thread_id < 0 || thread_id >= number_of_caches)
    thread_id = 0;
  return thread_id;
}

/*******************************************************************************
 * \brief Cache of the calling thread.
 * \return pointer to the thread-local cache entries
 * \author Frederick Stein
 ******************************************************************************/
static inline cache_entry *current_cache(void) {
  return &caches[current_cache_id() * FFTW_CACHE_SIZE];
}

/*******************************************************************************
 * \brief FIFO position of the calling thread's cache.
 * \return pointer to the thread-local FIFO index
 * \author Frederick Stein
 ******************************************************************************/
static inline int *current_cache_oldest_entry(void) {
  return &cache_oldest_entries[current_cache_id()];
}

/*******************************************************************************
 * \brief Fetches an fft plan from the cache. Returns NULL if not found.
 * \param key Encoded FFT key.
 * \author Ole Schuett, Frederick Stein
 ******************************************************************************/
static fftw_plan *lookup_plan_from_cache(const fft_key_t key) {
  cp_assert(is_initialized);
  omp_set_nest_lock(&cache_lock);
  cache_entry *cache = current_cache();
  fftw_plan *plan = NULL;
  for (int i = 0; i < FFTW_CACHE_SIZE; i++) {
    const int *x = cache[i].key;
    if (memcmp(key, x, sizeof(fft_key_t)) == 0) {
      plan = cache[i].plan;
      break;
    }
  }
  omp_unset_nest_lock(&cache_lock);
  return plan;
}

/*******************************************************************************
 * \brief Adds an fft plan to the cache. Assumes ownership of plan's memory.
 * \param key Encoded FFT key.
 * \param plan FFT plan handle.
 * \author Ole Schuett, Frederick Stein
 ******************************************************************************/
static void add_plan_to_cache(const fft_key_t key, fftw_plan *plan) {
  omp_set_nest_lock(&cache_lock);
  cache_entry *cache = current_cache();
  int *cache_oldest_entry = current_cache_oldest_entry();
  const int i = *cache_oldest_entry;
  *cache_oldest_entry = (*cache_oldest_entry + 1) % FFTW_CACHE_SIZE;
  if (cache[i].plan != NULL) {
    fftw_destroy_plan(*cache[i].plan);
    free(cache[i].plan);
  }
  memcpy(cache[i].key, key, sizeof(fft_key_t));
  cache[i].plan = plan;
  omp_unset_nest_lock(&cache_lock);
}
#endif

static bool has_guru_interface = false;

#if defined(__FFTW3)
/*******************************************************************************
 * \brief Check whether the FFTW guru interface is available.
 * \return True when a minimal guru plan can be created.
 * \author Frederick Stein
 ******************************************************************************/
bool is_guru_interface_available() {
  fftw_iodim dims[1];
  dims[0].n = 1;
  dims[0].is = 1;
  dims[0].os = 1;
  fftw_iodim howmany_dims[2];
  howmany_dims[0].n = 1;
  howmany_dims[0].is = 1;
  howmany_dims[0].os = 1;
  howmany_dims[1].n = 1;
  howmany_dims[1].is = 1;
  howmany_dims[1].os = 1;
  double complex buffer1;
  fftw_plan plan =
      fftw_plan_guru_dft(1, dims, 2, howmany_dims, &buffer1, &buffer1,
                         FFTW_FORWARD, fftw_planning_mode);

  if (plan != NULL) {
    fftw_destroy_plan(plan);
    return true;
  } else {
    return false;
  }
}
#endif

#if defined(__FFTW3)
/*******************************************************************************
 * \brief Decide whether a 3D FFT should be split into threaded 1D batches.
 * \param fft_size Size of the 3D transform.
 * \param transpose_in Whether the input layout is transposed.
 * \param transpose_out Whether the output layout is transposed.
 * \param in_place Whether the transform is in-place.
 * \return True when the threaded batched implementation should be used.
 * \author Frederick Stein
 ******************************************************************************/
static bool use_threaded_3d(const int fft_size[3], const bool transpose_in,
                            const bool transpose_out, const bool in_place) {
  return !transpose_in && !transpose_out && has_guru_interface &&
         ((((fft_size[0] >= 256 || fft_size[1] >= 256 || fft_size[2] >= 256) &&
            !in_place) ||
           ((fft_size[0] >= 128 || fft_size[1] >= 128 || fft_size[2] >= 128) &&
            in_place)) &&
          fft_team_size() > 1) &&
         (fftw_planning_mode == FFTW_ESTIMATE);
}
#endif

#if defined(__FFTW3)
fftw_plan *fft_fftw_create_3d_plan(const fft_key_t key, double complex *grid_in,
                                   double complex *grid_out);
fftw_plan *fft_fftw_create_3d_plan_r2c(const fft_key_t key, double *grid_in,
                                       double *grid_out);
fftw_plan *fft_fftw_create_guru_plan(const fft_key_t key,
                                     double complex *grid_in,
                                     double complex *grid_out);
fftw_plan *fft_fftw_create_guru_plan_r2c(const fft_key_t key, double *grid_in,
                                         double *grid_out);

/*******************************************************************************
 * \brief Register FFTW plans for a threaded batched complex 3D FFT.
 * \param direction FFT direction.
 * \param fft_size Size of the 3D transform.
 * \param number_of_ffts Number of batched transforms.
 * \param in_place Whether the transform is in-place.
 * \param grid_in Input complex grid buffer.
 * \param grid_out Output complex grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
static void
register_threaded_3d_c2c(const bool direction, const int fft_size[3],
                         const int number_of_ffts, const bool in_place,
                         double complex *grid_in, double complex *grid_out) {
  const int number_of_threads = fft_team_size();
  const int block_size =
      (number_of_ffts + number_of_threads - 1) / number_of_threads;
  const int number_of_full_blocks = number_of_ffts / block_size;
  fft_key_t key;
  get_key_3d(direction, fft_size, block_size, false, false, 1, in_place, key);
  if (lookup_plan_from_cache(key) == NULL)
    fft_fftw_create_3d_plan(key, grid_in, grid_out);
  if (block_size * number_of_full_blocks != number_of_ffts) {
    const int block_size_last_thread = number_of_ffts % block_size;
    fft_key_t key_last;
    get_key_3d(direction, fft_size, block_size_last_thread, false, false, 1,
               in_place, key_last);
    if (lookup_plan_from_cache(key_last) == NULL)
      fft_fftw_create_3d_plan(key_last, grid_in, grid_out);
  }
}

/*******************************************************************************
 * \brief Register FFTW plans for a threaded batched real-to-complex 3D FFT.
 * \param direction FFT direction.
 * \param fft_size Size of the 3D transform.
 * \param number_of_ffts Number of batched transforms.
 * \param in_place Whether the transform is in-place.
 * \param grid_in Input real grid buffer.
 * \param grid_out Output complex grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
static void register_threaded_3d_r2c(const bool direction,
                                     const int fft_size[3],
                                     const int number_of_ffts,
                                     const bool in_place, double *grid_in,
                                     double *grid_out) {
  const int number_of_threads = fft_team_size();
  const int block_size =
      (number_of_ffts + number_of_threads - 1) / number_of_threads;
  const int number_of_full_blocks = number_of_ffts / block_size;
  fft_key_t key;
  get_key_3d_r2c(direction, fft_size, block_size, false, false, 1, in_place,
                 key);
  if (lookup_plan_from_cache(key) == NULL)
    fft_fftw_create_3d_plan_r2c(key, grid_in, grid_out);
  if (block_size * number_of_full_blocks != number_of_ffts) {
    const int block_size_last_thread = number_of_ffts % block_size;
    fft_key_t key_last;
    get_key_3d_r2c(direction, fft_size, block_size_last_thread, false, false, 1,
                   in_place, key_last);
    if (lookup_plan_from_cache(key_last) == NULL)
      fft_fftw_create_3d_plan_r2c(key_last, grid_in, grid_out);
  }
}

/*******************************************************************************
 * \brief Execute a threaded batched complex-to-complex 3D FFT.
 * \param direction FFT direction.
 * \param fft_size Size of the 3D transform.
 * \param number_of_ffts Number of batched transforms.
 * \param in_place Whether the transform is in-place.
 * \param grid_in Input complex grid buffer.
 * \param grid_out Output complex grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
static void execute_threaded_3d_c2c(const bool direction, const int fft_size[3],
                                    const int number_of_ffts,
                                    const bool in_place,
                                    double complex *grid_in,
                                    double complex *grid_out) {
  const int number_of_threads = fft_team_size();
  const int block_size =
      (number_of_ffts + number_of_threads - 1) / number_of_threads;
  const int number_of_full_blocks = number_of_ffts / block_size;
  const int total_number_of_blocks =
      number_of_full_blocks + (number_of_ffts % block_size != 0);
  const int volume = fft_size[0] * fft_size[1] * fft_size[2];
#pragma omp parallel default(none) num_threads(number_of_threads) shared(      \
        block_size, number_of_full_blocks, total_number_of_blocks, direction,  \
            fft_size, grid_in, grid_out, in_place, number_of_ffts, volume)
  {
    const int thread_id = omp_get_thread_num();
    fft_key_t key;
    if (thread_id < number_of_full_blocks) {
      get_key_3d(direction, fft_size, block_size, false, false, 1, in_place,
                 key);
      fftw_plan *plan = lookup_plan_from_cache(key);
      if (plan == NULL)
        plan = fft_fftw_create_3d_plan(key, grid_in, grid_out);
      const int offset = thread_id * block_size * volume;
      fftw_execute_dft(*plan, grid_in + offset, grid_out + offset);
    } else if (thread_id + 1 == total_number_of_blocks) {
      const int block_size_last_thread = number_of_ffts % block_size;
      get_key_3d(direction, fft_size, block_size_last_thread, false, false, 1,
                 in_place, key);
      fftw_plan *plan = lookup_plan_from_cache(key);
      if (plan == NULL)
        plan = fft_fftw_create_3d_plan(key, grid_in, grid_out);
      const int offset = number_of_full_blocks * block_size * volume;
      fftw_execute_dft(*plan, grid_in + offset, grid_out + offset);
    }
  }
}

/*******************************************************************************
 * \brief Execute a threaded batched real-to-complex 3D FFT.
 * \param direction FFT direction.
 * \param fft_size Size of the 3D transform.
 * \param number_of_ffts Number of batched transforms.
 * \param in_place Whether the transform is in-place.
 * \param grid_in Input real grid buffer.
 * \param grid_out Output complex grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
static void execute_threaded_3d_r2c(const bool direction, const int fft_size[3],
                                    const int number_of_ffts,
                                    const bool in_place, double *grid_in,
                                    double *grid_out) {
  const int number_of_threads = fft_team_size();
  const int block_size =
      (number_of_ffts + number_of_threads - 1) / number_of_threads;
  const int number_of_full_blocks = number_of_ffts / block_size;
  const int total_number_of_blocks =
      number_of_full_blocks + (number_of_ffts % block_size != 0);
  const int real_volume = fft_size[0] * fft_size[1] * fft_size[2];
  const int complex_volume = fft_size[0] * fft_size[1] * (fft_size[2] / 2 + 1);
  const int real_offset = in_place ? 2 * complex_volume : real_volume;
#pragma omp parallel default(none) num_threads(number_of_threads)              \
    shared(block_size, number_of_full_blocks, total_number_of_blocks,          \
               direction, fft_size, grid_in, grid_out, in_place,               \
               number_of_ffts, complex_volume, real_offset)
  {
    const int thread_id = omp_get_thread_num();
    fft_key_t key;
    if (thread_id < number_of_full_blocks) {
      get_key_3d_r2c(direction, fft_size, block_size, false, false, 1, in_place,
                     key);
      fftw_plan *plan = lookup_plan_from_cache(key);
      if (plan == NULL)
        plan = fft_fftw_create_3d_plan_r2c(key, grid_in, grid_out);
      const int start = thread_id * block_size;
      if (direction) {
        fftw_execute_dft_r2c(*plan, grid_in + start * real_offset,
                             (double complex *)grid_out +
                                 start * complex_volume);
      } else {
        fftw_execute_dft_c2r(*plan,
                             (double complex *)grid_in + start * complex_volume,
                             grid_out + start * real_offset);
      }
    } else if (thread_id + 1 == total_number_of_blocks) {
      const int block_size_last_thread = number_of_ffts % block_size;
      get_key_3d_r2c(direction, fft_size, block_size_last_thread, false, false,
                     1, in_place, key);
      fftw_plan *plan = lookup_plan_from_cache(key);
      if (plan == NULL)
        plan = fft_fftw_create_3d_plan_r2c(key, grid_in, grid_out);
      const int start = number_of_full_blocks * block_size;
      if (direction) {
        fftw_execute_dft_r2c(*plan, grid_in + start * real_offset,
                             (double complex *)grid_out +
                                 start * complex_volume);
      } else {
        fftw_execute_dft_c2r(*plan,
                             (double complex *)grid_in + start * complex_volume,
                             grid_out + start * real_offset);
      }
    }
  }
}

/*******************************************************************************
 * \brief Register the thread-split 1D chain for a single real/complex 3D FFT.
 ******************************************************************************/
static void register_threaded_3d_r2c_1d(const bool direction,
                                        const int fft_size[3],
                                        const bool in_place, double *grid_real,
                                        double complex *grid_complex) {
  double complex *buffer = in_place
                               ? grid_complex
                               : fftw_alloc_complex(fft_size[0] * fft_size[1] *
                                                    (fft_size[2] / 2 + 1));
  double complex *input_c = direction ? buffer : grid_complex;
  double complex *output_c = direction ? grid_complex : buffer;
  const int number_of_threads = fft_team_size();
  int block_sizes[3];
  fftw_plan *plans[3], *plans_last_thread[3];
  for (int stage = 0; stage < 3; stage++) {
    plans[stage] = NULL;
    plans_last_thread[stage] = NULL;
    const int extent = stage == 0   ? fft_size[1] * (fft_size[2] / 2 + 1)
                       : stage == 1 ? fft_size[0]
                                    : fft_size[0] * fft_size[1];
    block_sizes[stage] = (extent + number_of_threads - 1) / number_of_threads;
    const bool needs_tail = block_sizes[stage] * number_of_threads != extent;
    const int tail = needs_tail ? extent % block_sizes[stage] : 0;
    if (stage == 0) {
      const int n = fft_size[1] * (fft_size[2] / 2 + 1);
      fft_iodim dim = {.n = fft_size[0], .is = n, .os = n};
      fft_iodim howmany = {.n = block_sizes[0], .is = 1, .os = 1};
      fft_key_t key;
      get_key_guru(direction, 1, &dim, 1, &howmany, 1, in_place, key);
      plans[0] = lookup_plan_from_cache(key);
      if (plans[0] == NULL)
        plans[0] = fft_fftw_create_guru_plan(key, input_c, output_c);
      if (needs_tail) {
        howmany.n = tail;
        get_key_guru(direction, 1, &dim, 1, &howmany, 1, in_place, key);
        plans_last_thread[0] = lookup_plan_from_cache(key);
        if (plans_last_thread[0] == NULL)
          plans_last_thread[0] =
              fft_fftw_create_guru_plan(key, input_c, output_c);
      }
    } else if (stage == 1) {
      const int half = fft_size[2] / 2 + 1;
      fft_iodim dim = {.n = fft_size[1], .is = half, .os = half};
      fft_iodim howmany[2] = {{.n = block_sizes[1],
                               .is = fft_size[1] * half,
                               .os = fft_size[1] * half},
                              {.n = half, .is = 1, .os = 1}};
      fft_key_t key;
      get_key_guru(direction, 1, &dim, 2, howmany, 1, in_place, key);
      plans[1] = lookup_plan_from_cache(key);
      if (plans[1] == NULL)
        plans[1] = fft_fftw_create_guru_plan(key, input_c, output_c);
      if (needs_tail) {
        howmany[0].n = tail;
        get_key_guru(direction, 1, &dim, 2, howmany, 1, in_place, key);
        plans_last_thread[1] = lookup_plan_from_cache(key);
        if (plans_last_thread[1] == NULL)
          plans_last_thread[1] =
              fft_fftw_create_guru_plan(key, input_c, output_c);
      }
    } else {
      const int half = fft_size[2] / 2 + 1;
      fft_iodim dim = {.n = fft_size[2], .is = 1, .os = 1};
      fft_iodim howmany = {
          .n = block_sizes[2],
          .is = direction ? (in_place ? 2 * half : fft_size[2]) : half,
          .os = direction ? half : (in_place ? 2 * half : fft_size[2])};
      fft_key_t key;
      get_key_guru_r2c(direction, 1, &dim, 1, &howmany, 1, in_place, key);
      plans[2] = lookup_plan_from_cache(key);
      if (plans[2] == NULL)
        plans[2] = fft_fftw_create_guru_plan_r2c(key, grid_real,
                                                 (double *)grid_complex);
      if (needs_tail) {
        howmany.n = tail;
        get_key_guru_r2c(direction, 1, &dim, 1, &howmany, 1, in_place, key);
        plans_last_thread[2] = lookup_plan_from_cache(key);
        if (plans_last_thread[2] == NULL)
          plans_last_thread[2] = fft_fftw_create_guru_plan_r2c(
              key, grid_real, (double *)grid_complex);
      }
    }
  }
  if (!in_place)
    fftw_free(buffer);
}

/*******************************************************************************
 * \brief Execute the thread-split 1D chain for a single real/complex 3D FFT.
 ******************************************************************************/
static void execute_threaded_3d_r2c_1d(const bool direction,
                                       const int fft_size[3],
                                       const bool in_place, double *grid_real,
                                       double complex *grid_complex) {
  const int half = fft_size[2] / 2 + 1;
  double complex *buffer =
      direction && in_place
          ? grid_complex
          : fftw_alloc_complex(fft_size[0] * fft_size[1] * half);
  double complex *input_c = direction ? buffer : grid_complex;
  double complex *output_c = direction ? grid_complex : buffer;
  const int number_of_threads = fft_team_size();
  int block_sizes[3];
  int number_of_full_blocks[3] = {0, 0, 0};
  int total_number_of_blocks[3] = {0, 0, 0};
  fftw_plan *plans[3], *plans_last_thread[3];
  for (int stage = 0; stage < 3; stage++) {
    plans[stage] = NULL;
    plans_last_thread[stage] = NULL;
    const int extent = stage == 0   ? fft_size[1] * half
                       : stage == 1 ? fft_size[0]
                                    : fft_size[0] * fft_size[1];
    block_sizes[stage] = (extent + number_of_threads - 1) / number_of_threads;
    number_of_full_blocks[stage] = extent / block_sizes[stage];
    total_number_of_blocks[stage] =
        (extent + block_sizes[stage] - 1) / block_sizes[stage];
    const bool needs_tail = block_sizes[stage] * number_of_threads != extent;
    const int tail = needs_tail ? extent % block_sizes[stage] : 0;
    if (stage == 0) {
      const int n = fft_size[1] * half;
      fft_iodim dim = {.n = fft_size[0], .is = n, .os = n};
      fft_iodim howmany = {.n = block_sizes[0], .is = 1, .os = 1};
      fft_key_t key;
      get_key_guru(direction, 1, &dim, 1, &howmany, 1, in_place, key);
      plans[0] = lookup_plan_from_cache(key);
      if (plans[0] == NULL) {
        double complex *tmp = fftw_alloc_complex(get_buffer_size_from_key(key));
        plans[0] = fft_fftw_create_guru_plan(key, tmp, output_c);
        fftw_free(tmp);
      }
      if (needs_tail) {
        howmany.n = tail;
        get_key_guru(direction, 1, &dim, 1, &howmany, 1, in_place, key);
        plans_last_thread[0] = lookup_plan_from_cache(key);
        if (plans_last_thread[0] == NULL) {
          double complex *tmp =
              fftw_alloc_complex(get_buffer_size_from_key(key));
          plans_last_thread[0] = fft_fftw_create_guru_plan(key, tmp, output_c);
          fftw_free(tmp);
        }
      }
    } else if (stage == 1) {
      fft_iodim dim = {.n = fft_size[1], .is = half, .os = half};
      fft_iodim howmany[2] = {{.n = block_sizes[1],
                               .is = fft_size[1] * half,
                               .os = fft_size[1] * half},
                              {.n = half, .is = 1, .os = 1}};
      fft_key_t key;
      get_key_guru(direction, 1, &dim, 2, howmany, 1, in_place, key);
      plans[1] = lookup_plan_from_cache(key);
      if (plans[1] == NULL) {
        double complex *tmp = fftw_alloc_complex(get_buffer_size_from_key(key));
        plans[1] = fft_fftw_create_guru_plan(key, tmp, output_c);
        fftw_free(tmp);
      }
      if (needs_tail) {
        howmany[0].n = tail;
        get_key_guru(direction, 1, &dim, 2, howmany, 1, in_place, key);
        plans_last_thread[1] = lookup_plan_from_cache(key);
        if (plans_last_thread[1] == NULL) {
          double complex *tmp =
              fftw_alloc_complex(get_buffer_size_from_key(key));
          plans_last_thread[1] = fft_fftw_create_guru_plan(key, tmp, output_c);
          fftw_free(tmp);
        }
      }
    } else {
      fft_iodim dim = {.n = fft_size[2], .is = 1, .os = 1};
      fft_iodim howmany = {
          .n = block_sizes[2],
          .is = direction ? (in_place ? 2 * half : fft_size[2]) : half,
          .os = direction ? half : (in_place ? 2 * half : fft_size[2])};
      fft_key_t key;
      get_key_guru_r2c(direction, 1, &dim, 1, &howmany, 1, in_place, key);
      plans[2] = lookup_plan_from_cache(key);
      if (plans[2] == NULL) {
        double *tmp = fftw_alloc_real(2 * get_buffer_size_from_key(key));
        plans[2] = fft_fftw_create_guru_plan_r2c(key, tmp, (double *)output_c);
        fftw_free(tmp);
      }
      if (needs_tail) {
        howmany.n = tail;
        get_key_guru_r2c(direction, 1, &dim, 1, &howmany, 1, in_place, key);
        plans_last_thread[2] = lookup_plan_from_cache(key);
        if (plans_last_thread[2] == NULL) {
          double *tmp = fftw_alloc_real(2 * get_buffer_size_from_key(key));
          plans_last_thread[2] =
              fft_fftw_create_guru_plan_r2c(key, tmp, (double *)output_c);
          fftw_free(tmp);
        }
      }
    }
  }
  const int real_stride =
      direction ? (in_place ? 2 * half : fft_size[2]) : half;
  const int complex_stride =
      direction ? half : (in_place ? 2 * half : fft_size[2]);
#pragma omp parallel default(none)                                             \
    shared(number_of_threads, number_of_full_blocks, total_number_of_blocks,   \
               plans, plans_last_thread, block_sizes, fft_size, grid_real,     \
               grid_complex, in_place, direction, input_c, output_c,           \
               real_stride, complex_stride)
  {
    const int thread_id = omp_get_thread_num();
    if (thread_id < number_of_full_blocks[0]) {
      fftw_execute_dft(*plans[0], input_c + block_sizes[0] * thread_id,
                       output_c + block_sizes[0] * thread_id);
    } else if (thread_id == total_number_of_blocks[0] - 1) {
      fftw_execute_dft(*plans_last_thread[0],
                       input_c + block_sizes[0] * thread_id,
                       output_c + block_sizes[0] * thread_id);
    }
#pragma omp barrier
    if (thread_id < number_of_full_blocks[1]) {
      const int offset =
          block_sizes[1] * fft_size[1] * (fft_size[2] / 2 + 1) * thread_id;
      fftw_execute_dft(*plans[1], output_c + offset, input_c + offset);
    } else if (thread_id == total_number_of_blocks[1] - 1) {
      const int offset =
          block_sizes[1] * fft_size[1] * (fft_size[2] / 2 + 1) * thread_id;
      fftw_execute_dft(*plans_last_thread[1], output_c + offset,
                       input_c + offset);
    }
#pragma omp barrier
    if (thread_id < number_of_full_blocks[2]) {
      const int start = block_sizes[2] * thread_id;
      if (direction) {
        fftw_execute_dft_r2c(*plans[2], grid_real + start * real_stride,
                             grid_complex + start * complex_stride);
      } else {
        fftw_execute_dft_c2r(*plans[2], grid_complex + start * complex_stride,
                             grid_real + start * real_stride);
      }
    } else if (thread_id == total_number_of_blocks[2] - 1) {
      const int start = block_sizes[2] * thread_id;
      if (direction) {
        fftw_execute_dft_r2c(*plans_last_thread[2],
                             grid_real + start * real_stride,
                             grid_complex + start * complex_stride);
      } else {
        fftw_execute_dft_c2r(*plans_last_thread[2],
                             grid_complex + start * complex_stride,
                             grid_real + start * real_stride);
      }
    }
  }
  if (!(direction && in_place))
    fftw_free(buffer);
}

/*******************************************************************************
 * \brief Acquire the plans for a thread-split single C2C 3D FFT.
 ******************************************************************************/
static void acquire_threaded_3d_c2c_1d(
    const bool direction, const int fft_size[3], const bool in_place,
    const bool for_execution, double complex *grid_in, double complex *grid_out,
    int block_sizes[3], fftw_plan *plans[3], fftw_plan *plans_last_thread[3],
    int number_of_full_blocks[3], int total_number_of_blocks[3]) {
  const int number_of_threads = fft_team_size();
  for (int stage = 0; stage < 3; stage++) {
    plans[stage] = NULL;
    plans_last_thread[stage] = NULL;
    const int extent = stage == 0   ? fft_size[1] * fft_size[2]
                       : stage == 1 ? fft_size[0]
                                    : fft_size[0] * fft_size[1];
    block_sizes[stage] = (extent + number_of_threads - 1) / number_of_threads;
    number_of_full_blocks[stage] = extent / block_sizes[stage];
    total_number_of_blocks[stage] =
        (extent + block_sizes[stage] - 1) / block_sizes[stage];
    const bool needs_tail = block_sizes[stage] * number_of_threads != extent;
    const int tail = needs_tail ? extent % block_sizes[stage] : 0;

    fft_iodim dim[1];
    fft_iodim howmany[2];
    int rank = 1;
    int howmany_rank = 1;
    if (stage == 0) {
      const int n = fft_size[1] * fft_size[2];
      dim[0] = (fft_iodim){.n = fft_size[0], .is = n, .os = n};
      howmany[0] = (fft_iodim){.n = block_sizes[0], .is = 1, .os = 1};
    } else if (stage == 1) {
      dim[0] =
          (fft_iodim){.n = fft_size[1], .is = fft_size[2], .os = fft_size[2]};
      howmany[0] = (fft_iodim){.n = block_sizes[1],
                               .is = fft_size[1] * fft_size[2],
                               .os = fft_size[1] * fft_size[2]};
      howmany[1] = (fft_iodim){.n = fft_size[2], .is = 1, .os = 1};
      howmany_rank = 2;
    } else {
      dim[0] = (fft_iodim){.n = fft_size[2], .is = 1, .os = 1};
      howmany[0] = (fft_iodim){
          .n = block_sizes[2], .is = fft_size[2], .os = fft_size[2]};
    }

    fft_key_t key;
    get_key_guru(direction, rank, dim, howmany_rank, howmany, 1, in_place, key);
    plans[stage] = lookup_plan_from_cache(key);
    if (plans[stage] == NULL) {
      if (for_execution) {
        double complex *buffer =
            fftw_alloc_complex(get_buffer_size_from_key(key));
        plans[stage] = fft_fftw_create_guru_plan(key, buffer,
                                                 in_place ? buffer : grid_out);
        fftw_free(buffer);
      } else {
        plans[stage] = fft_fftw_create_guru_plan(key, grid_in, grid_out);
      }
    }
    if (needs_tail) {
      howmany[0].n = tail;
      get_key_guru(direction, rank, dim, howmany_rank, howmany, 1, in_place,
                   key);
      plans_last_thread[stage] = lookup_plan_from_cache(key);
      if (plans_last_thread[stage] == NULL) {
        if (for_execution) {
          double complex *buffer =
              fftw_alloc_complex(get_buffer_size_from_key(key));
          plans_last_thread[stage] = fft_fftw_create_guru_plan(
              key, buffer, in_place ? buffer : grid_out);
          fftw_free(buffer);
        } else {
          plans_last_thread[stage] =
              fft_fftw_create_guru_plan(key, grid_in, grid_out);
        }
      }
    }
  }
}

/*******************************************************************************
 * \brief Register the thread-split 1D chain for a single C2C 3D FFT.
 ******************************************************************************/
static void register_threaded_3d_c2c_1d(const bool direction,
                                        const int fft_size[3],
                                        const bool in_place,
                                        double complex *grid_in,
                                        double complex *grid_out) {
  int block_sizes[3];
  int number_of_full_blocks[3];
  int total_number_of_blocks[3];
  fftw_plan *plans[3], *plans_last_thread[3];
  acquire_threaded_3d_c2c_1d(direction, fft_size, in_place, false, grid_in,
                             grid_out, block_sizes, plans, plans_last_thread,
                             number_of_full_blocks, total_number_of_blocks);
}

/*******************************************************************************
 * \brief Execute the thread-split 1D chain for a single C2C 3D FFT.
 ******************************************************************************/
static void execute_threaded_3d_c2c_1d(const bool direction,
                                       const int fft_size[3],
                                       const bool in_place,
                                       double complex *grid_in,
                                       double complex *grid_out) {
  int block_sizes[3];
  int number_of_full_blocks[3];
  int total_number_of_blocks[3];
  fftw_plan *plans[3], *plans_last_thread[3];
  acquire_threaded_3d_c2c_1d(direction, fft_size, in_place, true, grid_in,
                             grid_out, block_sizes, plans, plans_last_thread,
                             number_of_full_blocks, total_number_of_blocks);
#pragma omp parallel default(none)                                             \
    shared(number_of_full_blocks, total_number_of_blocks, plans,               \
               plans_last_thread, block_sizes, fft_size, grid_in, grid_out)
  {
    const int thread_id = omp_get_thread_num();
    if (thread_id < number_of_full_blocks[2]) {
      fftw_execute_dft(*plans[2],
                       grid_in + block_sizes[2] * fft_size[2] * thread_id,
                       grid_out + block_sizes[2] * fft_size[2] * thread_id);
    } else if (thread_id == total_number_of_blocks[2] - 1) {
      fftw_execute_dft(*plans_last_thread[2],
                       grid_in + block_sizes[2] * fft_size[2] * thread_id,
                       grid_out + block_sizes[2] * fft_size[2] * thread_id);
    }
#pragma omp barrier
    if (thread_id < number_of_full_blocks[1]) {
      const int offset = block_sizes[1] * fft_size[1] * fft_size[2] * thread_id;
      fftw_execute_dft(*plans[1], grid_out + offset, grid_in + offset);
    } else if (thread_id == total_number_of_blocks[1] - 1) {
      const int offset = block_sizes[1] * fft_size[1] * fft_size[2] * thread_id;
      fftw_execute_dft(*plans_last_thread[1], grid_out + offset,
                       grid_in + offset);
    }
#pragma omp barrier
    if (thread_id < number_of_full_blocks[0]) {
      fftw_execute_dft(*plans[0], grid_in + block_sizes[0] * thread_id,
                       grid_out + block_sizes[0] * thread_id);
    } else if (thread_id == total_number_of_blocks[0] - 1) {
      fftw_execute_dft(*plans_last_thread[0],
                       grid_in + block_sizes[0] * thread_id,
                       grid_out + block_sizes[0] * thread_id);
    }
  }
}

#endif

/*******************************************************************************
 * \\brief Initialize the FFT library (if not done externally).
 * \param fftw_planning_flag FFTW planning mode.
 * \param use_guru_interface Whether the FFTW Guru interface may be used.
 * \param wisdom_file Path to the FFTW wisdom file.
 * \author Frederick Stein
 ******************************************************************************/
int fft_fftw_init_lib(const fftw_plan_type fftw_planning_flag,
                      const bool use_guru_interface, const char *wisdom_file) {
#if defined(__FFTW3)
  cp_assert(omp_get_num_threads() == 1);
  if (is_initialized) {
    return 0;
  }
  number_of_caches =
      imax(1, omp_get_max_threads()); // NOLINT: per-thread caches
  caches =
      calloc((size_t)number_of_caches * FFTW_CACHE_SIZE, sizeof(cache_entry));
  cache_oldest_entries = calloc((size_t)number_of_caches, sizeof(int));
  cp_assert(caches != NULL && cache_oldest_entries != NULL);
  omp_init_nest_lock(&cache_lock);

  fftw_planning_mode = fftw_planning_flag;
  switch (fftw_planning_flag) {
  case FFT_ESTIMATE:
    fftw_planning_mode = FFTW_ESTIMATE;
    break;
  case FFT_MEASURE:
    fftw_planning_mode = FFTW_MEASURE;
    break;
  case FFT_PATIENT:
    fftw_planning_mode = FFTW_PATIENT;
    break;
  case FFT_EXHAUSTIVE:
    fftw_planning_mode = FFTW_EXHAUSTIVE;
    break;
  default:
    cp_assert(0 && "Unknown FFTW planning flag.");
  }

  // We need a threaded library!
  const int status = fftw_init_threads();
  if (status == 0)
    return 1;

  has_guru_interface = use_guru_interface;
  if (has_guru_interface) {
    has_guru_interface = is_guru_interface_available();
  }

  // Export wisdom after intializing the library to ensure correct threading
  // etc.
  if (wisdom_file != NULL) {
    const bool is_print_rank = cp_mpi_comm_rank(cp_mpi_get_comm_world()) == 0;
    const int error = fftw_import_wisdom_from_filename(wisdom_file);
    if (false && error != 0 && is_print_rank)
      fprintf(stderr,
              "Importing wisdom failed! Maybe the file does not exist.");
  }

  is_initialized = true;
  return 0;
#else
  cp_assert(false && "The grid library was not compiled with FFTW support.");
  (void)fftw_planning_flag;
  (void)use_guru_interface;
  (void)wisdom_file;
  return 0;
#endif
}

/*******************************************************************************
 * \brief Finalize the FFT library (if not done externally).
 * \param wisdom_file Path to the FFTW wisdom file.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_finalize_lib(const char *wisdom_file) {
#if defined(__FFTW3)
  cp_assert(omp_get_num_threads() == 1);
  if (!is_initialized) {
    return;
  }
  for (int t = 0; t < number_of_caches; t++) {
    cache_entry *cache = &caches[t * FFTW_CACHE_SIZE];
    for (int i = 0; i < FFTW_CACHE_SIZE; i++) {
      if (cache[i].plan != NULL) {
        fftw_destroy_plan(*cache[i].plan);
        free(cache[i].plan);
      }
    }
  }
  free(caches);
  caches = NULL;
  free(cache_oldest_entries);
  cache_oldest_entries = NULL;
  number_of_caches = 0;
  omp_destroy_nest_lock(&cache_lock);
  // Export wisdom before finalizing the library to ensure storing the correct
  // threading etc.
  if (wisdom_file != NULL) {
    const int error = fftw_export_wisdom_to_filename(wisdom_file);
    if (false && error != 0 && cp_mpi_comm_rank(cp_mpi_get_comm_world()))
      fprintf(stderr,
              "Exporting wisdom failed! Maybe writing access is missing.");
  }
  is_initialized = false;
  fftw_planning_mode = -1;
  fftw_cleanup_threads();
#else
  cp_assert(false && "The grid library was not compiled with FFTW support.");
  (void)wisdom_file;
#endif
}

/*******************************************************************************
 * \brief Whether a distributed FFT implementation is available.
 * \author Frederick Stein
 ******************************************************************************/
bool fft_fftw_lib_has_guru_interface() { return has_guru_interface; }

/*******************************************************************************
 * \brief Allocate buffer of type double.
 * \param length Number of elements.
 * \param buffer Buffer pointer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_allocate_double(const int length, double **buffer) {
#if defined(__FFTW3)
  cp_assert(is_initialized);
  cp_assert(buffer != NULL);
  cp_assert(*buffer == NULL);
  *buffer = fftw_alloc_real(imax(1, length));
#else
  (void)length;
  (void)buffer;
  cp_assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Allocate buffer of type double complex.
 * \param length Number of elements.
 * \param buffer Buffer pointer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_allocate_complex(const int length, double complex **buffer) {
#if defined(__FFTW3)
  cp_assert(is_initialized);
  cp_assert(buffer != NULL);
  cp_assert(*buffer == NULL);
  *buffer = fftw_alloc_complex(imax(1, length));
#else
  (void)length;
  (void)buffer;
  cp_assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Free buffer of type double.
 * \param buffer Buffer pointer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_free_double(double *buffer) {
#if defined(__FFTW3)
  cp_assert(is_initialized);
  fftw_free(buffer);
#else
  (void)buffer;
  cp_assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Free buffer of type double complex.
 * \param buffer Buffer pointer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_free_complex(double complex *buffer) {
#if defined(__FFTW3)
  cp_assert(is_initialized);
  fftw_free(buffer);
#else
  (void)buffer;
  cp_assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

#if defined(__FFTW3)

/*******************************************************************************
 * \brief Create plan of a local C2C 1D FFT.
 * \param key Encoded FFT key.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
fftw_plan *fft_fftw_create_1d_plan(const fft_key_t key, double complex *grid_in,
                                   double complex *grid_out) {
  omp_set_nest_lock(&cache_lock);
  fftw_plan *cached_plan = lookup_plan_from_cache(key);
  if (cached_plan != NULL) {
    omp_unset_nest_lock(&cache_lock);
    return cached_plan;
  }
  int rank, number_of_threads, number_of_ffts, istride, ostride, idist, odist;
  bool direction, inplace;
  int fft_size[3], inembed[3], onembed[3];
  fetch_data_from_key_nd(key, &rank, fft_size, &number_of_ffts,
                         &number_of_threads, &direction, &inplace, inembed,
                         onembed, &idist, &odist, &istride, &ostride);
  cp_assert(rank == 1);
  fftw_plan_with_nthreads(number_of_threads);
  fftw_plan *plan = malloc(sizeof(fftw_plan));
  if (grid_in == grid_out) {
    cp_assert(inplace);
    cp_assert(istride == ostride);
    cp_assert(idist == odist);
  } else {
    cp_assert(!inplace);
  }
  *plan = fftw_plan_many_dft(rank, fft_size, number_of_ffts, grid_in, inembed,
                             istride, idist, grid_out, onembed, ostride, odist,
                             direction ? FFTW_FORWARD : FFTW_BACKWARD,
                             fftw_planning_mode);
  cp_assert(plan != NULL);
  add_plan_to_cache(key, plan);
  omp_unset_nest_lock(&cache_lock);
  return plan;
}
/*******************************************************************************
 * \brief Create plan of a local R2C/C2R 1D FFT.
 * \param key Encoded FFT key.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
fftw_plan *fft_fftw_create_1d_plan_r2c(const fft_key_t key, double *grid_in,
                                       double *grid_out) {
  omp_set_nest_lock(&cache_lock);
  fftw_plan *cached_plan = lookup_plan_from_cache(key);
  if (cached_plan != NULL) {
    omp_unset_nest_lock(&cache_lock);
    return cached_plan;
  }
  int rank, number_of_threads, number_of_ffts, istride, ostride, idist, odist;
  bool direction, inplace;
  int fft_size[3], inembed[3], onembed[3];
  fetch_data_from_key_nd(key, &rank, fft_size, &number_of_ffts,
                         &number_of_threads, &direction, &inplace, inembed,
                         onembed, &idist, &odist, &istride, &ostride);
  cp_assert(rank == 1);
  fftw_plan_with_nthreads(number_of_threads);
  fftw_plan *plan = malloc(sizeof(fftw_plan));
  if (grid_in == grid_out) {
    if (istride == 1) {
      // Check whether the arrays are correctly padded
      if (direction) {
        cp_assert(idist == 2 * odist);
      } else {
        cp_assert(odist == 2 * idist);
      }
    }
    cp_assert(istride == ostride);
  }
  if (direction) {
    *plan =
        fftw_plan_many_dft_r2c(rank, fft_size, number_of_ffts, grid_in, inembed,
                               istride, idist, (double complex *)grid_out,
                               onembed, ostride, odist, fftw_planning_mode);
  } else {
    // We use the buffers the other way around to prevent
    // out-of-bounds-accesses of the planner if the output array has only the
    // minimum size
    *plan = fftw_plan_many_dft_c2r(
        rank, fft_size, number_of_ffts, (double complex *)grid_in, inembed,
        istride, idist, grid_out, onembed, ostride, odist, fftw_planning_mode);
  }
  cp_assert(plan != NULL);
  add_plan_to_cache(key, plan);
  omp_unset_nest_lock(&cache_lock);
  return plan;
}

/*******************************************************************************
 * \brief Create plan of a local C2C 2D FFT.
 * \param key Encoded FFT key.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
fftw_plan *fft_fftw_create_2d_plan(const fft_key_t key, double complex *grid_in,
                                   double complex *grid_out) {
  omp_set_nest_lock(&cache_lock);
  fftw_plan *cached_plan = lookup_plan_from_cache(key);
  if (cached_plan != NULL) {
    omp_unset_nest_lock(&cache_lock);
    return cached_plan;
  }
  int rank, number_of_threads, number_of_ffts, istride, ostride, idist, odist;
  bool direction, inplace;
  int fft_size[3], inembed[3], onembed[3];
  fetch_data_from_key_nd(key, &rank, fft_size, &number_of_ffts,
                         &number_of_threads, &direction, &inplace, inembed,
                         onembed, &idist, &odist, &istride, &ostride);
  cp_assert(rank == 2);
  fftw_plan_with_nthreads(number_of_threads);
  fftw_plan *plan = malloc(sizeof(fftw_plan));
  if (grid_in == grid_out) {
    // Array strides need to match for in-place FFTs
    cp_assert(istride == ostride);
    cp_assert(idist == odist);
  }
  *plan = fftw_plan_many_dft(rank, fft_size, number_of_ffts, grid_in, inembed,
                             istride, idist, grid_out, onembed, ostride, odist,
                             direction ? FFTW_FORWARD : FFTW_BACKWARD,
                             fftw_planning_mode);
  cp_assert(plan != NULL);
  add_plan_to_cache(key, plan);
  omp_unset_nest_lock(&cache_lock);
  return plan;
}

/*******************************************************************************
 * \brief Create plan of a local R2C/C2R 2D FFT.
 * \param key Encoded FFT key.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
fftw_plan *fft_fftw_create_2d_plan_r2c(const fft_key_t key, double *grid_in,
                                       double *grid_out) {
  omp_set_nest_lock(&cache_lock);
  fftw_plan *cached_plan = lookup_plan_from_cache(key);
  if (cached_plan != NULL) {
    omp_unset_nest_lock(&cache_lock);
    return cached_plan;
  }
  int rank, number_of_threads, number_of_ffts, istride, ostride, idist, odist;
  bool direction, inplace;
  int fft_size[3], inembed[3], onembed[3];
  fetch_data_from_key_nd(key, &rank, fft_size, &number_of_ffts,
                         &number_of_threads, &direction, &inplace, inembed,
                         onembed, &idist, &odist, &istride, &ostride);
  cp_assert(rank == 2);
  fftw_plan_with_nthreads(number_of_threads);
  // We need the guru interface here because cuts the last dimension in half
  // whereas we want the first dimension
  fftw_plan *plan = malloc(sizeof(fftw_plan));
  if (grid_in == grid_out) {
    if (istride == 1) {
      // Check whether the arrays are correctly padded
      if (direction) {
        cp_assert(idist == 2 * odist);
      } else {
        cp_assert(odist == 2 * idist);
      }
    }
    cp_assert(istride == ostride);
  }
  if (direction) {
    *plan =
        fftw_plan_many_dft_r2c(rank, fft_size, number_of_ffts, grid_in, inembed,
                               istride, idist, (double complex *)grid_out,
                               onembed, ostride, odist, fftw_planning_mode);
  } else {
    // We use the buffers the other way around to prevent
    // out-of-bounds-accesses of the planner if the output array has only the
    // minimum size
    *plan = fftw_plan_many_dft_c2r(
        rank, fft_size, number_of_ffts, (double complex *)grid_in, inembed,
        istride, idist, grid_out, onembed, ostride, odist, fftw_planning_mode);
  }
  cp_assert(plan != NULL);
  add_plan_to_cache(key, plan);
  omp_unset_nest_lock(&cache_lock);
  return plan;
}

/*******************************************************************************
 * \brief Create plan of a local C2C 3D FFT.
 * \param key Encoded FFT key.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
fftw_plan *fft_fftw_create_3d_plan(const fft_key_t key, double complex *grid_in,
                                   double complex *grid_out) {
  omp_set_nest_lock(&cache_lock);
  fftw_plan *cached_plan = lookup_plan_from_cache(key);
  if (cached_plan != NULL) {
    omp_unset_nest_lock(&cache_lock);
    return cached_plan;
  }
  int rank, number_of_threads, number_of_ffts, istride, ostride, idist, odist;
  bool direction, inplace;
  int fft_size[3], inembed[3], onembed[3];
  fetch_data_from_key_nd(key, &rank, fft_size, &number_of_ffts,
                         &number_of_threads, &direction, &inplace, inembed,
                         onembed, &idist, &odist, &istride, &ostride);
  cp_assert(rank == 3);
  fftw_plan_with_nthreads(number_of_threads);
  fftw_plan *plan = malloc(sizeof(fftw_plan));
  *plan = fftw_plan_many_dft(rank, fft_size, number_of_ffts, grid_in, inembed,
                             istride, idist, grid_out, onembed, ostride, odist,
                             direction ? FFTW_FORWARD : FFTW_BACKWARD,
                             fftw_planning_mode);
  add_plan_to_cache(key, plan);
  cp_assert(plan != NULL);
  omp_unset_nest_lock(&cache_lock);
  return plan;
}

/*******************************************************************************
 * \brief Create plan of a local R2C/C2R 3D FFT.
 * \param key Encoded FFT key.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
fftw_plan *fft_fftw_create_3d_plan_r2c(const fft_key_t key, double *grid_in,
                                       double *grid_out) {
  omp_set_nest_lock(&cache_lock);
  fftw_plan *cached_plan = lookup_plan_from_cache(key);
  if (cached_plan != NULL) {
    omp_unset_nest_lock(&cache_lock);
    return cached_plan;
  }
  int rank, number_of_threads, number_of_ffts, istride, ostride, idist, odist;
  bool direction, inplace;
  int fft_size[3], inembed[3], onembed[3];
  fetch_data_from_key_nd(key, &rank, fft_size, &number_of_ffts,
                         &number_of_threads, &direction, &inplace, inembed,
                         onembed, &idist, &odist, &istride, &ostride);
  cp_assert(rank == 3);
  fftw_plan_with_nthreads(number_of_threads);
  fftw_plan *plan = malloc(sizeof(fftw_plan));
  if (direction) {
    *plan =
        fftw_plan_many_dft_r2c(rank, fft_size, number_of_ffts, grid_in, inembed,
                               istride, idist, (double complex *)grid_out,
                               onembed, ostride, odist, fftw_planning_mode);
  } else {
    // We use the buffers the other way around to prevent
    // out-of-bounds-accesses of the planner if the output array has only the
    // minimum size
    *plan = fftw_plan_many_dft_c2r(
        rank, fft_size, number_of_ffts, (double complex *)grid_in, inembed,
        istride, idist, grid_out, onembed, ostride, odist, fftw_planning_mode);
  }
  add_plan_to_cache(key, plan);
  cp_assert(plan != NULL);
  omp_unset_nest_lock(&cache_lock);
  return plan;
}

/*******************************************************************************
 * \brief Create plan of a local C2C 3D FFT.
 * \param key Encoded FFT key.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
fftw_plan *fft_fftw_create_guru_plan(const fft_key_t key,
                                     double complex *grid_in,
                                     double complex *grid_out) {
  cp_assert(has_guru_interface);
  omp_set_nest_lock(&cache_lock);
  fftw_plan *cached_plan = lookup_plan_from_cache(key);
  if (cached_plan != NULL) {
    omp_unset_nest_lock(&cache_lock);
    return cached_plan;
  }
  bool direction, inplace;
  int rank, howmany_rank, number_of_threads;
  fft_iodim dims[3], howmany_dims[3];
  fetch_data_from_key_guru(key, &direction, &rank, dims, &howmany_rank,
                           howmany_dims, &number_of_threads, &inplace);

  fftw_plan_with_nthreads(number_of_threads);
  fftw_plan *plan = malloc(sizeof(fftw_plan));
  if ((double complex *)grid_in == grid_out) {
    cp_assert(inplace);
    for (int r = 0; r < rank; r++)
      cp_assert(dims[r].is == dims[r].os);
    for (int r = 0; r < howmany_rank; r++)
      cp_assert(howmany_dims[r].is == howmany_dims[r].os);
  } else {
    cp_assert(!inplace);
  }
  *plan = fftw_plan_guru_dft(rank, dims, howmany_rank, howmany_dims, grid_in,
                             grid_out, direction ? FFTW_FORWARD : FFTW_BACKWARD,
                             fftw_planning_mode);
  add_plan_to_cache(key, plan);
  cp_assert(plan != NULL);
  omp_unset_nest_lock(&cache_lock);
  return plan;
}

/*******************************************************************************
 * \brief Create plan of a local R2C/C2R 3D FFT.
 * \param key Encoded FFT key.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
fftw_plan *fft_fftw_create_guru_plan_r2c(const fft_key_t key, double *grid_in,
                                         double *grid_out) {
  cp_assert(has_guru_interface);
  omp_set_nest_lock(&cache_lock);
  fftw_plan *cached_plan = lookup_plan_from_cache(key);
  if (cached_plan != NULL) {
    omp_unset_nest_lock(&cache_lock);
    return cached_plan;
  }
  bool direction, inplace;
  int rank, howmany_rank, number_of_threads;
  fft_iodim dims[3], howmany_dims[3];
  fetch_data_from_key_guru(key, &direction, &rank, dims, &howmany_rank,
                           howmany_dims, &number_of_threads, &inplace);

  fftw_plan_with_nthreads(number_of_threads);
  fftw_plan *plan = malloc(sizeof(fftw_plan));
  if (direction) {
    if (grid_in == grid_out) {
      cp_assert(inplace);
      for (int r = 0; r < rank; r++)
        cp_assert(dims[r].is == 2 * dims[r].os ||
                  (dims[r].is == 1 && dims[r].os == 1));
      for (int r = 0; r < howmany_rank; r++)
        cp_assert(howmany_dims[r].is == 2 * howmany_dims[r].os);
    } else {
      cp_assert(!inplace);
    }
    *plan =
        fftw_plan_guru_dft_r2c(rank, dims, howmany_rank, howmany_dims, grid_in,
                               (double complex *)grid_out, fftw_planning_mode);
  } else {
    if (grid_in == grid_out) {
      cp_assert(inplace);
      for (int r = 0; r < rank; r++)
        cp_assert(2 * dims[r].is == dims[r].os ||
                  (dims[r].is == 1 && dims[r].os == 1));
      for (int r = 0; r < howmany_rank; r++)
        cp_assert(2 * howmany_dims[r].is == howmany_dims[r].os);
    } else {
      cp_assert(!inplace);
    }
    // We use the buffers the other way around to prevent
    // out-of-bounds-accesses of the planner if the output array has only the
    // minimum size
    *plan = fftw_plan_guru_dft_c2r(rank, dims, howmany_rank, howmany_dims,
                                   (double complex *)grid_in, grid_out,
                                   fftw_planning_mode);
  }
  add_plan_to_cache(key, plan);
  cp_assert(plan != NULL);
  omp_unset_nest_lock(&cache_lock);
  return plan;
}
#endif

/*******************************************************************************
 * \brief Performs a local C2C 1D FFT.
 * \param direction FFT direction (forward or backward).
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param leading_dimension_in Leading dimension of the input array.
 * \param leading_dimension_out Leading dimension of the output array.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_1d_local(
    const bool direction, const int fft_size, const int number_of_ffts,
    const bool transpose_in, const bool transpose_out,
    const int leading_dimension_in, const int leading_dimension_out,
    double complex *grid_in, double complex *grid_out) {
#if defined(__FFTW3)
  cp_assert(is_initialized);
  if (fft_size == 0 || number_of_ffts == 0)
    return;
  const bool in_place = grid_in == grid_out;
  const int number_of_threads = fft_team_size();
  if (fftw_planning_mode == FFTW_ESTIMATE) {
    fftw_plan *plan = NULL;
    const int block_size =
        (number_of_ffts + number_of_threads - 1) / number_of_threads;
    const int number_of_full_blocks = number_of_ffts / block_size;
    fft_key_t key;
    get_key_1d(direction, fft_size, block_size, transpose_in, transpose_out,
               leading_dimension_in, leading_dimension_out, 1, in_place, key);
    plan = lookup_plan_from_cache(key);
    if (plan == NULL)
      plan = fft_fftw_create_1d_plan(key, grid_in, grid_out);
    if (block_size * number_of_full_blocks != number_of_ffts) {
      const int block_size_last_thread = number_of_ffts % block_size;
      fft_key_t key;
      get_key_1d(direction, fft_size, block_size_last_thread, transpose_in,
                 transpose_out, leading_dimension_in, leading_dimension_out, 1,
                 in_place, key);
      plan = lookup_plan_from_cache(key);
      if (plan == NULL)
        plan = fft_fftw_create_1d_plan(key, grid_in, grid_out);
    }
  } else {
    fft_key_t key;
    get_key_1d(direction, fft_size, number_of_ffts, transpose_in, transpose_out,
               leading_dimension_in, leading_dimension_out, fft_team_size(),
               in_place, key);
    fftw_plan *plan = lookup_plan_from_cache(key);
    if (plan == NULL)
      plan = fft_fftw_create_1d_plan(key, grid_in, grid_out);
  }
#else
  (void)direction;
  (void)fft_size;
  (void)number_of_ffts;
  (void)grid_in;
  (void)grid_out;
  (void)transpose_in;
  (void)transpose_out;
  (void)leading_dimension_in;
  (void)leading_dimension_out;
  cp_assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local forward R2C FFT.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param leading_dimension_in Leading dimension of the input array.
 * \param leading_dimension_out Leading dimension of the output array.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_1d_r2c_local(const int fft_size,
                                    const int number_of_ffts,
                                    const bool transpose_in,
                                    const bool transpose_out,
                                    const int leading_dimension_in,
                                    const int leading_dimension_out,
                                    double *grid_in, double complex *grid_out) {
#if defined(__FFTW3)
  cp_assert(is_initialized);
  if (fft_size == 0 || number_of_ffts == 0)
    return;
  const bool in_place = grid_in == (double *)grid_out;
  const int number_of_threads = fft_team_size();
  // The self-written plans do not work with in-place transforms using the
  // transposed format
  if (fftw_planning_mode == FFTW_ESTIMATE && (!in_place || !transpose_in)) {
    fftw_plan *plan = NULL;
    const int block_size =
        (number_of_ffts + number_of_threads - 1) / number_of_threads;
    const int number_of_full_blocks = number_of_ffts / block_size;
    fft_key_t key;
    get_key_1d_r2c(true, fft_size, block_size, transpose_in, transpose_out,
                   leading_dimension_in, leading_dimension_out, 1, in_place,
                   key);
    plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      plan = fft_fftw_create_1d_plan_r2c(key, grid_in, (double *)grid_out);
    }
    if (block_size * number_of_full_blocks != number_of_ffts) {
      const int block_size_last_thread = number_of_ffts % block_size;
      fft_key_t key;
      get_key_1d_r2c(true, fft_size, block_size_last_thread, transpose_in,
                     transpose_out, leading_dimension_in, leading_dimension_out,
                     1, in_place, key);
      plan = lookup_plan_from_cache(key);
      if (plan == NULL) {
        plan = fft_fftw_create_1d_plan_r2c(key, grid_in, (double *)grid_out);
      }
    }
  } else {
    fft_key_t key;
    get_key_1d_r2c(true, fft_size, number_of_ffts, transpose_in, transpose_out,
                   leading_dimension_in, leading_dimension_out, fft_team_size(),
                   in_place, key);
    fftw_plan *plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      plan = fft_fftw_create_1d_plan_r2c(key, grid_in, (double *)grid_out);
    }
  }
#else
  (void)fft_size;
  (void)number_of_ffts;
  (void)grid_in;
  (void)grid_out;
  (void)transpose_in;
  (void)transpose_out;
  (void)leading_dimension_in;
  (void)leading_dimension_out;
  cp_assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local backwards C2R 1D FFT
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param leading_dimension_in Leading dimension of the input array.
 * \param leading_dimension_out Leading dimension of the output array.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_1d_c2r_local(const int fft_size,
                                    const int number_of_ffts,
                                    const bool transpose_in,
                                    const bool transpose_out,
                                    const int leading_dimension_in,
                                    const int leading_dimension_out,
                                    double complex *grid_in, double *grid_out) {
#if defined(__FFTW3)
  cp_assert(is_initialized);
  if (fft_size == 0 || number_of_ffts == 0)
    return;
  const bool in_place = (double *)grid_in == grid_out;
  const int number_of_threads = fft_team_size();
  // The self-written plans do not work with in-place transforms using the
  // transposed format
  if (fftw_planning_mode == FFTW_ESTIMATE && (!in_place || !transpose_in)) {
    fftw_plan *plan = NULL;
    const int block_size =
        (number_of_ffts + number_of_threads - 1) / number_of_threads;
    const int number_of_full_blocks = number_of_ffts / block_size;
    fft_key_t key;
    get_key_1d_r2c(false, fft_size, block_size, transpose_in, transpose_out,
                   leading_dimension_in, leading_dimension_out, 1, in_place,
                   key);
    plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      plan = fft_fftw_create_1d_plan_r2c(key, (double *)grid_in, grid_out);
    }
    if (block_size * number_of_full_blocks != number_of_ffts) {
      const int block_size_last_thread = number_of_ffts % block_size;
      fft_key_t key;
      get_key_1d_r2c(false, fft_size, block_size_last_thread, transpose_in,
                     transpose_out, leading_dimension_in, leading_dimension_out,
                     1, in_place, key);
      plan = lookup_plan_from_cache(key);
      if (plan == NULL) {
        plan = fft_fftw_create_1d_plan_r2c(key, (double *)grid_in, grid_out);
      }
    }
  } else {
    fft_key_t key;
    get_key_1d_r2c(false, fft_size, number_of_ffts, transpose_in, transpose_out,
                   leading_dimension_in, leading_dimension_out, fft_team_size(),
                   in_place, key);
    fftw_plan *plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      plan = fft_fftw_create_1d_plan_r2c(key, (double *)grid_in, grid_out);
    }
  }
#else
  (void)fft_size;
  (void)number_of_ffts;
  (void)grid_in;
  (void)grid_out;
  (void)transpose_in;
  (void)transpose_out;
  (void)leading_dimension_in;
  (void)leading_dimension_out;
  cp_assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local C2C 2D FFT
 * \param direction FFT direction (forward or backward).
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_2d_local(const bool direction, const int fft_size[2],
                                const int number_of_ffts,
                                const bool transpose_in,
                                const bool transpose_out,
                                double complex *grid_in,
                                double complex *grid_out) {
#if defined(__FFTW3)
  cp_assert(is_initialized);
  if (fft_size[0] == 0 || fft_size[1] == 0 || number_of_ffts == 0)
    return;
  const bool in_place = grid_in == grid_out;
  const int number_of_threads = fft_team_size();
  if (fftw_planning_mode == FFTW_ESTIMATE) {
    fftw_plan *plan = NULL, *plan_last_thread = NULL;
    const int block_size =
        (number_of_ffts + number_of_threads - 1) / number_of_threads;
    const int number_of_full_blocks = number_of_ffts / block_size;
    fft_key_t key;
    get_key_2d(direction, fft_size, block_size, transpose_in, transpose_out,
               transpose_in ? number_of_ffts : fft_size[0] * fft_size[1],
               transpose_out ? number_of_ffts : fft_size[0] * fft_size[1], 1,
               in_place, key);
    plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      plan = fft_fftw_create_2d_plan(key, grid_in, grid_out);
    }
    if (block_size * number_of_full_blocks != number_of_ffts) {
      const int block_size_last_thread = number_of_ffts % block_size;
      fft_key_t key;
      get_key_2d(direction, fft_size, block_size_last_thread, transpose_in,
                 transpose_out,
                 transpose_in ? number_of_ffts : fft_size[0] * fft_size[1],
                 transpose_out ? number_of_ffts : fft_size[0] * fft_size[1], 1,
                 in_place, key);
      plan_last_thread = lookup_plan_from_cache(key);
      if (plan_last_thread == NULL) {
        plan_last_thread = fft_fftw_create_2d_plan(key, grid_in, grid_out);
      }
    }
  } else {
    fft_key_t key;
    get_key_2d(direction, fft_size, number_of_ffts, transpose_in, transpose_out,
               transpose_in ? number_of_ffts : fft_size[0] * fft_size[1],
               transpose_out ? number_of_ffts : fft_size[0] * fft_size[1],
               fft_team_size(), grid_in == grid_out, key);
    fftw_plan *plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      plan = fft_fftw_create_2d_plan(key, grid_in, grid_out);
    }
  }
#else
  (void)direction;
  (void)fft_size;
  (void)number_of_ffts;
  (void)grid_in;
  (void)grid_out;
  (void)transpose_in;
  (void)transpose_out;
  cp_assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local forward R2C 2D FFT.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_2d_r2c_local(const int fft_size[2],
                                    const int number_of_ffts,
                                    const bool transpose_in,
                                    const bool transpose_out, double *grid_in,
                                    double complex *grid_out) {
#if defined(__FFTW3)
  cp_assert(is_initialized);
  if (fft_size[0] == 0 || fft_size[1] == 0 || number_of_ffts == 0)
    return;
  const bool in_place = grid_in == (double *)grid_out;
  const int number_of_threads = fft_team_size();
  // The self-written plans do not work with in-place transforms using the
  // transposed format
  if (fftw_planning_mode == FFTW_ESTIMATE && (!in_place || !transpose_in)) {
    fftw_plan *plan = NULL, *plan_last_thread = NULL;
    const int block_size =
        (number_of_ffts + number_of_threads - 1) / number_of_threads;
    const int number_of_full_blocks = number_of_ffts / block_size;
    fft_key_t key;
    get_key_2d_r2c(
        true, fft_size, block_size, transpose_in, transpose_out,
        transpose_in ? number_of_ffts
                     : fft_size[0] *
                           (in_place ? 2 * (fft_size[1] / 2 + 1) : fft_size[1]),
        transpose_out ? number_of_ffts : fft_size[0] * (fft_size[1] / 2 + 1), 1,
        in_place, key);
    plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      plan = fft_fftw_create_2d_plan_r2c(key, grid_in, (double *)grid_out);
    }
    if (block_size * number_of_full_blocks != number_of_ffts) {
      const int block_size_last_thread = number_of_ffts % block_size;
      fft_key_t key;
      get_key_2d_r2c(
          true, fft_size, block_size_last_thread, transpose_in, transpose_out,
          transpose_in ? number_of_ffts
                       : fft_size[0] * (in_place ? 2 * (fft_size[1] / 2 + 1)
                                                 : fft_size[1]),
          transpose_out ? number_of_ffts : fft_size[0] * (fft_size[1] / 2 + 1),
          1, in_place, key);
      plan_last_thread = lookup_plan_from_cache(key);
      if (plan_last_thread == NULL) {
        plan_last_thread =
            fft_fftw_create_2d_plan_r2c(key, grid_in, (double *)grid_out);
      }
    }
  } else {
    fft_key_t key;
    get_key_2d_r2c(
        true, fft_size, number_of_ffts, transpose_in, transpose_out,
        transpose_in ? number_of_ffts
                     : fft_size[0] *
                           (in_place ? 2 * (fft_size[1] / 2 + 1) : fft_size[1]),
        transpose_out ? number_of_ffts : fft_size[0] * (fft_size[1] / 2 + 1),
        fft_team_size(), grid_in == (double *)grid_out, key);
    fftw_plan *plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      plan = fft_fftw_create_2d_plan_r2c(key, grid_in, (double *)grid_out);
    }
  }
#else
  (void)fft_size;
  (void)number_of_ffts;
  (void)grid_in;
  (void)grid_out;
  (void)transpose_in;
  (void)transpose_out;
  cp_assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local backwards C2R 2D FFT.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_2d_c2r_local(const int fft_size[2],
                                    const int number_of_ffts,
                                    const bool transpose_in,
                                    const bool transpose_out,
                                    double complex *grid_in, double *grid_out) {
#if defined(__FFTW3)
  cp_assert(is_initialized);
  if (fft_size[0] == 0 || fft_size[1] == 0 || number_of_ffts == 0)
    return;
  const bool in_place = (double *)grid_in == grid_out;
  const int number_of_threads = fft_team_size();
  // The self-written plans do not work with in-place transforms using the
  // transposed format
  if (fftw_planning_mode == FFTW_ESTIMATE && (!in_place || !transpose_in)) {
    fftw_plan *plan = NULL, *plan_last_thread = NULL;
    const int block_size =
        (number_of_ffts + number_of_threads - 1) / number_of_threads;
    const int number_of_full_blocks = number_of_ffts / block_size;
    fft_key_t key;
    get_key_2d_r2c(
        false, fft_size, block_size, transpose_in, transpose_out,
        transpose_in ? number_of_ffts : fft_size[0] * (fft_size[1] / 2 + 1),
        transpose_out ? number_of_ffts
                      : fft_size[0] * (in_place ? 2 * (fft_size[1] / 2 + 1)
                                                : fft_size[1]),
        1, in_place, key);
    plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      plan = fft_fftw_create_2d_plan_r2c(key, (double *)grid_in, grid_out);
    }
    if (block_size * number_of_full_blocks != number_of_ffts) {
      const int block_size_last_thread = number_of_ffts % block_size;
      fft_key_t key;
      get_key_2d_r2c(
          false, fft_size, block_size_last_thread, transpose_in, transpose_out,
          transpose_in ? number_of_ffts : fft_size[0] * (fft_size[1] / 2 + 1),
          transpose_out ? number_of_ffts
                        : fft_size[0] * (in_place ? 2 * (fft_size[1] / 2 + 1)
                                                  : fft_size[1]),
          1, in_place, key);
      plan_last_thread = lookup_plan_from_cache(key);
      if (plan_last_thread == NULL) {
        plan_last_thread =
            fft_fftw_create_2d_plan_r2c(key, (double *)grid_in, grid_out);
      }
    }
  } else {
    fft_key_t key;
    get_key_2d_r2c(
        false, fft_size, number_of_ffts, transpose_in, transpose_out,
        transpose_in ? number_of_ffts : fft_size[0] * (fft_size[1] / 2 + 1),
        transpose_out ? number_of_ffts
                      : fft_size[0] * (in_place ? 2 * (fft_size[1] / 2 + 1)
                                                : fft_size[1]),
        fft_team_size(), grid_in == (double complex *)grid_out, key);
    fftw_plan *plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      plan = fft_fftw_create_2d_plan_r2c(key, (double *)grid_in, grid_out);
    }
  }
#else
  (void)fft_size;
  (void)number_of_ffts;
  (void)grid_in;
  (void)grid_out;
  (void)transpose_in;
  (void)transpose_out;
  cp_assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local C2C 3D FFT.
 * \param direction FFT direction (forward or backward).
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_3d_local(const bool direction, const int fft_size[3],
                                const int number_of_ffts,
                                const bool transpose_in,
                                const bool transpose_out,
                                double complex *grid_in,
                                double complex *grid_out) {
#if defined(__FFTW3)
  cp_assert(is_initialized);
  if (fft_size[0] == 0 || fft_size[1] == 0 || fft_size[2] == 0 ||
      number_of_ffts == 0)
    return;
  const bool in_place = grid_in == grid_out;
  if (number_of_ffts > 1 &&
      use_threaded_3d(fft_size, transpose_in, transpose_out, in_place)) {
    register_threaded_3d_c2c(direction, fft_size, number_of_ffts, in_place,
                             grid_in, grid_out);
  } else if (number_of_ffts == 1 &&
             use_threaded_3d(fft_size, transpose_in, transpose_out, in_place)) {
    register_threaded_3d_c2c_1d(direction, fft_size, in_place, grid_in,
                                grid_out);
  } else {
    fft_key_t key;
    get_key_3d(direction, fft_size, number_of_ffts, transpose_in, transpose_out,
               fft_team_size(), in_place, key);
    fftw_plan *plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      plan = fft_fftw_create_3d_plan(key, grid_in, grid_out);
    }
  }
#else
  (void)fft_size;
  (void)number_of_ffts;
  (void)transpose_in;
  (void)transpose_out;
  (void)grid_in;
  (void)grid_out;
  cp_assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local forward R2C 3D FFT.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_3d_r2c_local(const int fft_size[3],
                                    const int number_of_ffts,
                                    const bool transpose_in,
                                    const bool transpose_out, double *grid_in,
                                    double complex *grid_out) {
#if defined(__FFTW3)
  cp_assert(is_initialized);
  if (fft_size[0] == 0 || fft_size[1] == 0 || fft_size[2] == 0 ||
      number_of_ffts == 0)
    return;
  const bool in_place = grid_in == (double *)grid_out;
  if (number_of_ffts > 1 &&
      use_threaded_3d(fft_size, transpose_in, transpose_out, in_place)) {
    register_threaded_3d_r2c(true, fft_size, number_of_ffts, in_place, grid_in,
                             (double *)grid_out);
  } else if (number_of_ffts == 1 &&
             use_threaded_3d(fft_size, transpose_in, transpose_out, in_place)) {
    register_threaded_3d_r2c_1d(true, fft_size, in_place, grid_in, grid_out);
  } else {
    fft_key_t key;
    get_key_3d_r2c(true, fft_size, number_of_ffts, transpose_in, transpose_out,
                   fft_team_size(), in_place, key);
    fftw_plan *plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      plan = fft_fftw_create_3d_plan_r2c(key, grid_in, (double *)grid_out);
    }
  }
#else
  (void)fft_size;
  (void)number_of_ffts;
  (void)transpose_in;
  (void)transpose_out;
  (void)grid_in;
  (void)grid_out;
  cp_assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local backwards R2C 3D FFT.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_3d_c2r_local(const int fft_size[3],
                                    const int number_of_ffts,
                                    const bool transpose_in,
                                    const bool transpose_out,
                                    double complex *grid_in, double *grid_out) {
#if defined(__FFTW3)
  cp_assert(is_initialized);
  if (fft_size[0] == 0 || fft_size[1] == 0 || fft_size[2] == 0 ||
      number_of_ffts == 0)
    return;
  const bool in_place = (double *)grid_in == grid_out;
  if (number_of_ffts > 1 &&
      use_threaded_3d(fft_size, transpose_in, transpose_out, in_place)) {
    register_threaded_3d_r2c(false, fft_size, number_of_ffts, in_place,
                             (double *)grid_in, grid_out);
  } else if (number_of_ffts == 1 &&
             use_threaded_3d(fft_size, transpose_in, transpose_out, in_place)) {
    register_threaded_3d_r2c_1d(false, fft_size, in_place, grid_out, grid_in);
  } else {
    fft_key_t key;
    get_key_3d_r2c(false, fft_size, number_of_ffts, transpose_in, transpose_out,
                   fft_team_size(), in_place, key);
    fftw_plan *plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      plan = fft_fftw_create_3d_plan_r2c(key, (double *)grid_in, grid_out);
    }
  }
#else
  (void)fft_size;
  (void)number_of_ffts;
  (void)transpose_in;
  (void)transpose_out;
  (void)grid_in;
  (void)grid_out;
  cp_assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local forward C2C 1D FFT.
 * \param direction FFT direction (forward or backward).
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param leading_dimension_in Leading dimension of the input array.
 * \param leading_dimension_out Leading dimension of the output array.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_1d_local(const bool direction, const int fft_size,
                       const int number_of_ffts, const bool transpose_in,
                       const bool transpose_out, const int leading_dimension_in,
                       const int leading_dimension_out, double complex *grid_in,
                       double complex *grid_out) {
#if defined(__FFTW3)
  cp_assert(is_initialized);
  if (fft_size == 0 || number_of_ffts == 0)
    return;
  const bool in_place = grid_in == grid_out;
  const int number_of_threads = fft_team_size();
  if (number_of_threads == 1 && fftw_planning_mode == FFTW_ESTIMATE) {
    // Called from inside an OpenMP region: use the original batched strategy.
    fftw_plan *plan = NULL, *plan_last_thread = NULL;
    const int block_size =
        (number_of_ffts + number_of_threads - 1) / number_of_threads;
    const int number_of_full_blocks = number_of_ffts / block_size;
    const int total_number_of_blocks =
        number_of_full_blocks + (number_of_ffts % block_size != 0);
    fft_key_t key;
    get_key_1d(direction, fft_size, block_size, transpose_in, transpose_out,
               leading_dimension_in, leading_dimension_out, 1, in_place, key);
    plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      double complex *buffer =
          fftw_alloc_complex(get_buffer_size_from_key(key));
      plan = fft_fftw_create_1d_plan(key, buffer, in_place ? buffer : grid_out);
      fftw_free(buffer);
    }
    if (block_size * number_of_full_blocks != number_of_ffts) {
      const int block_size_last_thread = number_of_ffts % block_size;
      fft_key_t key_last;
      get_key_1d(direction, fft_size, block_size_last_thread, transpose_in,
                 transpose_out, leading_dimension_in, leading_dimension_out, 1,
                 in_place, key_last);
      plan_last_thread = lookup_plan_from_cache(key_last);
      if (plan_last_thread == NULL) {
        double complex *buffer =
            fftw_alloc_complex(get_buffer_size_from_key(key_last));
        plan_last_thread = fft_fftw_create_1d_plan(
            key_last, buffer, in_place ? buffer : grid_out);
        fftw_free(buffer);
      }
    }
    const int offset_in =
        transpose_in ? block_size : block_size * leading_dimension_in;
    const int offset_out =
        transpose_out ? block_size : block_size * leading_dimension_out;
#pragma omp parallel num_threads(number_of_threads) default(none)              \
    shared(number_of_full_blocks, plan, grid_in, grid_out, offset_in,          \
               offset_out, plan_last_thread, total_number_of_blocks)
    {
      const int thread_id = omp_get_thread_num();
      if (thread_id < number_of_full_blocks) {
        fftw_execute_dft(*plan, grid_in + thread_id * offset_in,
                         grid_out + thread_id * offset_out);
      } else if (thread_id + 1 == total_number_of_blocks) {
        fftw_execute_dft(*plan_last_thread, grid_in + thread_id * offset_in,
                         grid_out + thread_id * offset_out);
      }
    }
  } else {
    fft_key_t key;
    get_key_1d(direction, fft_size, number_of_ffts, transpose_in, transpose_out,
               leading_dimension_in, leading_dimension_out, number_of_threads,
               in_place, key);
    fftw_plan *plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      double complex *buffer =
          fftw_alloc_complex(get_buffer_size_from_key(key));
      plan = fft_fftw_create_1d_plan(key, buffer, in_place ? buffer : grid_out);
      fftw_free(buffer);
    }
    fftw_execute_dft(*plan, grid_in, grid_out);
  }
#else
  (void)direction;
  (void)fft_size;
  (void)number_of_ffts;
  (void)grid_in;
  (void)grid_out;
  (void)transpose_in;
  (void)transpose_out;
  (void)leading_dimension_in;
  (void)leading_dimension_out;
  cp_assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local forward R2C FFT.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param leading_dimension_in Leading dimension of the input array.
 * \param leading_dimension_out Leading dimension of the output array.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_1d_r2c_local(const int fft_size, const int number_of_ffts,
                           const bool transpose_in, const bool transpose_out,
                           const int leading_dimension_in,
                           const int leading_dimension_out, double *grid_in,
                           double complex *grid_out) {
#if defined(__FFTW3)
  cp_assert(is_initialized);
  if (fft_size == 0 || number_of_ffts == 0)
    return;
  const bool in_place = grid_in == (double *)grid_out;
  const int number_of_threads = fft_team_size();
  // The self-written plans do not work with in-place transforms using the
  // transposed format
  if (fftw_planning_mode == FFTW_ESTIMATE && (!in_place || !transpose_in)) {
    fftw_plan *plan = NULL, *plan_last_thread = NULL;
    const int block_size =
        (number_of_ffts + number_of_threads - 1) / number_of_threads;
    const int number_of_full_blocks = number_of_ffts / block_size;
    const int total_number_of_blocks =
        number_of_full_blocks + (number_of_ffts % block_size != 0);
    fft_key_t key;
    get_key_1d_r2c(true, fft_size, block_size, transpose_in, transpose_out,
                   leading_dimension_in, leading_dimension_out, 1, in_place,
                   key);
    plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      double *buffer = fftw_alloc_real(2 * get_buffer_size_from_key(key));
      plan = fft_fftw_create_1d_plan_r2c(
          key, buffer, in_place ? buffer : (double *)grid_out);
      fftw_free(buffer);
    }
    if (block_size * number_of_full_blocks != number_of_ffts) {
      const int block_size_last_thread = number_of_ffts % block_size;
      fft_key_t key;
      get_key_1d_r2c(true, fft_size, block_size_last_thread, transpose_in,
                     transpose_out, leading_dimension_in, leading_dimension_out,
                     1, in_place, key);
      plan_last_thread = lookup_plan_from_cache(key);
      if (plan_last_thread == NULL) {
        double *buffer = fftw_alloc_real(2 * get_buffer_size_from_key(key));
        plan_last_thread = fft_fftw_create_1d_plan_r2c(
            key, buffer, in_place ? buffer : (double *)grid_out);
        fftw_free(buffer);
      }
    }
    const int offset_in =
        transpose_in ? block_size : block_size * leading_dimension_in;
    const int offset_out =
        transpose_out ? block_size : block_size * leading_dimension_out;
#pragma omp parallel default(none)                                             \
    shared(number_of_full_blocks, plan, grid_in, grid_out, offset_in,          \
               offset_out, plan_last_thread, total_number_of_blocks)
    {
      const int thread_id = omp_get_thread_num();
      if (thread_id < number_of_full_blocks) {
        fftw_execute_dft_r2c(*plan, grid_in + thread_id * offset_in,
                             grid_out + thread_id * offset_out);
      } else if (thread_id + 1 == total_number_of_blocks) {
        fftw_execute_dft_r2c(*plan_last_thread, grid_in + thread_id * offset_in,
                             grid_out + thread_id * offset_out);
      }
    }
  } else {
    fft_key_t key;
    get_key_1d_r2c(true, fft_size, number_of_ffts, transpose_in, transpose_out,
                   leading_dimension_in, leading_dimension_out, fft_team_size(),
                   in_place, key);
    fftw_plan *plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      double *buffer = fftw_alloc_real(2 * get_buffer_size_from_key(key));
      plan = fft_fftw_create_1d_plan_r2c(
          key, buffer, in_place ? buffer : (double *)grid_out);
      fftw_free(buffer);
    }
    cp_assert(plan != NULL);
    fftw_execute_dft_r2c(*plan, grid_in, grid_out);
  }
#else
  (void)fft_size;
  (void)number_of_ffts;
  (void)grid_in;
  (void)grid_out;
  (void)transpose_in;
  (void)transpose_out;
  (void)leading_dimension_in;
  (void)leading_dimension_out;
  cp_assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local backwards C2R 1D FFT
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param leading_dimension_in Leading dimension of the input array.
 * \param leading_dimension_out Leading dimension of the output array.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_1d_c2r_local(const int fft_size, const int number_of_ffts,
                           const bool transpose_in, const bool transpose_out,
                           const int leading_dimension_in,
                           const int leading_dimension_out,
                           double complex *grid_in, double *grid_out) {
#if defined(__FFTW3)
  cp_assert(is_initialized);
  if (fft_size == 0 || number_of_ffts == 0)
    return;
  const bool in_place = (double *)grid_in == grid_out;
  const int number_of_threads = fft_team_size();
  // The self-written plans do not work with in-place transforms using the
  // transposed format
  if (fftw_planning_mode == FFTW_ESTIMATE && (!in_place || !transpose_in)) {
    fftw_plan *plan = NULL, *plan_last_thread = NULL;
    const int block_size =
        (number_of_ffts + number_of_threads - 1) / number_of_threads;
    const int number_of_full_blocks = number_of_ffts / block_size;
    const int total_number_of_blocks =
        number_of_full_blocks + (number_of_ffts % block_size != 0);
    fft_key_t key;
    get_key_1d_r2c(false, fft_size, block_size, transpose_in, transpose_out,
                   leading_dimension_in, leading_dimension_out, 1, in_place,
                   key);
    plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      double *buffer = fftw_alloc_real(2 * get_buffer_size_from_key(key));
      plan = fft_fftw_create_1d_plan_r2c(key, buffer,
                                         in_place ? buffer : grid_out);
      fftw_free(buffer);
    }
    if (block_size * number_of_full_blocks != number_of_ffts) {
      const int block_size_last_thread = number_of_ffts % block_size;
      fft_key_t key;
      get_key_1d_r2c(false, fft_size, block_size_last_thread, transpose_in,
                     transpose_out, leading_dimension_in, leading_dimension_out,
                     1, in_place, key);
      plan_last_thread = lookup_plan_from_cache(key);
      if (plan_last_thread == NULL) {
        double *buffer = fftw_alloc_real(2 * get_buffer_size_from_key(key));
        plan_last_thread = fft_fftw_create_1d_plan_r2c(
            key, buffer, in_place ? buffer : grid_out);
        fftw_free(buffer);
      }
    }
    const int offset_in =
        transpose_in ? block_size : block_size * leading_dimension_in;
    const int offset_out =
        transpose_out ? block_size : block_size * leading_dimension_out;
#pragma omp parallel default(none)                                             \
    shared(number_of_full_blocks, plan, grid_in, grid_out, offset_in,          \
               offset_out, plan_last_thread, total_number_of_blocks)
    {
      const int thread_id = omp_get_thread_num();
      if (thread_id < number_of_full_blocks) {
        fftw_execute_dft_c2r(*plan, grid_in + thread_id * offset_in,
                             grid_out + thread_id * offset_out);
      } else if (thread_id + 1 == total_number_of_blocks) {
        fftw_execute_dft_c2r(*plan_last_thread, grid_in + thread_id * offset_in,
                             grid_out + thread_id * offset_out);
      }
    }
  } else {
    fft_key_t key;
    get_key_1d_r2c(false, fft_size, number_of_ffts, transpose_in, transpose_out,
                   leading_dimension_in, leading_dimension_out, fft_team_size(),
                   grid_in == (double complex *)grid_out, key);
    fftw_plan *plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      double *buffer = fftw_alloc_real(2 * get_buffer_size_from_key(key));
      plan = fft_fftw_create_1d_plan_r2c(key, buffer,
                                         in_place ? buffer : grid_out);
      fftw_free(buffer);
    }
    cp_assert(plan != NULL);
    fftw_execute_dft_c2r(*plan, grid_in, grid_out);
  }
#else
  (void)fft_size;
  (void)number_of_ffts;
  (void)grid_in;
  (void)grid_out;
  (void)transpose_in;
  (void)transpose_out;
  (void)leading_dimension_in;
  (void)leading_dimension_out;
  cp_assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local C2C 2D FFT
 * \param direction FFT direction (forward or backward).
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_2d_local(const bool direction, const int fft_size[2],
                       const int number_of_ffts, const bool transpose_in,
                       const bool transpose_out, double complex *grid_in,
                       double complex *grid_out) {
#if defined(__FFTW3)
  cp_assert(is_initialized);
  if (fft_size[0] == 0 || fft_size[1] == 0 || number_of_ffts == 0)
    return;
  const bool in_place = grid_in == grid_out;
  const int number_of_threads = fft_team_size();
  if (fftw_planning_mode == FFTW_ESTIMATE) {
    fftw_plan *plan = NULL, *plan_last_thread = NULL;
    const int block_size =
        (number_of_ffts + number_of_threads - 1) / number_of_threads;
    const int number_of_full_blocks = number_of_ffts / block_size;
    const int total_number_of_blocks =
        number_of_full_blocks + (number_of_ffts % block_size != 0);
    fft_key_t key;
    get_key_2d(direction, fft_size, block_size, transpose_in, transpose_out,
               transpose_in ? number_of_ffts : fft_size[0] * fft_size[1],
               transpose_out ? number_of_ffts : fft_size[0] * fft_size[1], 1,
               in_place, key);
    plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      double complex *buffer =
          fftw_alloc_complex(get_buffer_size_from_key(key));
      plan = fft_fftw_create_2d_plan(key, buffer, in_place ? buffer : grid_out);
      fftw_free(buffer);
    }
    if (block_size * number_of_full_blocks != number_of_ffts) {
      const int block_size_last_thread = number_of_ffts % block_size;
      fft_key_t key;
      get_key_2d(direction, fft_size, block_size_last_thread, transpose_in,
                 transpose_out,
                 transpose_in ? number_of_ffts : fft_size[0] * fft_size[1],
                 transpose_out ? number_of_ffts : fft_size[0] * fft_size[1], 1,
                 in_place, key);
      plan_last_thread = lookup_plan_from_cache(key);
      if (plan_last_thread == NULL) {
        double complex *buffer =
            fftw_alloc_complex(get_buffer_size_from_key(key));
        plan_last_thread =
            fft_fftw_create_2d_plan(key, buffer, in_place ? buffer : grid_out);
        fftw_free(buffer);
      }
    }
    const int offset_in =
        transpose_in ? block_size : block_size * fft_size[0] * fft_size[1];
    const int offset_out =
        transpose_out ? block_size : block_size * fft_size[0] * fft_size[1];
#pragma omp parallel default(none)                                             \
    shared(number_of_full_blocks, plan, grid_in, grid_out, offset_in,          \
               offset_out, plan_last_thread, total_number_of_blocks)
    {
      const int thread_id = omp_get_thread_num();
      if (thread_id < number_of_full_blocks) {
        fftw_execute_dft(*plan, grid_in + thread_id * offset_in,
                         grid_out + thread_id * offset_out);
      } else if (thread_id + 1 == total_number_of_blocks) {
        fftw_execute_dft(*plan_last_thread, grid_in + thread_id * offset_in,
                         grid_out + thread_id * offset_out);
      }
    }
  } else {
    fft_key_t key;
    get_key_2d(direction, fft_size, number_of_ffts, transpose_in, transpose_out,
               transpose_in ? number_of_ffts : fft_size[0] * fft_size[1],
               transpose_out ? number_of_ffts : fft_size[0] * fft_size[1],
               fft_team_size(), grid_in == grid_out, key);
    fftw_plan *plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      double complex *buffer =
          fftw_alloc_complex(get_buffer_size_from_key(key));
      plan = fft_fftw_create_2d_plan(key, buffer,
                                     grid_in == grid_out ? buffer : grid_out);
      fftw_free(buffer);
    }
    fftw_execute_dft(*plan, grid_in, grid_out);
  }
#else
  (void)direction;
  (void)fft_size;
  (void)number_of_ffts;
  (void)grid_in;
  (void)grid_out;
  (void)transpose_in;
  (void)transpose_out;
  cp_assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local forward R2C 2D FFT.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_2d_r2c_local(const int fft_size[2], const int number_of_ffts,
                           const bool transpose_in, const bool transpose_out,
                           double *grid_in, double complex *grid_out) {
#if defined(__FFTW3)
  cp_assert(is_initialized);
  if (fft_size[0] == 0 || fft_size[1] == 0 || number_of_ffts == 0)
    return;
  const bool in_place = grid_in == (double *)grid_out;
  const int number_of_threads = fft_team_size();
  // The self-written plans do not work with in-place transforms using the
  // transposed format
  if (fftw_planning_mode == FFTW_ESTIMATE && (!in_place || !transpose_in)) {
    fftw_plan *plan = NULL, *plan_last_thread = NULL;
    const int block_size =
        (number_of_ffts + number_of_threads - 1) / number_of_threads;
    const int number_of_full_blocks = number_of_ffts / block_size;
    const int total_number_of_blocks =
        number_of_full_blocks + (number_of_ffts % block_size != 0);
    fft_key_t key;
    get_key_2d_r2c(
        true, fft_size, block_size, transpose_in, transpose_out,
        transpose_in ? number_of_ffts
                     : fft_size[0] *
                           (in_place ? 2 * (fft_size[1] / 2 + 1) : fft_size[1]),
        transpose_out ? number_of_ffts : fft_size[0] * (fft_size[1] / 2 + 1), 1,
        in_place, key);
    plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      double *buffer = fftw_alloc_real(2 * get_buffer_size_from_key(key));
      plan = fft_fftw_create_2d_plan_r2c(
          key, buffer, in_place ? buffer : (double *)grid_out);
      fftw_free(buffer);
    }
    if (block_size * number_of_full_blocks != number_of_ffts) {
      const int block_size_last_thread = number_of_ffts % block_size;
      fft_key_t key;
      get_key_2d_r2c(
          true, fft_size, block_size_last_thread, transpose_in, transpose_out,
          transpose_in ? number_of_ffts
                       : fft_size[0] * (in_place ? 2 * (fft_size[1] / 2 + 1)
                                                 : fft_size[1]),
          transpose_out ? number_of_ffts : fft_size[0] * (fft_size[1] / 2 + 1),
          1, in_place, key);
      plan_last_thread = lookup_plan_from_cache(key);
      if (plan_last_thread == NULL) {
        double *buffer = fftw_alloc_real(2 * get_buffer_size_from_key(key));
        plan_last_thread = fft_fftw_create_2d_plan_r2c(
            key, buffer, in_place ? buffer : (double *)grid_out);
        fftw_free(buffer);
      }
    }
    const int offset_in =
        transpose_in ? block_size
                     : block_size * fft_size[0] *
                           (in_place ? 2 * (fft_size[1] / 2 + 1) : fft_size[1]);
    const int offset_out =
        transpose_out ? block_size
                      : block_size * fft_size[0] * (fft_size[1] / 2 + 1);
#pragma omp parallel default(none)                                             \
    shared(number_of_full_blocks, plan, grid_in, grid_out, offset_in,          \
               offset_out, plan_last_thread, total_number_of_blocks)
    {
      const int thread_id = omp_get_thread_num();
      if (thread_id < number_of_full_blocks) {
        fftw_execute_dft_r2c(*plan, grid_in + thread_id * offset_in,
                             grid_out + thread_id * offset_out);
      } else if (thread_id + 1 == total_number_of_blocks) {
        fftw_execute_dft_r2c(*plan_last_thread, grid_in + thread_id * offset_in,
                             grid_out + thread_id * offset_out);
      }
    }
  } else {
    fft_key_t key;
    get_key_2d_r2c(
        true, fft_size, number_of_ffts, transpose_in, transpose_out,
        transpose_in ? number_of_ffts
                     : fft_size[0] *
                           (in_place ? 2 * (fft_size[1] / 2 + 1) : fft_size[1]),
        transpose_out ? number_of_ffts : fft_size[0] * (fft_size[1] / 2 + 1),
        fft_team_size(), grid_in == (double *)grid_out, key);
    fftw_plan *plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      double *buffer = fftw_alloc_real(2 * get_buffer_size_from_key(key));
      plan = fft_fftw_create_2d_plan_r2c(
          key, buffer,
          grid_in == (double *)grid_out ? buffer : (double *)grid_out);
      fftw_free(buffer);
    }
    fftw_execute_dft_r2c(*plan, grid_in, grid_out);
  }
#else
  (void)fft_size;
  (void)number_of_ffts;
  (void)grid_in;
  (void)grid_out;
  (void)transpose_in;
  (void)transpose_out;
  cp_assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local backwards C2R 2D FFT.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_2d_c2r_local(const int fft_size[2], const int number_of_ffts,
                           const bool transpose_in, const bool transpose_out,
                           double complex *grid_in, double *grid_out) {
#if defined(__FFTW3)
  cp_assert(is_initialized);
  if (fft_size[0] == 0 || fft_size[1] == 0 || number_of_ffts == 0)
    return;
  const bool in_place = (double *)grid_in == grid_out;
  const int number_of_threads = fft_team_size();
  // The self-written plans do not work with in-place transforms using the
  // transposed format
  if (fftw_planning_mode == FFTW_ESTIMATE && (!in_place || !transpose_in)) {
    fftw_plan *plan = NULL, *plan_last_thread = NULL;
    const int block_size =
        (number_of_ffts + number_of_threads - 1) / number_of_threads;
    const int number_of_full_blocks = number_of_ffts / block_size;
    const int total_number_of_blocks =
        number_of_full_blocks + (number_of_ffts % block_size != 0);
    fft_key_t key;
    get_key_2d_r2c(
        false, fft_size, block_size, transpose_in, transpose_out,
        transpose_in ? number_of_ffts : fft_size[0] * (fft_size[1] / 2 + 1),
        transpose_out ? number_of_ffts
                      : fft_size[0] * (in_place ? 2 * (fft_size[1] / 2 + 1)
                                                : fft_size[1]),
        1, in_place, key);
    plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      double *buffer = fftw_alloc_real(2 * get_buffer_size_from_key(key));
      plan = fft_fftw_create_2d_plan_r2c(key, buffer,
                                         in_place ? buffer : grid_out);
      fftw_free(buffer);
    }
    if (block_size * number_of_full_blocks != number_of_ffts) {
      const int block_size_last_thread = number_of_ffts % block_size;
      fft_key_t key;
      get_key_2d_r2c(
          false, fft_size, block_size_last_thread, transpose_in, transpose_out,
          transpose_in ? number_of_ffts : fft_size[0] * (fft_size[1] / 2 + 1),
          transpose_out ? number_of_ffts
                        : fft_size[0] * (in_place ? 2 * (fft_size[1] / 2 + 1)
                                                  : fft_size[1]),
          1, in_place, key);
      plan_last_thread = lookup_plan_from_cache(key);
      if (plan_last_thread == NULL) {
        double *buffer = fftw_alloc_real(2 * get_buffer_size_from_key(key));
        plan_last_thread = fft_fftw_create_2d_plan_r2c(
            key, buffer, in_place ? buffer : grid_out);
        fftw_free(buffer);
      }
    }
    const int offset_in =
        transpose_in ? block_size
                     : block_size * fft_size[0] * (fft_size[1] / 2 + 1);
    const int offset_out =
        transpose_out
            ? block_size
            : block_size * fft_size[0] *
                  (in_place ? 2 * (fft_size[1] / 2 + 1) : fft_size[1]);
#pragma omp parallel default(none)                                             \
    shared(number_of_full_blocks, plan, grid_in, grid_out, offset_in,          \
               offset_out, plan_last_thread, total_number_of_blocks)
    {
      const int thread_id = omp_get_thread_num();
      if (thread_id < number_of_full_blocks) {
        fftw_execute_dft_c2r(*plan, grid_in + thread_id * offset_in,
                             grid_out + thread_id * offset_out);
      } else if (thread_id + 1 == total_number_of_blocks) {
        fftw_execute_dft_c2r(*plan_last_thread, grid_in + thread_id * offset_in,
                             grid_out + thread_id * offset_out);
      }
    }
  } else {
    fft_key_t key;
    get_key_2d_r2c(
        false, fft_size, number_of_ffts, transpose_in, transpose_out,
        transpose_in ? number_of_ffts : fft_size[0] * (fft_size[1] / 2 + 1),
        transpose_out ? number_of_ffts
                      : fft_size[0] * (in_place ? 2 * (fft_size[1] / 2 + 1)
                                                : fft_size[1]),
        fft_team_size(), grid_in == (double complex *)grid_out, key);
    fftw_plan *plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      double complex *buffer =
          fftw_alloc_complex(get_buffer_size_from_key(key));
      plan = fft_fftw_create_2d_plan_r2c(
          key, (double *)buffer,
          (double *)grid_in == grid_out ? (double *)buffer : grid_out);
      fftw_free(buffer);
    }
    fftw_execute_dft_c2r(*plan, grid_in, grid_out);
  }
#else
  (void)fft_size;
  (void)number_of_ffts;
  (void)grid_in;
  (void)grid_out;
  (void)transpose_in;
  (void)transpose_out;
  cp_assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local C2C 3D FFT.
 * \param direction FFT direction (forward or backward).
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_3d_local(const bool direction, const int fft_size[3],
                       const int number_of_ffts, const bool transpose_in,
                       const bool transpose_out, double complex *grid_in,
                       double complex *grid_out) {
#if defined(__FFTW3)
  cp_assert(is_initialized);
  if (fft_size[0] == 0 || fft_size[1] == 0 || fft_size[2] == 0 ||
      number_of_ffts == 0)
    return;
  const bool in_place = grid_in == grid_out;
  if (number_of_ffts > 1 &&
      use_threaded_3d(fft_size, transpose_in, transpose_out, in_place)) {
    execute_threaded_3d_c2c(direction, fft_size, number_of_ffts, in_place,
                            grid_in, grid_out);
  } else if (number_of_ffts == 1 &&
             use_threaded_3d(fft_size, transpose_in, transpose_out, in_place)) {
    execute_threaded_3d_c2c_1d(direction, fft_size, in_place, grid_in,
                               grid_out);
  } else {
    fft_key_t key;
    get_key_3d(direction, fft_size, number_of_ffts, transpose_in, transpose_out,
               fft_team_size(), in_place, key);
    fftw_plan *plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      double complex *buffer =
          fftw_alloc_complex(get_buffer_size_from_key(key));
      plan = fft_fftw_create_3d_plan(key, buffer, in_place ? buffer : grid_out);
      fftw_free(buffer);
    }
    cp_assert(plan != NULL);
    cp_assert(grid_in != NULL);
    cp_assert(grid_out != NULL);
    fftw_execute_dft(*plan, grid_in, grid_out);
  }
#else
  (void)fft_size;
  (void)number_of_ffts;
  (void)transpose_in;
  (void)transpose_out;
  (void)grid_in;
  (void)grid_out;
  cp_assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local forward R2C 3D FFT.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_3d_r2c_local(const int fft_size[3], const int number_of_ffts,
                           const bool transpose_in, const bool transpose_out,
                           double *grid_in, double complex *grid_out) {
#if defined(__FFTW3)
  cp_assert(is_initialized);
  if (fft_size[0] == 0 || fft_size[1] == 0 || fft_size[2] == 0 ||
      number_of_ffts == 0)
    return;
  const bool in_place = grid_in == (double *)grid_out;
  if (number_of_ffts > 1 &&
      use_threaded_3d(fft_size, transpose_in, transpose_out, in_place)) {
    execute_threaded_3d_r2c(true, fft_size, number_of_ffts, in_place, grid_in,
                            (double *)grid_out);
  } else if (number_of_ffts == 1 &&
             use_threaded_3d(fft_size, transpose_in, transpose_out, in_place)) {
    execute_threaded_3d_r2c_1d(true, fft_size, in_place, grid_in, grid_out);
  } else {
    fft_key_t key;
    get_key_3d_r2c(true, fft_size, number_of_ffts, transpose_in, transpose_out,
                   fft_team_size(), in_place, key);
    fftw_plan *plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      double *buffer = fftw_alloc_real(2 * get_buffer_size_from_key(key));
      plan = fft_fftw_create_3d_plan_r2c(
          key, buffer, in_place ? buffer : (double *)grid_out);
      fftw_free(buffer);
    }
    cp_assert(plan != NULL);
    cp_assert(grid_in != NULL);
    cp_assert(grid_out != NULL);
    fftw_execute_dft_r2c(*plan, grid_in, grid_out);
  }
#else
  (void)fft_size;
  (void)number_of_ffts;
  (void)transpose_in;
  (void)transpose_out;
  (void)grid_in;
  (void)grid_out;
  cp_assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local backwards R2C 3D FFT.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_3d_c2r_local(const int fft_size[3], const int number_of_ffts,
                           const bool transpose_in, const bool transpose_out,
                           double complex *grid_in, double *grid_out) {
#if defined(__FFTW3)
  cp_assert(is_initialized);
  if (fft_size[0] == 0 || fft_size[1] == 0 || fft_size[2] == 0 ||
      number_of_ffts == 0)
    return;
  const bool in_place = (double *)grid_in == grid_out;
  if (number_of_ffts > 1 &&
      use_threaded_3d(fft_size, transpose_in, transpose_out, in_place)) {
    execute_threaded_3d_r2c(false, fft_size, number_of_ffts, in_place,
                            (double *)grid_in, grid_out);
  } else if (number_of_ffts == 1 &&
             use_threaded_3d(fft_size, transpose_in, transpose_out, in_place)) {
    execute_threaded_3d_r2c_1d(false, fft_size, in_place, grid_out, grid_in);
  } else {
    fft_key_t key;
    get_key_3d_r2c(false, fft_size, number_of_ffts, transpose_in, transpose_out,
                   fft_team_size(), in_place, key);
    fftw_plan *plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      double complex *buffer =
          fftw_alloc_complex(get_buffer_size_from_key(key));
      plan = fft_fftw_create_3d_plan_r2c(
          key, (double *)buffer, in_place ? (double *)buffer : grid_out);
      fftw_free(buffer);
    }
    cp_assert(plan != NULL);
    cp_assert(grid_in != NULL);
    cp_assert(grid_out != NULL);
    fftw_execute_dft_c2r(*plan, grid_in, grid_out);
  }
#else
  (void)fft_size;
  (void)number_of_ffts;
  (void)transpose_in;
  (void)transpose_out;
  (void)grid_in;
  (void)grid_out;
  cp_assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

// EOF
