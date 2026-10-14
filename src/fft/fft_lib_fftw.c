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

#include <assert.h>
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
// (local/distributed, C2C/R2C) This works to run all tests
#define FFTW_CACHE_SIZE 256
static cache_entry cache[FFTW_CACHE_SIZE];
static int cache_oldest_entry = 0; // used for LRU eviction

static bool is_initialized = false;

static int fftw_planning_mode = -1;

/*******************************************************************************
 * \brief Fetches an fft plan from the cache. Returns NULL if not found.
 * \author Ole Schuett, Frederick Stein
 ******************************************************************************/
static fftw_plan *lookup_plan_from_cache(const fft_key_t key) {
  assert(is_initialized);
  for (int i = 0; i < FFTW_CACHE_SIZE; i++) {
    const int *x = cache[i].key;
    if (memcmp(key, x, sizeof(fft_key_t)) == 0) {
      return cache[i].plan;
    }
  }
  return NULL;
}

/*******************************************************************************
 * \brief Adds an fft plan to the cache. Assumes ownership of plan's memory.
 * \author Ole Schuett, Frederick Stein
 ******************************************************************************/
static void add_plan_to_cache(const fft_key_t key, fftw_plan *plan) {
  const int i = cache_oldest_entry;
  cache_oldest_entry = (cache_oldest_entry + 1) % FFTW_CACHE_SIZE;
  if (cache[i].plan != NULL) {
    fftw_destroy_plan(*cache[i].plan);
    free(cache[i].plan);
  }
  memcpy(cache[i].key, key, sizeof(fft_key_t));
  cache[i].plan = plan;
}
#endif

static bool has_guru_interface = false;

#if defined(__FFTW3)
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

/*******************************************************************************
 * \brief Initialize the FFT library (if not done externally).
 * \author Frederick Stein
 ******************************************************************************/
int fft_fftw_init_lib(const fftw_plan_type fftw_planning_flag,
                      const bool use_guru_interface, const char *wisdom_file) {
#if defined(__FFTW3)
  assert(omp_get_num_threads() == 1);
  if (is_initialized) {
    return 0;
  }
  memset(cache, 0, sizeof(cache_entry) * FFTW_CACHE_SIZE);
  cache_oldest_entry = 0;

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
    assert(0 && "Unknown FFTW planning flag.");
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
  assert(false && "The grid library was not compiled with FFTW support.");
  (void)fftw_planning_flag;
  (void)use_guru_interface;
  (void)wisdom_file;
  return 0;
#endif
}

/*******************************************************************************
 * \brief Finalize the FFT library (if not done externally).
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_finalize_lib(const char *wisdom_file) {
#if defined(__FFTW3)
  assert(omp_get_num_threads() == 1);
  if (!is_initialized) {
    return;
  }
  for (int i = 0; i < FFTW_CACHE_SIZE; i++) {
    if (cache[i].plan != NULL) {
      fftw_destroy_plan(*cache[i].plan);
      free(cache[i].plan);
    }
  }
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
  assert(false && "The grid library was not compiled with FFTW support.");
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
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_allocate_double(const int length, double **buffer) {
#if defined(__FFTW3)
  assert(is_initialized);
  assert(buffer != NULL);
  assert(*buffer == NULL);
  *buffer = fftw_alloc_real(imax(1, length));
#else
  (void)length;
  (void)buffer;
  assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Allocate buffer of type double complex.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_allocate_complex(const int length, double complex **buffer) {
#if defined(__FFTW3)
  assert(is_initialized);
  assert(buffer != NULL);
  assert(*buffer == NULL);
  *buffer = fftw_alloc_complex(imax(1, length));
#else
  (void)length;
  (void)buffer;
  assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Free buffer of type double.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_free_double(double *buffer) {
#if defined(__FFTW3)
  assert(is_initialized);
  fftw_free(buffer);
#else
  (void)buffer;
  assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Free buffer of type double complex.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_free_complex(double complex *buffer) {
#if defined(__FFTW3)
  assert(is_initialized);
  fftw_free(buffer);
#else
  (void)buffer;
  assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

#if defined(__FFTW3)

/*******************************************************************************
 * \brief Create plan of a local C2C 1D FFT.
 * \author Frederick Stein
 ******************************************************************************/
fftw_plan *fft_fftw_create_1d_plan(const fft_key_t key, double complex *grid_in,
                                   double complex *grid_out) {
  int rank, number_of_threads, number_of_ffts, istride, ostride, idist, odist;
  bool direction, inplace;
  int fft_size[3], inembed[3], onembed[3];
  fetch_data_from_key_nd(key, &rank, fft_size, &number_of_ffts,
                         &number_of_threads, &direction, &inplace, inembed,
                         onembed, &idist, &odist, &istride, &ostride);
  assert(rank == 1);
  fftw_plan_with_nthreads(number_of_threads);
  fftw_plan *plan = malloc(sizeof(fftw_plan));
  if (grid_in == grid_out) {
    assert(inplace);
    assert(istride == ostride);
    assert(idist == odist);
  } else {
    assert(!inplace);
  }
  *plan = fftw_plan_many_dft(rank, fft_size, number_of_ffts, grid_in, inembed,
                             istride, idist, grid_out, onembed, ostride, odist,
                             direction ? FFTW_FORWARD : FFTW_BACKWARD,
                             fftw_planning_mode);
  assert(plan != NULL);
  add_plan_to_cache(key, plan);
  return plan;
}
/*******************************************************************************
 * \brief Create plan of a local R2C/C2R 1D FFT.
 * \author Frederick Stein
 ******************************************************************************/
fftw_plan *fft_fftw_create_1d_plan_r2c(const fft_key_t key, double *grid_in,
                                       double *grid_out) {
  int rank, number_of_threads, number_of_ffts, istride, ostride, idist, odist;
  bool direction, inplace;
  int fft_size[3], inembed[3], onembed[3];
  fetch_data_from_key_nd(key, &rank, fft_size, &number_of_ffts,
                         &number_of_threads, &direction, &inplace, inembed,
                         onembed, &idist, &odist, &istride, &ostride);
  assert(rank == 1);
  fftw_plan_with_nthreads(number_of_threads);
  fftw_plan *plan = malloc(sizeof(fftw_plan));
  if (grid_in == grid_out) {
    if (istride == 1) {
      // Check whether the arrays are correctly padded
      if (direction) {
        assert(idist == 2 * odist);
      } else {
        assert(odist == 2 * idist);
      }
    }
    assert(istride == ostride);
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
  assert(plan != NULL);
  add_plan_to_cache(key, plan);
  return plan;
}

/*******************************************************************************
 * \brief Create plan of a local C2C 2D FFT.
 * \author Frederick Stein
 ******************************************************************************/
fftw_plan *fft_fftw_create_2d_plan(const fft_key_t key, double complex *grid_in,
                                   double complex *grid_out) {
  int rank, number_of_threads, number_of_ffts, istride, ostride, idist, odist;
  bool direction, inplace;
  int fft_size[3], inembed[3], onembed[3];
  fetch_data_from_key_nd(key, &rank, fft_size, &number_of_ffts,
                         &number_of_threads, &direction, &inplace, inembed,
                         onembed, &idist, &odist, &istride, &ostride);
  assert(rank == 2);
  fftw_plan_with_nthreads(number_of_threads);
  fftw_plan *plan = malloc(sizeof(fftw_plan));
  if (grid_in == grid_out) {
    // Array strides need to match for in-place FFTs
    assert(istride == ostride);
    assert(idist == odist);
  }
  *plan = fftw_plan_many_dft(rank, fft_size, number_of_ffts, grid_in, inembed,
                             istride, idist, grid_out, onembed, ostride, odist,
                             direction ? FFTW_FORWARD : FFTW_BACKWARD,
                             fftw_planning_mode);
  assert(plan != NULL);
  add_plan_to_cache(key, plan);
  return plan;
}

/*******************************************************************************
 * \brief Create plan of a local R2C/C2R 2D FFT.
 * \author Frederick Stein
 ******************************************************************************/
fftw_plan *fft_fftw_create_2d_plan_r2c(const fft_key_t key, double *grid_in,
                                       double *grid_out) {
  int rank, number_of_threads, number_of_ffts, istride, ostride, idist, odist;
  bool direction, inplace;
  int fft_size[3], inembed[3], onembed[3];
  fetch_data_from_key_nd(key, &rank, fft_size, &number_of_ffts,
                         &number_of_threads, &direction, &inplace, inembed,
                         onembed, &idist, &odist, &istride, &ostride);
  assert(rank == 2);
  fftw_plan_with_nthreads(number_of_threads);
  // We need the guru interface here because cuts the last dimension in half
  // whereas we want the first dimension
  fftw_plan *plan = malloc(sizeof(fftw_plan));
  if (grid_in == grid_out) {
    if (istride == 1) {
      // Check whether the arrays are correctly padded
      if (direction) {
        assert(idist == 2 * odist);
      } else {
        assert(odist == 2 * idist);
      }
    }
    assert(istride == ostride);
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
  assert(plan != NULL);
  add_plan_to_cache(key, plan);
  return plan;
}

/*******************************************************************************
 * \brief Create plan of a local C2C 3D FFT.
 * \author Frederick Stein
 ******************************************************************************/
fftw_plan *fft_fftw_create_3d_plan(const fft_key_t key, double complex *grid_in,
                                   double complex *grid_out) {
  int rank, number_of_threads, number_of_ffts, istride, ostride, idist, odist;
  bool direction, inplace;
  int fft_size[3], inembed[3], onembed[3];
  fetch_data_from_key_nd(key, &rank, fft_size, &number_of_ffts,
                         &number_of_threads, &direction, &inplace, inembed,
                         onembed, &idist, &odist, &istride, &ostride);
  assert(rank == 3);
  assert(istride == 1);
  assert(ostride == 1);
  assert(number_of_ffts == 1);
  fftw_plan_with_nthreads(number_of_threads);
  fftw_plan *plan = malloc(sizeof(fftw_plan));
  *plan = fftw_plan_dft_3d(fft_size[0], fft_size[1], fft_size[2], grid_in,
                           grid_out, direction ? FFTW_FORWARD : FFTW_BACKWARD,
                           fftw_planning_mode);
  add_plan_to_cache(key, plan);
  assert(plan != NULL);
  return plan;
}

/*******************************************************************************
 * \brief Create plan of a local R2C/C2R 3D FFT.
 * \author Frederick Stein
 ******************************************************************************/
fftw_plan *fft_fftw_create_3d_plan_r2c(const fft_key_t key, double *grid_in,
                                       double *grid_out) {
  int rank, number_of_threads, number_of_ffts, istride, ostride, idist, odist;
  bool direction, inplace;
  int fft_size[3], inembed[3], onembed[3];
  fetch_data_from_key_nd(key, &rank, fft_size, &number_of_ffts,
                         &number_of_threads, &direction, &inplace, inembed,
                         onembed, &idist, &odist, &istride, &ostride);
  assert(rank == 3);
  assert(istride == 1);
  assert(ostride == 1);
  assert(number_of_ffts == 1);
  fftw_plan_with_nthreads(number_of_threads);
  fftw_plan *plan = malloc(sizeof(fftw_plan));
  if (direction) {
    *plan =
        fftw_plan_dft_r2c_3d(fft_size[0], fft_size[1], fft_size[2], grid_in,
                             (double complex *)grid_out, fftw_planning_mode);
  } else {
    // We use the buffers the other way around to prevent
    // out-of-bounds-accesses of the planner if the output array has only the
    // minimum size
    *plan = fftw_plan_dft_c2r_3d(fft_size[0], fft_size[1], fft_size[2],
                                 (double complex *)grid_in, grid_out,
                                 fftw_planning_mode);
  }
  add_plan_to_cache(key, plan);
  assert(plan != NULL);
  return plan;
}

/*******************************************************************************
 * \brief Create plan of a local C2C 3D FFT.
 * \author Frederick Stein
 ******************************************************************************/
fftw_plan *fft_fftw_create_guru_plan(const fft_key_t key,
                                     double complex *grid_in,
                                     double complex *grid_out) {
  assert(has_guru_interface);
  bool direction, inplace;
  int rank, howmany_rank, number_of_threads;
  fft_iodim dims[3], howmany_dims[3];
  fetch_data_from_key_guru(key, &direction, &rank, dims, &howmany_rank,
                           howmany_dims, &number_of_threads, &inplace);

  fftw_plan_with_nthreads(number_of_threads);
  fftw_plan *plan = malloc(sizeof(fftw_plan));
  if ((double complex *)grid_in == grid_out) {
    assert(inplace);
    for (int r = 0; r < rank; r++)
      assert(dims[r].is == dims[r].os);
    for (int r = 0; r < howmany_rank; r++)
      assert(howmany_dims[r].is == howmany_dims[r].os);
  } else {
    assert(!inplace);
  }
  *plan = fftw_plan_guru_dft(rank, dims, howmany_rank, howmany_dims, grid_in,
                             grid_out, direction ? FFTW_FORWARD : FFTW_BACKWARD,
                             fftw_planning_mode);
  add_plan_to_cache(key, plan);
  assert(plan != NULL);
  return plan;
}

/*******************************************************************************
 * \brief Create plan of a local R2C/C2R 3D FFT.
 * \author Frederick Stein
 ******************************************************************************/
fftw_plan *fft_fftw_create_guru_plan_r2c(const fft_key_t key, double *grid_in,
                                         double *grid_out) {
  assert(has_guru_interface);
  bool direction, inplace;
  int rank, howmany_rank, number_of_threads;
  fft_iodim dims[3], howmany_dims[3];
  fetch_data_from_key_guru(key, &direction, &rank, dims, &howmany_rank,
                           howmany_dims, &number_of_threads, &inplace);

  fftw_plan_with_nthreads(number_of_threads);
  fftw_plan *plan = malloc(sizeof(fftw_plan));
  if (direction) {
    if (grid_in == grid_out) {
      assert(inplace);
      for (int r = 0; r < rank; r++)
        assert(dims[r].is == 2 * dims[r].os ||
               (dims[r].is == 1 && dims[r].os == 1));
      for (int r = 0; r < howmany_rank; r++)
        assert(howmany_dims[r].is == 2 * howmany_dims[r].os);
    } else {
      assert(!inplace);
    }
    *plan =
        fftw_plan_guru_dft_r2c(rank, dims, howmany_rank, howmany_dims, grid_in,
                               (double complex *)grid_out, fftw_planning_mode);
  } else {
    if (grid_in == grid_out) {
      assert(inplace);
      for (int r = 0; r < rank; r++)
        assert(2 * dims[r].is == dims[r].os ||
               (dims[r].is == 1 && dims[r].os == 1));
      for (int r = 0; r < howmany_rank; r++)
        assert(2 * howmany_dims[r].is == howmany_dims[r].os);
    } else {
      assert(!inplace);
    }
    // We use the buffers the other way around to prevent
    // out-of-bounds-accesses of the planner if the output array has only the
    // minimum size
    *plan = fftw_plan_guru_dft_c2r(rank, dims, howmany_rank, howmany_dims,
                                   (double complex *)grid_in, grid_out,
                                   fftw_planning_mode);
  }
  add_plan_to_cache(key, plan);
  assert(plan != NULL);
  return plan;
}
#endif

/*******************************************************************************
 * \brief Performs a local C2C 1D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_1d_local(
    const bool direction, const int fft_size, const int number_of_ffts,
    const bool transpose_in, const bool transpose_out,
    const int leading_dimension_in, const int leading_dimension_out,
    double complex *grid_in, double complex *grid_out) {
#if defined(__FFTW3)
  assert(omp_get_num_threads() == 1);
  assert(is_initialized);
  if (fft_size == 0 || number_of_ffts == 0)
    return;
  const bool in_place = grid_in == grid_out;
  int number_of_threads = 1;
#pragma omp parallel default(none) shared(number_of_threads)
  {
#pragma omp single
    {
      number_of_threads = omp_get_num_threads();
    }
  }
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
               leading_dimension_in, leading_dimension_out,
               omp_get_max_threads(), in_place, key);
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
  assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local forward R2C FFT.
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
  assert(omp_get_num_threads() == 1);
  assert(is_initialized);
  if (fft_size == 0 || number_of_ffts == 0)
    return;
  const bool in_place = grid_in == (double *)grid_out;
  int number_of_threads = 1;
#pragma omp parallel default(none) shared(number_of_threads)
  {
#pragma omp single
    {
      number_of_threads = omp_get_num_threads();
    }
  }
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
                   leading_dimension_in, leading_dimension_out,
                   omp_get_max_threads(), in_place, key);
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
  assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local backwards C2R 1D FFT
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
  assert(omp_get_num_threads() == 1);
  assert(is_initialized);
  if (fft_size == 0 || number_of_ffts == 0)
    return;
  const bool in_place = (double *)grid_in == grid_out;
  int number_of_threads = 1;
#pragma omp parallel default(none) shared(number_of_threads)
  {
#pragma omp single
    {
      number_of_threads = omp_get_num_threads();
    }
  }
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
                   leading_dimension_in, leading_dimension_out,
                   omp_get_max_threads(), in_place, key);
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
  assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local C2C 2D FFT
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_2d_local(const bool direction, const int fft_size[2],
                                const int number_of_ffts,
                                const bool transpose_in,
                                const bool transpose_out,
                                double complex *grid_in,
                                double complex *grid_out) {
#if defined(__FFTW3)
  assert(omp_get_num_threads() == 1);
  assert(is_initialized);
  if (fft_size[0] == 0 || fft_size[1] == 0 || number_of_ffts == 0)
    return;
  const bool in_place = grid_in == grid_out;
  int number_of_threads = 1;
#pragma omp parallel default(none) shared(number_of_threads)
  {
#pragma omp single
    {
      number_of_threads = omp_get_num_threads();
    }
  }
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
               omp_get_max_threads(), grid_in == grid_out, key);
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
  assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local forward R2C 2D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_2d_r2c_local(const int fft_size[2],
                                    const int number_of_ffts,
                                    const bool transpose_in,
                                    const bool transpose_out, double *grid_in,
                                    double complex *grid_out) {
#if defined(__FFTW3)
  assert(omp_get_num_threads() == 1);
  assert(is_initialized);
  if (fft_size[0] == 0 || fft_size[1] == 0 || number_of_ffts == 0)
    return;
  const bool in_place = grid_in == (double *)grid_out;
  int number_of_threads = 1;
#pragma omp parallel default(none) shared(number_of_threads)
  {
#pragma omp single
    {
      number_of_threads = omp_get_num_threads();
    }
  }
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
        omp_get_max_threads(), grid_in == (double *)grid_out, key);
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
  assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local backwards C2R 2D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_2d_c2r_local(const int fft_size[2],
                                    const int number_of_ffts,
                                    const bool transpose_in,
                                    const bool transpose_out,
                                    double complex *grid_in, double *grid_out) {
#if defined(__FFTW3)
  assert(omp_get_num_threads() == 1);
  assert(is_initialized);
  if (fft_size[0] == 0 || fft_size[1] == 0 || number_of_ffts == 0)
    return;
  const bool in_place = (double *)grid_in == grid_out;
  int number_of_threads = 1;
#pragma omp parallel default(none) shared(number_of_threads)
  {
#pragma omp single
    {
      number_of_threads = omp_get_num_threads();
    }
  }
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
        omp_get_max_threads(), grid_in == (double complex *)grid_out, key);
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
  assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local C2C 3D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_3d_local(const bool direction, const int fft_size[3],
                                double complex *grid_in,
                                double complex *grid_out) {
#if defined(__FFTW3)
  assert(omp_get_num_threads() == 1);
  assert(is_initialized);
  if (fft_size[0] == 0 || fft_size[1] == 0 || fft_size[2] == 0)
    return;
  const bool in_place = grid_in == grid_out;
  if (has_guru_interface &&
      ((((fft_size[0] >= 256 || fft_size[1] >= 256 || fft_size[2] >= 256) &&
         !in_place) ||
        ((fft_size[0] >= 128 || fft_size[1] >= 128 || fft_size[2] >= 128) &&
         in_place)) &&
       omp_get_max_threads() > 1) &&
      (fftw_planning_mode == FFTW_ESTIMATE)) {
    // The 3D FFT is not efficient with threading and estimate planning mode
    // So, we decompose it in a sequence of 1D FFTs
    int number_of_threads = 1;
#pragma omp parallel default(none) shared(number_of_threads)
    {
#pragma omp single
      {
        number_of_threads = omp_get_num_threads();
      }
    }
    int block_sizes[3];
    fftw_plan *plans[3], *plans_last_thread[3];
    {
      const int number_of_ffts = fft_size[1] * fft_size[2];
      block_sizes[0] =
          (number_of_ffts + number_of_threads - 1) / number_of_threads;
      fft_iodim dim = {
          .n = fft_size[0], .is = number_of_ffts, .os = number_of_ffts};
      fft_iodim howmany_dim = {.n = block_sizes[0], .is = 1, .os = 1};
      fft_key_t key;
      get_key_guru(direction, 1, &dim, 1, &howmany_dim, 1, in_place, key);

      plans[0] = lookup_plan_from_cache(key);
      if (plans[0] == NULL) {
        plans[0] = fft_fftw_create_guru_plan(key, grid_in, grid_out);
      }
      if (block_sizes[0] * number_of_threads != number_of_ffts) {
        const int block_size_last_thread = number_of_ffts % block_sizes[0];
        howmany_dim.n = block_size_last_thread;
        fft_key_t key;
        get_key_guru(direction, 1, &dim, 1, &howmany_dim, 1, in_place, key);
        plans_last_thread[0] = lookup_plan_from_cache(key);
        if (plans_last_thread[0] == NULL) {
          plans_last_thread[0] =
              fft_fftw_create_guru_plan(key, grid_in, grid_out);
        }
      }
    }
    {
      const int number_of_ffts = fft_size[0];
      block_sizes[1] =
          (number_of_ffts + number_of_threads - 1) / number_of_threads;
      fft_iodim dim = {.n = fft_size[1], .is = fft_size[2], .os = fft_size[2]};
      fft_iodim howmany_dims[2] = {{.n = block_sizes[1],
                                    .is = fft_size[1] * fft_size[2],
                                    .os = fft_size[1] * fft_size[2]},
                                   {.n = fft_size[2], .is = 1, .os = 1}};
      fft_key_t key;
      get_key_guru(direction, 1, &dim, 2, howmany_dims, 1, in_place, key);
      plans[1] = lookup_plan_from_cache(key);
      if (plans[1] == NULL) {
        plans[1] = fft_fftw_create_guru_plan(key, grid_in, grid_out);
      }
      if (block_sizes[1] * number_of_threads != number_of_ffts) {
        const int block_size_last_thread = number_of_ffts % block_sizes[1];
        howmany_dims[0].n = block_size_last_thread;
        fft_key_t key;
        get_key_guru(direction, 1, &dim, 2, howmany_dims, 1, in_place, key);
        plans_last_thread[1] = lookup_plan_from_cache(key);
        if (plans_last_thread[1] == NULL) {
          plans_last_thread[1] =
              fft_fftw_create_guru_plan(key, grid_in, grid_out);
        }
      }
    }
    {
      const int number_of_ffts = fft_size[0] * fft_size[1];
      block_sizes[2] =
          (number_of_ffts + number_of_threads - 1) / number_of_threads;
      fft_iodim dim = {.n = fft_size[2], .is = 1, .os = 1};
      fft_iodim howmany_dim = {
          .n = block_sizes[2], .is = fft_size[2], .os = fft_size[2]};
      fft_key_t key;
      get_key_guru(direction, 1, &dim, 1, &howmany_dim, 1, in_place, key);
      plans[2] = lookup_plan_from_cache(key);
      if (plans[2] == NULL) {
        plans[2] = fft_fftw_create_guru_plan(key, grid_in, grid_out);
      }
      if (block_sizes[2] * number_of_threads != number_of_ffts) {
        const int block_size_last_thread = number_of_ffts % block_sizes[2];
        howmany_dim.n = block_size_last_thread;
        fft_key_t key;
        get_key_guru(direction, 1, &dim, 1, &howmany_dim, 1, in_place, key);

        plans_last_thread[2] = lookup_plan_from_cache(key);
        if (plans_last_thread[2] == NULL) {
          plans_last_thread[2] =
              fft_fftw_create_guru_plan(key, grid_in, grid_out);
        }
      }
    }
  } else {
    fft_key_t key;
    get_key_3d(direction, fft_size, omp_get_max_threads(), in_place, key);
    fftw_plan *plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      plan = fft_fftw_create_3d_plan(key, grid_in, grid_out);
    }
  }
#else
  (void)fft_size;
  (void)grid_in;
  (void)grid_out;
  assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local forward R2C 3D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_3d_r2c_local(const int fft_size[3], double *grid_in,
                                    double complex *grid_out) {
#if defined(__FFTW3)
  assert(omp_get_num_threads() == 1);
  assert(is_initialized);
  if (fft_size[0] == 0 || fft_size[1] == 0 || fft_size[2] == 0)
    return;
  const bool in_place = grid_in == (double *)grid_out;
  if (has_guru_interface &&
      ((((fft_size[0] >= 256 || fft_size[1] >= 256 || fft_size[2] >= 256) &&
         !in_place) ||
        ((fft_size[0] >= 128 || fft_size[1] >= 128 || fft_size[2] >= 128) &&
         in_place)) &&
       omp_get_max_threads() > 1) &&
      (fftw_planning_mode == FFTW_ESTIMATE)) {
    // We cannot use the input buffer for planning because it may be too small
    // So, we allocate a new one
    double complex *buffer =
        in_place ? (double complex *)grid_in
                 : fftw_alloc_complex(fft_size[0] * fft_size[1] *
                                      (fft_size[2] / 2 + 1));
    // The 3D FFT is not efficient with threading and estimate planning mode
    // So, we decompose it in a sequence of 1D FFTs
    int number_of_threads = 1;
#pragma omp parallel default(none) shared(number_of_threads)
    {
#pragma omp single
      {
        number_of_threads = omp_get_num_threads();
      }
    }
    int block_sizes[3];
    fftw_plan *plans[3], *plans_last_thread[3];
    {
      const int number_of_ffts = fft_size[1] * (fft_size[2] / 2 + 1);
      block_sizes[0] =
          (number_of_ffts + number_of_threads - 1) / number_of_threads;
      fft_iodim dim = {
          .n = fft_size[0], .is = number_of_ffts, .os = number_of_ffts};
      fft_iodim howmany_dim = {.n = block_sizes[0], .is = 1, .os = 1};
      fft_key_t key;
      get_key_guru(true, 1, &dim, 1, &howmany_dim, 1, in_place, key);

      plans[0] = lookup_plan_from_cache(key);
      if (plans[0] == NULL) {
        plans[0] = fft_fftw_create_guru_plan(key, buffer, grid_out);
      }
      if (block_sizes[0] * number_of_threads != number_of_ffts) {
        const int block_size_last_thread = number_of_ffts % block_sizes[0];
        howmany_dim.n = block_size_last_thread;
        fft_key_t key;
        get_key_guru(true, 1, &dim, 1, &howmany_dim, 1, in_place, key);

        plans_last_thread[0] = lookup_plan_from_cache(key);
        if (plans_last_thread[0] == NULL) {
          plans_last_thread[0] =
              fft_fftw_create_guru_plan(key, buffer, grid_out);
        }
      }
    }
    {
      const int number_of_ffts = fft_size[0];
      block_sizes[1] =
          (number_of_ffts + number_of_threads - 1) / number_of_threads;
      fft_iodim dim = {.n = fft_size[1],
                       .is = fft_size[2] / 2 + 1,
                       .os = fft_size[2] / 2 + 1};
      fft_iodim howmany_dims[2] = {
          {.n = block_sizes[1],
           .is = fft_size[1] * (fft_size[2] / 2 + 1),
           .os = fft_size[1] * (fft_size[2] / 2 + 1)},
          {.n = (fft_size[2] / 2 + 1), .is = 1, .os = 1}};
      fft_key_t key;
      get_key_guru(true, 1, &dim, 2, howmany_dims, 1, in_place, key);

      plans[1] = lookup_plan_from_cache(key);
      if (plans[1] == NULL) {
        plans[1] = fft_fftw_create_guru_plan(key, buffer, grid_out);
      }
      if (block_sizes[1] * number_of_threads != number_of_ffts) {
        const int block_size_last_thread = number_of_ffts % block_sizes[1];
        howmany_dims[0].n = block_size_last_thread;
        fft_key_t key;
        get_key_guru(true, 1, &dim, 2, howmany_dims, 1, in_place, key);

        plans_last_thread[1] = lookup_plan_from_cache(key);
        if (plans_last_thread[1] == NULL) {
          plans_last_thread[1] =
              fft_fftw_create_guru_plan(key, buffer, grid_out);
        }
      }
    }
    {
      const int number_of_ffts = fft_size[0] * fft_size[1];
      block_sizes[2] =
          (number_of_ffts + number_of_threads - 1) / number_of_threads;
      fft_iodim dim = {.n = fft_size[2], .is = 1, .os = 1};
      fft_iodim howmany_dim = {.n = block_sizes[2],
                               .is = in_place ? 2 * (fft_size[2] / 2 + 1)
                                              : fft_size[2],
                               .os = fft_size[2] / 2 + 1};
      fft_key_t key;
      get_key_guru_r2c(true, 1, &dim, 1, &howmany_dim, 1, in_place, key);

      plans[2] = lookup_plan_from_cache(key);
      if (plans[2] == NULL) {
        plans[2] =
            fft_fftw_create_guru_plan_r2c(key, grid_in, (double *)grid_out);
      }
      if (block_sizes[2] * number_of_threads != number_of_ffts) {
        const int block_size_last_thread = number_of_ffts % block_sizes[2];
        howmany_dim.n = block_size_last_thread;
        fft_key_t key;
        get_key_guru_r2c(true, 1, &dim, 1, &howmany_dim, 1, in_place, key);

        plans_last_thread[2] = lookup_plan_from_cache(key);
        if (plans_last_thread[2] == NULL) {
          plans_last_thread[2] =
              fft_fftw_create_guru_plan_r2c(key, grid_in, (double *)grid_out);
        }
      }
    }
    if (!in_place)
      fftw_free(buffer);
  } else {
    fft_key_t key;
    get_key_3d_r2c(true, fft_size, omp_get_max_threads(), in_place, key);
    fftw_plan *plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      plan = fft_fftw_create_3d_plan_r2c(key, grid_in, (double *)grid_out);
    }
  }
#else
  (void)fft_size;
  (void)grid_in;
  (void)grid_out;
  assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local backwards R2C 3D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_3d_c2r_local(const int fft_size[3],
                                    double complex *grid_in, double *grid_out) {
#if defined(__FFTW3)
  assert(omp_get_num_threads() == 1);
  assert(is_initialized);
  if (fft_size[0] == 0 || fft_size[1] == 0 || fft_size[2] == 0)
    return;
  const bool in_place = (double *)grid_in == grid_out;
  if (has_guru_interface &&
      ((((fft_size[0] >= 256 || fft_size[1] >= 256 || fft_size[2] >= 256) &&
         !in_place) ||
        ((fft_size[0] >= 128 || fft_size[1] >= 128 || fft_size[2] >= 128) &&
         in_place)) &&
       omp_get_max_threads() > 1) &&
      (fftw_planning_mode == FFTW_ESTIMATE)) {
    // We cannot use the output buffer for planning because it may be too small
    // So, we allocate a new one
    double complex *buffer =
        in_place ? grid_in
                 : fftw_alloc_complex(fft_size[0] * fft_size[1] *
                                      (fft_size[2] / 2 + 1));
    // The 3D FFT is not efficient with threading and estimate planning mode
    // So, we decompose it in a sequence of 1D FFTs
    int number_of_threads = 1;
#pragma omp parallel default(none) shared(number_of_threads)
    {
#pragma omp single
      {
        number_of_threads = omp_get_num_threads();
      }
    }
    int block_sizes[3];
    fftw_plan *plans[3], *plans_last_thread[3];
    {
      const int number_of_ffts = fft_size[1] * (fft_size[2] / 2 + 1);
      block_sizes[0] =
          (number_of_ffts + number_of_threads - 1) / number_of_threads;
      fft_iodim dim = {
          .n = fft_size[0], .is = number_of_ffts, .os = number_of_ffts};
      fft_iodim howmany_dim = {.n = block_sizes[0], .is = 1, .os = 1};
      fft_key_t key;
      get_key_guru(false, 1, &dim, 1, &howmany_dim, 1, in_place, key);

      plans[0] = lookup_plan_from_cache(key);
      if (plans[0] == NULL) {
        plans[0] = fft_fftw_create_guru_plan(key, grid_in, buffer);
      }
      if (block_sizes[0] * number_of_threads != number_of_ffts) {
        const int block_size_last_thread = number_of_ffts % block_sizes[0];
        howmany_dim.n = block_size_last_thread;
        fft_key_t key;
        get_key_guru(false, 1, &dim, 1, &howmany_dim, 1, in_place, key);

        plans_last_thread[0] = lookup_plan_from_cache(key);
        if (plans_last_thread[0] == NULL) {
          plans_last_thread[0] =
              fft_fftw_create_guru_plan(key, grid_in, buffer);
        }
      }
    }
    {
      const int number_of_ffts = fft_size[0];
      block_sizes[1] =
          (number_of_ffts + number_of_threads - 1) / number_of_threads;
      fft_iodim dim = {.n = fft_size[1],
                       .is = fft_size[2] / 2 + 1,
                       .os = fft_size[2] / 2 + 1};
      fft_iodim howmany_dims[2] = {
          {.n = block_sizes[1],
           .is = fft_size[1] * (fft_size[2] / 2 + 1),
           .os = fft_size[1] * (fft_size[2] / 2 + 1)},
          {.n = (fft_size[2] / 2 + 1), .is = 1, .os = 1}};
      fft_key_t key;
      get_key_guru(false, 1, &dim, 2, howmany_dims, 1, in_place, key);

      plans[1] = lookup_plan_from_cache(key);
      if (plans[1] == NULL) {
        plans[1] = fft_fftw_create_guru_plan(key, grid_in, buffer);
      }
      if (block_sizes[1] * number_of_threads != number_of_ffts) {
        const int block_size_last_thread = number_of_ffts % block_sizes[1];
        howmany_dims[0].n = block_size_last_thread;
        fft_key_t key;
        get_key_guru(false, 1, &dim, 2, howmany_dims, 1, in_place, key);

        plans_last_thread[1] = lookup_plan_from_cache(key);
        if (plans_last_thread[1] == NULL) {
          plans_last_thread[1] =
              fft_fftw_create_guru_plan(key, grid_in, buffer);
        }
      }
    }
    {
      const int number_of_ffts = fft_size[0] * fft_size[1];
      block_sizes[2] =
          (number_of_ffts + number_of_threads - 1) / number_of_threads;
      fft_iodim dim = {.n = fft_size[2], .is = 1, .os = 1};
      fft_iodim howmany_dim = {.n = block_sizes[2],
                               .is = fft_size[2] / 2 + 1,
                               .os = in_place ? 2 * (fft_size[2] / 2 + 1)
                                              : fft_size[2]};
      fft_key_t key;
      get_key_guru_r2c(false, 1, &dim, 1, &howmany_dim, 1, in_place, key);

      plans[2] = lookup_plan_from_cache(key);
      if (plans[2] == NULL) {
        plans[2] =
            fft_fftw_create_guru_plan_r2c(key, (double *)grid_in, grid_out);
      }
      if (block_sizes[2] * number_of_threads != number_of_ffts) {
        const int block_size_last_thread = number_of_ffts % block_sizes[2];
        howmany_dim.n = block_size_last_thread;
        fft_key_t key;
        get_key_guru_r2c(false, 1, &dim, 1, &howmany_dim, 1, in_place, key);

        plans_last_thread[2] = lookup_plan_from_cache(key);
        if (plans_last_thread[2] == NULL) {
          plans_last_thread[2] =
              fft_fftw_create_guru_plan_r2c(key, (double *)grid_in, grid_out);
        }
      }
    }
    if (!in_place)
      fftw_free(buffer);
  } else {
    fft_key_t key;
    get_key_3d_r2c(false, fft_size, omp_get_max_threads(), in_place, key);
    fftw_plan *plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      plan = fft_fftw_create_3d_plan_r2c(key, (double *)grid_in, grid_out);
    }
  }
#else
  (void)fft_size;
  (void)grid_in;
  (void)grid_out;
  assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local forward C2C 1D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_1d_local(const bool direction, const int fft_size,
                       const int number_of_ffts, const bool transpose_in,
                       const bool transpose_out, const int leading_dimension_in,
                       const int leading_dimension_out, double complex *grid_in,
                       double complex *grid_out) {
#if defined(__FFTW3)
  assert(omp_get_num_threads() == 1);
  assert(is_initialized);
  if (fft_size == 0 || number_of_ffts == 0)
    return;
  const bool in_place = grid_in == grid_out;
  int number_of_threads = 1;
#pragma omp parallel default(none) shared(number_of_threads)
  {
#pragma omp single
    {
      number_of_threads = omp_get_num_threads();
    }
  }
  if (fftw_planning_mode == FFTW_ESTIMATE) {
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
      fft_key_t key;
      get_key_1d(direction, fft_size, block_size_last_thread, transpose_in,
                 transpose_out, leading_dimension_in, leading_dimension_out, 1,
                 in_place, key);
      plan_last_thread = lookup_plan_from_cache(key);
      if (plan_last_thread == NULL) {
        double complex *buffer =
            fftw_alloc_complex(get_buffer_size_from_key(key));
        plan_last_thread =
            fft_fftw_create_1d_plan(key, buffer, in_place ? buffer : grid_out);
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
               leading_dimension_in, leading_dimension_out,
               omp_get_max_threads(), in_place, key);
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
  assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local forward R2C FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_1d_r2c_local(const int fft_size, const int number_of_ffts,
                           const bool transpose_in, const bool transpose_out,
                           const int leading_dimension_in,
                           const int leading_dimension_out, double *grid_in,
                           double complex *grid_out) {
#if defined(__FFTW3)
  assert(omp_get_num_threads() == 1);
  assert(is_initialized);
  if (fft_size == 0 || number_of_ffts == 0)
    return;
  const bool in_place = grid_in == (double *)grid_out;
  int number_of_threads = 1;
#pragma omp parallel default(none) shared(number_of_threads)
  {
#pragma omp single
    {
      number_of_threads = omp_get_num_threads();
    }
  }
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
                   leading_dimension_in, leading_dimension_out,
                   omp_get_max_threads(), in_place, key);
    fftw_plan *plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      double *buffer = fftw_alloc_real(2 * get_buffer_size_from_key(key));
      plan = fft_fftw_create_1d_plan_r2c(
          key, buffer, in_place ? buffer : (double *)grid_out);
      fftw_free(buffer);
    }
    assert(plan != NULL);
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
  assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local backwards C2R 1D FFT
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_1d_c2r_local(const int fft_size, const int number_of_ffts,
                           const bool transpose_in, const bool transpose_out,
                           const int leading_dimension_in,
                           const int leading_dimension_out,
                           double complex *grid_in, double *grid_out) {
#if defined(__FFTW3)
  assert(omp_get_num_threads() == 1);
  assert(is_initialized);
  if (fft_size == 0 || number_of_ffts == 0)
    return;
  const bool in_place = (double *)grid_in == grid_out;
  int number_of_threads = 1;
#pragma omp parallel default(none) shared(number_of_threads)
  {
#pragma omp single
    {
      number_of_threads = omp_get_num_threads();
    }
  }
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
                   leading_dimension_in, leading_dimension_out,
                   omp_get_max_threads(), grid_in == (double complex *)grid_out,
                   key);
    fftw_plan *plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      double *buffer = fftw_alloc_real(2 * get_buffer_size_from_key(key));
      plan = fft_fftw_create_1d_plan_r2c(key, buffer,
                                         in_place ? buffer : grid_out);
      fftw_free(buffer);
    }
    assert(plan != NULL);
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
  assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local C2C 2D FFT
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_2d_local(const bool direction, const int fft_size[2],
                       const int number_of_ffts, const bool transpose_in,
                       const bool transpose_out, double complex *grid_in,
                       double complex *grid_out) {
#if defined(__FFTW3)
  assert(omp_get_num_threads() == 1);
  assert(is_initialized);
  if (fft_size[0] == 0 || fft_size[1] == 0 || number_of_ffts == 0)
    return;
  const bool in_place = grid_in == grid_out;
  int number_of_threads = 1;
#pragma omp parallel default(none) shared(number_of_threads)
  {
#pragma omp single
    {
      number_of_threads = omp_get_num_threads();
    }
  }
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
               omp_get_max_threads(), grid_in == grid_out, key);
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
  assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local forward R2C 2D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_2d_r2c_local(const int fft_size[2], const int number_of_ffts,
                           const bool transpose_in, const bool transpose_out,
                           double *grid_in, double complex *grid_out) {
#if defined(__FFTW3)
  assert(omp_get_num_threads() == 1);
  assert(is_initialized);
  if (fft_size[0] == 0 || fft_size[1] == 0 || number_of_ffts == 0)
    return;
  const bool in_place = grid_in == (double *)grid_out;
  int number_of_threads = 1;
#pragma omp parallel default(none) shared(number_of_threads)
  {
#pragma omp single
    {
      number_of_threads = omp_get_num_threads();
    }
  }
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
        omp_get_max_threads(), grid_in == (double *)grid_out, key);
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
  assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local backwards C2R 2D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_2d_c2r_local(const int fft_size[2], const int number_of_ffts,
                           const bool transpose_in, const bool transpose_out,
                           double complex *grid_in, double *grid_out) {
#if defined(__FFTW3)
  assert(omp_get_num_threads() == 1);
  assert(is_initialized);
  if (fft_size[0] == 0 || fft_size[1] == 0 || number_of_ffts == 0)
    return;
  const bool in_place = (double *)grid_in == grid_out;
  int number_of_threads = 1;
#pragma omp parallel default(none) shared(number_of_threads)
  {
#pragma omp single
    {
      number_of_threads = omp_get_num_threads();
    }
  }
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
        omp_get_max_threads(), grid_in == (double complex *)grid_out, key);
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
  assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local C2C 3D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_3d_local(const bool direction, const int fft_size[3],
                       double complex *grid_in, double complex *grid_out) {
#if defined(__FFTW3)
  assert(omp_get_num_threads() == 1);
  assert(is_initialized);
  if (fft_size[0] == 0 || fft_size[1] == 0 || fft_size[2] == 0)
    return;
  const bool in_place = grid_in == grid_out;
  if (has_guru_interface &&
      ((((fft_size[0] >= 256 || fft_size[1] >= 256 || fft_size[2] >= 256) &&
         !in_place) ||
        ((fft_size[0] >= 128 || fft_size[1] >= 128 || fft_size[2] >= 128) &&
         in_place)) &&
       omp_get_max_threads() > 1) &&
      (fftw_planning_mode == FFTW_ESTIMATE)) {
    // The 3D FFT is not efficient with threading and estimate planning mode
    // So, we decompose it in a sequence of 1D FFTs
    int number_of_threads = 1;
#pragma omp parallel default(none) shared(number_of_threads)
    {
#pragma omp single
      {
        number_of_threads = omp_get_num_threads();
      }
    }
    int block_sizes[3];
    fftw_plan *plans[3], *plans_last_thread[3];
    int number_of_full_blocks[3] = {0, 0, 0};
    int total_number_of_blocks[3] = {0, 0, 0};
    {
      const int number_of_ffts = fft_size[1] * fft_size[2];
      block_sizes[0] =
          (number_of_ffts + number_of_threads - 1) / number_of_threads;
      number_of_full_blocks[0] = number_of_ffts / block_sizes[0];
      total_number_of_blocks[0] =
          (number_of_ffts + block_sizes[0] - 1) / block_sizes[0];
      fft_iodim dim = {
          .n = fft_size[0], .is = number_of_ffts, .os = number_of_ffts};
      fft_iodim howmany_dim = {.n = block_sizes[0], .is = 1, .os = 1};
      fft_key_t key;
      get_key_guru(direction, 1, &dim, 1, &howmany_dim, 1, in_place, key);

      plans[0] = lookup_plan_from_cache(key);
      if (plans[0] == NULL) {
        double complex *buffer =
            fftw_alloc_complex(get_buffer_size_from_key(key));
        plans[0] = fft_fftw_create_guru_plan(key, buffer,
                                             in_place ? buffer : grid_out);
        fftw_free(buffer);
      }
      if (block_sizes[0] * number_of_threads != number_of_ffts) {
        const int block_size_last_thread = number_of_ffts % block_sizes[0];
        howmany_dim.n = block_size_last_thread;
        fft_key_t key;
        get_key_guru(direction, 1, &dim, 1, &howmany_dim, 1, in_place, key);
        plans_last_thread[0] = lookup_plan_from_cache(key);
        if (plans_last_thread[0] == NULL) {
          double complex *buffer =
              fftw_alloc_complex(get_buffer_size_from_key(key));
          plans_last_thread[0] = fft_fftw_create_guru_plan(
              key, buffer, in_place ? buffer : grid_out);
          fftw_free(buffer);
        }
      }
    }
    {
      const int number_of_ffts = fft_size[0];
      block_sizes[1] =
          (number_of_ffts + number_of_threads - 1) / number_of_threads;
      number_of_full_blocks[1] = number_of_ffts / block_sizes[1];
      total_number_of_blocks[1] =
          (number_of_ffts + block_sizes[1] - 1) / block_sizes[1];
      fft_iodim dim = {.n = fft_size[1], .is = fft_size[2], .os = fft_size[2]};
      fft_iodim howmany_dims[2] = {{.n = block_sizes[1],
                                    .is = fft_size[1] * fft_size[2],
                                    .os = fft_size[1] * fft_size[2]},
                                   {.n = fft_size[2], .is = 1, .os = 1}};
      fft_key_t key;
      get_key_guru(direction, 1, &dim, 2, howmany_dims, 1, in_place, key);
      plans[1] = lookup_plan_from_cache(key);
      if (plans[1] == NULL) {
        double complex *buffer =
            fftw_alloc_complex(get_buffer_size_from_key(key));
        plans[1] = fft_fftw_create_guru_plan(key, buffer,
                                             in_place ? buffer : grid_out);
        fftw_free(buffer);
      }
      if (block_sizes[1] * number_of_threads != number_of_ffts) {
        const int block_size_last_thread = number_of_ffts % block_sizes[1];
        howmany_dims[0].n = block_size_last_thread;
        fft_key_t key;
        get_key_guru(direction, 1, &dim, 2, howmany_dims, 1, in_place, key);
        plans_last_thread[1] = lookup_plan_from_cache(key);
        if (plans_last_thread[1] == NULL) {
          double complex *buffer =
              fftw_alloc_complex(get_buffer_size_from_key(key));
          plans_last_thread[1] = fft_fftw_create_guru_plan(
              key, buffer, in_place ? buffer : grid_out);
          fftw_free(buffer);
        }
      }
    }
    {
      const int number_of_ffts = fft_size[0] * fft_size[1];
      block_sizes[2] =
          (number_of_ffts + number_of_threads - 1) / number_of_threads;
      number_of_full_blocks[2] = number_of_ffts / block_sizes[2];
      total_number_of_blocks[2] =
          (number_of_ffts + block_sizes[2] - 1) / block_sizes[2];
      fft_iodim dim = {.n = fft_size[2], .is = 1, .os = 1};
      fft_iodim howmany_dim = {
          .n = block_sizes[2], .is = fft_size[2], .os = fft_size[2]};
      fft_key_t key;
      get_key_guru(direction, 1, &dim, 1, &howmany_dim, 1, in_place, key);
      plans[2] = lookup_plan_from_cache(key);
      if (plans[2] == NULL) {
        double complex *buffer =
            fftw_alloc_complex(get_buffer_size_from_key(key));
        plans[2] = fft_fftw_create_guru_plan(key, buffer,
                                             in_place ? buffer : grid_out);
        fftw_free(buffer);
      }
      if (block_sizes[2] * number_of_threads != number_of_ffts) {
        const int block_size_last_thread = number_of_ffts % block_sizes[2];
        howmany_dim.n = block_size_last_thread;
        fft_key_t key;
        get_key_guru(direction, 1, &dim, 1, &howmany_dim, 1, in_place, key);

        plans_last_thread[2] = lookup_plan_from_cache(key);
        if (plans_last_thread[2] == NULL) {
          double complex *buffer =
              fftw_alloc_complex(get_buffer_size_from_key(key));
          plans_last_thread[2] = fft_fftw_create_guru_plan(
              key, buffer, in_place ? buffer : grid_out);
          fftw_free(buffer);
        }
      }
    }
#pragma omp parallel default(none)                                             \
    shared(number_of_threads, number_of_full_blocks, total_number_of_blocks,   \
               plans, plans_last_thread, block_sizes, fft_size, grid_in,       \
               grid_out)
    {
      const int thread_id = omp_get_thread_num();
      if (thread_id < number_of_full_blocks[2]) {
        assert(*plans[2] != NULL);
        fftw_execute_dft(*plans[2],
                         grid_in + block_sizes[2] * fft_size[2] * thread_id,
                         grid_out + block_sizes[2] * fft_size[2] * thread_id);
      } else if (thread_id == total_number_of_blocks[2] - 1) {
        assert(*plans_last_thread[2] != NULL);
        fftw_execute_dft(*plans_last_thread[2],
                         grid_in + block_sizes[2] * fft_size[2] * thread_id,
                         grid_out + block_sizes[2] * fft_size[2] * thread_id);
      }
#pragma omp barrier
      if (thread_id < number_of_full_blocks[1]) {
        assert(*plans[1] != NULL);
        fftw_execute_dft(
            *plans[1],
            grid_out + block_sizes[1] * fft_size[1] * fft_size[2] * thread_id,
            grid_in + block_sizes[1] * fft_size[1] * fft_size[2] * thread_id);
      } else if (thread_id == total_number_of_blocks[1] - 1) {
        assert(*plans_last_thread[1] != NULL);
        fftw_execute_dft(
            *plans_last_thread[1],
            grid_out + block_sizes[1] * fft_size[1] * fft_size[2] * thread_id,
            grid_in + block_sizes[1] * fft_size[1] * fft_size[2] * thread_id);
      }
#pragma omp barrier
      if (thread_id < number_of_full_blocks[0]) {
        assert(*plans[0] != NULL);
        fftw_execute_dft(*plans[0], grid_in + block_sizes[0] * thread_id,
                         grid_out + block_sizes[0] * thread_id);
      } else if (thread_id == total_number_of_blocks[0] - 1) {
        assert(*plans_last_thread[0] != NULL);
        fftw_execute_dft(*plans_last_thread[0],
                         grid_in + block_sizes[0] * thread_id,
                         grid_out + block_sizes[0] * thread_id);
      }
    }
  } else {
    fft_key_t key;
    get_key_3d(direction, fft_size, omp_get_max_threads(), in_place, key);
    fftw_plan *plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      double complex *buffer =
          fftw_alloc_complex(get_buffer_size_from_key(key));
      plan = fft_fftw_create_3d_plan(key, buffer, in_place ? buffer : grid_out);
      fftw_free(buffer);
    }
    assert(plan != NULL);
    assert(grid_in != NULL);
    assert(grid_out != NULL);
    fftw_execute_dft(*plan, grid_in, grid_out);
  }
#else
  (void)fft_size;
  (void)grid_in;
  (void)grid_out;
  assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local forward R2C 3D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_3d_r2c_local(const int fft_size[3], double *grid_in,
                           double complex *grid_out) {
#if defined(__FFTW3)
  assert(omp_get_num_threads() == 1);
  assert(is_initialized);
  if (fft_size[0] == 0 || fft_size[1] == 0 || fft_size[2] == 0)
    return;
  const bool in_place = grid_in == (double *)grid_out;
  if (has_guru_interface &&
      ((((fft_size[0] >= 256 || fft_size[1] >= 256 || fft_size[2] >= 256) &&
         !in_place) ||
        ((fft_size[0] >= 128 || fft_size[1] >= 128 || fft_size[2] >= 128) &&
         in_place)) &&
       omp_get_max_threads() > 1) &&
      (fftw_planning_mode == FFTW_ESTIMATE)) {
    // The 3D FFT is not efficient with threading and estimate planning mode
    // So, we decompose it in a sequence of 1D FFTs
    int number_of_threads = 1;
#pragma omp parallel default(none) shared(number_of_threads)
    {
#pragma omp single
      {
        number_of_threads = omp_get_num_threads();
      }
    }
    int block_sizes[3];
    fftw_plan *plans[3], *plans_last_thread[3];
    int number_of_full_blocks[3] = {0, 0, 0};
    int total_number_of_blocks[3] = {0, 0, 0};
    {
      const int number_of_ffts = fft_size[1] * (fft_size[2] / 2 + 1);
      block_sizes[0] =
          (number_of_ffts + number_of_threads - 1) / number_of_threads;
      number_of_full_blocks[0] = number_of_ffts / block_sizes[0];
      total_number_of_blocks[0] =
          (number_of_ffts + block_sizes[0] - 1) / block_sizes[0];
      fft_iodim dim = {
          .n = fft_size[0], .is = number_of_ffts, .os = number_of_ffts};
      fft_iodim howmany_dim = {.n = block_sizes[0], .is = 1, .os = 1};
      fft_key_t key;
      get_key_guru(true, 1, &dim, 1, &howmany_dim, 1, in_place, key);

      plans[0] = lookup_plan_from_cache(key);
      if (plans[0] == NULL) {
        double complex *buffer =
            fftw_alloc_complex(get_buffer_size_from_key(key));
        plans[0] = fft_fftw_create_guru_plan(key, buffer, grid_out);
        fftw_free(buffer);
      }
      if (block_sizes[0] * number_of_threads != number_of_ffts) {
        const int block_size_last_thread = number_of_ffts % block_sizes[0];
        howmany_dim.n = block_size_last_thread;
        fft_key_t key;
        get_key_guru(true, 1, &dim, 1, &howmany_dim, 1, in_place, key);

        plans_last_thread[0] = lookup_plan_from_cache(key);
        if (plans_last_thread[0] == NULL) {
          double complex *buffer =
              fftw_alloc_complex(get_buffer_size_from_key(key));
          plans_last_thread[0] =
              fft_fftw_create_guru_plan(key, buffer, grid_out);
          fftw_free(buffer);
        }
      }
    }
    {
      const int number_of_ffts = fft_size[0];
      block_sizes[1] =
          (number_of_ffts + number_of_threads - 1) / number_of_threads;
      number_of_full_blocks[1] = number_of_ffts / block_sizes[1];
      total_number_of_blocks[1] =
          (number_of_ffts + block_sizes[1] - 1) / block_sizes[1];
      fft_iodim dim = {.n = fft_size[1],
                       .is = fft_size[2] / 2 + 1,
                       .os = fft_size[2] / 2 + 1};
      fft_iodim howmany_dims[2] = {
          {.n = block_sizes[1],
           .is = fft_size[1] * (fft_size[2] / 2 + 1),
           .os = fft_size[1] * (fft_size[2] / 2 + 1)},
          {.n = (fft_size[2] / 2 + 1), .is = 1, .os = 1}};
      fft_key_t key;
      get_key_guru(true, 1, &dim, 2, howmany_dims, 1, in_place, key);

      plans[1] = lookup_plan_from_cache(key);
      if (plans[1] == NULL) {
        double complex *buffer =
            fftw_alloc_complex(get_buffer_size_from_key(key));
        plans[1] = fft_fftw_create_guru_plan(key, buffer, grid_out);
        fftw_free(buffer);
      }
      if (block_sizes[1] * number_of_threads != number_of_ffts) {
        const int block_size_last_thread = number_of_ffts % block_sizes[1];
        howmany_dims[0].n = block_size_last_thread;
        fft_key_t key;
        get_key_guru(true, 1, &dim, 2, howmany_dims, 1, in_place, key);

        plans_last_thread[1] = lookup_plan_from_cache(key);
        if (plans_last_thread[1] == NULL) {
          double complex *buffer =
              fftw_alloc_complex(get_buffer_size_from_key(key));
          plans_last_thread[1] =
              fft_fftw_create_guru_plan(key, buffer, grid_out);
          fftw_free(buffer);
        }
      }
    }
    {
      const int number_of_ffts = fft_size[0] * fft_size[1];
      block_sizes[2] =
          (number_of_ffts + number_of_threads - 1) / number_of_threads;
      number_of_full_blocks[2] = number_of_ffts / block_sizes[2];
      total_number_of_blocks[2] =
          (number_of_ffts + block_sizes[2] - 1) / block_sizes[2];
      fft_iodim dim = {.n = fft_size[2], .is = 1, .os = 1};
      fft_iodim howmany_dim = {.n = block_sizes[2],
                               .is = in_place ? 2 * (fft_size[2] / 2 + 1)
                                              : fft_size[2],
                               .os = fft_size[2] / 2 + 1};
      fft_key_t key;
      get_key_guru_r2c(true, 1, &dim, 1, &howmany_dim, 1, in_place, key);

      plans[2] = lookup_plan_from_cache(key);
      if (plans[2] == NULL) {
        double *buffer = fftw_alloc_real(2 * get_buffer_size_from_key(key));
        plans[2] =
            fft_fftw_create_guru_plan_r2c(key, buffer, (double *)grid_out);
        fftw_free(buffer);
      }
      if (block_sizes[2] * number_of_threads != number_of_ffts) {
        const int block_size_last_thread = number_of_ffts % block_sizes[2];
        howmany_dim.n = block_size_last_thread;
        fft_key_t key;
        get_key_guru_r2c(true, 1, &dim, 1, &howmany_dim, 1, in_place, key);

        plans_last_thread[2] = lookup_plan_from_cache(key);
        if (plans_last_thread[2] == NULL) {
          double *buffer = fftw_alloc_real(2 * get_buffer_size_from_key(key));
          plans_last_thread[2] =
              fft_fftw_create_guru_plan_r2c(key, buffer, (double *)grid_out);
          fftw_free(buffer);
        }
      }
    }
#pragma omp parallel default(none)                                             \
    shared(number_of_threads, total_number_of_blocks, number_of_full_blocks,   \
               plans, plans_last_thread, block_sizes, fft_size, grid_in,       \
               grid_out, in_place)
    {
      const int thread_id = omp_get_thread_num();
      if (thread_id < number_of_full_blocks[2]) {
        assert(*plans[2] != NULL);
        fftw_execute_dft_r2c(
            *plans[2],
            grid_in + block_sizes[2] *
                          (in_place ? 2 * (fft_size[2] / 2 + 1) : fft_size[2]) *
                          thread_id,
            grid_out + block_sizes[2] * (fft_size[2] / 2 + 1) * thread_id);
      } else if (thread_id == total_number_of_blocks[2] - 1) {
        assert(*plans_last_thread[2] != NULL);
        fftw_execute_dft_r2c(
            *plans_last_thread[2],
            grid_in + block_sizes[2] *
                          (in_place ? 2 * (fft_size[2] / 2 + 1) : fft_size[2]) *
                          thread_id,
            grid_out + block_sizes[2] * (fft_size[2] / 2 + 1) * thread_id);
      }
#pragma omp barrier
      if (thread_id < number_of_full_blocks[1]) {
        assert(*plans[1] != NULL);
        fftw_execute_dft(*plans[1],
                         grid_out + block_sizes[1] * fft_size[1] *
                                        (fft_size[2] / 2 + 1) * thread_id,
                         ((double complex *)grid_in) +
                             block_sizes[1] * fft_size[1] *
                                 (fft_size[2] / 2 + 1) * thread_id);
      } else if (thread_id == total_number_of_blocks[1] - 1) {
        assert(*plans_last_thread[1] != NULL);
        fftw_execute_dft(*plans_last_thread[1],
                         grid_out + block_sizes[1] * fft_size[1] *
                                        (fft_size[2] / 2 + 1) * thread_id,
                         ((double complex *)grid_in) +
                             block_sizes[1] * fft_size[1] *
                                 (fft_size[2] / 2 + 1) * thread_id);
      }
#pragma omp barrier
      if (thread_id < number_of_full_blocks[0]) {
        assert(*plans[0] != NULL);
        fftw_execute_dft(*plans[0],
                         (double complex *)grid_in + block_sizes[0] * thread_id,
                         grid_out + block_sizes[0] * thread_id);
      } else if (thread_id == total_number_of_blocks[0] - 1) {
        assert(*plans_last_thread[0] != NULL);
        fftw_execute_dft(*plans_last_thread[0],
                         (double complex *)grid_in + block_sizes[0] * thread_id,
                         grid_out + block_sizes[0] * thread_id);
      }
    }
  } else {
    fft_key_t key;
    get_key_3d_r2c(true, fft_size, omp_get_max_threads(), in_place, key);
    fftw_plan *plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      double *buffer = fftw_alloc_real(2 * get_buffer_size_from_key(key));
      plan = fft_fftw_create_3d_plan_r2c(
          key, buffer, in_place ? buffer : (double *)grid_out);
      fftw_free(buffer);
    }
    assert(plan != NULL);
    assert(grid_in != NULL);
    assert(grid_out != NULL);
    fftw_execute_dft_r2c(*plan, grid_in, grid_out);
  }
#else
  (void)fft_size;
  (void)grid_in;
  (void)grid_out;
  assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

/*******************************************************************************
 * \brief Performs a local backwards R2C 3D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_3d_c2r_local(const int fft_size[3], double complex *grid_in,
                           double *grid_out) {
#if defined(__FFTW3)
  assert(omp_get_num_threads() == 1);
  assert(is_initialized);
  if (fft_size[0] == 0 || fft_size[1] == 0 || fft_size[2] == 0)
    return;
  const bool in_place = (double *)grid_in == grid_out;
  if (has_guru_interface &&
      ((((fft_size[0] >= 256 || fft_size[1] >= 256 || fft_size[2] >= 256) &&
         !in_place) ||
        ((fft_size[0] >= 128 || fft_size[1] >= 128 || fft_size[2] >= 128) &&
         in_place)) &&
       omp_get_max_threads() > 1) &&
      (fftw_planning_mode == FFTW_ESTIMATE)) {
    // We cannot use the output buffer for planning because it may be too small
    // So, we allocate a new one
    double complex *buffer =
        fftw_alloc_complex(fft_size[0] * fft_size[1] * (fft_size[2] / 2 + 1));
    // The 3D FFT is not efficient with threading and estimate planning mode
    // So, we decompose it in a sequence of 1D FFTs
    int number_of_threads = 1;
#pragma omp parallel default(none) shared(number_of_threads)
    {
#pragma omp single
      {
        number_of_threads = omp_get_num_threads();
      }
    }
    int block_sizes[3];
    fftw_plan *plans[3], *plans_last_thread[3];
    int number_of_full_blocks[3] = {0, 0, 0};
    int total_number_of_blocks[3] = {0, 0, 0};
    {
      const int number_of_ffts = fft_size[1] * (fft_size[2] / 2 + 1);
      block_sizes[0] =
          (number_of_ffts + number_of_threads - 1) / number_of_threads;
      number_of_full_blocks[0] = number_of_ffts / block_sizes[0];
      total_number_of_blocks[0] =
          (number_of_ffts + block_sizes[0] - 1) / block_sizes[0];
      fft_iodim dim = {
          .n = fft_size[0], .is = number_of_ffts, .os = number_of_ffts};
      fft_iodim howmany_dim = {.n = block_sizes[0], .is = 1, .os = 1};
      fft_key_t key;
      get_key_guru(false, 1, &dim, 1, &howmany_dim, 1, in_place, key);

      plans[0] = lookup_plan_from_cache(key);
      if (plans[0] == NULL) {
        plans[0] =
            fft_fftw_create_guru_plan(key, (double complex *)grid_out, buffer);
      }
      if (block_sizes[0] * number_of_threads != number_of_ffts) {
        const int block_size_last_thread = number_of_ffts % block_sizes[0];
        howmany_dim.n = block_size_last_thread;
        fft_key_t key;
        get_key_guru(false, 1, &dim, 1, &howmany_dim, 1, in_place, key);

        plans_last_thread[0] = lookup_plan_from_cache(key);
        if (plans_last_thread[0] == NULL) {
          double complex *buffer2 =
              fftw_alloc_complex(get_buffer_size_from_key(key));
          plans_last_thread[0] =
              fft_fftw_create_guru_plan(key, buffer2, buffer);
          fftw_free(buffer2);
        }
      }
    }
    {
      const int number_of_ffts = fft_size[0];
      block_sizes[1] =
          (number_of_ffts + number_of_threads - 1) / number_of_threads;
      number_of_full_blocks[1] = number_of_ffts / block_sizes[1];
      total_number_of_blocks[1] =
          (number_of_ffts + block_sizes[1] - 1) / block_sizes[1];
      fft_iodim dim = {.n = fft_size[1],
                       .is = fft_size[2] / 2 + 1,
                       .os = fft_size[2] / 2 + 1};
      fft_iodim howmany_dims[2] = {
          {.n = block_sizes[1],
           .is = fft_size[1] * (fft_size[2] / 2 + 1),
           .os = fft_size[1] * (fft_size[2] / 2 + 1)},
          {.n = (fft_size[2] / 2 + 1), .is = 1, .os = 1}};
      fft_key_t key;
      get_key_guru(false, 1, &dim, 2, howmany_dims, 1, in_place, key);

      plans[1] = lookup_plan_from_cache(key);
      if (plans[1] == NULL) {
        double complex *buffer2 =
            fftw_alloc_complex(get_buffer_size_from_key(key));
        plans[1] = fft_fftw_create_guru_plan(key, buffer2, buffer);
        fftw_free(buffer2);
      }
      if (block_sizes[1] * number_of_threads != number_of_ffts) {
        const int block_size_last_thread = number_of_ffts % block_sizes[1];
        howmany_dims[0].n = block_size_last_thread;
        fft_key_t key;
        get_key_guru(false, 1, &dim, 2, howmany_dims, 1, in_place, key);

        plans_last_thread[1] = lookup_plan_from_cache(key);
        if (plans_last_thread[1] == NULL) {
          double complex *buffer2 =
              fftw_alloc_complex(get_buffer_size_from_key(key));
          plans_last_thread[1] =
              fft_fftw_create_guru_plan(key, buffer2, buffer);
          fftw_free(buffer2);
        }
      }
    }
    {
      const int number_of_ffts = fft_size[0] * fft_size[1];
      block_sizes[2] =
          (number_of_ffts + number_of_threads - 1) / number_of_threads;
      number_of_full_blocks[2] = number_of_ffts / block_sizes[2];
      total_number_of_blocks[2] =
          (number_of_ffts + block_sizes[2] - 1) / block_sizes[2];
      fft_iodim dim = {.n = fft_size[2], .is = 1, .os = 1};
      fft_iodim howmany_dim = {.n = block_sizes[2],
                               .is = fft_size[2] / 2 + 1,
                               .os = in_place ? 2 * (fft_size[2] / 2 + 1)
                                              : fft_size[2]};
      fft_key_t key;
      get_key_guru_r2c(false, 1, &dim, 1, &howmany_dim, 1, in_place, key);

      plans[2] = lookup_plan_from_cache(key);
      if (plans[2] == NULL) {
        double complex *buffer2 =
            fftw_alloc_complex(get_buffer_size_from_key(key));
        plans[2] = fft_fftw_create_guru_plan_r2c(key, (double *)buffer,
                                                 (double *)buffer2);
        fftw_free(buffer2);
      }
      if (block_sizes[2] * number_of_threads != number_of_ffts) {
        const int block_size_last_thread = number_of_ffts % block_sizes[2];
        howmany_dim.n = block_size_last_thread;
        fft_key_t key;
        get_key_guru_r2c(false, 1, &dim, 1, &howmany_dim, 1, in_place, key);

        plans_last_thread[2] = lookup_plan_from_cache(key);
        if (plans_last_thread[2] == NULL) {
          double complex *buffer2 =
              fftw_alloc_complex(get_buffer_size_from_key(key));
          plans_last_thread[2] = fft_fftw_create_guru_plan_r2c(
              key, (double *)buffer, (double *)buffer2);
          fftw_free(buffer2);
        }
      }
    }
#pragma omp parallel default(none)                                             \
    shared(number_of_threads, total_number_of_blocks, number_of_full_blocks,   \
               plans, plans_last_thread, block_sizes, fft_size, grid_in,       \
               grid_out, buffer, in_place)
    {
      const int thread_id = omp_get_thread_num();
      if (thread_id < number_of_full_blocks[0]) {
        assert(*plans[0] != NULL);
        fftw_execute_dft(*plans[0], grid_in + block_sizes[0] * thread_id,
                         buffer + block_sizes[0] * thread_id);
      } else if (thread_id == total_number_of_blocks[0] - 1) {
        assert(*plans_last_thread[0] != NULL);
        fftw_execute_dft(*plans_last_thread[0],
                         grid_in + block_sizes[0] * thread_id,
                         buffer + block_sizes[0] * thread_id);
      }
#pragma omp barrier
      if (thread_id < number_of_full_blocks[1]) {
        assert(*plans[1] != NULL);
        fftw_execute_dft(*plans[1],
                         buffer + block_sizes[1] * fft_size[1] *
                                      (fft_size[2] / 2 + 1) * thread_id,
                         grid_in + block_sizes[1] * fft_size[1] *
                                       (fft_size[2] / 2 + 1) * thread_id);
      } else if (thread_id == total_number_of_blocks[1] - 1) {
        assert(*plans_last_thread[1] != NULL);
        fftw_execute_dft(*plans_last_thread[1],
                         buffer + block_sizes[1] * fft_size[1] *
                                      (fft_size[2] / 2 + 1) * thread_id,
                         grid_in + block_sizes[1] * fft_size[1] *
                                       (fft_size[2] / 2 + 1) * thread_id);
      }
#pragma omp barrier
      if (thread_id < number_of_full_blocks[2]) {
        assert(*plans[2] != NULL);
        fftw_execute_dft_c2r(
            *plans[2],
            grid_in + block_sizes[2] * (fft_size[2] / 2 + 1) * thread_id,
            grid_out +
                block_sizes[2] *
                    (in_place ? 2 * (fft_size[2] / 2 + 1) : fft_size[2]) *
                    thread_id);
      } else if (thread_id == total_number_of_blocks[2] - 1) {
        assert(*plans_last_thread[2] != NULL);
        fftw_execute_dft_c2r(
            *plans_last_thread[2],
            grid_in + block_sizes[2] * (fft_size[2] / 2 + 1) * thread_id,
            grid_out +
                block_sizes[2] *
                    (in_place ? 2 * (fft_size[2] / 2 + 1) : fft_size[2]) *
                    thread_id);
      }
    }
    fftw_free(buffer);
  } else {
    fft_key_t key;
    get_key_3d_r2c(false, fft_size, omp_get_max_threads(), in_place, key);
    fftw_plan *plan = lookup_plan_from_cache(key);
    if (plan == NULL) {
      double complex *buffer =
          fftw_alloc_complex(get_buffer_size_from_key(key));
      plan = fft_fftw_create_3d_plan_r2c(
          key, (double *)buffer, in_place ? (double *)buffer : grid_out);
      fftw_free(buffer);
    }
    assert(plan != NULL);
    assert(grid_in != NULL);
    assert(grid_out != NULL);
    fftw_execute_dft_c2r(*plan, grid_in, grid_out);
  }
#else
  (void)fft_size;
  (void)grid_in;
  (void)grid_out;
  assert(0 && "The grid library was not compiled with FFTW support.");
#endif
}

// EOF
