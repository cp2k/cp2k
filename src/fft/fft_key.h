/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#ifndef FFT_KEY_H
#define FFT_KEY_H

// Does not work with fft_gpu, so each file including fft_key.h must include
// cp_mpi.h on this own
// #include "../mpiwrap/cp_mpi.h"
#include "fft_utils.h"

#if defined(__FFTW3)
#include <fftw3.h>

#include <assert.h>

#define KEY_SIZE 12
typedef int fft_key_t[KEY_SIZE];
// Meaning of the elements:
// 0: rank, r2c, inplace, forward (see below)
// 1: communicator (currently only the null-communicator -> local FFTs, to be
// used later for distributed FFTs of FFTW3) 2: number of threads 3-5: FFT sizes
// (or number of FFTs) 6-8: input strides 9-11: output strides

// These constants are used to encode general information on a FFT into bits
// Modulo 4 encodes the rank (1, 2, 3)
// The third bit encodes whether a R2C/C2R-FFT is used
// 4 == 2^2
#define FFT_KEY_R2C 4
// The fourth bit encodes inplace FFTs
// 8 == 2^3
#define FFT_KEY_INPLACE 8
// The fifth bit encodes the direction (1 is forward, 0 is backward FFT)
// 16 == 2^4
#define FFT_KEY_FORWARD 16

/*******************************************************************************
 * \brief Get key from FFT input parameters (1D, C2C-case)
 * \author Frederick Stein
 ******************************************************************************/
static inline void
get_key_1d(const bool direction, const int fft_size, const int number_of_ffts,
           const bool transpose_in, const bool transpose_out,
           const int leading_dimension_in, const int leading_dimension_out,
           const int number_of_threads, const bool inplace, int *key) {
  assert((!inplace || (transpose_in == transpose_out)) &&
         "Inplace plans need the same strides for input and output arrays!");
  key[0] = 1 + FFT_KEY_INPLACE * inplace + direction * FFT_KEY_FORWARD;
  key[1] = cp_mpi_comm_c2f(cp_mpi_get_comm_null());
  key[2] = number_of_threads;
  key[3] = fft_size;
  key[4] = number_of_ffts;
  key[5] = 0;
  key[6] = transpose_in ? leading_dimension_in : 1;
  key[7] = transpose_in ? 1 : leading_dimension_in;
  key[8] = 0;
  key[9] = transpose_out ? leading_dimension_out : 1;
  key[10] = transpose_out ? 1 : leading_dimension_out;
  key[11] = 0;
}

/*******************************************************************************
 * \brief Get key from FFT input parameters (1D, R2C-case)
 * \author Frederick Stein
 ******************************************************************************/
static inline void
get_key_1d_r2c(const bool direction, const int fft_size,
               const int number_of_ffts, const bool transpose_in,
               const bool transpose_out, const int leading_dimension_in,
               const int leading_dimension_out, const int number_of_threads,
               const bool inplace, int *key) {
  assert((!inplace || (transpose_in == transpose_out)) &&
         "Inplace plans need the same strides for input and output arrays!");
  key[0] =
      1 + FFT_KEY_R2C + FFT_KEY_INPLACE * inplace + direction * FFT_KEY_FORWARD;
  key[1] = cp_mpi_comm_c2f(cp_mpi_get_comm_null());
  key[2] = number_of_threads;
  key[3] = fft_size;
  key[4] = number_of_ffts;
  key[5] = 0;
  key[6] = transpose_in ? leading_dimension_in : 1;
  key[7] = transpose_in ? 1 : leading_dimension_in;
  key[8] = 0;
  key[9] = transpose_out ? leading_dimension_out : 1;
  key[10] = transpose_out ? 1 : leading_dimension_out;
  key[11] = 0;
}

/*******************************************************************************
 * \brief Get key from FFT input parameters (2D, C2C-case)
 * \author Frederick Stein
 ******************************************************************************/
static inline void get_key_2d(const bool direction, const int fft_size[2],
                              const int number_of_ffts, const bool transpose_in,
                              const bool transpose_out,
                              const int leading_dimension_in,
                              const int leading_dimension_out,
                              const int number_of_threads, const bool inplace,
                              int *key) {
  assert((!inplace || (transpose_in == transpose_out)) &&
         "Inplace plans need the same strides for input and output arrays!");
  key[0] = 2 + FFT_KEY_INPLACE * inplace + direction * FFT_KEY_FORWARD;
  key[1] = cp_mpi_comm_c2f(cp_mpi_get_comm_null());
  key[2] = number_of_threads;
  key[3] = fft_size[0];
  key[4] = fft_size[1];
  key[5] = number_of_ffts;
  key[6] = (transpose_in ? leading_dimension_in : 1) * fft_size[1];
  key[7] = transpose_in ? leading_dimension_in : 1;
  key[8] = transpose_in ? 1 : leading_dimension_in;
  key[9] = (transpose_out ? leading_dimension_out : 1) * fft_size[1];
  key[10] = transpose_out ? leading_dimension_out : 1;
  key[11] = transpose_out ? 1 : leading_dimension_out;
}

/*******************************************************************************
 * \brief Get key from FFT input parameters (2D, R2C-case)
 * \author Frederick Stein
 ******************************************************************************/
static inline void
get_key_2d_r2c(const bool direction, const int fft_size[2],
               const int number_of_ffts, const bool transpose_in,
               const bool transpose_out, const int leading_dimension_in,
               const int leading_dimension_out, const int number_of_threads,
               const bool inplace, int *key) {
  assert((!inplace || (transpose_in == transpose_out)) &&
         "Inplace plans need the same strides for input and output arrays!");
  key[0] =
      2 + FFT_KEY_R2C + FFT_KEY_INPLACE * inplace + direction * FFT_KEY_FORWARD;
  key[1] = cp_mpi_comm_c2f(cp_mpi_get_comm_null());
  key[2] = number_of_threads;
  key[3] = fft_size[0];
  key[4] = fft_size[1];
  key[5] = number_of_ffts;
  key[6] = (transpose_in ? leading_dimension_in : 1) *
           (direction ? (!inplace ? fft_size[1] : 2 * (fft_size[1] / 2 + 1))
                      : fft_size[1] / 2 + 1);
  key[7] = transpose_in ? leading_dimension_in : 1;
  key[8] = transpose_in ? 1 : leading_dimension_in;
  key[9] = (transpose_out ? leading_dimension_out : 1) *
           (!direction ? (!inplace ? fft_size[1] : 2 * (fft_size[1] / 2 + 1))
                       : fft_size[1] / 2 + 1);
  key[10] = transpose_out ? leading_dimension_out : 1;
  key[11] = transpose_out ? 1 : leading_dimension_out;
}

/*******************************************************************************
 * \brief Get key from FFT input parameters (3D, C2C-case)
 * \author Frederick Stein
 ******************************************************************************/
static inline void get_key_3d(const bool direction, const int fft_size[3],
                              const int number_of_threads, const bool inplace,
                              int *key) {
  key[0] = 3 + FFT_KEY_INPLACE * inplace + direction * FFT_KEY_FORWARD;
  key[1] = cp_mpi_comm_c2f(cp_mpi_get_comm_null());
  key[2] = number_of_threads;
  key[3] = fft_size[0];
  key[4] = fft_size[1];
  key[5] = fft_size[2];
  key[6] = fft_size[1] * fft_size[2];
  key[7] = fft_size[2];
  key[8] = 1;
  key[9] = fft_size[1] * fft_size[2];
  key[10] = fft_size[2];
  key[11] = 1;
}

/*******************************************************************************
 * \brief Get key from FFT input parameters (3D, R2C-case)
 * \author Frederick Stein
 ******************************************************************************/
static inline void get_key_3d_r2c(const bool direction, const int fft_size[3],
                                  const int number_of_threads,
                                  const bool inplace, int *key) {
  key[0] =
      3 + FFT_KEY_R2C + FFT_KEY_INPLACE * inplace + direction * FFT_KEY_FORWARD;
  key[1] = cp_mpi_comm_c2f(cp_mpi_get_comm_null());
  key[2] = number_of_threads;
  key[3] = fft_size[0];
  key[4] = fft_size[1];
  key[5] = fft_size[2];
  key[6] = fft_size[1] * (direction ? fft_size[2] : fft_size[2] / 2 + 1);
  key[7] = direction ? fft_size[2] : fft_size[2] / 2 + 1;
  key[8] = 1;
  key[9] = fft_size[1] * (!direction ? fft_size[2] : fft_size[2] / 2 + 1);
  key[10] = !direction ? fft_size[2] : fft_size[2] / 2 + 1;
  key[11] = 1;
}

/*******************************************************************************
 * \brief Get key from FFT input parameters (Guru, C2C-case)
 * \author Frederick Stein
 ******************************************************************************/
static inline void get_key_guru(const bool direction, int rank,
                                const fftw_iodim *dims, int howmany_rank,
                                const fftw_iodim *howmany_dims,
                                const int number_of_threads, const bool inplace,
                                int *key) {
  assert(rank + howmany_rank <= 3 &&
         "Larger combined ranks than 3 are not implemented\n");

  key[0] = rank + FFT_KEY_INPLACE * inplace + direction * FFT_KEY_FORWARD;
  key[1] = cp_mpi_comm_c2f(cp_mpi_get_comm_null());
  key[2] = number_of_threads;
  key[3] =
      rank > 0 ? dims[0].n : (rank + howmany_rank > 0 ? howmany_dims[0].n : 0);
  key[4] = rank > 1 ? dims[1].n
                    : (rank + howmany_rank > 1 ? howmany_dims[1 - rank].n : 0);
  key[5] = rank > 2 ? dims[2].n
                    : (rank + howmany_rank > 2 ? howmany_dims[2 - rank].n : 0);
  key[6] = rank > 0 ? dims[0].is
                    : (rank + howmany_rank > 0 ? howmany_dims[0].is : 0);
  key[7] = rank > 1 ? dims[1].is
                    : (rank + howmany_rank > 1 ? howmany_dims[1 - rank].is : 0);
  key[8] = rank > 2 ? dims[2].is
                    : (rank + howmany_rank > 2 ? howmany_dims[2 - rank].is : 0);
  key[9] = rank > 0 ? dims[0].os
                    : (rank + howmany_rank > 0 ? howmany_dims[0].os : 0);
  key[10] = rank > 1
                ? dims[1].os
                : (rank + howmany_rank > 1 ? howmany_dims[1 - rank].os : 0);
  key[11] = rank > 2
                ? dims[2].os
                : (rank + howmany_rank > 2 ? howmany_dims[2 - rank].os : 0);
}

/*******************************************************************************
 * \brief Get key from FFT input parameters (Guru, R2C-case)
 * \author Frederick Stein
 ******************************************************************************/
static inline void get_key_guru_r2c(const bool direction, int rank,
                                    const fftw_iodim *dims, int howmany_rank,
                                    const fftw_iodim *howmany_dims,
                                    const int number_of_threads,
                                    const bool inplace, int *key) {
  assert(rank + howmany_rank <= 3 &&
         "Larger combined ranks than 3 are not implemented\n");

  key[0] = rank + FFT_KEY_R2C + FFT_KEY_INPLACE * inplace +
           direction * FFT_KEY_FORWARD;
  key[1] = cp_mpi_comm_c2f(cp_mpi_get_comm_null());
  key[2] = number_of_threads;
  key[3] =
      rank > 0 ? dims[0].n : (rank + howmany_rank > 0 ? howmany_dims[0].n : 0);
  key[4] = rank > 1 ? dims[1].n
                    : (rank + howmany_rank > 1 ? howmany_dims[1 - rank].n : 0);
  key[5] = rank > 2 ? dims[2].n
                    : (rank + howmany_rank > 2 ? howmany_dims[2 - rank].n : 0);
  key[6] = rank > 0 ? dims[0].is
                    : (rank + howmany_rank > 0 ? howmany_dims[0].is : 0);
  key[7] = rank > 1 ? dims[1].is
                    : (rank + howmany_rank > 1 ? howmany_dims[1 - rank].is : 0);
  key[8] = rank > 2 ? dims[2].is
                    : (rank + howmany_rank > 2 ? howmany_dims[2 - rank].is : 0);
  key[9] = rank > 0 ? dims[0].os
                    : (rank + howmany_rank > 0 ? howmany_dims[0].os : 0);
  key[10] = rank > 1
                ? dims[1].os
                : (rank + howmany_rank > 1 ? howmany_dims[1 - rank].os : 0);
  key[11] = rank > 2
                ? dims[2].os
                : (rank + howmany_rank > 2 ? howmany_dims[2 - rank].os : 0);
}

/*******************************************************************************
 * \brief Determine buffer size for a local FFT from a key
 * \author Frederick Stein
 ******************************************************************************/
static inline int get_buffer_size_from_key(const fft_key_t key) {
  int buffer_size = 0;
  for (int r = 0; r < 3; r++) {
    if (r != key[0] % 4 - 1 || ((key[0] & FFT_KEY_R2C) != FFT_KEY_R2C)) {
      buffer_size += key[3 + r] * key[6 + r];
    } else {
      buffer_size += (key[3 + r] / 2 + 1) * key[4 + r];
    }
  }
  int buffer_size_out = 0;
  for (int r = 0; r < 3; r++) {
    if (r != key[0] % 4 - 1 || ((key[0] & FFT_KEY_R2C) != FFT_KEY_R2C)) {
      buffer_size_out += key[3 + r] * key[9 + r];
    } else {
      buffer_size_out += (key[3 + r] / 2 + 1) * key[9 + r];
    }
  }
  return imax(buffer_size, buffer_size_out);
}

static inline void fetch_data_from_key_nd(
    const fft_key_t key, int *rank, int *fft_size, int *number_of_ffts,
    int *number_of_threads, bool *direction, bool *inplace, int *inembed,
    int *onembed, int *idist, int *odist, int *istride, int *ostride) {
  assert(key[1] == cp_mpi_comm_c2f(cp_mpi_get_comm_null()) &&
         "Distributed FFTs are not supported in this function!");
  *rank = key[0] % 4;
  *inplace = (key[0] & FFT_KEY_INPLACE) == FFT_KEY_INPLACE;
  *number_of_threads = key[2];
  *direction = (key[0] & FFT_KEY_FORWARD) == FFT_KEY_FORWARD;
  const int is_r2c = (key[0] & FFT_KEY_R2C) == FFT_KEY_R2C;
  for (int r = 0; r < *rank; r++) {
    fft_size[r] = key[3 + r];
    if (r < *rank - 1) {
      inembed[r] = key[3 + r];
      onembed[r] = key[3 + r];
    } else {
      if (*direction) {
        inembed[r] = is_r2c && *inplace ? 2 * (key[3 + r] / 2 + 1) : key[3 + r];
        onembed[r] = is_r2c ? key[3 + r] / 2 + 1 : key[3 + r];
      } else {
        inembed[r] = is_r2c ? key[3 + r] / 2 + 1 : key[3 + r];
        onembed[r] = is_r2c && *inplace ? 2 * (key[3 + r] / 2 + 1) : key[3 + r];
      }
    }
  }
  *number_of_threads = key[2];
  *number_of_ffts = *rank < 3 ? key[3 + *rank] : 1;
  if (*rank < 3) {
    *idist = key[6 + *rank];
    *odist = key[9 + *rank];
    *istride = key[5 + *rank];
    *ostride = key[8 + *rank];
  } else {
    *idist = 1;
    *odist = 1;
    *istride = 1;
    *ostride = 1;
  }
}

static inline void fetch_data_from_key_guru(const fft_key_t key,
                                            bool *direction, int *rank,
                                            fftw_iodim *dims, int *howmany_rank,
                                            fftw_iodim *howmany_dims,
                                            int *number_of_threads,
                                            bool *inplace) {
  assert(key[1] == cp_mpi_comm_c2f(cp_mpi_get_comm_null()) &&
         "Distributed FFTs are not supported in this function!");
  *direction = (key[0] & FFT_KEY_FORWARD) == FFT_KEY_FORWARD;
  *inplace = (key[0] & FFT_KEY_INPLACE) == FFT_KEY_INPLACE;
  *rank = key[0] % 4;
  *howmany_rank = 3 - *rank - (key[5] == 0) - (key[4] == 0) - (key[3] == 0);
  for (int r = 0; r < *rank; r++) {
    dims[r].n = key[3 + r];
    dims[r].is = key[6 + r];
    dims[r].os = key[9 + r];
  }
  for (int r = 0; r < *howmany_rank; r++) {
    howmany_dims[r].n = key[3 + *rank + r];
    howmany_dims[r].is = key[6 + *rank + r];
    howmany_dims[r].os = key[9 + *rank + r];
  }
  *number_of_threads = key[2];
}
#endif

#endif /* FFT_KEY_H */