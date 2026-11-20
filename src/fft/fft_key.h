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

#include "../base/cp_assert.h"

#define KEY_SIZE 15
typedef int fft_key_t[KEY_SIZE];
// Meaning of the elements:
// 0: rank, r2c, inplace, forward (see below)
// 1: communicator (currently only the null-communicator -> local FFTs, to be
// used later for distributed FFTs of FFTW3)
// 2: number of threads
// 3-6: FFT sizes for the transform dimensions followed by the batch size
// (unused trailing slots are zero, so rank + 1 dimensions are stored)
// 7-10: input strides in the same order (transform dimensions, then batch)
// 11-14: output strides in the same order (transform dimensions, then batch)

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
 * \brief Fill the size and stride fields of a key for a non-Guru FFT.
 * \param key Encoded FFT key.
 * \param rank Number of transform dimensions.
 * \param fft_size Transform sizes (up to three dimensions).
 * \param number_of_ffts Number of batched transforms.
 * \param direction FFT direction (forward or backward).
 * \param inplace Whether the transform is in-place.
 * \param is_r2c Whether this is a real-to-complex or complex-to-real FFT.
 * \param istride Stride of the first input dimension.
 * \param idist Distance between consecutive input transforms.
 * \param ostride Stride of the first output dimension.
 * \param odist Distance between consecutive output transforms.
 ******************************************************************************/
static inline void set_nd_key(int *key, const int rank, const int fft_size[3],
                              const int number_of_ffts, const bool direction,
                              const bool inplace, const bool is_r2c,
                              const int istride, const int idist,
                              const int ostride, const int odist) {
  int inembed[3] = {0, 0, 0};
  int onembed[3] = {0, 0, 0};
  for (int r = 0; r < rank; r++) {
    inembed[r] = fft_size[r];
    onembed[r] = fft_size[r];
  }
  if (is_r2c && rank > 0) {
    const int r = rank - 1;
    if (direction) {
      onembed[r] = fft_size[r] / 2 + 1;
      if (inplace)
        inembed[r] = 2 * (fft_size[r] / 2 + 1);
    } else {
      inembed[r] = fft_size[r] / 2 + 1;
      if (inplace)
        onembed[r] = 2 * (fft_size[r] / 2 + 1);
    }
  }
  for (int d = 0; d < 4; d++) {
    key[3 + d] = 0;
    key[7 + d] = 0;
    key[11 + d] = 0;
  }
  for (int d = 0; d < rank; d++)
    key[3 + d] = fft_size[d];
  key[3 + rank] = number_of_ffts;
  key[7] = istride;
  for (int d = 1; d < rank; d++)
    key[7 + d] = key[7 + d - 1] * inembed[d - 1];
  key[7 + rank] = idist;
  key[11] = ostride;
  for (int d = 1; d < rank; d++)
    key[11 + d] = key[11 + d - 1] * onembed[d - 1];
  key[11 + rank] = odist;
}

/*******************************************************************************
 * \brief Get key from FFT input parameters (1D, C2C-case)
 * \author Frederick Stein
 ******************************************************************************/
static inline void
get_key_1d(const bool direction, const int fft_size, const int number_of_ffts,
           const bool transpose_in, const bool transpose_out,
           const int leading_dimension_in, const int leading_dimension_out,
           const int number_of_threads, const bool inplace, int *key) {
  cp_assert((!inplace || (transpose_in == transpose_out)) &&
            "Inplace plans need the same strides for input and output arrays!");
  const int fft_size_3d[3] = {fft_size, 0, 0};
  const int istride = transpose_in ? leading_dimension_in : 1;
  const int idist = transpose_in ? 1 : leading_dimension_in;
  const int ostride = transpose_out ? leading_dimension_out : 1;
  const int odist = transpose_out ? 1 : leading_dimension_out;
  key[0] = 1 + FFT_KEY_INPLACE * inplace + direction * FFT_KEY_FORWARD;
  key[1] = cp_mpi_comm_c2f(cp_mpi_get_comm_null());
  key[2] = number_of_threads;
  set_nd_key(key, 1, fft_size_3d, number_of_ffts, direction, inplace, false,
             istride, idist, ostride, odist);
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
  cp_assert((!inplace || (transpose_in == transpose_out)) &&
            "Inplace plans need the same strides for input and output arrays!");
  const int fft_size_3d[3] = {fft_size, 0, 0};
  const int istride = transpose_in ? leading_dimension_in : 1;
  const int idist = transpose_in ? 1 : leading_dimension_in;
  const int ostride = transpose_out ? leading_dimension_out : 1;
  const int odist = transpose_out ? 1 : leading_dimension_out;
  key[0] =
      1 + FFT_KEY_R2C + FFT_KEY_INPLACE * inplace + direction * FFT_KEY_FORWARD;
  key[1] = cp_mpi_comm_c2f(cp_mpi_get_comm_null());
  key[2] = number_of_threads;
  set_nd_key(key, 1, fft_size_3d, number_of_ffts, direction, inplace, true,
             istride, idist, ostride, odist);
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
  cp_assert((!inplace || (transpose_in == transpose_out)) &&
            "Inplace plans need the same strides for input and output arrays!");
  const int fft_size_3d[3] = {fft_size[0], fft_size[1], 0};
  const int istride = transpose_in ? leading_dimension_in : 1;
  const int idist = transpose_in ? 1 : leading_dimension_in;
  const int ostride = transpose_out ? leading_dimension_out : 1;
  const int odist = transpose_out ? 1 : leading_dimension_out;
  key[0] = 2 + FFT_KEY_INPLACE * inplace + direction * FFT_KEY_FORWARD;
  key[1] = cp_mpi_comm_c2f(cp_mpi_get_comm_null());
  key[2] = number_of_threads;
  set_nd_key(key, 2, fft_size_3d, number_of_ffts, direction, inplace, false,
             istride, idist, ostride, odist);
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
  cp_assert((!inplace || (transpose_in == transpose_out)) &&
            "Inplace plans need the same strides for input and output arrays!");
  const int fft_size_3d[3] = {fft_size[0], fft_size[1], 0};
  const int istride = transpose_in ? leading_dimension_in : 1;
  const int idist = transpose_in ? 1 : leading_dimension_in;
  const int ostride = transpose_out ? leading_dimension_out : 1;
  const int odist = transpose_out ? 1 : leading_dimension_out;
  key[0] =
      2 + FFT_KEY_R2C + FFT_KEY_INPLACE * inplace + direction * FFT_KEY_FORWARD;
  key[1] = cp_mpi_comm_c2f(cp_mpi_get_comm_null());
  key[2] = number_of_threads;
  set_nd_key(key, 2, fft_size_3d, number_of_ffts, direction, inplace, true,
             istride, idist, ostride, odist);
}

/*******************************************************************************
 * \brief Get key from FFT input parameters (3D, C2C-case)
 * \author Frederick Stein
 ******************************************************************************/
static inline void get_key_3d(const bool direction, const int fft_size[3],
                              const int number_of_ffts, const bool transpose_in,
                              const bool transpose_out,
                              const int number_of_threads, const bool inplace,
                              int *key) {
  cp_assert((!inplace || (transpose_in == transpose_out)) &&
            "Inplace plans need the same strides for input and output arrays!");
  const int volume = fft_size[0] * fft_size[1] * fft_size[2];
  const int istride = transpose_in ? number_of_ffts : 1;
  const int idist = transpose_in ? 1 : volume;
  const int ostride = transpose_out ? number_of_ffts : 1;
  const int odist = transpose_out ? 1 : volume;
  key[0] = 3 + FFT_KEY_INPLACE * inplace + direction * FFT_KEY_FORWARD;
  key[1] = cp_mpi_comm_c2f(cp_mpi_get_comm_null());
  key[2] = number_of_threads;
  set_nd_key(key, 3, fft_size, number_of_ffts, direction, inplace, false,
             istride, idist, ostride, odist);
}

/*******************************************************************************
 * \brief Get key from FFT input parameters (3D, R2C-case)
 * \author Frederick Stein
 ******************************************************************************/
static inline void get_key_3d_r2c(const bool direction, const int fft_size[3],
                                  const int number_of_ffts,
                                  const bool transpose_in,
                                  const bool transpose_out,
                                  const int number_of_threads,
                                  const bool inplace, int *key) {
  cp_assert((!inplace || (transpose_in == transpose_out)) &&
            "Inplace plans need the same strides for input and output arrays!");
  const int real_volume = fft_size[0] * fft_size[1] * fft_size[2];
  const int complex_volume = fft_size[0] * fft_size[1] * (fft_size[2] / 2 + 1);
  const int input_volume =
      direction ? (inplace ? 2 * complex_volume : real_volume) : complex_volume;
  const int output_volume =
      direction ? complex_volume : (inplace ? 2 * complex_volume : real_volume);
  const int istride = transpose_in ? number_of_ffts : 1;
  const int idist = transpose_in ? 1 : input_volume;
  const int ostride = transpose_out ? number_of_ffts : 1;
  const int odist = transpose_out ? 1 : output_volume;
  key[0] =
      3 + FFT_KEY_R2C + FFT_KEY_INPLACE * inplace + direction * FFT_KEY_FORWARD;
  key[1] = cp_mpi_comm_c2f(cp_mpi_get_comm_null());
  key[2] = number_of_threads;
  set_nd_key(key, 3, fft_size, number_of_ffts, direction, inplace, true,
             istride, idist, ostride, odist);
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
  key[0] = rank + FFT_KEY_INPLACE * inplace + direction * FFT_KEY_FORWARD;
  cp_assert(rank + howmany_rank <= 4 &&
            "Larger combined ranks than 4 are not implemented\n");

  key[1] = cp_mpi_comm_c2f(cp_mpi_get_comm_null());
  key[2] = number_of_threads;
  for (int d = 0; d < 4; d++) {
    key[3 + d] = 0;
    key[7 + d] = 0;
    key[11 + d] = 0;
  }
  for (int d = 0; d < rank; d++) {
    key[3 + d] = dims[d].n;
    key[7 + d] = dims[d].is;
    key[11 + d] = dims[d].os;
  }
  for (int d = 0; d < howmany_rank; d++) {
    key[3 + rank + d] = howmany_dims[d].n;
    key[7 + rank + d] = howmany_dims[d].is;
    key[11 + rank + d] = howmany_dims[d].os;
  }
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
  key[0] = rank + FFT_KEY_R2C + FFT_KEY_INPLACE * inplace +
           direction * FFT_KEY_FORWARD;
  cp_assert(rank + howmany_rank <= 4 &&
            "Larger combined ranks than 4 are not implemented\n");

  key[1] = cp_mpi_comm_c2f(cp_mpi_get_comm_null());
  key[2] = number_of_threads;
  for (int d = 0; d < 4; d++) {
    key[3 + d] = 0;
    key[7 + d] = 0;
    key[11 + d] = 0;
  }
  for (int d = 0; d < rank; d++) {
    key[3 + d] = dims[d].n;
    key[7 + d] = dims[d].is;
    key[11 + d] = dims[d].os;
  }
  for (int d = 0; d < howmany_rank; d++) {
    key[3 + rank + d] = howmany_dims[d].n;
    key[7 + rank + d] = howmany_dims[d].is;
    key[11 + rank + d] = howmany_dims[d].os;
  }
}

/*******************************************************************************
 * \brief Determine buffer size for a local FFT from a key
 * \author Frederick Stein
 ******************************************************************************/
static inline int get_buffer_size_from_key(const fft_key_t key) {
  const int rank = key[0] % 4;
  const bool is_r2c = (key[0] & FFT_KEY_R2C) == FFT_KEY_R2C;
  const bool direction = (key[0] & FFT_KEY_FORWARD) == FFT_KEY_FORWARD;
  const bool inplace = (key[0] & FFT_KEY_INPLACE) == FFT_KEY_INPLACE;
  int used = 4;
  while (used > 0 && key[3 + used - 1] == 0)
    used--;
  if (used == 0)
    return 0;
  int max_span = 1;
  for (int side = 0; side < 2; side++) {
    const int *strides = (side == 0) ? &key[7] : &key[11];
    int span = 1;
    for (int d = 0; d < used; d++) {
      int n = key[3 + d];
      const bool complex_side =
          (side == 0 && !direction) || (side == 1 && direction);
      if (is_r2c && d == rank - 1) {
        const int half = n / 2 + 1;
        if (complex_side)
          n = inplace ? 2 * half : half;
        else
          n = inplace ? 2 * half : n;
      }
      if (n > 1)
        span += (n - 1) * strides[d];
    }
    if (span > max_span)
      max_span = span;
  }
  return max_span;
}

/*******************************************************************************
 * \brief Decode a key for a local (non-Guru) FFT into its plan parameters.
 * \param key Encoded FFT key.
 * \param rank Number of transform dimensions.
 * \param fft_size Transform sizes (up to three dimensions).
 * \param number_of_ffts Number of batched transforms.
 * \param number_of_threads Number of FFTW threads used by the plan.
 * \param direction FFT direction (forward or backward).
 * \param inplace Whether the transform is in-place.
 * \param inembed Input array dimensions.
 * \param onembed Output array dimensions.
 * \param idist Distance between consecutive input transforms.
 * \param odist Distance between consecutive output transforms.
 * \param istride Stride between consecutive input elements.
 * \param ostride Stride between consecutive output elements.
 ******************************************************************************/
static inline void fetch_data_from_key_nd(
    const fft_key_t key, int *rank, int *fft_size, int *number_of_ffts,
    int *number_of_threads, bool *direction, bool *inplace, int *inembed,
    int *onembed, int *idist, int *odist, int *istride, int *ostride) {
  cp_assert(key[1] == cp_mpi_comm_c2f(cp_mpi_get_comm_null()) &&
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
  *number_of_ffts = key[3 + *rank];
  *idist = key[7 + *rank];
  *odist = key[11 + *rank];
  *istride = key[7];
  *ostride = key[11];
}

/*******************************************************************************
 * \brief Decode a key for a local Guru FFT into its plan parameters.
 * \param key Encoded FFT key.
 * \param direction FFT direction (forward or backward).
 * \param rank Number of transform dimensions.
 * \param dims Transform dimensions with their input and output strides.
 * \param howmany_rank Number of howmany dimensions.
 * \param howmany_dims Howmany dimensions with their in/output strides.
 * \param number_of_threads Number of FFTW threads used by the plan.
 * \param inplace Whether the transform is in-place.
 ******************************************************************************/
static inline void fetch_data_from_key_guru(const fft_key_t key,
                                            bool *direction, int *rank,
                                            fftw_iodim *dims, int *howmany_rank,
                                            fftw_iodim *howmany_dims,
                                            int *number_of_threads,
                                            bool *inplace) {
  cp_assert(key[1] == cp_mpi_comm_c2f(cp_mpi_get_comm_null()) &&
            "Distributed FFTs are not supported in this function!");
  *direction = (key[0] & FFT_KEY_FORWARD) == FFT_KEY_FORWARD;
  *inplace = (key[0] & FFT_KEY_INPLACE) == FFT_KEY_INPLACE;
  *rank = key[0] % 4;
  int used = 4;
  while (used > 0 && key[3 + used - 1] == 0)
    used--;
  *howmany_rank = used > *rank ? used - *rank : 0;
  for (int r = 0; r < *rank; r++) {
    dims[r].n = key[3 + r];
    dims[r].is = key[7 + r];
    dims[r].os = key[11 + r];
  }
  for (int r = 0; r < *howmany_rank; r++) {
    howmany_dims[r].n = key[3 + *rank + r];
    howmany_dims[r].is = key[7 + *rank + r];
    howmany_dims[r].os = key[11 + *rank + r];
  }
  *number_of_threads = key[2];
}
#endif

#endif /* FFT_KEY_H */