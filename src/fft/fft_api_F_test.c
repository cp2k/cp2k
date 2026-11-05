/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include "../mpiwrap/cp_mpi.h"
#include "fft_api_F.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


static inline int fft3d_index(const int fft_size[3], const int x, const int y,
                              const int z, const int batch,
                              const int number_of_ffts, const bool transpose) {
  const int spatial = (x * fft_size[1] + y) * fft_size[2] + z;
  return transpose ? spatial * number_of_ffts + batch
                   : batch * (fft_size[0] * fft_size[1] * fft_size[2]) + spatial;
}

static inline int fft3d_index_complex(const int fft_size[3], const int x,
                                      const int y, const int z, const int batch,
                                      const int number_of_ffts,
                                      const bool transpose) {
  const int spatial = (x * fft_size[1] + y) * (fft_size[2] / 2 + 1) + z;
  return transpose
             ? spatial * number_of_ffts + batch
             : batch * (fft_size[0] * fft_size[1] * (fft_size[2] / 2 + 1)) +
                   spatial;
}

static inline int fft3d_index_padded(const int fft_size[3], const int x,
                                     const int y, const int z, const int batch,
                                     const int number_of_ffts,
                                     const bool transpose) {
  const int spatial =
      (x * fft_size[1] + y) * 2 * (fft_size[2] / 2 + 1) + z;
  return transpose ? spatial * number_of_ffts + batch
                   : batch * (2 * fft_size[0] * fft_size[1] *
                              (fft_size[2] / 2 + 1)) +
                         spatial;
}

/*******************************************************************************
 * \brief Function to test the local FFT backend.
 * \author Frederick Stein
 ******************************************************************************/
int fft_test_F_1d_local_low(const int fft_size, const int number_of_ffts,
                            const int transpose_rs, const int transpose_gs) {
  const int my_process = cp_mpi_comm_rank(cp_mpi_get_comm_world());

  int errors = 0;

  const double pi = acos(-1);

  double complex *input_array = NULL, *output_array = NULL;
  fft_allocate_complex_F(fft_size * number_of_ffts, &input_array);
  fft_allocate_complex_F(fft_size * number_of_ffts, &output_array);

  fft_register_1d_local_F(
      true, fft_size, number_of_ffts, transpose_rs, transpose_gs,
      transpose_rs ? number_of_ffts : fft_size,
      transpose_gs ? number_of_ffts : fft_size, input_array, output_array);
  fft_register_1d_local_F(
      false, fft_size, number_of_ffts, transpose_gs, transpose_rs,
      transpose_gs ? number_of_ffts : fft_size,
      transpose_rs ? number_of_ffts : fft_size, output_array, input_array);

  memset(input_array, 0, fft_size * number_of_ffts * sizeof(double complex));

  // Check the forward FFT
  if (transpose_rs) {
#pragma omp parallel for default(none)                                         \
    shared(input_array, fft_size, number_of_ffts)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      input_array[(number_of_fft % fft_size) * number_of_ffts + number_of_fft] =
          1.0;
    }
  } else {
#pragma omp parallel for default(none)                                         \
    shared(input_array, fft_size, number_of_ffts)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      input_array[(number_of_fft % fft_size) + number_of_fft * fft_size] = 1.0;
    }
  }

  fft_1d_local_F(true, fft_size, number_of_ffts, transpose_rs, transpose_gs,
                 transpose_rs ? number_of_ffts : fft_size,
                 transpose_gs ? number_of_ffts : fft_size, input_array,
                 output_array);

  double max_error = 0.0;
  if (transpose_gs) {
#pragma omp parallel for default(none)                                         \
    shared(output_array, fft_size, number_of_ffts, pi, my_process)             \
    reduction(max : max_error) collapse(2)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index = 0; index < fft_size; index++) {
        const double complex my_value =
            output_array[number_of_fft + index * number_of_ffts];
        const double complex ref_value =
            cexp(-2.0 * I * pi * (number_of_fft % fft_size) * index / fft_size);
        const double current_error = cabs(my_value - ref_value);
        if (my_process == 0 && current_error > 1.0e-4)
          printf("ERROR %i %i/%i %i: (%f %f) (%f %f)\n", index, number_of_fft,
                 fft_size, number_of_ffts, creal(my_value), cimag(my_value),
                 creal(ref_value), cimag(ref_value));
        max_error = fmax(max_error, current_error);
      }
    }
  } else {
#pragma omp parallel for default(none)                                         \
    shared(output_array, fft_size, number_of_ffts, pi, my_process)             \
    reduction(max : max_error) collapse(2)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index = 0; index < fft_size; index++) {
        const double complex my_value =
            output_array[number_of_fft * fft_size + index];
        const double complex ref_value =
            cexp(-2.0 * I * pi * (number_of_fft % fft_size) * index / fft_size);
        const double current_error = cabs(my_value - ref_value);
        if (my_process == 0 && current_error > 1.0e-4)
          printf("ERROR %i %i/%i %i: (%f %f) (%f %f)\n", index, number_of_fft,
                 fft_size, number_of_ffts, creal(my_value), cimag(my_value),
                 creal(ref_value), cimag(ref_value));
        max_error = fmax(max_error, current_error);
      }
    }
  }
  fflush(stdout);

  if (max_error > 1.0e-4) {
    if (my_process == 0) {
      printf("The fw 1D-FFT does not work correctly (%i %i): %f!\n", fft_size,
             number_of_ffts, max_error);
      fflush(stdout);
    }
    errors++;
  }

  // Check the backward FFT
  memset(output_array, 0, fft_size * number_of_ffts * sizeof(double complex));

  if (transpose_gs) {
#pragma omp parallel for default(none)                                         \
    shared(output_array, fft_size, number_of_ffts)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      output_array[number_of_fft + number_of_fft % fft_size * number_of_ffts] =
          1.0;
    }
  } else {
#pragma omp parallel for default(none)                                         \
    shared(output_array, fft_size, number_of_ffts)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      output_array[number_of_fft * fft_size + number_of_fft % fft_size] = 1.0;
    }
  }

  fft_1d_local_F(false, fft_size, number_of_ffts, transpose_gs, transpose_rs,
                 transpose_gs ? number_of_ffts : fft_size,
                 transpose_rs ? number_of_ffts : fft_size, output_array,
                 input_array);

  max_error = 0.0;
  if (transpose_rs) {
#pragma omp parallel for default(none)                                         \
    shared(input_array, fft_size, number_of_ffts, pi, my_process)              \
    reduction(max : max_error) collapse(2)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index = 0; index < fft_size; index++) {
        const double complex my_value =
            input_array[index * number_of_ffts + number_of_fft];
        const double complex ref_value =
            cexp(2.0 * I * pi * (number_of_fft % fft_size) * index / fft_size);
        const double current_error = cabs(my_value - ref_value);
        if (my_process == 0 && current_error > 1.0e-12)
          printf("ERROR %i %i/%i %i: (%f %f) (%f %f)\n", index, number_of_fft,
                 fft_size, number_of_ffts, creal(my_value), cimag(my_value),
                 creal(ref_value), cimag(ref_value));
        max_error = fmax(max_error, current_error);
      }
    }
  } else {
#pragma omp parallel for default(none)                                         \
    shared(input_array, fft_size, number_of_ffts, pi, my_process)              \
    reduction(max : max_error) collapse(2)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index = 0; index < fft_size; index++) {
        const double complex my_value =
            input_array[index + number_of_fft * fft_size];
        const double complex ref_value =
            cexp(2.0 * I * pi * (number_of_fft % fft_size) * index / fft_size);
        const double current_error = cabs(my_value - ref_value);
        if (my_process == 0 && current_error > 1.0e-12)
          printf("ERROR %i %i/%i %i: (%f %f) (%f %f)\n", index, number_of_fft,
                 fft_size, number_of_ffts, creal(my_value), cimag(my_value),
                 creal(ref_value), cimag(ref_value));
        max_error = fmax(max_error, current_error);
      }
    }
  }
  fflush(stdout);

  fft_free_complex_F(input_array);
  fft_free_complex_F(output_array);

  if (max_error > 1e-12) {
    if (my_process == 0) {
      printf("The bw 1D FFT does not work correctly (%i %i): %f!\n", fft_size,
             number_of_ffts, max_error);
      fflush(stdout);
    }
    errors++;
  }

  if (errors == 0 && my_process == 0)
    printf("The 1D FFT does work correctly (%i %i)!\n", fft_size,
           number_of_ffts);
  return errors;
}

/*******************************************************************************
 * \brief Function to test the local FFT backend (inplace).
 * \author Frederick Stein
 ******************************************************************************/
int fft_test_F_1d_local_inplace_low(const int fft_size,
                                    const int number_of_ffts,
                                    const int transpose) {
  const int my_process = cp_mpi_comm_rank(cp_mpi_get_comm_world());

  int errors = 0;

  const double pi = acos(-1);

  double complex *input_array = NULL;
  fft_allocate_complex_F(fft_size * number_of_ffts, &input_array);

  fft_register_1d_local_inplace_F(true, fft_size, number_of_ffts, transpose,
                                  transpose ? number_of_ffts : fft_size,
                                  input_array);
  fft_register_1d_local_inplace_F(false, fft_size, number_of_ffts, transpose,
                                  transpose ? number_of_ffts : fft_size,
                                  input_array);

  memset(input_array, 0, fft_size * number_of_ffts * sizeof(double complex));

  // Check the forward FFT
  if (transpose) {
#pragma omp parallel for default(none)                                         \
    shared(input_array, fft_size, number_of_ffts)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      input_array[(number_of_fft % fft_size) * number_of_ffts + number_of_fft] =
          1.0;
    }
  } else {
#pragma omp parallel for default(none)                                         \
    shared(input_array, fft_size, number_of_ffts)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      input_array[(number_of_fft % fft_size) + number_of_fft * fft_size] = 1.0;
    }
  }

  fft_1d_local_inplace_F(true, fft_size, number_of_ffts, transpose,
                         transpose ? number_of_ffts : fft_size, input_array);

  double max_error = 0.0;
  if (transpose) {
#pragma omp parallel for default(none)                                         \
    shared(input_array, fft_size, number_of_ffts, pi, my_process)              \
    reduction(max : max_error) collapse(2)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index = 0; index < fft_size; index++) {
        const double complex my_value =
            input_array[number_of_fft + index * number_of_ffts];
        const double complex ref_value =
            cexp(-2.0 * I * pi * (number_of_fft % fft_size) * index / fft_size);
        const double current_error = cabs(my_value - ref_value);
        if (my_process == 0 && current_error > 1.0e-4)
          printf("ERROR %i %i/%i %i: (%f %f) (%f %f)\n", index, number_of_fft,
                 fft_size, number_of_ffts, creal(my_value), cimag(my_value),
                 creal(ref_value), cimag(ref_value));
        max_error = fmax(max_error, current_error);
      }
    }
  } else {
#pragma omp parallel for default(none)                                         \
    shared(input_array, fft_size, number_of_ffts, pi, my_process)              \
    reduction(max : max_error) collapse(2)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index = 0; index < fft_size; index++) {
        const double complex my_value =
            input_array[number_of_fft * fft_size + index];
        const double complex ref_value =
            cexp(-2.0 * I * pi * (number_of_fft % fft_size) * index / fft_size);
        const double current_error = cabs(my_value - ref_value);
        if (my_process == 0 && current_error > 1.0e-4)
          printf("ERROR %i %i/%i %i: (%f %f) (%f %f)\n", index, number_of_fft,
                 fft_size, number_of_ffts, creal(my_value), cimag(my_value),
                 creal(ref_value), cimag(ref_value));
        max_error = fmax(max_error, current_error);
      }
    }
  }
  fflush(stdout);

  if (max_error > 1.0e-4) {
    if (my_process == 0) {
      printf("The fw 1D-FFT (inplace) does not work correctly (%i %i): %f!\n",
             fft_size, number_of_ffts, max_error);
      fflush(stdout);
    }
    errors++;
  }

  // Check the backward FFT
  memset(input_array, 0, fft_size * number_of_ffts * sizeof(double complex));

  if (transpose) {
#pragma omp parallel for default(none)                                         \
    shared(input_array, fft_size, number_of_ffts)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      input_array[number_of_fft + number_of_fft % fft_size * number_of_ffts] =
          1.0;
    }
  } else {
#pragma omp parallel for default(none)                                         \
    shared(input_array, fft_size, number_of_ffts)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      input_array[number_of_fft * fft_size + number_of_fft % fft_size] = 1.0;
    }
  }

  fft_1d_local_inplace_F(false, fft_size, number_of_ffts, transpose,
                         transpose ? number_of_ffts : fft_size, input_array);

  max_error = 0.0;
  if (transpose) {
#pragma omp parallel for default(none)                                         \
    shared(input_array, fft_size, number_of_ffts, pi, my_process)              \
    reduction(max : max_error) collapse(2)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index = 0; index < fft_size; index++) {
        const double complex my_value =
            input_array[index * number_of_ffts + number_of_fft];
        const double complex ref_value =
            cexp(2.0 * I * pi * (number_of_fft % fft_size) * index / fft_size);
        const double current_error = cabs(my_value - ref_value);
        if (my_process == 0 && current_error > 1.0e-12)
          printf("ERROR %i %i/%i %i: (%f %f) (%f %f)\n", index, number_of_fft,
                 fft_size, number_of_ffts, creal(my_value), cimag(my_value),
                 creal(ref_value), cimag(ref_value));
        max_error = fmax(max_error, current_error);
      }
    }
  } else {
#pragma omp parallel for default(none)                                         \
    shared(input_array, fft_size, number_of_ffts, pi, my_process)              \
    reduction(max : max_error) collapse(2)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index = 0; index < fft_size; index++) {
        const double complex my_value =
            input_array[index + number_of_fft * fft_size];
        const double complex ref_value =
            cexp(2.0 * I * pi * (number_of_fft % fft_size) * index / fft_size);
        const double current_error = cabs(my_value - ref_value);
        if (my_process == 0 && current_error > 1.0e-12)
          printf("ERROR %i %i/%i %i: (%f %f) (%f %f)\n", index, number_of_fft,
                 fft_size, number_of_ffts, creal(my_value), cimag(my_value),
                 creal(ref_value), cimag(ref_value));
        max_error = fmax(max_error, current_error);
      }
    }
  }
  fflush(stdout);

  fft_free_complex_F(input_array);

  if (max_error > 1e-12) {
    if (my_process == 0) {
      printf("The bw 1D FFT (inplace) does not work correctly (%i %i): %f!\n",
             fft_size, number_of_ffts, max_error);
      fflush(stdout);
    }
    errors++;
  }

  if (errors == 0 && my_process == 0)
    printf("The 1D FFT (inplace) does work correctly (%i %i)!\n", fft_size,
           number_of_ffts);
  return errors;
}

/*******************************************************************************
 * \brief Function to test the local FFT backend.
 * \author Frederick Stein
 ******************************************************************************/
int fft_test_F_2d_local_low(const int fft_size[2], const int number_of_ffts,
                            const int transpose_rs, const int transpose_gs) {
  const int my_process = cp_mpi_comm_rank(cp_mpi_get_comm_world());

  int errors = 0;

  const double pi = acos(-1);

  double complex *input_array = NULL, *output_array = NULL;
  const int elements_per_fft = fft_size[0] * fft_size[1];
  fft_allocate_complex_F(elements_per_fft * number_of_ffts, &input_array);
  fft_allocate_complex_F(elements_per_fft * number_of_ffts, &output_array);

  const int fft_size_F[] = {fft_size[1], fft_size[0]};

  fft_register_2d_local_F(true, fft_size_F, number_of_ffts, transpose_rs,
                          transpose_gs, input_array, output_array);
  fft_register_2d_local_F(false, fft_size_F, number_of_ffts, transpose_gs,
                          transpose_rs, output_array, input_array);

  memset(input_array, 0,
         elements_per_fft * number_of_ffts * sizeof(double complex));
  double max_error = 0.0;
  // Check the forward FFT
  if (transpose_rs) {
#pragma omp parallel for default(none)                                         \
    shared(input_array, elements_per_fft, number_of_ffts)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      input_array[number_of_fft % elements_per_fft * number_of_ffts +
                  number_of_fft] = 1.0;
    }
  } else {
#pragma omp parallel for default(none)                                         \
    shared(input_array, elements_per_fft, number_of_ffts)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      input_array[number_of_fft % elements_per_fft +
                  number_of_fft * elements_per_fft] = 1.0;
    }
  }

  fft_2d_local_F(true, fft_size_F, number_of_ffts, transpose_rs, transpose_gs,
                 input_array, output_array);

  if (transpose_gs) {
#pragma omp parallel for default(none)                                         \
    shared(output_array, fft_size, number_of_ffts, pi, my_process)             \
    reduction(max : max_error) collapse(3)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index_0 = 0; index_0 < fft_size[0]; index_0++) {
        for (int index_1 = 0; index_1 < fft_size[1]; index_1++) {
          const double complex my_value =
              output_array[number_of_fft +
                           (index_0 * fft_size[1] + index_1) * number_of_ffts];
          const double complex ref_value = cexp(
              -2.0 * I * pi *
              ((double)(number_of_fft / fft_size[1]) * index_0 / fft_size[0] +
               (double)(number_of_fft % fft_size[1]) * index_1 / fft_size[1]));
          double current_error = cabs(my_value - ref_value);
          if (my_process == 0 && current_error > 1.0e-4)
            printf("ERROR %i %i %i/%i %i %i: (%f %f) (%f %f)\n", index_0,
                   index_1, number_of_fft, fft_size[0], fft_size[1],
                   number_of_ffts, creal(my_value), cimag(my_value),
                   creal(ref_value), cimag(ref_value));
          max_error = fmax(max_error, current_error);
        }
      }
    }
  } else {
#pragma omp parallel for default(none)                                         \
    shared(output_array, fft_size, number_of_ffts, pi, my_process)             \
    reduction(max : max_error) collapse(3)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index_0 = 0; index_0 < fft_size[0]; index_0++) {
        for (int index_1 = 0; index_1 < fft_size[1]; index_1++) {
          const double complex my_value =
              output_array[(number_of_fft * fft_size[0] + index_0) *
                               fft_size[1] +
                           index_1];
          const double complex ref_value = cexp(
              -2.0 * I * pi *
              ((double)(number_of_fft / fft_size[1]) * index_0 / fft_size[0] +
               (double)(number_of_fft % fft_size[1]) * index_1 / fft_size[1]));
          double current_error = cabs(my_value - ref_value);
          if (my_process == 0 && current_error > 1.0e-3)
            printf("ERROR %i %i %i/%i %i %i: (%f %f) (%f %f)\n", index_0,
                   index_1, number_of_fft, fft_size[0], fft_size[1],
                   number_of_ffts, creal(my_value), cimag(my_value),
                   creal(ref_value), cimag(ref_value));
          max_error = fmax(max_error, current_error);
        }
      }
    }
  }
  fflush(stdout);

  if (max_error > 1.0e-3) {
    if (my_process == 0) {
      printf("The fw 2D-FFT does not work correctly (%i %i/%i): %f!\n",
             fft_size[0], fft_size[1], number_of_ffts, max_error);
      fflush(stdout);
    }
    errors++;
  }

  // Check the backward FFT
  memset(output_array, 0,
         elements_per_fft * number_of_ffts * sizeof(double complex));

  if (transpose_gs) {
#pragma omp parallel for default(none)                                         \
    shared(output_array, elements_per_fft, number_of_ffts)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      output_array[number_of_fft +
                   (number_of_fft % elements_per_fft) * number_of_ffts] = 1.0;
    }
  } else {
#pragma omp parallel for default(none)                                         \
    shared(output_array, elements_per_fft, number_of_ffts)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      output_array[number_of_fft * elements_per_fft +
                   number_of_fft % elements_per_fft] = 1.0;
    }
  }

  fft_2d_local_F(false, fft_size_F, number_of_ffts, transpose_gs, transpose_rs,
                 output_array, input_array);

  max_error = 0.0;
  if (transpose_rs) {
#pragma omp parallel for default(none)                                         \
    shared(input_array, fft_size, number_of_ffts, pi, my_process)              \
    reduction(max : max_error) collapse(3)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index_0 = 0; index_0 < fft_size[0]; index_0++) {
        for (int index_1 = 0; index_1 < fft_size[1]; index_1++) {
          const double complex my_value =
              input_array[(index_0 * fft_size[1] + index_1) * number_of_ffts +
                          number_of_fft];
          const double complex ref_value = cexp(
              2.0 * I * pi *
              ((double)(number_of_fft / fft_size[1]) * index_0 / fft_size[0] +
               (double)(number_of_fft % fft_size[1]) * index_1 / fft_size[1]));
          double current_error = cabs(my_value - ref_value);
          if (my_process == 0 && current_error > 1e-12)
            printf("Error %i %i %i: (%f %f) (%f %f)\n", index_0, index_1,
                   number_of_fft, creal(my_value), cimag(my_value),
                   creal(ref_value), cimag(ref_value));
          max_error = fmax(max_error, current_error);
        }
      }
    }
  } else {
#pragma omp parallel for default(none)                                         \
    shared(input_array, fft_size, number_of_ffts, pi, my_process)              \
    reduction(max : max_error) collapse(3)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index_0 = 0; index_0 < fft_size[0]; index_0++) {
        for (int index_1 = 0; index_1 < fft_size[1]; index_1++) {
          const double complex my_value =
              input_array[index_0 * fft_size[1] + index_1 +
                          number_of_fft * fft_size[0] * fft_size[1]];
          const double complex ref_value = cexp(
              2.0 * I * pi *
              ((double)(number_of_fft / fft_size[1]) * index_0 / fft_size[0] +
               (double)(number_of_fft % fft_size[1]) * index_1 / fft_size[1]));
          double current_error = cabs(my_value - ref_value);
          if (my_process == 0 && current_error > 1e-12)
            printf("Error %i %i %i: (%f %f) (%f %f)\n", index_0, index_1,
                   number_of_fft, creal(my_value), cimag(my_value),
                   creal(ref_value), cimag(ref_value));
          max_error = fmax(max_error, current_error);
        }
      }
    }
  }
  fflush(stdout);

  fft_free_complex_F(input_array);
  fft_free_complex_F(output_array);

  if (max_error > 1e-12) {
    if (my_process == 0) {
      printf("The bw 2D-FFT does not work correctly (%i %i/%i): %f!\n",
             fft_size[0], fft_size[1], number_of_ffts, max_error);
      fflush(stdout);
    }
    errors++;
  }

  if (errors == 0 && my_process == 0)
    printf("The 2D FFT does work correctly (%i %i/%i)!\n", fft_size[0],
           fft_size[1], number_of_ffts);
  return errors;
}

/*******************************************************************************
 * \brief Function to test the local FFT backend.
 * \author Frederick Stein
 ******************************************************************************/
int fft_test_F_2d_local_inplace_low(const int fft_size[2],
                                    const int number_of_ffts,
                                    const int transpose) {
  const int my_process = cp_mpi_comm_rank(cp_mpi_get_comm_world());

  int errors = 0;

  const double pi = acos(-1);

  double complex *input_array = NULL;
  const int elements_per_fft = fft_size[0] * fft_size[1];
  fft_allocate_complex_F(elements_per_fft * number_of_ffts, &input_array);
  double complex *output_array = input_array;

  const int fft_size_F[] = {fft_size[1], fft_size[0]};

  fft_register_2d_local_inplace_F(true, fft_size_F, number_of_ffts, transpose,
                                  input_array);
  fft_register_2d_local_inplace_F(false, fft_size_F, number_of_ffts, transpose,
                                  output_array);

  memset(input_array, 0,
         elements_per_fft * number_of_ffts * sizeof(double complex));
  double max_error = 0.0;
  // Check the forward FFT
  if (transpose) {
#pragma omp parallel for default(none)                                         \
    shared(input_array, elements_per_fft, number_of_ffts)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      input_array[number_of_fft % elements_per_fft * number_of_ffts +
                  number_of_fft] = 1.0;
    }
  } else {
#pragma omp parallel for default(none)                                         \
    shared(input_array, elements_per_fft, number_of_ffts)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      input_array[number_of_fft % elements_per_fft +
                  number_of_fft * elements_per_fft] = 1.0;
    }
  }

  fft_2d_local_inplace_F(true, fft_size_F, number_of_ffts, transpose,
                         input_array);

  if (transpose) {
#pragma omp parallel for default(none)                                         \
    shared(output_array, fft_size, number_of_ffts, pi, my_process)             \
    reduction(max : max_error) collapse(3)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index_0 = 0; index_0 < fft_size[0]; index_0++) {
        for (int index_1 = 0; index_1 < fft_size[1]; index_1++) {
          const double complex my_value =
              output_array[number_of_fft +
                           (index_0 * fft_size[1] + index_1) * number_of_ffts];
          const double complex ref_value = cexp(
              -2.0 * I * pi *
              ((double)(number_of_fft / fft_size[1]) * index_0 / fft_size[0] +
               (double)(number_of_fft % fft_size[1]) * index_1 / fft_size[1]));
          double current_error = cabs(my_value - ref_value);
          if (my_process == 0 && current_error > 1.0e-4)
            printf("ERROR %i %i %i/%i %i %i: (%f %f) (%f %f)\n", index_0,
                   index_1, number_of_fft, fft_size[0], fft_size[1],
                   number_of_ffts, creal(my_value), cimag(my_value),
                   creal(ref_value), cimag(ref_value));
          max_error = fmax(max_error, current_error);
        }
      }
    }
  } else {
#pragma omp parallel for default(none)                                         \
    shared(output_array, fft_size, number_of_ffts, pi, my_process)             \
    reduction(max : max_error) collapse(3)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index_0 = 0; index_0 < fft_size[0]; index_0++) {
        for (int index_1 = 0; index_1 < fft_size[1]; index_1++) {
          const double complex my_value =
              output_array[(number_of_fft * fft_size[0] + index_0) *
                               fft_size[1] +
                           index_1];
          const double complex ref_value = cexp(
              -2.0 * I * pi *
              ((double)(number_of_fft / fft_size[1]) * index_0 / fft_size[0] +
               (double)(number_of_fft % fft_size[1]) * index_1 / fft_size[1]));
          double current_error = cabs(my_value - ref_value);
          if (my_process == 0 && current_error > 1.0e-3)
            printf("ERROR %i %i %i/%i %i %i: (%f %f) (%f %f)\n", index_0,
                   index_1, number_of_fft, fft_size[0], fft_size[1],
                   number_of_ffts, creal(my_value), cimag(my_value),
                   creal(ref_value), cimag(ref_value));
          max_error = fmax(max_error, current_error);
        }
      }
    }
  }
  fflush(stdout);

  if (max_error > 1.0e-3) {
    if (my_process == 0) {
      printf(
          "The fw 2D-FFT (inplace) does not work correctly (%i %i/%i): %f!\n",
          fft_size[0], fft_size[1], number_of_ffts, max_error);
      fflush(stdout);
    }
    errors++;
  }

  // Check the backward FFT
  memset(output_array, 0,
         elements_per_fft * number_of_ffts * sizeof(double complex));

  if (transpose) {
#pragma omp parallel for default(none)                                         \
    shared(output_array, elements_per_fft, number_of_ffts)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      output_array[number_of_fft +
                   (number_of_fft % elements_per_fft) * number_of_ffts] = 1.0;
    }
  } else {
#pragma omp parallel for default(none)                                         \
    shared(output_array, elements_per_fft, number_of_ffts)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      output_array[number_of_fft * elements_per_fft +
                   number_of_fft % elements_per_fft] = 1.0;
    }
  }

  fft_2d_local_inplace_F(false, fft_size_F, number_of_ffts, transpose,
                         output_array);

  max_error = 0.0;
  if (transpose) {
#pragma omp parallel for default(none)                                         \
    shared(input_array, fft_size, number_of_ffts, pi, my_process)              \
    reduction(max : max_error) collapse(3)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index_0 = 0; index_0 < fft_size[0]; index_0++) {
        for (int index_1 = 0; index_1 < fft_size[1]; index_1++) {
          const double complex my_value =
              input_array[(index_0 * fft_size[1] + index_1) * number_of_ffts +
                          number_of_fft];
          const double complex ref_value = cexp(
              2.0 * I * pi *
              ((double)(number_of_fft / fft_size[1]) * index_0 / fft_size[0] +
               (double)(number_of_fft % fft_size[1]) * index_1 / fft_size[1]));
          double current_error = cabs(my_value - ref_value);
          if (my_process == 0 && current_error > 1e-12)
            printf("Error %i %i %i: (%f %f) (%f %f)\n", index_0, index_1,
                   number_of_fft, creal(my_value), cimag(my_value),
                   creal(ref_value), cimag(ref_value));
          max_error = fmax(max_error, current_error);
        }
      }
    }
  } else {
#pragma omp parallel for default(none)                                         \
    shared(input_array, fft_size, number_of_ffts, pi, my_process)              \
    reduction(max : max_error) collapse(3)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index_0 = 0; index_0 < fft_size[0]; index_0++) {
        for (int index_1 = 0; index_1 < fft_size[1]; index_1++) {
          const double complex my_value =
              input_array[index_0 * fft_size[1] + index_1 +
                          number_of_fft * fft_size[0] * fft_size[1]];
          const double complex ref_value = cexp(
              2.0 * I * pi *
              ((double)(number_of_fft / fft_size[1]) * index_0 / fft_size[0] +
               (double)(number_of_fft % fft_size[1]) * index_1 / fft_size[1]));
          double current_error = cabs(my_value - ref_value);
          if (my_process == 0 && current_error > 1e-12)
            printf("Error %i %i %i: (%f %f) (%f %f)\n", index_0, index_1,
                   number_of_fft, creal(my_value), cimag(my_value),
                   creal(ref_value), cimag(ref_value));
          max_error = fmax(max_error, current_error);
        }
      }
    }
  }
  fflush(stdout);

  fft_free_complex_F(input_array);

  if (max_error > 1e-12) {
    if (my_process == 0) {
      printf(
          "The bw 2D-FFT (inplace) does not work correctly (%i %i/%i): %f!\n",
          fft_size[0], fft_size[1], number_of_ffts, max_error);
      fflush(stdout);
    }
    errors++;
  }

  if (errors == 0 && my_process == 0)
    printf("The 2D FFT (inplace) does work correctly (%i %i/%i)!\n",
           fft_size[0], fft_size[1], number_of_ffts);
  return errors;
}

/*******************************************************************************
 * \brief Function to test the local FFT backend.
 * \author Frederick Stein
 ******************************************************************************/
int fft_test_F_3d_local_low(const int fft_size[3], const int test_every,
                            const bool transpose_in, const bool transpose_out) {
  (void)test_every;
  const int my_process = cp_mpi_comm_rank(cp_mpi_get_comm_world());

  int errors = 0;

  const double pi = acos(-1);
  const int volume = fft_size[0] * fft_size[1] * fft_size[2];
  const int number_of_ffts =
      test_every > 0 ? (volume + test_every - 1) / test_every : volume;
  const int fft_size_F[] = {fft_size[2], fft_size[1], fft_size[0]};

  double complex *input_array = NULL, *output_array = NULL;
  fft_allocate_complex_F(volume * number_of_ffts, &input_array);
  fft_allocate_complex_F(volume * number_of_ffts, &output_array);

  fft_register_3d_local_F(true, fft_size_F, number_of_ffts, transpose_in, transpose_out,
                        input_array, output_array);
  fft_register_3d_local_F(false, fft_size_F, number_of_ffts, transpose_out, transpose_in,
                        output_array, input_array);

  double max_error = 0.0;
  memset(input_array, 0, volume * number_of_ffts * sizeof(double complex));
#pragma omp parallel for collapse(3)
  for (int mx = 0; mx < fft_size[0]; mx++) {
    for (int my = 0; my < fft_size[1]; my++) {
      for (int mz = 0; mz < fft_size[2]; mz++) {
        const int index = (mx * fft_size[1] + my) * fft_size[2] + mz;
        if (test_every > 0 && index % test_every != 0)
          continue;
        const int batch = test_every > 0 ? index / test_every : index;
        input_array[fft3d_index(fft_size, mx, my, mz, batch, number_of_ffts, transpose_in)] = 1.0;
      }
    }
  }
  fft_3d_local_F(true, fft_size_F, number_of_ffts, transpose_in, transpose_out, input_array,
               output_array);

#pragma omp parallel for default(none) firstprivate(test_every) firstprivate(transpose_in, transpose_out)                                         \
    shared(output_array, fft_size, pi, my_process, volume, number_of_ffts)     \
    reduction(max : max_error) collapse(3)
  for (int mx = 0; mx < fft_size[0]; mx++) {
    for (int my = 0; my < fft_size[1]; my++) {
      for (int mz = 0; mz < fft_size[2]; mz++) {
        const int basis = (mx * fft_size[1] + my) * fft_size[2] + mz;
        if (test_every > 0 && basis % test_every != 0)
          continue;
        const int batch = test_every > 0 ? basis / test_every : basis;
        for (int nx = 0; nx < fft_size[0]; nx++) {
          for (int ny = 0; ny < fft_size[1]; ny++) {
            for (int nz = 0; nz < fft_size[2]; nz++) {
              const int index = (nx * fft_size[1] + ny) * fft_size[2] + nz;
              const double complex my_value =
                  output_array[fft3d_index(fft_size, nx, ny, nz, batch, number_of_ffts, transpose_out)];
              const double complex ref_value =
                  cexp(-2.0 * I * pi *
                       (((double)mx) * nx / fft_size[0] +
                        ((double)my) * ny / fft_size[1] +
                        ((double)mz) * nz / fft_size[2]));
              const double current_error = cabs(my_value - ref_value);
              if (my_process == 0 && current_error > 1e-6)
                printf("Error %i %i %i/%i %i %i: (%f %f) (%f %f)\n", nx, ny, nz,
                       mx, my, mz, creal(my_value), cimag(my_value),
                       creal(ref_value), cimag(ref_value));
              max_error = fmax(max_error, current_error);
            }
          }
        }
      }
    }
  }
  fflush(stdout);

  if (max_error > 1.0e-12) {
    if (my_process == 0) {
      printf("The fw 3D-FFT does not work correctly (%i %i %i): %f!\n",
             fft_size[0], fft_size[1], fft_size[2], max_error);
      fflush(stdout);
    }
    errors++;
  }

  max_error = 0.0;
  memset(output_array, 0, volume * number_of_ffts * sizeof(double complex));
#pragma omp parallel for collapse(3)
  for (int mx = 0; mx < fft_size[0]; mx++) {
    for (int my = 0; my < fft_size[1]; my++) {
      for (int mz = 0; mz < fft_size[2]; mz++) {
        const int index = (mx * fft_size[1] + my) * fft_size[2] + mz;
        if (test_every > 0 && index % test_every != 0)
          continue;
        const int batch = test_every > 0 ? index / test_every : index;
        output_array[fft3d_index(fft_size, mx, my, mz, batch, number_of_ffts, transpose_out)] = 1.0;
      }
    }
  }
  fft_3d_local_F(false, fft_size_F, number_of_ffts, transpose_out, transpose_in, output_array,
               input_array);

#pragma omp parallel for default(none) firstprivate(test_every) firstprivate(transpose_in, transpose_out)                                         \
    shared(input_array, fft_size, pi, my_process, volume, number_of_ffts)      \
    reduction(max : max_error) collapse(3)
  for (int mx = 0; mx < fft_size[0]; mx++) {
    for (int my = 0; my < fft_size[1]; my++) {
      for (int mz = 0; mz < fft_size[2]; mz++) {
        const int basis = (mx * fft_size[1] + my) * fft_size[2] + mz;
        if (test_every > 0 && basis % test_every != 0)
          continue;
        const int batch = test_every > 0 ? basis / test_every : basis;
        for (int nx = 0; nx < fft_size[0]; nx++) {
          for (int ny = 0; ny < fft_size[1]; ny++) {
            for (int nz = 0; nz < fft_size[2]; nz++) {
              const int index = (nx * fft_size[1] + ny) * fft_size[2] + nz;
              const double complex my_value =
                  input_array[fft3d_index(fft_size, nx, ny, nz, batch, number_of_ffts, transpose_in)];
              const double complex ref_value =
                  cexp(2.0 * I * pi *
                       (((double)mx) * nx / fft_size[0] +
                        ((double)my) * ny / fft_size[1] +
                        ((double)mz) * nz / fft_size[2]));
              const double current_error = cabs(my_value - ref_value);
              if (my_process == 0 && current_error > 1e-6)
                printf("Error %i %i %i/%i %i %i: (%f %f) (%f %f)\n", nx, ny, nz,
                       mx, my, mz, creal(my_value), cimag(my_value),
                       creal(ref_value), cimag(ref_value));
              max_error = fmax(max_error, current_error);
            }
          }
        }
      }
    }
  }
  fflush(stdout);

  fft_free_complex_F(input_array);
  fft_free_complex_F(output_array);

  if (max_error > 1e-12) {
    if (my_process == 0) {
      printf("The bw 3D-FFT does not work correctly (%i %i %i): %f!\n",
             fft_size[0], fft_size[1], fft_size[2], max_error);
      fflush(stdout);
    }
    errors++;
  }

  if (errors == 0 && my_process == 0)
    printf("The 3D FFT does work correctly (%i %i %i)!\n", fft_size[0],
           fft_size[1], fft_size[2]);
  return errors;
}

/*******************************************************************************
 * \brief Function to test the local FFT backend.
 * \author Frederick Stein
 ******************************************************************************/
int fft_test_F_3d_local_inplace_low(const int fft_size[3], const int test_every,
                                    const bool transpose) {
  (void)test_every;
  const int my_process = cp_mpi_comm_rank(cp_mpi_get_comm_world());

  int errors = 0;

  const double pi = acos(-1);
  const int volume = fft_size[0] * fft_size[1] * fft_size[2];
  const int number_of_ffts =
      test_every > 0 ? (volume + test_every - 1) / test_every : volume;
  const int fft_size_F[] = {fft_size[2], fft_size[1], fft_size[0]};

  double complex *input_array = NULL;
  fft_allocate_complex_F(volume * number_of_ffts, &input_array);
  double complex *output_array = input_array;

  fft_register_3d_local_inplace_F(true, fft_size_F, number_of_ffts, transpose,
                                   input_array);
  fft_register_3d_local_inplace_F(false, fft_size_F, number_of_ffts, transpose,
                                   output_array);

  double max_error = 0.0;
  memset(input_array, 0, volume * number_of_ffts * sizeof(double complex));
#pragma omp parallel for collapse(3)
  for (int mx = 0; mx < fft_size[0]; mx++) {
    for (int my = 0; my < fft_size[1]; my++) {
      for (int mz = 0; mz < fft_size[2]; mz++) {
        const int index = (mx * fft_size[1] + my) * fft_size[2] + mz;
        if (test_every > 0 && index % test_every != 0)
          continue;
        const int batch = test_every > 0 ? index / test_every : index;
        input_array[fft3d_index(fft_size, mx, my, mz, batch, number_of_ffts, transpose)] = 1.0;
      }
    }
  }
  fft_3d_local_inplace_F(true, fft_size_F, number_of_ffts, transpose,
                          input_array);

#pragma omp parallel for default(none) firstprivate(test_every) firstprivate(transpose)                                         \
    shared(output_array, fft_size, pi, volume, number_of_ffts)                  \
    reduction(max : max_error) collapse(3)
  for (int mx = 0; mx < fft_size[0]; mx++) {
    for (int my = 0; my < fft_size[1]; my++) {
      for (int mz = 0; mz < fft_size[2]; mz++) {
        const int basis = (mx * fft_size[1] + my) * fft_size[2] + mz;
        if (test_every > 0 && basis % test_every != 0)
          continue;
        const int batch = test_every > 0 ? basis / test_every : basis;
        for (int nx = 0; nx < fft_size[0]; nx++) {
          for (int ny = 0; ny < fft_size[1]; ny++) {
            for (int nz = 0; nz < fft_size[2]; nz++) {
              const int index = (nx * fft_size[1] + ny) * fft_size[2] + nz;
              const double complex my_value =
                  output_array[fft3d_index(fft_size, nx, ny, nz, batch, number_of_ffts, transpose)];
              const double complex ref_value =
                  cexp(-2.0 * I * pi *
                       (((double)mx) * nx / fft_size[0] +
                        ((double)my) * ny / fft_size[1] +
                        ((double)mz) * nz / fft_size[2]));
              max_error = fmax(max_error, cabs(my_value - ref_value));
            }
          }
        }
      }
    }
  }
  fflush(stdout);

  if (max_error > 1.0e-12) {
    if (my_process == 0) {
      printf(
          "The fw 3D-FFT (inplace) does not work correctly (%i %i %i): %f!\n",
          fft_size[0], fft_size[1], fft_size[2], max_error);
      fflush(stdout);
    }
    errors++;
  }

  max_error = 0.0;
  memset(output_array, 0, volume * number_of_ffts * sizeof(double complex));
#pragma omp parallel for collapse(3)
  for (int mx = 0; mx < fft_size[0]; mx++) {
    for (int my = 0; my < fft_size[1]; my++) {
      for (int mz = 0; mz < fft_size[2]; mz++) {
        const int index = (mx * fft_size[1] + my) * fft_size[2] + mz;
        if (test_every > 0 && index % test_every != 0)
          continue;
        const int batch = test_every > 0 ? index / test_every : index;
        output_array[fft3d_index(fft_size, mx, my, mz, batch, number_of_ffts, transpose)] = 1.0;
      }
    }
  }
  fft_3d_local_inplace_F(false, fft_size_F, number_of_ffts, transpose,
                           output_array);

#pragma omp parallel for default(none) firstprivate(test_every) firstprivate(transpose)                                         \
    shared(input_array, fft_size, pi, volume, number_of_ffts)                   \
    reduction(max : max_error) collapse(3)
  for (int mx = 0; mx < fft_size[0]; mx++) {
    for (int my = 0; my < fft_size[1]; my++) {
      for (int mz = 0; mz < fft_size[2]; mz++) {
        const int basis = (mx * fft_size[1] + my) * fft_size[2] + mz;
        if (test_every > 0 && basis % test_every != 0)
          continue;
        const int batch = test_every > 0 ? basis / test_every : basis;
        for (int nx = 0; nx < fft_size[0]; nx++) {
          for (int ny = 0; ny < fft_size[1]; ny++) {
            for (int nz = 0; nz < fft_size[2]; nz++) {
              const int index = (nx * fft_size[1] + ny) * fft_size[2] + nz;
              const double complex my_value =
                  input_array[fft3d_index(fft_size, nx, ny, nz, batch, number_of_ffts, transpose)];
              const double complex ref_value =
                  cexp(2.0 * I * pi *
                       (((double)mx) * nx / fft_size[0] +
                        ((double)my) * ny / fft_size[1] +
                        ((double)mz) * nz / fft_size[2]));
              max_error = fmax(max_error, cabs(my_value - ref_value));
            }
          }
        }
      }
    }
  }
  fflush(stdout);

  fft_free_complex_F(input_array);

  if (max_error > 1e-12) {
    if (my_process == 0) {
      printf(
          "The bw 3D-FFT (inplace) does not work correctly (%i %i %i): %f!\n",
          fft_size[0], fft_size[1], fft_size[2], max_error);
      fflush(stdout);
    }
    errors++;
  }

  if (errors == 0 && my_process == 0)
    printf("The 3D FFT (inplace) does work correctly (%i %i %i)!\n",
           fft_size[0], fft_size[1], fft_size[2]);
  return errors;
}

/*******************************************************************************
 * \brief Function to test the Fortran.
 * \author Frederick Stein
 ******************************************************************************/
int fft_test_api_local() {
  int errors = 0;

  fft_library_init_F(FFT_LIBRARY_BACKEND_DEFAULT,
                     FFT_LIBRARY_FFTW_MODE_ESTIMATE, true, NULL);

  clock_t begin = clock();
  errors += fft_test_F_1d_local_low(15, 13, true, true);
  errors += fft_test_F_1d_local_low(9, 11, true, false);
  errors += fft_test_F_1d_local_low(10, 14, false, true);
  errors += fft_test_F_1d_local_low(14, 13, false, false);

  errors += fft_test_F_1d_local_inplace_low(15, 13, true);
  errors += fft_test_F_1d_local_inplace_low(14, 13, false);

  errors += fft_test_F_2d_local_low((const int[2]){10, 10}, 10, true, true);
  errors += fft_test_F_2d_local_low((const int[2]){8, 9}, 10, true, false);
  errors += fft_test_F_2d_local_low((const int[2]){7, 10}, 11, false, true);
  errors += fft_test_F_2d_local_low((const int[2]){6, 7}, 10, false, false);

  errors += fft_test_F_2d_local_inplace_low((const int[2]){8, 9}, 10, true);
  errors += fft_test_F_2d_local_inplace_low((const int[2]){7, 10}, 11, false);

  // Reduce tests to ca 10 per set
  errors += fft_test_F_3d_local_low((const int[3]){8, 8, 8}, 10, false, false);
  errors += fft_test_F_3d_local_low((const int[3]){7, 5, 3}, 10, true, true);

  errors += fft_test_F_3d_local_inplace_low((const int[3]){8, 8, 8}, 10, false);
  errors += fft_test_F_3d_local_inplace_low((const int[3]){7, 5, 3}, 10, true);

  clock_t end = clock();
  if (cp_mpi_comm_rank(cp_mpi_get_comm_world()) == 0)
    printf("Time to test local FFTs: %f\n",
           (double)(end - begin) / CLOCKS_PER_SEC);

  fft_print_timing_report_F(0.01);

  fft_library_finalize_F(NULL);

  return errors;
}

// EOF
