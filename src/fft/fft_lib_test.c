/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include "fft_lib_test.h"
#include "fft_utils.h"

#include "../mpiwrap/cp_mpi.h"
#include "fft_lib.h"
#include "fft_utils.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/*******************************************************************************
 * \brief Function to test the local FFT backend.
 * \author Frederick Stein
 ******************************************************************************/
int fft_test_1d_local_low(const int fft_size, const int number_of_ffts,
                          const int transpose_rs, const int transpose_gs) {
  const int my_process = cp_mpi_comm_rank(cp_mpi_get_comm_world());

  int errors = 0;

  const double pi = acos(-1);

  double complex *input_array = NULL, *output_array = NULL;
  fft_allocate_complex(fft_size * number_of_ffts, &input_array);
  fft_allocate_complex(fft_size * number_of_ffts, &output_array);

  fft_register_1d_local(true, fft_size, number_of_ffts, transpose_rs,
                        transpose_gs, transpose_rs ? number_of_ffts : fft_size,
                        transpose_gs ? number_of_ffts : fft_size, input_array,
                        output_array);
  fft_register_1d_local(false, fft_size, number_of_ffts, transpose_gs,
                        transpose_rs, transpose_gs ? number_of_ffts : fft_size,
                        transpose_rs ? number_of_ffts : fft_size, output_array,
                        input_array);

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

  fft_1d_local(true, fft_size, number_of_ffts, transpose_rs, transpose_gs,
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

  fft_1d_local(false, fft_size, number_of_ffts, transpose_gs, transpose_rs,
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

  fft_free_complex(input_array);
  fft_free_complex(output_array);

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
int fft_test_1d_local_inplace_low(const int fft_size, const int number_of_ffts,
                                  const int transpose) {
  const int my_process = cp_mpi_comm_rank(cp_mpi_get_comm_world());

  int errors = 0;

  const double pi = acos(-1);

  double complex *input_array = NULL;
  fft_allocate_complex(fft_size * number_of_ffts, &input_array);

  fft_register_1d_local(true, fft_size, number_of_ffts, transpose, transpose,
                        transpose ? number_of_ffts : fft_size,
                        transpose ? number_of_ffts : fft_size, input_array,
                        input_array);
  fft_register_1d_local(false, fft_size, number_of_ffts, transpose, transpose,
                        transpose ? number_of_ffts : fft_size,
                        transpose ? number_of_ffts : fft_size, input_array,
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

  fft_1d_local(true, fft_size, number_of_ffts, transpose, transpose,
               transpose ? number_of_ffts : fft_size,
               transpose ? number_of_ffts : fft_size, input_array, input_array);

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

  fft_1d_local(false, fft_size, number_of_ffts, transpose, transpose,
               transpose ? number_of_ffts : fft_size,
               transpose ? number_of_ffts : fft_size, input_array, input_array);

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

  fft_free_complex(input_array);

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
int fft_test_1d_local_r2c_low(const int fft_size, const int number_of_ffts,
                              const int transpose_rs, const int transpose_gs) {
  const int my_process = cp_mpi_comm_rank(cp_mpi_get_comm_world());

  int errors = 0;

  const double pi = acos(-1);

  double *input_array = NULL;
  double complex *output_array = NULL;
  fft_allocate_double(2 * (fft_size / 2 + 1) * number_of_ffts, &input_array);
  fft_allocate_complex((fft_size / 2 + 1) * number_of_ffts + 4, &output_array);

  fft_register_1d_r2c_local(fft_size, number_of_ffts, transpose_rs,
                            transpose_gs,
                            transpose_rs ? number_of_ffts : fft_size,
                            transpose_gs ? number_of_ffts : fft_size / 2 + 1,
                            input_array, output_array);
  fft_register_1d_c2r_local(
      fft_size, number_of_ffts, transpose_gs, transpose_rs,
      transpose_gs ? number_of_ffts : fft_size / 2 + 1,
      transpose_rs ? number_of_ffts : fft_size, output_array, input_array);

  memset(input_array, 0, fft_size * number_of_ffts * sizeof(double));
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

  fft_1d_r2c_local(fft_size, number_of_ffts, transpose_rs, transpose_gs,
                   transpose_rs ? number_of_ffts : fft_size,
                   transpose_gs ? number_of_ffts : fft_size / 2 + 1,
                   input_array, output_array);

  double max_error = 0.0;
  if (transpose_gs) {
#pragma omp parallel for default(none)                                         \
    shared(output_array, fft_size, number_of_ffts, pi, my_process)             \
    reduction(max : max_error) collapse(2)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index = 0; index < fft_size / 2 + 1; index++) {
        const double complex my_value =
            output_array[number_of_fft + index * number_of_ffts];
        const double complex ref_value =
            cexp(-2.0 * I * pi * (number_of_fft % fft_size) * index / fft_size);
        const double current_error = cabs(my_value - ref_value);
        if (my_process == 0 && current_error > 1e-12)
          printf("Error %i %i / %i %i: (%f %f) (%f %f)\n", index, number_of_fft,
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
      for (int index = 0; index < fft_size / 2 + 1; index++) {
        const double complex my_value =
            output_array[number_of_fft * (fft_size / 2 + 1) + index];
        const double complex ref_value =
            cexp(-2.0 * I * pi * (number_of_fft % fft_size) * index / fft_size);
        const double current_error = cabs(my_value - ref_value);
        if (my_process == 0 && current_error > 1e-12)
          printf("Error %i %i / %i %i: (%f %f) (%f %f)\n", index, number_of_fft,
                 fft_size, number_of_ffts, creal(my_value), cimag(my_value),
                 creal(ref_value), cimag(ref_value));
        max_error = fmax(max_error, current_error);
      }
    }
  }
  fflush(stdout);

  if (max_error > 1.0e-8) {
    if (my_process == 0) {
      printf("The fw R2C 1D-FFT does not work correctly (%i %i): %f!\n",
             fft_size, number_of_ffts, max_error);
      fflush(stdout);
    }
    errors++;
  }

  // Check the backward FFT
  memset(output_array, 0,
         (fft_size / 2 + 1) * number_of_ffts * sizeof(double complex));

  if (transpose_gs) {
#pragma omp parallel for default(none) collapse(2)                             \
    shared(output_array, fft_size, number_of_ffts)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index = 0; index < fft_size / 2 + 1; index++) {
        output_array[number_of_fft + index * number_of_ffts] =
            cexp(-2.0 * I * acos(-1) * index * number_of_fft / fft_size);
      }
    }
  } else {
#pragma omp parallel for default(none)                                         \
    shared(output_array, fft_size, number_of_ffts)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index = 0; index < fft_size / 2 + 1; index++) {
        output_array[number_of_fft * (fft_size / 2 + 1) + index] =
            cexp(-2.0 * I * acos(-1) * index * number_of_fft / fft_size);
      }
    }
  }

  fft_1d_c2r_local(fft_size, number_of_ffts, transpose_gs, transpose_rs,
                   transpose_gs ? number_of_ffts : fft_size / 2 + 1,
                   transpose_rs ? number_of_ffts : fft_size, output_array,
                   input_array);

  max_error = 0.0;
  if (transpose_rs) {
#pragma omp parallel for default(none)                                         \
    shared(input_array, fft_size, number_of_ffts, my_process)                  \
    reduction(max : max_error) collapse(2)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index = 0; index < fft_size; index++) {
        const double my_value =
            input_array[index * number_of_ffts + number_of_fft];
        const double ref_value =
            (number_of_fft % fft_size == index ? (double)fft_size : 0.0);
        const double current_error = fabs(my_value - ref_value);
        if (my_process == 0 && current_error > 1e-4)
          printf("ERROR %i %i / %i %i : %f %f\n", number_of_fft, index,
                 number_of_ffts, fft_size, my_value, ref_value);
        max_error = fmax(max_error, current_error);
      }
    }
  } else {
#pragma omp parallel for default(none)                                         \
    shared(input_array, fft_size, number_of_ffts, my_process)                  \
    reduction(max : max_error) collapse(2)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index = 0; index < fft_size; index++) {
        const double my_value = input_array[index + number_of_fft * fft_size];
        const double ref_value =
            (number_of_fft % fft_size == index ? (double)fft_size : 0.0);
        const double current_error = fabs(my_value - ref_value);
        if (my_process == 0 && current_error > 1e-4)
          printf("ERROR %i %i / %i %i : %f %f\n", number_of_fft, index,
                 number_of_ffts, fft_size, my_value, ref_value);
        max_error = fmax(max_error, current_error);
      }
    }
  }
  fflush(stdout);

  if (max_error > 1e-8) {
    if (my_process == 0) {
      printf("The bw C2R-1D FFT does not work correctly (%i %i): %f!\n",
             fft_size, number_of_ffts, max_error);
      fflush(stdout);
    }
    errors++;
  }

  fft_free_double(input_array);
  fft_free_complex(output_array);

  if (errors == 0 && my_process == 0)
    printf("The 1D R/C FFT does work correctly (%i %i)!\n", fft_size,
           number_of_ffts);
  return errors;
}

/*******************************************************************************
 * \brief Function to test the local FFT backend.
 * \author Frederick Stein
 ******************************************************************************/
int fft_test_1d_local_r2c_inplace_low(const int fft_size,
                                      const int number_of_ffts,
                                      const int transpose) {
  const int my_process = cp_mpi_comm_rank(cp_mpi_get_comm_world());

  int errors = 0;

  const double pi = acos(-1);

  double *input_array = NULL;
  fft_allocate_double(2 * (fft_size / 2 + 1) * number_of_ffts, &input_array);
  double complex *output_array = (double complex *)input_array;

  fft_register_1d_r2c_local(fft_size, number_of_ffts, transpose, transpose,
                            transpose ? number_of_ffts : 2 * (fft_size / 2 + 1),
                            transpose ? number_of_ffts : fft_size / 2 + 1,
                            input_array, output_array);
  fft_register_1d_c2r_local(fft_size, number_of_ffts, transpose, transpose,
                            transpose ? number_of_ffts : fft_size / 2 + 1,
                            transpose ? number_of_ffts : 2 * (fft_size / 2 + 1),
                            output_array, input_array);

  memset(input_array, 0,
         2 * (fft_size / 2 + 1) * number_of_ffts * sizeof(double));
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
      input_array[(number_of_fft % fft_size) +
                  number_of_fft * 2 * (fft_size / 2 + 1)] = 1.0;
    }
  }

  fft_1d_r2c_local(fft_size, number_of_ffts, transpose, transpose,
                   transpose ? number_of_ffts : 2 * (fft_size / 2 + 1),
                   transpose ? number_of_ffts : fft_size / 2 + 1, input_array,
                   output_array);

  double max_error = 0.0;
  if (transpose) {
#pragma omp parallel for default(none)                                         \
    shared(output_array, fft_size, number_of_ffts, pi, my_process)             \
    reduction(max : max_error) collapse(2)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index = 0; index < fft_size / 2 + 1; index++) {
        const double complex my_value =
            output_array[number_of_fft + index * number_of_ffts];
        const double complex ref_value =
            cexp(-2.0 * I * pi * (number_of_fft % fft_size) * index / fft_size);
        const double current_error = cabs(my_value - ref_value);
        if (my_process == 0 && current_error > 1e-12)
          printf("Error %i %i / %i %i: (%f %f) (%f %f)\n", index, number_of_fft,
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
      for (int index = 0; index < fft_size / 2 + 1; index++) {
        const double complex my_value =
            output_array[number_of_fft * (fft_size / 2 + 1) + index];
        const double complex ref_value =
            cexp(-2.0 * I * pi * (number_of_fft % fft_size) * index / fft_size);
        const double current_error = cabs(my_value - ref_value);
        if (my_process == 0 && current_error > 1e-12)
          printf("Error %i %i / %i %i: (%f %f) (%f %f)\n", index, number_of_fft,
                 fft_size, number_of_ffts, creal(my_value), cimag(my_value),
                 creal(ref_value), cimag(ref_value));
        max_error = fmax(max_error, current_error);
      }
    }
  }
  fflush(stdout);

  if (max_error > 1.0e-8) {
    if (my_process == 0) {
      printf(
          "The fw R2C 1D-FFT (inplace) does not work correctly (%i %i): %f!\n",
          fft_size, number_of_ffts, max_error);
      fflush(stdout);
    }
    errors++;
  }

  // Check the backward FFT
  memset(output_array, 0,
         (fft_size / 2 + 1) * number_of_ffts * sizeof(double complex));

  if (transpose) {
#pragma omp parallel for default(none) collapse(2)                             \
    shared(output_array, fft_size, number_of_ffts)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index = 0; index < fft_size / 2 + 1; index++) {
        output_array[number_of_fft + index * number_of_ffts] =
            cexp(-2.0 * I * acos(-1) * index * number_of_fft / fft_size);
      }
    }
  } else {
#pragma omp parallel for default(none)                                         \
    shared(output_array, fft_size, number_of_ffts)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index = 0; index < fft_size / 2 + 1; index++) {
        output_array[number_of_fft * (fft_size / 2 + 1) + index] =
            cexp(-2.0 * I * acos(-1) * index * number_of_fft / fft_size);
      }
    }
  }

  fft_1d_c2r_local(fft_size, number_of_ffts, transpose, transpose,
                   transpose ? number_of_ffts : fft_size / 2 + 1,
                   transpose ? number_of_ffts : 2 * (fft_size / 2 + 1),
                   output_array, input_array);

  max_error = 0.0;
  if (transpose) {
#pragma omp parallel for default(none)                                         \
    shared(input_array, fft_size, number_of_ffts, my_process)                  \
    reduction(max : max_error) collapse(2)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index = 0; index < fft_size; index++) {
        const double my_value =
            input_array[index * number_of_ffts + number_of_fft];
        const double ref_value =
            (number_of_fft % fft_size == index ? (double)fft_size : 0.0);
        const double current_error = fabs(my_value - ref_value);
        if (my_process == 0 && current_error > 1e-4)
          printf("ERROR %i %i / %i %i : %f %f\n", number_of_fft, index,
                 number_of_ffts, fft_size, my_value, ref_value);
        max_error = fmax(max_error, current_error);
      }
    }
  } else {
#pragma omp parallel for default(none)                                         \
    shared(input_array, fft_size, number_of_ffts, my_process)                  \
    reduction(max : max_error) collapse(2)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index = 0; index < fft_size; index++) {
        const double my_value =
            input_array[index + number_of_fft * 2 * (fft_size / 2 + 1)];
        const double ref_value =
            (number_of_fft % fft_size == index ? (double)fft_size : 0.0);
        const double current_error = fabs(my_value - ref_value);
        if (my_process == 0 && current_error > 1e-4)
          printf("ERROR %i %i / %i %i : %f %f\n", number_of_fft, index,
                 number_of_ffts, fft_size, my_value, ref_value);
        max_error = fmax(max_error, current_error);
      }
    }
  }
  fflush(stdout);

  if (max_error > 1e-8) {
    if (my_process == 0) {
      printf(
          "The bw C2R-1D FFT (inplace) does not work correctly (%i %i): %f!\n",
          fft_size, number_of_ffts, max_error);
      fflush(stdout);
    }
    errors++;
  }

  fft_free_double(input_array);

  if (errors == 0 && my_process == 0)
    printf("The 1D R2C FFT (inplace) does work correctly (%i %i)!\n", fft_size,
           number_of_ffts);
  return errors;
}

/*******************************************************************************
 * \brief Function to test the local FFT backend.
 * \author Frederick Stein
 ******************************************************************************/
int fft_test_2d_local_low(const int fft_size[2], const int number_of_ffts,
                          const int transpose_rs, const int transpose_gs) {
  const int my_process = cp_mpi_comm_rank(cp_mpi_get_comm_world());

  int errors = 0;

  const double pi = acos(-1);

  double complex *input_array = NULL, *output_array = NULL;
  const int elements_per_fft = fft_size[0] * fft_size[1];
  fft_allocate_complex(elements_per_fft * number_of_ffts, &input_array);
  fft_allocate_complex(elements_per_fft * number_of_ffts, &output_array);

  fft_register_2d_local(true, fft_size, number_of_ffts, transpose_rs,
                        transpose_gs, input_array, output_array);
  fft_register_2d_local(false, fft_size, number_of_ffts, transpose_gs,
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

  fft_2d_local(true, fft_size, number_of_ffts, transpose_rs, transpose_gs,
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

  fft_2d_local(false, fft_size, number_of_ffts, transpose_gs, transpose_rs,
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

  fft_free_complex(input_array);
  fft_free_complex(output_array);

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
int fft_test_2d_local_inplace_low(const int fft_size[2],
                                  const int number_of_ffts,
                                  const int transpose) {
  const int my_process = cp_mpi_comm_rank(cp_mpi_get_comm_world());

  int errors = 0;

  const double pi = acos(-1);

  double complex *input_array = NULL;
  const int elements_per_fft = fft_size[0] * fft_size[1];
  fft_allocate_complex(elements_per_fft * number_of_ffts, &input_array);
  double complex *output_array = input_array;

  fft_register_2d_local(true, fft_size, number_of_ffts, transpose, transpose,
                        input_array, output_array);
  fft_register_2d_local(false, fft_size, number_of_ffts, transpose, transpose,
                        output_array, input_array);

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

  fft_2d_local(true, fft_size, number_of_ffts, transpose, transpose,
               input_array, output_array);

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

  fft_2d_local(false, fft_size, number_of_ffts, transpose, transpose,
               output_array, input_array);

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

  fft_free_complex(input_array);

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
int fft_test_2d_local_r2c_low(const int fft_size[2], const int number_of_ffts,
                              const int transpose_rs, const int transpose_gs) {
  const int my_process = cp_mpi_comm_rank(cp_mpi_get_comm_world());

  int errors = 0;

  const double pi = acos(-1);

  double *real_buffer = NULL;
  double complex *complex_buffer = NULL;
  fft_allocate_double(2 * (fft_size[1] / 2 + 1) * fft_size[0] * number_of_ffts,
                      &real_buffer);
  fft_allocate_complex((fft_size[1] / 2 + 1) * fft_size[0] * number_of_ffts,
                       &complex_buffer);

  fft_register_2d_r2c_local(fft_size, number_of_ffts, transpose_rs,
                            transpose_gs, real_buffer, complex_buffer);
  fft_register_2d_c2r_local(fft_size, number_of_ffts, transpose_gs,
                            transpose_rs, complex_buffer, real_buffer);

  memset(real_buffer, 0,
         2 * (fft_size[1] / 2 + 1) * fft_size[0] * number_of_ffts *
             sizeof(double));

  double max_error = 0.0;
  // Check the forward FFT
  if (transpose_rs) {
#pragma omp parallel for default(none)                                         \
    shared(real_buffer, fft_size, number_of_ffts)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      real_buffer[(number_of_fft / fft_size[1] % fft_size[0] * fft_size[1] +
                   number_of_fft % fft_size[1]) *
                      number_of_ffts +
                  number_of_fft] = 1.0;
    }
  } else {
#pragma omp parallel for default(none)                                         \
    shared(real_buffer, fft_size, number_of_ffts)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      real_buffer[(number_of_fft / fft_size[1] % fft_size[0] * fft_size[1] +
                   number_of_fft % fft_size[1]) +
                  number_of_fft * (fft_size[0] * fft_size[1])] = 1.0;
    }
  }

  fft_2d_r2c_local(fft_size, number_of_ffts, transpose_rs, transpose_gs,
                   real_buffer, complex_buffer);

  if (transpose_gs) {
#pragma omp parallel for default(none)                                         \
    shared(complex_buffer, fft_size, number_of_ffts, pi, my_process)           \
    reduction(max : max_error) collapse(3)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index_0 = 0; index_0 < fft_size[0]; index_0++) {
        for (int index_1 = 0; index_1 < fft_size[1] / 2 + 1; index_1++) {
          const double complex my_value =
              complex_buffer[number_of_fft +
                             (index_0 * (fft_size[1] / 2 + 1) + index_1) *
                                 number_of_ffts];
          const double complex ref_value = cexp(
              -2.0 * I * pi *
              ((double)(number_of_fft / fft_size[1]) * index_0 / fft_size[0] +
               (double)(number_of_fft % fft_size[1]) * index_1 / fft_size[1]));
          double current_error = cabs(my_value - ref_value);
          if (my_process == 0 && current_error > 1e-6)
            printf("Error %i %i %i/%i %i %i: (%f %f) (%f %f)\n", index_0,
                   index_1, number_of_fft, fft_size[0], fft_size[1],
                   number_of_ffts, creal(my_value), cimag(my_value),
                   creal(ref_value), cimag(ref_value));
          max_error = fmax(max_error, current_error);
        }
      }
    }
  } else {
#pragma omp parallel for default(none)                                         \
    shared(complex_buffer, fft_size, number_of_ffts, pi, my_process)           \
    reduction(max : max_error) collapse(3)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index_0 = 0; index_0 < fft_size[0]; index_0++) {
        for (int index_1 = 0; index_1 < fft_size[1] / 2 + 1; index_1++) {
          const double complex my_value =
              complex_buffer[number_of_fft * fft_size[0] *
                                 (fft_size[1] / 2 + 1) +
                             index_0 * (fft_size[1] / 2 + 1) + index_1];
          const double complex ref_value = cexp(
              -2.0 * I * pi *
              ((double)(number_of_fft / fft_size[1]) * index_0 / fft_size[0] +
               (double)(number_of_fft % fft_size[1]) * index_1 / fft_size[1]));
          double current_error = cabs(my_value - ref_value);
          if (my_process == 0 && current_error > 1e-6)
            printf("Error %i %i %i/%i %i %i: (%f %f) (%f %f)\n", index_0,
                   index_1, number_of_fft, fft_size[0], fft_size[1],
                   number_of_ffts, creal(my_value), cimag(my_value),
                   creal(ref_value), cimag(ref_value));
          max_error = fmax(max_error, current_error);
        }
      }
    }
  }
  fflush(stdout);

  if (max_error > 1.0e-12) {
    if (my_process == 0) {
      printf("The fw R2C 2D-FFT does not work correctly (%i %i/%i): %f!\n",
             fft_size[0], fft_size[1], number_of_ffts, max_error);
      fflush(stdout);
    }
    errors++;
  }

  // Check the backward FFT
  memset(complex_buffer, 0,
         fft_size[0] * (fft_size[1] / 2 + 1) * number_of_ffts *
             sizeof(double complex));

  if (transpose_gs) {
#pragma omp parallel for default(none)                                         \
    shared(complex_buffer, fft_size, number_of_ffts, pi)                       \
    reduction(max : max_error) collapse(3)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index_0 = 0; index_0 < fft_size[0]; index_0++) {
        for (int index_1 = 0; index_1 < fft_size[1] / 2 + 1; index_1++) {
          complex_buffer[number_of_fft +
                         (index_0 * (fft_size[1] / 2 + 1) + index_1) *
                             number_of_ffts] =
              cexp(-2.0 * I * pi *
                   ((double)(number_of_fft / fft_size[1]) * index_0 /
                        fft_size[0] +
                    (double)(number_of_fft % fft_size[1]) * index_1 /
                        fft_size[1]));
        }
      }
    }
  } else {
#pragma omp parallel for default(none)                                         \
    shared(complex_buffer, fft_size, number_of_ffts, pi)                       \
    reduction(max : max_error) collapse(3)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index_0 = 0; index_0 < fft_size[0]; index_0++) {
        for (int index_1 = 0; index_1 < fft_size[1] / 2 + 1; index_1++) {
          complex_buffer[number_of_fft * fft_size[0] * (fft_size[1] / 2 + 1) +
                         index_0 * (fft_size[1] / 2 + 1) + index_1] =
              cexp(-2.0 * I * pi *
                   ((double)(number_of_fft / fft_size[1]) * index_0 /
                        fft_size[0] +
                    (double)(number_of_fft % fft_size[1]) * index_1 /
                        fft_size[1]));
        }
      }
    }
  }

  fft_2d_c2r_local(fft_size, number_of_ffts, transpose_gs, transpose_rs,
                   complex_buffer, real_buffer);

  max_error = 0.0;
  if (transpose_rs) {
#pragma omp parallel for default(none)                                         \
    shared(real_buffer, fft_size, number_of_ffts, pi, my_process)              \
    reduction(max : max_error) collapse(3)
    for (int index_0 = 0; index_0 < fft_size[0]; index_0++) {
      for (int index_1 = 0; index_1 < fft_size[1]; index_1++) {
        for (int number_of_fft = 0; number_of_fft < number_of_ffts;
             number_of_fft++) {
          const double my_value =
              real_buffer[(index_0 * fft_size[1] + index_1) * number_of_ffts +
                          number_of_fft];
          const double ref_value =
              index_0 == number_of_fft / fft_size[1] % fft_size[0] &&
                      index_1 == number_of_fft % fft_size[1]
                  ? (double)(fft_size[0] * fft_size[1])
                  : 0.0;
          const double current_error = fabs(my_value - ref_value);
          if (my_process == 0 && current_error > 1e-6)
            printf("Error %i %i %i/%i %i %i: %f %f\n", index_0, index_1,
                   number_of_fft, fft_size[0], fft_size[1], number_of_ffts,
                   my_value, ref_value);
          max_error = fmax(max_error, current_error);
        }
      }
    }
  } else {
#pragma omp parallel for default(none)                                         \
    shared(real_buffer, fft_size, number_of_ffts, pi, my_process)              \
    reduction(max : max_error) collapse(3)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index_0 = 0; index_0 < fft_size[0]; index_0++) {
        for (int index_1 = 0; index_1 < fft_size[1]; index_1++) {
          const double my_value =
              real_buffer[index_0 * fft_size[1] + index_1 +
                          number_of_fft * fft_size[0] * fft_size[1]];
          const double ref_value =
              index_0 == number_of_fft / fft_size[1] % fft_size[0] &&
                      index_1 == number_of_fft % fft_size[1]
                  ? (double)(fft_size[0] * fft_size[1])
                  : 0.0;
          double current_error = fabs(my_value - ref_value);
          if (my_process == 0 && current_error > 1e-6)
            printf("Error %i %i %i/%i %i %i: %f %f\n", index_0, index_1,
                   number_of_fft, fft_size[0], fft_size[1], number_of_ffts,
                   my_value, ref_value);
          max_error = fmax(max_error, current_error);
        }
      }
    }
  }
  fflush(stdout);

  fft_free_double(real_buffer);
  fft_free_complex(complex_buffer);

  if (max_error > 1e-8) {
    if (my_process == 0) {
      printf("The bw C2R 2D-FFT does not work correctly (%i %i/%i): %f!\n",
             fft_size[0], fft_size[1], number_of_ffts, max_error);
      fflush(stdout);
    }
    errors++;
  }

  if (errors == 0 && my_process == 0)
    printf("The 2D R2C/C2R FFT does work correctly (%i %i/%i)!\n", fft_size[0],
           fft_size[1], number_of_ffts);
  return errors;
}

/*******************************************************************************
 * \brief Function to test the local FFT backend.
 * \author Frederick Stein
 ******************************************************************************/
int fft_test_2d_local_r2c_inplace_low(const int fft_size[2],
                                      const int number_of_ffts,
                                      const int transpose) {
  const int my_process = cp_mpi_comm_rank(cp_mpi_get_comm_world());

  int errors = 0;

  const double pi = acos(-1);

  double *real_buffer = NULL;
  fft_allocate_double(2 * (fft_size[1] / 2 + 1) * fft_size[0] * number_of_ffts,
                      &real_buffer);
  double complex *complex_buffer = (double complex *)real_buffer;

  fft_register_2d_r2c_local(fft_size, number_of_ffts, transpose, transpose,
                            real_buffer, complex_buffer);
  fft_register_2d_c2r_local(fft_size, number_of_ffts, transpose, transpose,
                            complex_buffer, real_buffer);

  memset(real_buffer, 0,
         2 * (fft_size[1] / 2 + 1) * fft_size[0] * number_of_ffts *
             sizeof(double));

  double max_error = 0.0;
  // Check the forward FFT
  if (transpose) {
#pragma omp parallel for default(none)                                         \
    shared(real_buffer, fft_size, number_of_ffts)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      real_buffer[(number_of_fft / fft_size[1] % fft_size[0] * 2 *
                       (fft_size[1] / 2 + 1) +
                   number_of_fft % fft_size[1]) *
                      number_of_ffts +
                  number_of_fft] = 1.0;
    }
  } else {
#pragma omp parallel for default(none)                                         \
    shared(real_buffer, fft_size, number_of_ffts)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      real_buffer[(number_of_fft / fft_size[1] % fft_size[0] * 2 *
                       (fft_size[1] / 2 + 1) +
                   number_of_fft % fft_size[1]) +
                  number_of_fft * (fft_size[0] * 2 * (fft_size[1] / 2 + 1))] =
          1.0;
    }
  }

  fft_2d_r2c_local(fft_size, number_of_ffts, transpose, transpose, real_buffer,
                   complex_buffer);

  if (transpose) {
#pragma omp parallel for default(none)                                         \
    shared(complex_buffer, fft_size, number_of_ffts, pi, my_process)           \
    reduction(max : max_error) collapse(3)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index_0 = 0; index_0 < fft_size[0]; index_0++) {
        for (int index_1 = 0; index_1 < fft_size[1] / 2 + 1; index_1++) {
          const double complex my_value =
              complex_buffer[number_of_fft +
                             (index_0 * (fft_size[1] / 2 + 1) + index_1) *
                                 number_of_ffts];
          const double complex ref_value = cexp(
              -2.0 * I * pi *
              ((double)(number_of_fft / fft_size[1]) * index_0 / fft_size[0] +
               (double)(number_of_fft % fft_size[1]) * index_1 / fft_size[1]));
          double current_error = cabs(my_value - ref_value);
          if (my_process == 0 && current_error > 1e-6)
            printf("Error %i %i %i/%i %i %i: (%f %f) (%f %f)\n", index_0,
                   index_1, number_of_fft, fft_size[0], fft_size[1],
                   number_of_ffts, creal(my_value), cimag(my_value),
                   creal(ref_value), cimag(ref_value));
          max_error = fmax(max_error, current_error);
        }
      }
    }
  } else {
#pragma omp parallel for default(none)                                         \
    shared(complex_buffer, fft_size, number_of_ffts, pi, my_process)           \
    reduction(max : max_error) collapse(3)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index_0 = 0; index_0 < fft_size[0]; index_0++) {
        for (int index_1 = 0; index_1 < fft_size[1] / 2 + 1; index_1++) {
          const double complex my_value =
              complex_buffer[number_of_fft * fft_size[0] *
                                 (fft_size[1] / 2 + 1) +
                             index_0 * (fft_size[1] / 2 + 1) + index_1];
          const double complex ref_value = cexp(
              -2.0 * I * pi *
              ((double)(number_of_fft / fft_size[1]) * index_0 / fft_size[0] +
               (double)(number_of_fft % fft_size[1]) * index_1 / fft_size[1]));
          double current_error = cabs(my_value - ref_value);
          if (my_process == 0 && current_error > 1e-6)
            printf("Error %i %i %i/%i %i %i: (%f %f) (%f %f)\n", index_0,
                   index_1, number_of_fft, fft_size[0], fft_size[1],
                   number_of_ffts, creal(my_value), cimag(my_value),
                   creal(ref_value), cimag(ref_value));
          max_error = fmax(max_error, current_error);
        }
      }
    }
  }
  fflush(stdout);

  if (max_error > 1.0e-12) {
    if (my_process == 0) {
      printf("The fw R2C 2D-FFT (inplace) does not work correctly (%i %i/%i): "
             "%f!\n",
             fft_size[0], fft_size[1], number_of_ffts, max_error);
      fflush(stdout);
    }
    errors++;
  }

  // Check the backward FFT
  memset(complex_buffer, 0,
         fft_size[0] * (fft_size[1] / 2 + 1) * number_of_ffts *
             sizeof(double complex));

  if (transpose) {
#pragma omp parallel for default(none)                                         \
    shared(complex_buffer, fft_size, number_of_ffts, pi)                       \
    reduction(max : max_error) collapse(3)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index_0 = 0; index_0 < fft_size[0]; index_0++) {
        for (int index_1 = 0; index_1 < fft_size[1] / 2 + 1; index_1++) {
          complex_buffer[number_of_fft +
                         (index_0 * (fft_size[1] / 2 + 1) + index_1) *
                             number_of_ffts] =
              cexp(-2.0 * I * pi *
                   ((double)(number_of_fft / fft_size[1]) * index_0 /
                        fft_size[0] +
                    (double)(number_of_fft % fft_size[1]) * index_1 /
                        fft_size[1]));
        }
      }
    }
  } else {
#pragma omp parallel for default(none)                                         \
    shared(complex_buffer, fft_size, number_of_ffts, pi)                       \
    reduction(max : max_error) collapse(3)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index_0 = 0; index_0 < fft_size[0]; index_0++) {
        for (int index_1 = 0; index_1 < fft_size[1] / 2 + 1; index_1++) {
          complex_buffer[number_of_fft * fft_size[0] * (fft_size[1] / 2 + 1) +
                         index_0 * (fft_size[1] / 2 + 1) + index_1] =
              cexp(-2.0 * I * pi *
                   ((double)(number_of_fft / fft_size[1]) * index_0 /
                        fft_size[0] +
                    (double)(number_of_fft % fft_size[1]) * index_1 /
                        fft_size[1]));
        }
      }
    }
  }

  fft_2d_c2r_local(fft_size, number_of_ffts, transpose, transpose,
                   complex_buffer, real_buffer);

  max_error = 0.0;
  if (transpose) {
#pragma omp parallel for default(none)                                         \
    shared(real_buffer, fft_size, number_of_ffts, pi, my_process)              \
    reduction(max : max_error) collapse(3)
    for (int index_0 = 0; index_0 < fft_size[0]; index_0++) {
      for (int index_1 = 0; index_1 < fft_size[1]; index_1++) {
        for (int number_of_fft = 0; number_of_fft < number_of_ffts;
             number_of_fft++) {
          const double my_value =
              real_buffer[(index_0 * 2 * (fft_size[1] / 2 + 1) + index_1) *
                              number_of_ffts +
                          number_of_fft];
          const double ref_value =
              index_0 == number_of_fft / fft_size[1] % fft_size[0] &&
                      index_1 == number_of_fft % fft_size[1]
                  ? (double)(fft_size[0] * fft_size[1])
                  : 0.0;
          const double current_error = fabs(my_value - ref_value);
          if (my_process == 0 && current_error > 1e-6)
            printf("Error %i %i %i/%i %i %i: %f %f\n", index_0, index_1,
                   number_of_fft, fft_size[0], fft_size[1], number_of_ffts,
                   my_value, ref_value);
          max_error = fmax(max_error, current_error);
        }
      }
    }
  } else {
#pragma omp parallel for default(none)                                         \
    shared(real_buffer, fft_size, number_of_ffts, pi, my_process)              \
    reduction(max : max_error) collapse(3)
    for (int number_of_fft = 0; number_of_fft < number_of_ffts;
         number_of_fft++) {
      for (int index_0 = 0; index_0 < fft_size[0]; index_0++) {
        for (int index_1 = 0; index_1 < fft_size[1]; index_1++) {
          const double my_value =
              real_buffer[index_0 * 2 * (fft_size[1] / 2 + 1) + index_1 +
                          number_of_fft * fft_size[0] * 2 *
                              (fft_size[1] / 2 + 1)];
          const double ref_value =
              index_0 == number_of_fft / fft_size[1] % fft_size[0] &&
                      index_1 == number_of_fft % fft_size[1]
                  ? (double)(fft_size[0] * fft_size[1])
                  : 0.0;
          double current_error = fabs(my_value - ref_value);
          if (my_process == 0 && current_error > 1e-6)
            printf("Error %i %i %i/%i %i %i: %f %f\n", index_0, index_1,
                   number_of_fft, fft_size[0], fft_size[1], number_of_ffts,
                   my_value, ref_value);
          max_error = fmax(max_error, current_error);
        }
      }
    }
  }
  fflush(stdout);

  fft_free_double(real_buffer);

  if (max_error > 1e-8) {
    if (my_process == 0) {
      printf("The bw C2R 2D-FFT (inplace) does not work correctly (%i %i/%i): "
             "%f!\n",
             fft_size[0], fft_size[1], number_of_ffts, max_error);
      fflush(stdout);
    }
    errors++;
  }

  if (errors == 0 && my_process == 0)
    printf("The 2D R2C/C2R FFT (inplace) does work correctly (%i %i/%i)!\n",
           fft_size[0], fft_size[1], number_of_ffts);
  return errors;
}

/*******************************************************************************
 * \brief Function to test the local FFT backend.
 * \author Frederick Stein
 ******************************************************************************/
int fft_test_3d_local_low(const int fft_size[3], const int test_every) {
  const int my_process = cp_mpi_comm_rank(cp_mpi_get_comm_world());

  int errors = 0;

  const double pi = acos(-1);

  double complex *input_array = NULL, *output_array = NULL;
  fft_allocate_complex(fft_size[0] * fft_size[1] * fft_size[2], &input_array);
  fft_allocate_complex(fft_size[0] * fft_size[1] * fft_size[2], &output_array);

  fft_register_3d_local(true, fft_size, input_array, output_array);
  fft_register_3d_local(false, fft_size, output_array, input_array);

  double max_error = 0.0;
  int number_of_tests = 0;
  for (int mx = 0; mx < fft_size[0]; mx++) {
    for (int my = 0; my < fft_size[1]; my++) {
      for (int mz = 0; mz < fft_size[2]; mz++) {
        if (test_every > 0 && number_of_tests % test_every != 0) {
          number_of_tests++;
          continue;
        }
        number_of_tests++;
        memset(input_array, 0,
               fft_size[0] * fft_size[1] * fft_size[2] *
                   sizeof(double complex));
        input_array[(mx * fft_size[1] + my) * fft_size[2] + mz] = 1.0;
        fft_3d_local(true, fft_size, input_array, output_array);

#pragma omp parallel for default(none)                                         \
    shared(output_array, fft_size, pi, mx, my, mz, my_process)                 \
    reduction(max : max_error) collapse(3)
        for (int nx = 0; nx < fft_size[0]; nx++) {
          for (int ny = 0; ny < fft_size[1]; ny++) {
            for (int nz = 0; nz < fft_size[2]; nz++) {
              const double complex my_value =
                  output_array[(nx * fft_size[1] + ny) * fft_size[2] + nz];
              const double complex ref_value =
                  cexp(-2.0 * I * pi *
                       (((double)mx) * nx / fft_size[0] +
                        ((double)my) * ny / fft_size[1] +
                        ((double)mz) * nz / fft_size[2]));
              double current_error = cabs(my_value - ref_value);
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
  number_of_tests = 0;
  for (int mx = 0; mx < fft_size[0]; mx++) {
    for (int my = 0; my < fft_size[1]; my++) {
      for (int mz = 0; mz < fft_size[2]; mz++) {
        if (test_every > 0 && number_of_tests % test_every != 0) {
          number_of_tests++;
          continue;
        }
        number_of_tests++;
        memset(output_array, 0, product3(fft_size) * sizeof(double complex));
        output_array[(mx * fft_size[1] + my) * fft_size[2] + mz] = 1.0;
        fft_3d_local(false, fft_size, output_array, input_array);

#pragma omp parallel for default(none)                                         \
    shared(input_array, fft_size, pi, mx, my, mz, my_process)                  \
    reduction(max : max_error) collapse(3)
        for (int nx = 0; nx < fft_size[0]; nx++) {
          for (int ny = 0; ny < fft_size[1]; ny++) {
            for (int nz = 0; nz < fft_size[2]; nz++) {
              const double complex my_value =
                  input_array[(nx * fft_size[1] + ny) * fft_size[2] + nz];
              const double complex ref_value =
                  cexp(2.0 * I * pi *
                       (((double)mx) * nx / fft_size[0] +
                        ((double)my) * ny / fft_size[1] +
                        ((double)mz) * nz / fft_size[2]));
              double current_error = cabs(my_value - ref_value);
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

  fft_free_complex(input_array);
  fft_free_complex(output_array);

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
int fft_test_3d_local_inplace_low(const int fft_size[3], const int test_every) {
  const int my_process = cp_mpi_comm_rank(cp_mpi_get_comm_world());

  int errors = 0;

  const double pi = acos(-1);

  double complex *input_array = NULL;
  fft_allocate_complex(fft_size[0] * fft_size[1] * fft_size[2], &input_array);
  double complex *output_array = input_array;

  fft_register_3d_local(true, fft_size, input_array, output_array);
  fft_register_3d_local(false, fft_size, output_array, input_array);

  double max_error = 0.0;
  int number_of_tests = 0;
  for (int mx = 0; mx < fft_size[0]; mx++) {
    for (int my = 0; my < fft_size[1]; my++) {
      for (int mz = 0; mz < fft_size[2]; mz++) {
        if (test_every > 0 && number_of_tests % test_every != 0) {
          number_of_tests++;
          continue;
        }
        number_of_tests++;
        memset(input_array, 0,
               fft_size[0] * fft_size[1] * fft_size[2] *
                   sizeof(double complex));
        input_array[(mx * fft_size[1] + my) * fft_size[2] + mz] = 1.0;
        fft_3d_local(true, fft_size, input_array, output_array);

#pragma omp parallel for default(none)                                         \
    shared(output_array, fft_size, pi, mx, my, mz) reduction(max : max_error)  \
    collapse(3)
        for (int nx = 0; nx < fft_size[0]; nx++) {
          for (int ny = 0; ny < fft_size[1]; ny++) {
            for (int nz = 0; nz < fft_size[2]; nz++) {
              const double complex my_value =
                  output_array[(nx * fft_size[1] + ny) * fft_size[2] + nz];
              const double complex ref_value =
                  cexp(-2.0 * I * pi *
                       (((double)mx) * nx / fft_size[0] +
                        ((double)my) * ny / fft_size[1] +
                        ((double)mz) * nz / fft_size[2]));
              double current_error = cabs(my_value - ref_value);
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
      printf(
          "The fw 3D-FFT (inplace) does not work correctly (%i %i %i): %f!\n",
          fft_size[0], fft_size[1], fft_size[2], max_error);
      fflush(stdout);
    }
    errors++;
  }

  max_error = 0.0;
  number_of_tests = 0;
  for (int mx = 0; mx < fft_size[0]; mx++) {
    for (int my = 0; my < fft_size[1]; my++) {
      for (int mz = 0; mz < fft_size[2]; mz++) {
        if (test_every > 0 && number_of_tests % test_every != 0) {
          number_of_tests++;
          continue;
        }
        number_of_tests++;
        memset(output_array, 0, product3(fft_size) * sizeof(double complex));
        output_array[(mx * fft_size[1] + my) * fft_size[2] + mz] = 1.0;
        fft_3d_local(false, fft_size, output_array, input_array);

#pragma omp parallel for default(none)                                         \
    shared(input_array, fft_size, pi, mx, my, mz) reduction(max : max_error)   \
    collapse(3)
        for (int nx = 0; nx < fft_size[0]; nx++) {
          for (int ny = 0; ny < fft_size[1]; ny++) {
            for (int nz = 0; nz < fft_size[2]; nz++) {
              const double complex my_value =
                  input_array[(nx * fft_size[1] + ny) * fft_size[2] + nz];
              const double complex ref_value =
                  cexp(2.0 * I * pi *
                       (((double)mx) * nx / fft_size[0] +
                        ((double)my) * ny / fft_size[1] +
                        ((double)mz) * nz / fft_size[2]));
              double current_error = cabs(my_value - ref_value);
              max_error = fmax(max_error, current_error);
            }
          }
        }
      }
    }
  }
  fflush(stdout);

  fft_free_complex(input_array);

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
 * \brief Function to test the local FFT backend.
 * \author Frederick Stein
 ******************************************************************************/
int fft_test_3d_local_r2c_low(const int fft_size[3], const int test_every) {
  const int my_process = cp_mpi_comm_rank(cp_mpi_get_comm_world());

  int errors = 0;

  const double pi = acos(-1);

  double *double_buffer = NULL;
  double complex *complex_buffer = NULL;
  fft_allocate_double(2 * fft_size[0] * fft_size[1] * (fft_size[2] / 2 + 1),
                      &double_buffer);
  fft_allocate_complex(fft_size[0] * fft_size[1] * (fft_size[2] / 2 + 1),
                       &complex_buffer);

  fft_register_3d_r2c_local(fft_size, double_buffer, complex_buffer);
  fft_register_3d_c2r_local(fft_size, complex_buffer, double_buffer);

  double max_error = 0.0;
  int number_of_tests = 0;
  for (int mx = 0; mx < fft_size[0]; mx++) {
    for (int my = 0; my < fft_size[1]; my++) {
      for (int mz = 0; mz < fft_size[2]; mz++) {
        if (test_every > 0 && number_of_tests % test_every != 0) {
          number_of_tests++;
          continue;
        }
        number_of_tests++;
        memset(double_buffer, 0, product3(fft_size) * sizeof(double));
        double_buffer[(mx * fft_size[1] + my) * fft_size[2] + mz] = 1.0;
        fft_3d_r2c_local(fft_size, double_buffer, complex_buffer);

#pragma omp parallel for default(none)                                         \
    shared(complex_buffer, fft_size, pi, mx, my, mz, my_process)               \
    reduction(max : max_error) collapse(3)
        for (int nx = 0; nx < fft_size[0]; nx++) {
          for (int ny = 0; ny < fft_size[1]; ny++) {
            for (int nz = 0; nz < fft_size[2] / 2 + 1; nz++) {
              const double complex my_value =
                  complex_buffer[(nx * fft_size[1] + ny) *
                                     (fft_size[2] / 2 + 1) +
                                 nz];
              const double complex ref_value =
                  cexp(-2.0 * I * pi *
                       (((double)mx) * nx / fft_size[0] +
                        ((double)my) * ny / fft_size[1] +
                        ((double)mz) * nz / fft_size[2]));
              double current_error = cabs(my_value - ref_value);
              if (my_process == 0 && current_error > 1e-6) {
                printf("ERROR %i %i %i/%i %i %i: (%f %f) (%f %f)\n", nx, ny, nz,
                       mx, my, mz, creal(my_value), cimag(my_value),
                       creal(ref_value), cimag(ref_value));
              }
              max_error = fmax(max_error, current_error);
            }
          }
        }
      }
    }
  }
  fflush(stdout);

  if (max_error > 1.0e-6) {
    if (my_process == 0) {
      printf("The fw R2C 3D-FFT does not work correctly (%i %i %i): %f!\n",
             fft_size[0], fft_size[1], fft_size[2], max_error);
      fflush(stdout);
    }
    errors++;
  }

  max_error = 0.0;
  number_of_tests = 0;
  for (int mx = 0; mx < fft_size[0]; mx++) {
    for (int my = 0; my < fft_size[1]; my++) {
      for (int mz = 0; mz < fft_size[2]; mz++) {
        if (test_every > 0 && number_of_tests % test_every != 0) {
          number_of_tests++;
          continue;
        }
        number_of_tests++;

#pragma omp parallel for default(none)                                         \
    shared(complex_buffer, fft_size, pi, mx, my, mz, my_process) collapse(3)
        for (int nz = 0; nz < fft_size[2] / 2 + 1; nz++) {
          for (int ny = 0; ny < fft_size[1]; ny++) {
            for (int nx = 0; nx < fft_size[0]; nx++) {
              complex_buffer[(nx * fft_size[1] + ny) * (fft_size[2] / 2 + 1) +
                             nz] = cexp(-2.0 * I * pi *
                                        (((double)mx) * nx / fft_size[0] +
                                         ((double)my) * ny / fft_size[1] +
                                         ((double)mz) * nz / fft_size[2]));
            }
          }
        }
        fft_3d_c2r_local(fft_size, complex_buffer, double_buffer);

#pragma omp parallel for default(none)                                         \
    shared(double_buffer, fft_size, pi, mx, my, mz, my_process)                \
    reduction(max : max_error) collapse(3)
        for (int nx = 0; nx < fft_size[0]; nx++) {
          for (int ny = 0; ny < fft_size[1]; ny++) {
            for (int nz = 0; nz < fft_size[2]; nz++) {
              const double my_value =
                  double_buffer[(nx * fft_size[1] + ny) * fft_size[2] + nz];
              const double ref_value = (nx == mx && ny == my && nz == mz)
                                           ? (double)product3(fft_size)
                                           : 0.0;
              double current_error = fabs(my_value - ref_value);
              if (my_process == 0 && current_error > 1e-6) {
                printf("ERROR %i %i %i/%i %i %i: %f %f\n", nx, ny, nz, mx, my,
                       mz, my_value, ref_value);
              }
              max_error = fmax(max_error, current_error);
            }
          }
        }
      }
    }
  }
  fflush(stdout);

  fft_free_double(double_buffer);
  fft_free_complex(complex_buffer);

  if (max_error > 1e-6) {
    if (my_process == 0) {
      printf("The bw 3D C2R FFT does not work correctly (%i %i %i): %f!\n",
             fft_size[0], fft_size[1], fft_size[2], max_error);
      fflush(stdout);
    }
    errors++;
  }

  if (errors == 0 && my_process == 0)
    printf("The 3D R2C/C2R FFT does work correctly (%i %i %i)!\n", fft_size[0],
           fft_size[1], fft_size[2]);
  return errors;
}

/*******************************************************************************
 * \brief Function to test the local FFT backend.
 * \author Frederick Stein
 ******************************************************************************/
int fft_test_3d_local_r2c_inplace_low(const int fft_size[3],
                                      const int test_every) {
  const int my_process = cp_mpi_comm_rank(cp_mpi_get_comm_world());

  int errors = 0;

  const double pi = acos(-1);

  const int buffer_size = fft_size[0] * fft_size[1] * (fft_size[2] / 2 + 1);
  double *double_buffer = NULL;
  fft_allocate_double(2 * buffer_size, &double_buffer);
  double complex *complex_buffer = (double complex *)double_buffer;

  fft_register_3d_r2c_local(fft_size, double_buffer, complex_buffer);
  fft_register_3d_c2r_local(fft_size, complex_buffer, double_buffer);

  double max_error = 0.0;
  int number_of_tests = 0;
  for (int mx = 0; mx < fft_size[0]; mx++) {
    for (int my = 0; my < fft_size[1]; my++) {
      for (int mz = 0; mz < fft_size[2]; mz++) {
        if (test_every > 0 && number_of_tests % test_every != 0) {
          number_of_tests++;
          continue;
        }
        number_of_tests++;
        memset(double_buffer, 0, 2 * buffer_size * sizeof(double));
        double_buffer[(mx * fft_size[1] + my) * 2 * (fft_size[2] / 2 + 1) +
                      mz] = 1.0;
        fft_3d_r2c_local(fft_size, double_buffer, complex_buffer);

#pragma omp parallel for default(none)                                         \
    shared(complex_buffer, fft_size, pi, mx, my, mz, my_process)               \
    reduction(max : max_error) collapse(3)
        for (int nx = 0; nx < fft_size[0]; nx++) {
          for (int ny = 0; ny < fft_size[1]; ny++) {
            for (int nz = 0; nz < fft_size[2] / 2 + 1; nz++) {
              const double complex my_value =
                  complex_buffer[(nx * fft_size[1] + ny) *
                                     (fft_size[2] / 2 + 1) +
                                 nz];
              const double complex ref_value =
                  cexp(-2.0 * I * pi *
                       (((double)mx) * nx / fft_size[0] +
                        ((double)my) * ny / fft_size[1] +
                        ((double)mz) * nz / fft_size[2]));
              double current_error = cabs(my_value - ref_value);
              if (my_process == 0 && current_error > 1e-6) {
                printf("ERROR %i %i %i/%i %i %i: (%f %f) (%f %f)\n", nx, ny, nz,
                       mx, my, mz, creal(my_value), cimag(my_value),
                       creal(ref_value), cimag(ref_value));
              }
              max_error = fmax(max_error, current_error);
            }
          }
        }
      }
    }
  }
  fflush(stdout);

  if (max_error > 1.0e-6) {
    if (my_process == 0) {
      printf("The fw R2C 3D-FFT (inplace) does not work correctly (%i %i %i): "
             "%f!\n",
             fft_size[0], fft_size[1], fft_size[2], max_error);
      fflush(stdout);
    }
    errors++;
  }

  max_error = 0.0;
  number_of_tests = 0;
  for (int mx = 0; mx < fft_size[0]; mx++) {
    for (int my = 0; my < fft_size[1]; my++) {
      for (int mz = 0; mz < fft_size[2]; mz++) {
        if (test_every > 0 && number_of_tests % test_every != 0) {
          number_of_tests++;
          continue;
        }
        number_of_tests++;

#pragma omp parallel for default(none)                                         \
    shared(complex_buffer, fft_size, pi, mx, my, mz, my_process) collapse(3)
        for (int nz = 0; nz < fft_size[2] / 2 + 1; nz++) {
          for (int ny = 0; ny < fft_size[1]; ny++) {
            for (int nx = 0; nx < fft_size[0]; nx++) {
              complex_buffer[(nx * fft_size[1] + ny) * (fft_size[2] / 2 + 1) +
                             nz] = cexp(-2.0 * I * pi *
                                        (((double)mx) * nx / fft_size[0] +
                                         ((double)my) * ny / fft_size[1] +
                                         ((double)mz) * nz / fft_size[2]));
            }
          }
        }
        fft_3d_c2r_local(fft_size, complex_buffer, double_buffer);

#pragma omp parallel for default(none)                                         \
    shared(double_buffer, fft_size, pi, mx, my, mz, my_process)                \
    reduction(max : max_error) collapse(3)
        for (int nx = 0; nx < fft_size[0]; nx++) {
          for (int ny = 0; ny < fft_size[1]; ny++) {
            for (int nz = 0; nz < fft_size[2]; nz++) {
              const double my_value =
                  double_buffer[(nx * fft_size[1] + ny) * 2 *
                                    (fft_size[2] / 2 + 1) +
                                nz];
              const double ref_value = (nx == mx && ny == my && nz == mz)
                                           ? (double)product3(fft_size)
                                           : 0.0;
              double current_error = fabs(my_value - ref_value);
              if (my_process == 0 && current_error > 1e-6) {
                printf("ERROR %i %i %i/%i %i %i: %f %f\n", nx, ny, nz, mx, my,
                       mz, my_value, ref_value);
              }
              max_error = fmax(max_error, current_error);
            }
          }
        }
      }
    }
  }
  fflush(stdout);

  fft_free_double(double_buffer);

  if (max_error > 1e-6) {
    if (my_process == 0) {
      printf("The bw 3D C2R FFT (inplace) does not work correctly (%i %i %i): "
             "%f!\n",
             fft_size[0], fft_size[1], fft_size[2], max_error);
      fflush(stdout);
    }
    errors++;
  }

  if (errors == 0 && my_process == 0)
    printf("The 3D R2C/C2R FFT (inplace) does work correctly (%i %i %i)!\n",
           fft_size[0], fft_size[1], fft_size[2]);
  return errors;
}

/*******************************************************************************
 * \brief Function to test the local FFT backend (1-3D).
 * \author Frederick Stein
 ******************************************************************************/
int fft_test_local() {
  int errors = 0;

  clock_t begin = clock();
  errors += fft_test_1d_local_low(15, 13, true, true);
  errors += fft_test_1d_local_low(9, 11, true, false);
  errors += fft_test_1d_local_low(10, 14, false, true);
  errors += fft_test_1d_local_low(14, 13, false, false);

  errors += fft_test_1d_local_inplace_low(15, 13, true);
  errors += fft_test_1d_local_inplace_low(14, 13, false);

  errors += fft_test_1d_local_r2c_low(15, 13, true, false);
  errors += fft_test_1d_local_r2c_low(9, 11, false, false);
  errors += fft_test_1d_local_r2c_low(10, 14, false, true);
  errors += fft_test_1d_local_r2c_low(14, 13, true, true);

  errors += fft_test_1d_local_r2c_inplace_low(15, 11, false);
  errors += fft_test_1d_local_r2c_inplace_low(14, 13, true);

  errors += fft_test_2d_local_low((const int[2]){10, 10}, 10, true, true);
  errors += fft_test_2d_local_low((const int[2]){8, 9}, 10, true, false);
  errors += fft_test_2d_local_low((const int[2]){7, 10}, 11, false, true);
  errors += fft_test_2d_local_low((const int[2]){6, 7}, 10, false, false);

  errors += fft_test_2d_local_inplace_low((const int[2]){8, 9}, 10, true);
  errors += fft_test_2d_local_inplace_low((const int[2]){7, 10}, 11, false);

  errors += fft_test_2d_local_r2c_low((const int[2]){10, 10}, 10, true, true);
  errors += fft_test_2d_local_r2c_low((const int[2]){8, 9}, 10, true, false);
  errors += fft_test_2d_local_r2c_low((const int[2]){7, 10}, 10, false, true);
  errors += fft_test_2d_local_r2c_low((const int[2]){6, 7}, 11, false, false);

  errors += fft_test_2d_local_r2c_inplace_low((const int[2]){8, 9}, 10, true);
  errors += fft_test_2d_local_r2c_inplace_low((const int[2]){7, 10}, 10, false);

  // Reduce tests to ca 10 per set
  errors += fft_test_3d_local_low((const int[3]){8, 8, 8}, 10);
  errors += fft_test_3d_local_low((const int[3]){7, 5, 3}, 10);

  errors += fft_test_3d_local_inplace_low((const int[3]){8, 8, 8}, 10);
  errors += fft_test_3d_local_inplace_low((const int[3]){7, 5, 3}, 10);

  errors += fft_test_3d_local_r2c_low((const int[3]){8, 8, 8}, 10);
  errors += fft_test_3d_local_r2c_low((const int[3]){7, 5, 3}, 10);

  errors += fft_test_3d_local_r2c_inplace_low((const int[3]){8, 8, 8}, 10);
  errors += fft_test_3d_local_r2c_inplace_low((const int[3]){7, 5, 3}, 10);

  clock_t end = clock();
  if (cp_mpi_comm_rank(cp_mpi_get_comm_world()) == 0)
    printf("Time to test local FFTs: %f\n",
           (double)(end - begin) / CLOCKS_PER_SEC);

  return errors;
}

// EOF
