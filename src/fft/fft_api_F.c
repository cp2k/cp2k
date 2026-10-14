/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include "fft_api_F.h"
#include "fft_lib.h"
#include "fft_timer.h"

#include <assert.h>
#include <stdio.h>
#include <string.h>

int fft_library_init_F(const int backend_F, const int fftw_plan,
                       const bool use_guru_interface, const char *wisdom_file) {
  fft_lib backend;
  switch (backend_F) {
  case FFT_LIBRARY_BACKEND_DEFAULT:
    backend = FFT_LIB_DEFAULT;
    break;
  case FFT_LIBRARY_BACKEND_FFTW:
    backend = FFT_LIB_FFTW;
    break;
  case FFT_LIBRARY_BACKEND_GPU:
    backend = FFT_LIB_GPU;
    break;
  default:
    return 1;
  }
  fftw_plan_type plan_type;
  switch (fftw_plan) {
  case FFT_LIBRARY_FFTW_MODE_ESTIMATE:
    plan_type = FFT_ESTIMATE;
    break;
  case FFT_LIBRARY_FFTW_MODE_MEASURE:
    plan_type = FFT_MEASURE;
    break;
  case FFT_LIBRARY_FFTW_MODE_PATIENT:
    plan_type = FFT_PATIENT;
    break;
  case FFT_LIBRARY_FFTW_MODE_EXHAUSTIVE:
    plan_type = FFT_EXHAUSTIVE;
    break;
  default:
    return 1;
  }
  const int status =
      fft_init_lib(backend, plan_type, use_guru_interface, wisdom_file);
  fft_init_timer(false);
  return status;
}

void fft_library_finalize_F(const char *wisdom_file) {
  fft_finalize_lib(wisdom_file);
}

void fft_print_timing_report_F(const double threshold) {
  fft_print_timing_report(threshold);
}

/*******************************************************************************
 * \brief Allocate buffer of type double.
 * \author Frederick Stein
 ******************************************************************************/
void fft_allocate_double_F(const int length, double **buffer) {
  fft_allocate_double(length, buffer);
}

/*******************************************************************************
 * \brief Allocate buffer of type double complex.
 * \author Frederick Stein
 ******************************************************************************/
void fft_allocate_complex_F(const int length, double complex **buffer) {
  fft_allocate_complex(length, buffer);
}

/*******************************************************************************
 * \brief Allocate buffer of type double.
 * \author Frederick Stein
 ******************************************************************************/
void fft_free_double_F(double *buffer) { fft_free_double(buffer); }

/*******************************************************************************
 * \brief Allocate buffer of type double complex.
 * \author Frederick Stein
 ******************************************************************************/
void fft_free_complex_F(double complex *buffer) { fft_free_complex(buffer); }

int fft_backend_in_use_F() { return fft_lib_backend_in_use(); }

bool fft_has_guru_interface_F() { return fft_lib_has_guru_interface(); }

void fft_register_1d_local_F(const bool dir, const int fft_size,
                             const int number_of_ffts, const bool transpose_in,
                             const bool transpose_out,
                             const int leading_dimension_in,
                             const int leading_dimension_out,
                             double complex *grid_in,
                             double complex *grid_out) {
  fft_register_1d_local(dir, fft_size, number_of_ffts, transpose_in,
                        transpose_out, leading_dimension_in,
                        leading_dimension_out, grid_in, grid_out);
}

void fft_register_2d_local_F(const bool dir, const int fft_size[2],
                             const int number_of_ffts, const bool transpose_in,
                             const bool transpose_out, double complex *grid_in,
                             double complex *grid_out) {
  fft_register_2d_local(dir, (const int[]){fft_size[1], fft_size[0]},
                        number_of_ffts, transpose_in, transpose_out, grid_in,
                        grid_out);
}

void fft_register_3d_local_F(const bool dir, const int fft_size[3],
                             double complex *grid_in,
                             double complex *grid_out) {
  fft_register_3d_local(dir,
                        (const int[]){fft_size[2], fft_size[1], fft_size[0]},
                        grid_in, grid_out);
}

void fft_register_1d_local_inplace_F(const bool dir, const int fft_size,
                                     const int number_of_ffts,
                                     const bool transpose,
                                     const int leading_dimension,
                                     double complex *grid) {
  fft_register_1d_local(dir, fft_size, number_of_ffts, transpose, transpose,
                        leading_dimension, leading_dimension, grid, grid);
}

void fft_register_2d_local_inplace_F(const bool dir, const int fft_size[2],
                                     const int number_of_ffts,
                                     const bool transpose,
                                     double complex *grid) {
  fft_register_2d_local(dir, (const int[]){fft_size[1], fft_size[0]},
                        number_of_ffts, transpose, transpose, grid, grid);
}

void fft_register_3d_local_inplace_F(const bool dir, const int fft_size[3],
                                     double complex *grid) {
  fft_register_3d_local(
      dir, (const int[]){fft_size[2], fft_size[1], fft_size[0]}, grid, grid);
}

void fft_1d_local_F(const bool dir, const int fft_size,
                    const int number_of_ffts, const bool transpose_in,
                    const bool transpose_out, const int leading_dimension_in,
                    const int leading_dimension_out, double complex *grid_in,
                    double complex *grid_out) {
  fft_1d_local(dir, fft_size, number_of_ffts, transpose_in, transpose_out,
               leading_dimension_in, leading_dimension_out, grid_in, grid_out);
}

void fft_2d_local_F(const bool dir, const int fft_size[2],
                    const int number_of_ffts, const bool transpose_in,
                    const bool transpose_out, double complex *grid_in,
                    double complex *grid_out) {
  fft_2d_local(dir, (const int[]){fft_size[1], fft_size[0]}, number_of_ffts,
               transpose_in, transpose_out, grid_in, grid_out);
}

void fft_3d_local_F(const bool dir, const int fft_size[3],
                    double complex *grid_in, double complex *grid_out) {
  fft_3d_local(dir, (const int[]){fft_size[2], fft_size[1], fft_size[0]},
               grid_in, grid_out);
}

void fft_1d_local_inplace_F(const bool dir, const int fft_size,
                            const int number_of_ffts, const bool transpose,
                            const int leading_dimension, double complex *grid) {
  fft_1d_local(dir, fft_size, number_of_ffts, transpose, transpose,
               leading_dimension, leading_dimension, grid, grid);
}

void fft_2d_local_inplace_F(const bool dir, const int fft_size[2],
                            const int number_of_ffts, const bool transpose,
                            double complex *grid) {
  fft_2d_local(dir, (const int[]){fft_size[1], fft_size[0]}, number_of_ffts,
               transpose, transpose, grid, grid);
}

void fft_3d_local_inplace_F(const bool dir, const int fft_size[3],
                            double complex *grid) {
  fft_3d_local(dir, (const int[]){fft_size[2], fft_size[1], fft_size[0]}, grid,
               grid);
}

// EOF
