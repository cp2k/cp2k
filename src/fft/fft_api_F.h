/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#ifndef FFT_API_F_H
#define FFT_API_F_H

#include <complex.h>
#include <stdbool.h>

// Keep in accordance to fft_api.F
#define FFT_LIBRARY_BACKEND_DEFAULT 1
#define FFT_LIBRARY_BACKEND_FFTW 2
#define FFT_LIBRARY_BACKEND_GPU 3
#define FFT_LIBRARY_FFTW_MODE_ESTIMATE 11
#define FFT_LIBRARY_FFTW_MODE_MEASURE 12
#define FFT_LIBRARY_FFTW_MODE_PATIENT 13
#define FFT_LIBRARY_FFTW_MODE_EXHAUSTIVE 14

int fft_library_init_F(const int backend_F, const int fftw_plan,
                       const bool use_guru_interface, const char *wisdom_file);

void fft_library_finalize_F(const char *wisdom_file);

void fft_print_timing_report_F(const double threshold);

int fft_backend_in_use_F();

bool fft_has_guru_interface_F();

/*******************************************************************************
 * \brief Allocate buffer of type double.
 * \author Frederick Stein
 ******************************************************************************/
void fft_allocate_double_F(const int length, double **buffer);

/*******************************************************************************
 * \brief Allocate buffer of type double complex.
 * \author Frederick Stein
 ******************************************************************************/
void fft_allocate_complex_F(const int length, double complex **buffer);

/*******************************************************************************
 * \brief Allocate buffer of type double.
 * \author Frederick Stein
 ******************************************************************************/
void fft_free_double_F(double *buffer);

/*******************************************************************************
 * \brief Allocate buffer of type double complex.
 * \author Frederick Stein
 ******************************************************************************/
void fft_free_complex_F(double complex *buffer);

void fft_register_1d_local_F(const bool dir, const int fft_size,
                             const int number_of_ffts, const bool transpose_in,
                             const bool transpose_out,
                             const int leading_dimension_in,
                             const int leading_dimension_out,
                             double complex *grid_in, double complex *grid_out);

void fft_register_2d_local_F(const bool dir, const int fft_size[2],
                             const int number_of_ffts, const bool transpose_in,
                             const bool transpose_out, double complex *grid_in,
                             double complex *grid_out);

void fft_register_3d_local_F(const bool dir, const int fft_size[3],
                             double complex *grid_in, double complex *grid_out);

void fft_register_1d_local_inplace_F(const bool dir, const int fft_size,
                                     const int number_of_ffts,
                                     const bool transpose,
                                     const int leading_dimension,
                                     double complex *grid);

void fft_register_2d_local_inplace_F(const bool dir, const int fft_size[2],
                                     const int number_of_ffts,
                                     const bool transpose,
                                     double complex *grid);

void fft_register_3d_local_inplace_F(const bool dir, const int fft_size[3],
                                     double complex *grid);

void fft_1d_local_F(const bool dir, const int fft_size,
                    const int number_of_ffts, const bool transpose_in,
                    const bool transpose_out, const int leading_dimension_in,
                    const int leading_dimension_out, double complex *grid_in,
                    double complex *grid_out);

void fft_2d_local_F(const bool dir, const int fft_size[2],
                    const int number_of_ffts, const bool transpose_in,
                    const bool transpose_out, double complex *grid_in,
                    double complex *grid_out);

void fft_3d_local_F(const bool dir, const int fft_size[3],
                    double complex *grid_in, double complex *grid_out);

void fft_1d_local_inplace_F(const bool dir, const int fft_size,
                            const int number_of_ffts, const bool transpose,
                            const int leading_dimension, double complex *grid);

void fft_2d_local_inplace_F(const bool dir, const int fft_size[2],
                            const int number_of_ffts, const bool transpose,
                            double complex *grid);

void fft_3d_local_inplace_F(const bool dir, const int fft_size[3],
                            double complex *grid);

#endif
// EOF
