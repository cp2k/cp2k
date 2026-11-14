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

/*******************************************************************************
 * \brief Fortran interface for fft_library_init.
 *
 * \param backend_F Fortran-side FFT backend selector.
 * \param fftw_plan FFTW planning mode.
 * \param use_guru_interface Whether the FFTW guru interface may be used.
 * \param wisdom_file Path to the FFTW wisdom file.
 * \author Frederick Stein
 ******************************************************************************/
int fft_library_init_F(const int backend_F, const int fftw_plan,
                       const bool use_guru_interface, const char *wisdom_file);

/*******************************************************************************
 * \brief Fortran interface for fft_library_finalize.
 *
 * \param wisdom_file Path to the FFTW wisdom file.
 * \author Frederick Stein
 ******************************************************************************/
void fft_library_finalize_F(const char *wisdom_file);

typedef void (*fft_print_func_t)(const char *, int, int);
/*******************************************************************************
 * \brief Fortran interface for fft_print_timing_report.
 *
 * \param threshold Minimum relative timing value to report.
 * \param fortran_comm Fortran MPI communicator handle.
 * \author Frederick Stein
 ******************************************************************************/
void fft_print_timing_report_F(const double threshold, const int fortran_comm);
/*******************************************************************************
 * \brief Fortran interface for fft_print_timing_report_unit.
 *
 * \param threshold Minimum relative timing value to report.
 * \param print_func Callback used to write the timing report.
 * \param output_unit Fortran output unit for the timing report.
 * \param fortran_comm Fortran MPI communicator handle.
 * \author Frederick Stein
 ******************************************************************************/
void fft_print_timing_report_unit_F(const double threshold,
                                    fft_print_func_t print_func,
                                    const int output_unit,
                                    const int fortran_comm);

/*******************************************************************************
 * \brief Fortran interface for fft_backend_in_use.
 *
 * \author Frederick Stein
 ******************************************************************************/
int fft_backend_in_use_F();

/*******************************************************************************
 * \brief Fortran interface for fft_has_guru_interface.
 *
 * \author Frederick Stein
 ******************************************************************************/
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

/*******************************************************************************
 * \brief Fortran interface for fft_register_1d_local.
 *
 * \param dir FFT direction: true for forward, false for backward.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose_in Whether the input has the transposed layout.
 * \param transpose_out Whether the output has the transposed layout.
 * \param leading_dimension_in Leading dimension of the input buffer.
 * \param leading_dimension_out Leading dimension of the output buffer.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_1d_local_F(const bool dir, const int fft_size,
                             const int number_of_ffts, const bool transpose_in,
                             const bool transpose_out,
                             const int leading_dimension_in,
                             const int leading_dimension_out,
                             double complex *grid_in, double complex *grid_out);

/*******************************************************************************
 * \brief Fortran interface for fft_register_2d_local.
 *
 * \param dir FFT direction: true for forward, false for backward.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose_in Whether the input has the transposed layout.
 * \param transpose_out Whether the output has the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_2d_local_F(const bool dir, const int fft_size[2],
                             const int number_of_ffts, const bool transpose_in,
                             const bool transpose_out, double complex *grid_in,
                             double complex *grid_out);

/*******************************************************************************
 * \brief Fortran interface for fft_register_3d_local.
 *
 * \param dir FFT direction: true for forward, false for backward.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose_in Whether the input has the transposed layout.
 * \param transpose_out Whether the output has the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_3d_local_F(const bool dir, const int fft_size[3],
                             const int number_of_ffts, const bool transpose_in,
                             const bool transpose_out, double complex *grid_in,
                             double complex *grid_out);

/*******************************************************************************
 * \brief Fortran interface for fft_register_1d_local_inplace.
 *
 * \param dir FFT direction: true for forward, false for backward.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose Whether the in-place buffer has the transposed layout.
 * \param leading_dimension Leading dimension of the in-place buffer.
 * \param grid In-place grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_1d_local_inplace_F(const bool dir, const int fft_size,
                                     const int number_of_ffts,
                                     const bool transpose,
                                     const int leading_dimension,
                                     double complex *grid);

/*******************************************************************************
 * \brief Fortran interface for fft_register_2d_local_inplace.
 *
 * \param dir FFT direction: true for forward, false for backward.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose Whether the in-place buffer has the transposed layout.
 * \param grid In-place grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_2d_local_inplace_F(const bool dir, const int fft_size[2],
                                     const int number_of_ffts,
                                     const bool transpose,
                                     double complex *grid);

/*******************************************************************************
 * \brief Fortran interface for fft_register_3d_local_inplace.
 *
 * \param dir FFT direction: true for forward, false for backward.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose Whether the in-place buffer has the transposed layout.
 * \param grid In-place grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_3d_local_inplace_F(const bool dir, const int fft_size[3],
                                     const int number_of_ffts,
                                     const bool transpose,
                                     double complex *grid);

/*******************************************************************************
 * \brief Fortran interface for fft_1d_local.
 *
 * \param dir FFT direction: true for forward, false for backward.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose_in Whether the input has the transposed layout.
 * \param transpose_out Whether the output has the transposed layout.
 * \param leading_dimension_in Leading dimension of the input buffer.
 * \param leading_dimension_out Leading dimension of the output buffer.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_1d_local_F(const bool dir, const int fft_size,
                    const int number_of_ffts, const bool transpose_in,
                    const bool transpose_out, const int leading_dimension_in,
                    const int leading_dimension_out, double complex *grid_in,
                    double complex *grid_out);

/*******************************************************************************
 * \brief Fortran interface for fft_2d_local.
 *
 * \param dir FFT direction: true for forward, false for backward.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose_in Whether the input has the transposed layout.
 * \param transpose_out Whether the output has the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_2d_local_F(const bool dir, const int fft_size[2],
                    const int number_of_ffts, const bool transpose_in,
                    const bool transpose_out, double complex *grid_in,
                    double complex *grid_out);

/*******************************************************************************
 * \brief Fortran interface for fft_3d_local.
 *
 * \param dir FFT direction: true for forward, false for backward.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose_in Whether the input has the transposed layout.
 * \param transpose_out Whether the output has the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_3d_local_F(const bool dir, const int fft_size[3],
                    const int number_of_ffts, const bool transpose_in,
                    const bool transpose_out, double complex *grid_in,
                    double complex *grid_out);

/*******************************************************************************
 * \brief Fortran interface for fft_1d_local_inplace.
 *
 * \param dir FFT direction: true for forward, false for backward.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose Whether the in-place buffer has the transposed layout.
 * \param leading_dimension Leading dimension of the in-place buffer.
 * \param grid In-place grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_1d_local_inplace_F(const bool dir, const int fft_size,
                            const int number_of_ffts, const bool transpose,
                            const int leading_dimension, double complex *grid);

/*******************************************************************************
 * \brief Fortran interface for fft_2d_local_inplace.
 *
 * \param dir FFT direction: true for forward, false for backward.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose Whether the in-place buffer has the transposed layout.
 * \param grid In-place grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_2d_local_inplace_F(const bool dir, const int fft_size[2],
                            const int number_of_ffts, const bool transpose,
                            double complex *grid);

/*******************************************************************************
 * \brief Fortran interface for fft_3d_local_inplace.
 *
 * \param dir FFT direction: true for forward, false for backward.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose Whether the in-place buffer has the transposed layout.
 * \param grid In-place grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_3d_local_inplace_F(const bool dir, const int fft_size[3],
                            const int number_of_ffts, const bool transpose,
                            double complex *grid);

/*******************************************************************************
 * \brief Fortran interface for fft_register_1d_r2c_local.
 *
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose_in Whether the input has the transposed layout.
 * \param transpose_out Whether the output has the transposed layout.
 * \param leading_dimension_in Leading dimension of the input buffer.
 * \param leading_dimension_out Leading dimension of the output buffer.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_1d_r2c_local_F(const int fft_size, const int number_of_ffts,
                                 const bool transpose_in,
                                 const bool transpose_out,
                                 const int leading_dimension_in,
                                 const int leading_dimension_out,
                                 double *grid_in, double complex *grid_out);
/*******************************************************************************
 * \brief Fortran interface for fft_register_1d_c2r_local.
 *
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose_in Whether the input has the transposed layout.
 * \param transpose_out Whether the output has the transposed layout.
 * \param leading_dimension_in Leading dimension of the input buffer.
 * \param leading_dimension_out Leading dimension of the output buffer.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_1d_c2r_local_F(const int fft_size, const int number_of_ffts,
                                 const bool transpose_in,
                                 const bool transpose_out,
                                 const int leading_dimension_in,
                                 const int leading_dimension_out,
                                 double complex *grid_in, double *grid_out);
/*******************************************************************************
 * \brief Fortran interface for fft_register_2d_r2c_local.
 *
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose_in Whether the input has the transposed layout.
 * \param transpose_out Whether the output has the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_2d_r2c_local_F(const int fft_size[2],
                                 const int number_of_ffts,
                                 const bool transpose_in,
                                 const bool transpose_out, double *grid_in,
                                 double complex *grid_out);
/*******************************************************************************
 * \brief Fortran interface for fft_register_2d_c2r_local.
 *
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose_in Whether the input has the transposed layout.
 * \param transpose_out Whether the output has the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_2d_c2r_local_F(const int fft_size[2],
                                 const int number_of_ffts,
                                 const bool transpose_in,
                                 const bool transpose_out,
                                 double complex *grid_in, double *grid_out);
/*******************************************************************************
 * \brief Fortran interface for fft_register_3d_r2c_local.
 *
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose_in Whether the input has the transposed layout.
 * \param transpose_out Whether the output has the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_3d_r2c_local_F(const int fft_size[3],
                                 const int number_of_ffts,
                                 const bool transpose_in,
                                 const bool transpose_out, double *grid_in,
                                 double complex *grid_out);
/*******************************************************************************
 * \brief Fortran interface for fft_register_3d_c2r_local.
 *
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose_in Whether the input has the transposed layout.
 * \param transpose_out Whether the output has the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_3d_c2r_local_F(const int fft_size[3],
                                 const int number_of_ffts,
                                 const bool transpose_in,
                                 const bool transpose_out,
                                 double complex *grid_in, double *grid_out);

/*******************************************************************************
 * \brief Fortran interface for fft_1d_r2c_local.
 *
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose_in Whether the input has the transposed layout.
 * \param transpose_out Whether the output has the transposed layout.
 * \param leading_dimension_in Leading dimension of the input buffer.
 * \param leading_dimension_out Leading dimension of the output buffer.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_1d_r2c_local_F(const int fft_size, const int number_of_ffts,
                        const bool transpose_in, const bool transpose_out,
                        const int leading_dimension_in,
                        const int leading_dimension_out, double *grid_in,
                        double complex *grid_out);
/*******************************************************************************
 * \brief Fortran interface for fft_1d_c2r_local.
 *
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose_in Whether the input has the transposed layout.
 * \param transpose_out Whether the output has the transposed layout.
 * \param leading_dimension_in Leading dimension of the input buffer.
 * \param leading_dimension_out Leading dimension of the output buffer.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_1d_c2r_local_F(const int fft_size, const int number_of_ffts,
                        const bool transpose_in, const bool transpose_out,
                        const int leading_dimension_in,
                        const int leading_dimension_out,
                        double complex *grid_in, double *grid_out);
/*******************************************************************************
 * \brief Fortran interface for fft_2d_r2c_local.
 *
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose_in Whether the input has the transposed layout.
 * \param transpose_out Whether the output has the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_2d_r2c_local_F(const int fft_size[2], const int number_of_ffts,
                        const bool transpose_in, const bool transpose_out,
                        double *grid_in, double complex *grid_out);
/*******************************************************************************
 * \brief Fortran interface for fft_2d_c2r_local.
 *
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose_in Whether the input has the transposed layout.
 * \param transpose_out Whether the output has the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_2d_c2r_local_F(const int fft_size[2], const int number_of_ffts,
                        const bool transpose_in, const bool transpose_out,
                        double complex *grid_in, double *grid_out);
/*******************************************************************************
 * \brief Fortran interface for fft_3d_r2c_local.
 *
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose_in Whether the input has the transposed layout.
 * \param transpose_out Whether the output has the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_3d_r2c_local_F(const int fft_size[3], const int number_of_ffts,
                        const bool transpose_in, const bool transpose_out,
                        double *grid_in, double complex *grid_out);
/*******************************************************************************
 * \brief Fortran interface for fft_3d_c2r_local.
 *
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose_in Whether the input has the transposed layout.
 * \param transpose_out Whether the output has the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_3d_c2r_local_F(const int fft_size[3], const int number_of_ffts,
                        const bool transpose_in, const bool transpose_out,
                        double complex *grid_in, double *grid_out);

/*******************************************************************************
 * \brief Fortran interface for fft_register_1d_r2c_local_inplace.
 *
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose Whether the in-place buffer has the transposed layout.
 * \param leading_dimension Leading dimension of the in-place buffer.
 * \param grid In-place grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_1d_r2c_local_inplace_F(const int fft_size,
                                         const int number_of_ffts,
                                         const bool transpose,
                                         const int leading_dimension,
                                         double *grid);
/*******************************************************************************
 * \brief Fortran interface for fft_register_1d_c2r_local_inplace.
 *
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose Whether the in-place buffer has the transposed layout.
 * \param leading_dimension Leading dimension of the in-place buffer.
 * \param grid In-place grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_1d_c2r_local_inplace_F(const int fft_size,
                                         const int number_of_ffts,
                                         const bool transpose,
                                         const int leading_dimension,
                                         double *grid);
/*******************************************************************************
 * \brief Fortran interface for fft_register_2d_r2c_local_inplace.
 *
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose Whether the in-place buffer has the transposed layout.
 * \param grid In-place grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_2d_r2c_local_inplace_F(const int fft_size[2],
                                         const int number_of_ffts,
                                         const bool transpose, double *grid);
/*******************************************************************************
 * \brief Fortran interface for fft_register_2d_c2r_local_inplace.
 *
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose Whether the in-place buffer has the transposed layout.
 * \param grid In-place grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_2d_c2r_local_inplace_F(const int fft_size[2],
                                         const int number_of_ffts,
                                         const bool transpose, double *grid);
/*******************************************************************************
 * \brief Fortran interface for fft_register_3d_r2c_local_inplace.
 *
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose Whether the in-place buffer has the transposed layout.
 * \param grid In-place grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_3d_r2c_local_inplace_F(const int fft_size[3],
                                         const int number_of_ffts,
                                         const bool transpose, double *grid);
/*******************************************************************************
 * \brief Fortran interface for fft_register_3d_c2r_local_inplace.
 *
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose Whether the in-place buffer has the transposed layout.
 * \param grid In-place grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_3d_c2r_local_inplace_F(const int fft_size[3],
                                         const int number_of_ffts,
                                         const bool transpose, double *grid);

/*******************************************************************************
 * \brief Fortran interface for fft_1d_r2c_local_inplace.
 *
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose Whether the in-place buffer has the transposed layout.
 * \param leading_dimension Leading dimension of the in-place buffer.
 * \param grid In-place grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_1d_r2c_local_inplace_F(const int fft_size, const int number_of_ffts,
                                const bool transpose,
                                const int leading_dimension, double *grid);
/*******************************************************************************
 * \brief Fortran interface for fft_1d_c2r_local_inplace.
 *
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose Whether the in-place buffer has the transposed layout.
 * \param leading_dimension Leading dimension of the in-place buffer.
 * \param grid In-place grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_1d_c2r_local_inplace_F(const int fft_size, const int number_of_ffts,
                                const bool transpose,
                                const int leading_dimension, double *grid);
/*******************************************************************************
 * \brief Fortran interface for fft_2d_r2c_local_inplace.
 *
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose Whether the in-place buffer has the transposed layout.
 * \param grid In-place grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_2d_r2c_local_inplace_F(const int fft_size[2], const int number_of_ffts,
                                const bool transpose, double *grid);
/*******************************************************************************
 * \brief Fortran interface for fft_2d_c2r_local_inplace.
 *
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose Whether the in-place buffer has the transposed layout.
 * \param grid In-place grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_2d_c2r_local_inplace_F(const int fft_size[2], const int number_of_ffts,
                                const bool transpose, double *grid);
/*******************************************************************************
 * \brief Fortran interface for fft_3d_r2c_local_inplace.
 *
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose Whether the in-place buffer has the transposed layout.
 * \param grid In-place grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_3d_r2c_local_inplace_F(const int fft_size[3], const int number_of_ffts,
                                const bool transpose, double *grid);
/*******************************************************************************
 * \brief Fortran interface for fft_3d_c2r_local_inplace.
 *
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent transforms (batch size).
 * \param transpose Whether the in-place buffer has the transposed layout.
 * \param grid In-place grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_3d_c2r_local_inplace_F(const int fft_size[3], const int number_of_ffts,
                                const bool transpose, double *grid);

#endif
// EOF
