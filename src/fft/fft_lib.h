/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#ifndef FFT_LIB_H
#define FFT_LIB_H

#include "../mpiwrap/cp_mpi.h"
#include "fft_lib_fftw.h"

#include <complex.h>
#include <stdbool.h>

typedef enum { FFT_LIB_FFTW, FFT_LIB_GPU, FFT_LIB_DEFAULT } fft_lib;

/*******************************************************************************
 * \brief Initialize the FFT library (if not done externally).
 * \author Frederick Stein
 ******************************************************************************/
int fft_init_lib(const fft_lib lib, const int fftw_planning_flag,
                 const bool use_guru_interface, const char *wisdom_file);

/*******************************************************************************
 * \brief Finalize the FFT library (if not done externally).
 * \author Frederick Stein
 ******************************************************************************/
void fft_finalize_lib(const char *wisdom_file);

/*******************************************************************************
 * \brief Get the default library (GPU if offloading was enabled, else FFTW3).
 * \author Frederick Stein
 ******************************************************************************/
int fft_lib_default_library();

/*******************************************************************************
 * \brief Inquire the library backend in use.
 * \author Frederick Stein
 ******************************************************************************/
int fft_lib_backend_in_use();

/*******************************************************************************
 * \brief Whether compound MPI implementations are available.
 * \author Frederick Stein
 ******************************************************************************/
bool fft_lib_has_guru_interface();

/*******************************************************************************
 * \brief Allocate buffer of type double.
 * \author Frederick Stein
 ******************************************************************************/
void fft_allocate_double(const int length, double **buffer);

/*******************************************************************************
 * \brief Allocate buffer of type double complex.
 * \author Frederick Stein
 ******************************************************************************/
void fft_allocate_complex(const int length, double complex **buffer);

/*******************************************************************************
 * \brief Allocate buffer of type double.
 * \author Frederick Stein
 ******************************************************************************/
void fft_free_double(double *buffer);

/*******************************************************************************
 * \brief Allocate buffer of type double complex.
 * \author Frederick Stein
 ******************************************************************************/
void fft_free_complex(double complex *buffer);

/*******************************************************************************
 * \brief Register a local 1D C2C FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_1d_local(const bool dir, const int fft_size,
                           const int number_of_ffts, const bool transpose_rs,
                           const bool transpose_gs,
                           const int leading_dimension_rs,
                           const int leading_dimension_gs,
                           double complex *grid_in, double complex *grid_out);

/*******************************************************************************
 * \brief Register a local 1D R2C FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_1d_r2c_local(const int fft_size, const int number_of_ffts,
                               const bool transpose_rs, const bool transpose_gs,
                               const int leading_dimension_rs,
                               const int leading_dimension_gs, double *grid_in,
                               double complex *grid_out);

/*******************************************************************************
 * \brief Register a local 1D C2R FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_1d_c2r_local(const int fft_size, const int number_of_ffts,
                               const bool transpose_rs, const bool transpose_gs,
                               const int leading_dimension_rs,
                               const int leading_dimension_gs,
                               double complex *grid_in, double *grid_out);

/*******************************************************************************
 * \brief Register a local 2D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_2d_local(const bool dir, const int fft_size[2],
                           const int number_of_ffts, const bool transpose_rs,
                           const bool transpose_gs, double complex *grid_in,
                           double complex *grid_out);

/*******************************************************************************
 * \brief Register a local 2D R2C FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_2d_r2c_local(const int fft_size[2], const int number_of_ffts,
                               const bool transpose_rs, const bool transpose_gs,
                               double *grid_in, double complex *grid_out);

/*******************************************************************************
 * \brief Register a local 2D C2R FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_2d_c2r_local(const int fft_size[2], const int number_of_ffts,
                               const bool transpose_rs, const bool transpose_gs,
                               double complex *grid_in, double *grid_out);

/*******************************************************************************
 * \brief Register a local 3D FFT.
 * \note fft_3d_bw_local(grid_gs, grid_rs, n) is the reverse to
 * fft_3d_rw_local(grid_rs, grid_gs, n) (ignoring normalization).
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_3d_local(const bool dir, const int fft_size[3],
                           double complex *grid_in, double complex *grid_out);

/*******************************************************************************
 * \brief Register a local 3D R2C FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_3d_r2c_local(const int fft_size[3], double *grid_in,
                               double complex *grid_out);

/*******************************************************************************
 * \brief Register a local 3D C2R FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_3d_c2r_local(const int fft_size[3], double complex *grid_in,
                               double *grid_out);

/*******************************************************************************
 * \brief Perform a local 1D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_1d_local(const bool dir, const int fft_size, const int number_of_ffts,
                  const bool transpose_rs, const bool transpose_gs,
                  const int leading_dimension_rs,
                  const int leading_dimension_gs, double complex *grid_in,
                  double complex *grid_out);

/*******************************************************************************
 * \brief Perform a local 1D R2C FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_1d_r2c_local(const int fft_size, const int number_of_ffts,
                      const bool transpose_rs, const bool transpose_gs,
                      const int leading_dimension_rs,
                      const int leading_dimension_gs, double *grid_in,
                      double complex *grid_out);

/*******************************************************************************
 * \brief Perform a local 1D C2R FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_1d_c2r_local(const int fft_size, const int number_of_ffts,
                      const bool transpose_rs, const bool transpose_gs,
                      const int leading_dimension_rs,
                      const int leading_dimension_gs, double complex *grid_in,
                      double *grid_out);

/*******************************************************************************
 * \brief Perform a local 2D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_2d_local(const bool dir, const int fft_size[2],
                  const int number_of_ffts, const bool transpose_rs,
                  const bool transpose_gs, double complex *grid_in,
                  double complex *grid_out);

/*******************************************************************************
 * \brief Perform a local 2D R2C FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_2d_r2c_local(const int fft_size[2], const int number_of_ffts,
                      const bool transpose_rs, const bool transpose_gs,
                      double *grid_in, double complex *grid_out);

/*******************************************************************************
 * \brief Perform a local 2D C2R FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_2d_c2r_local(const int fft_size[2], const int number_of_ffts,
                      const bool transpose_rs, const bool transpose_gs,
                      double complex *grid_in, double *grid_out);

/*******************************************************************************
 * \brief Perform a local 3D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_3d_local(const bool dir, const int fft_size[3],
                  double complex *grid_in, double complex *grid_out);

/*******************************************************************************
 * \brief Perform a local 3D R2C FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_3d_r2c_local(const int fft_size[3], double *grid_in,
                      double complex *grid_out);

/*******************************************************************************
 * \brief Perform a local 3D C2R FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_3d_c2r_local(const int fft_size[3], double complex *grid_in,
                      double *grid_out);

#endif /* FFT_LIB_H */

// EOF
