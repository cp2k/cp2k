/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#ifndef FFT_LIB_FFTW_H
#define FFT_LIB_FFTW_H

#include "../mpiwrap/cp_mpi.h"

// We include this first to require FFTW to use C complex numbers
#include <complex.h>

#if defined(__FFTW3)
#include <fftw3.h>
typedef fftw_iodim fft_iodim;
#else
typedef struct {
  int n;
  int is;
  int os;
} fft_iodim;
#endif
#include <stdbool.h>

typedef enum {
  FFT_ESTIMATE,
  FFT_MEASURE,
  FFT_PATIENT,
  FFT_EXHAUSTIVE
} fftw_plan_type;

/*******************************************************************************
 * \brief Initialize the FFT library (if not done externally).
 * \author Frederick Stein, Ole Schuett
 ******************************************************************************/
int fft_fftw_init_lib(const fftw_plan_type fftw_planning_flag,
                      const bool use_guru_interface, const char *wisdom_file);

/*******************************************************************************
 * \brief Finalize the FFT library (if not done externally).
 * \author Frederick Stein, Ole Schuett
 ******************************************************************************/
void fft_fftw_finalize_lib(const char *wisdom_file);

/*******************************************************************************
 * \brief Whether a distributed FFT implementation is available.
 * \author Frederick Stein
 ******************************************************************************/
bool fft_fftw_lib_use_mpi();

/*******************************************************************************
 * \brief Whether the guru interface is available.
 * \author Frederick Stein
 ******************************************************************************/
bool fft_fftw_lib_has_guru_interface();

/*******************************************************************************
 * \brief Allocate buffer of type double.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_allocate_double(const int length, double **buffer);

/*******************************************************************************
 * \brief Allocate buffer of type double complex.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_allocate_complex(const int length, double complex **buffer);

/*******************************************************************************
 * \brief Free buffer of type double.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_free_double(double *buffer);

/*******************************************************************************
 * \brief Free buffer of type double complex.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_free_complex(double complex *buffer);

/*******************************************************************************
 * \brief Performs a local C2C 1D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_1d_local(
    const bool direction, const int fft_size, const int number_of_ffts,
    const bool transpose_in, const bool transpose_out,
    const int leading_dimension_in, const int leading_dimension_out,
    double complex *grid_in, double complex *grid_out);

/*******************************************************************************
 * \brief Performs a local forward R2C FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_1d_r2c_local(
    const int fft_size, const int number_of_ffts, const bool transpose_in,
    const bool transpose_out, const int leading_dimension_in,
    const int leading_dimension_out, double *grid_in, double complex *grid_out);

/*******************************************************************************
 * \brief Performs a local backwards C2R 1D FFT
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_1d_c2r_local(
    const int fft_size, const int number_of_ffts, const bool transpose_in,
    const bool transpose_out, const int leading_dimension_in,
    const int leading_dimension_out, double complex *grid_in, double *grid_out);

/*******************************************************************************
 * \brief Performs a local C2C 2D FFT
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_2d_local(const bool direction, const int fft_size[2],
                                const int number_of_ffts,
                                const bool transpose_in,
                                const bool transpose_out,
                                double complex *grid_in,
                                double complex *grid_out);

/*******************************************************************************
 * \brief Performs a local forward R2C 2D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_2d_r2c_local(const int fft_size[2],
                                    const int number_of_ffts,
                                    const bool transpose_in,
                                    const bool transpose_out, double *grid_in,
                                    double complex *grid_out);

/*******************************************************************************
 * \brief Performs a local backwards C2R 2D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_2d_c2r_local(const int fft_size[2],
                                    const int number_of_ffts,
                                    const bool transpose_in,
                                    const bool transpose_out,
                                    double complex *grid_in, double *grid_out);

/*******************************************************************************
 * \brief Performs a local C2C 3D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_guru(const bool direction, int rank,
                            const fft_iodim *dims, int howmany_rank,
                            const fft_iodim *howmany_dims,
                            const int number_of_threads,
                            double complex *grid_in, double complex *grid_out);

/*******************************************************************************
 * \brief Performs a local forward R2C 3D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_r2c_guru(int rank, const fft_iodim *dims,
                                int howmany_rank, const fft_iodim *howmany_dims,
                                const int number_of_threads, double *grid_in,
                                double complex *grid_out);

/*******************************************************************************
 * \brief Performs a local backwards R2C 3D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_c2r_guru(int rank, const fft_iodim *dims,
                                int howmany_rank, const fft_iodim *howmany_dims,
                                const int number_of_threads,
                                double complex *grid_in, double *grid_out);

/*******************************************************************************
 * \brief Performs a local C2C 3D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_3d_local(const bool direction, const int fft_size[3],
                                double complex *grid_in,
                                double complex *grid_out);

/*******************************************************************************
 * \brief Performs a local forward R2C 3D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_3d_r2c_local(const int fft_size[3], double *grid_in,
                                    double complex *grid_out);

/*******************************************************************************
 * \brief Performs a local backwards R2C 3D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_3d_c2r_local(const int fft_size[3],
                                    double complex *grid_in, double *grid_out);

/*******************************************************************************
 * \brief Performs a distributed C2C 2D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_2d_distributed(const bool direction,
                                      const int npts_global[2],
                                      const int number_of_ffts,
                                      const cp_mpi_comm_t comm,
                                      double complex *grid_in,
                                      double complex *grid_out);

/*******************************************************************************
 * \brief Performs a distributed forward R2C 2D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_2d_r2c_distributed(const int npts_global[2],
                                          const int number_of_ffts,
                                          const cp_mpi_comm_t comm,
                                          double *grid_in,
                                          double complex *grid_out);

/*******************************************************************************
 * \brief Performs a distributed backwards C2R 2D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_2d_c2r_distributed(const int npts_global[2],
                                          const int number_of_ffts,
                                          const cp_mpi_comm_t comm,
                                          double complex *grid_in,
                                          double *grid_out);

/*******************************************************************************
 * \brief Performs a distributed C2C 3D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_3d_distributed(const bool direction,
                                      const int npts_global[3],
                                      const cp_mpi_comm_t comm,
                                      double complex *grid_in,
                                      double complex *grid_out);

/*******************************************************************************
 * \brief Performs a distributed forward R2C 3D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_3d_r2c_distributed(const int npts_global[3],
                                          const cp_mpi_comm_t comm,
                                          double *grid_in,
                                          double complex *grid_out);

/*******************************************************************************
 * \brief Performs a distributed backwards C2R 3D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_register_3d_c2r_distributed(const int npts_global[3],
                                          const cp_mpi_comm_t comm,
                                          double complex *grid_in,
                                          double *grid_out);

/*******************************************************************************
 * \brief Performs a local forward C2C 1D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_1d_local(const bool direction, const int fft_size,
                       const int number_of_ffts, const bool transpose_in,
                       const bool transpose_out, const int leading_dimension_in,
                       const int leading_dimension_out, double complex *grid_in,
                       double complex *grid_out);

/*******************************************************************************
 * \brief Performs a local forward R2C FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_1d_r2c_local(const int fft_size, const int number_of_ffts,
                           const bool transpose_in, const bool transpose_out,
                           const int leading_dimension_in,
                           const int leading_dimension_out, double *grid_in,
                           double complex *grid_out);

/*******************************************************************************
 * \brief Performs a local backwards C2R 1D FFT
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_1d_c2r_local(const int fft_size, const int number_of_ffts,
                           const bool transpose_in, const bool transpose_out,
                           const int leading_dimension_in,
                           const int leading_dimension_out,
                           double complex *grid_in, double *grid_out);

/*******************************************************************************
 * \brief Performs a local C2C 2D FFT
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_2d_local(const bool direction, const int fft_size[2],
                       const int number_of_ffts, const bool transpose_in,
                       const bool transpose_out, double complex *grid_in,
                       double complex *grid_out);

/*******************************************************************************
 * \brief Performs a local forward R2C 2D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_2d_r2c_local(const int fft_size[2], const int number_of_ffts,
                           const bool transpose_in, const bool transpose_out,
                           double *grid_in, double complex *grid_out);

/*******************************************************************************
 * \brief Performs a local backwards C2R 2D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_2d_c2r_local(const int fft_size[2], const int number_of_ffts,
                           const bool transpose_in, const bool transpose_out,
                           double complex *grid_in, double *grid_out);

/*******************************************************************************
 * \brief Performs a local C2C 3D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_guru(const bool direction, const int rank, const fft_iodim *dims,
                   const int howmany_rank, const fft_iodim *howmany_dims,
                   const int number_of_threads, double complex *grid_in,
                   double complex *grid_out);

/*******************************************************************************
 * \brief Performs a local forward R2C 3D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_guru_r2c(int rank, const fft_iodim *dims, int howmany_rank,
                       const fft_iodim *howmany_dims,
                       const int number_of_threads, double *grid_in,
                       double complex *grid_out);

/*******************************************************************************
 * \brief Performs a local backwards R2C 3D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_guru_c2r(int rank, const fft_iodim *dims, int howmany_rank,
                       const fft_iodim *howmany_dims,
                       const int number_of_threads, double complex *grid_in,
                       double *grid_out);

/*******************************************************************************
 * \brief Performs a local C2C 3D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_3d_local(const bool direction, const int fft_size[3],
                       double complex *grid_in, double complex *grid_out);

/*******************************************************************************
 * \brief Performs a local forward R2C 3D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_3d_r2c_local(const int fft_size[3], double *grid_in,
                           double complex *grid_out);

/*******************************************************************************
 * \brief Performs a local backwards R2C 3D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_3d_c2r_local(const int fft_size[3], double complex *grid_in,
                           double *grid_out);

/*******************************************************************************
 * \brief Returns sizes and starts of distributed C2C 2D FFTs.
 * \author Frederick Stein
 ******************************************************************************/
int fft_fftw_2d_distributed_sizes(const int npts_global[2],
                                  const int number_of_ffts,
                                  const cp_mpi_comm_t comm, int *local_n0,
                                  int *local_n0_start, int *local_n1,
                                  int *local_n1_start);

/*******************************************************************************
 * \brief Returns sizes and starts of distributed R2C/C2R 2D FFTs.
 * \author Frederick Stein
 ******************************************************************************/
int fft_fftw_2d_distributed_sizes_r2c(const int npts_global[2],
                                      const int number_of_ffts,
                                      const cp_mpi_comm_t comm, int *local_n0,
                                      int *local_n0_start, int *local_n1,
                                      int *local_n1_start);

/*******************************************************************************
 * \brief Returns sizes and starts of distributed C2C 3D FFTs.
 * \author Frederick Stein
 ******************************************************************************/
int fft_fftw_3d_distributed_sizes(const int npts_global[3],
                                  const cp_mpi_comm_t comm, int *local_n0,
                                  int *local_n0_start, int *local_n1,
                                  int *local_n1_start);

/*******************************************************************************
 * \brief Returns sizes and starts of distributed R2C/C2R 3D FFTs.
 * \author Frederick Stein
 ******************************************************************************/
int fft_fftw_3d_distributed_sizes_r2c(const int npts_global[3],
                                      const cp_mpi_comm_t comm, int *local_n0,
                                      int *local_n0_start, int *local_n1,
                                      int *local_n1_start);

/*******************************************************************************
 * \brief Performs a distributed forward C2C 2D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_2d_distributed(const bool direction, const int npts_global[2],
                             const int number_of_ffts, const cp_mpi_comm_t comm,
                             double complex *restrict grid_in,
                             double complex *restrict grid_out);

/*******************************************************************************
 * \brief Performs a distributed forward R2C 2D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_2d_r2c_distributed(const int npts_global[2],
                                 const int number_of_ffts,
                                 const cp_mpi_comm_t comm,
                                 double *restrict grid_in,
                                 double complex *restrict grid_out);

/*******************************************************************************
 * \brief Performs a distributed backwards C2R 2D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_2d_c2r_distributed(const int npts_global[2],
                                 const int number_of_ffts,
                                 const cp_mpi_comm_t comm,
                                 double complex *restrict grid_in,
                                 double *restrict grid_out);

/*******************************************************************************
 * \brief Performs a distributed C2C 3D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_3d_distributed(const bool direction, const int npts_global[3],
                             const cp_mpi_comm_t comm,
                             double complex *restrict grid_in,
                             double complex *restrict grid_out);

/*******************************************************************************
 * \brief Performs a distributed forward R2C 3D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_3d_r2c_distributed(const int npts_global[3],
                                 const cp_mpi_comm_t comm,
                                 double *restrict grid_in,
                                 double complex *restrict grid_out);

/*******************************************************************************
 * \brief Performs a distributed backwards C2R 3D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_fftw_3d_c2r_distributed(const int npts_global[3],
                                 const cp_mpi_comm_t comm,
                                 double complex *restrict grid_in,
                                 double *restrict grid_out);

#endif /* FFT_LIB_FFTW_H */

// EOF
