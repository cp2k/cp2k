/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#include <complex.h>
#include <stdbool.h>
#include <stddef.h>

/*******************************************************************************
 * \brief Initializes the fft_gpu library.
 * \author Ole Schuett
 ******************************************************************************/
void fft_gpu_init(void);

/*******************************************************************************
 * \brief Releases resources held by the fft_gpu library.
 * \author Ole Schuett
 ******************************************************************************/
void fft_gpu_finalize(void);

/*******************************************************************************
 * \brief Checks size of device buffers and re-allocates them if necessary.
 * \author Ole Schuett
 ******************************************************************************/
void ensure_memory_sizes(const size_t requested_buffer_size,
                         const size_t requested_map_size);

/*******************************************************************************
 * \brief Allocate buffer of type double.
 * \author Frederick Stein
 ******************************************************************************/
void fft_gpu_allocate_double(const int length, double **buffer);

/*******************************************************************************
 * \brief Allocate buffer of type double complex.
 * \author Frederick Stein
 ******************************************************************************/
void fft_gpu_allocate_complex(const int length, double complex **buffer);

/*******************************************************************************
 * \brief Allocate buffer of type double.
 * \author Frederick Stein
 ******************************************************************************/
void fft_gpu_free_double(double *buffer);

/*******************************************************************************
 * \brief Allocate buffer of type double complex.
 * \author Frederick Stein
 ******************************************************************************/
void fft_gpu_free_complex(double complex *buffer);

/*******************************************************************************
 * \brief   Performs a (double precision complex) 1D-FFT on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void fft_register_gpu_f(const bool dir, const int n, const int m,
                        const bool transpose_in, const bool transpose_out,
                        const int leading_dimension_in,
                        const int leading_dimension_out, const bool inplace);

/*******************************************************************************
 * \brief   Performs a (double precision complex) R2C 1D-FFT on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void fft_register_r2c_gpu_f(const bool dir, const int n, const int m,
                            const bool transpose_in, const bool transpose_out,
                            const int leading_dimension_in,
                            const int leading_dimension_out,
                            const bool inplace);

/*******************************************************************************
 * \brief   Performs a (double precision complex) 1D-FFT on the GPU.
 * \author  Frederick Stein
 ******************************************************************************/
void fft_register_gpu_ff(const bool dir, const int n[2], const int m,
                         const bool transpose_in, const bool transpose_out,
                         const bool inplace);

/*******************************************************************************
 * \brief   Performs a (double precision complex) R2C 1D-FFT on the GPU.
 * \author  Frederick Stein
 ******************************************************************************/
void fft_register_r2c_gpu_ff(const bool dir, const int n[2], const int m,
                             const bool transpose_in, const bool transpose_out,
                             const bool inplace);

/*******************************************************************************
 * \brief   Performs a (double precision complex) 3D-FFT on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void fft_register_gpu_fff(const bool dir, const int *npts, const bool inplace);

/*******************************************************************************
 * \brief   Performs a 3D-R2C/C2R-FFT, on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void fft_register_r2c_gpu_fff(const bool dir, const int *npts,
                              const bool inplace);

/*******************************************************************************
 * \brief   Performs a (double precision complex) 1D-FFT on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void fft_gpu_f(const double *zin, double *zout, const bool dir, const int n,
               const int m, const bool transpose_in, const bool transpose_out,
               const int leading_dimension_in, const int leading_dimension_out);

/*******************************************************************************
 * \brief   Performs a (double precision complex) R2C 1D-FFT on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void fft_r2c_gpu_f(const double *zin, double *zout, const bool dir, const int n,
                   const int m, const bool transpose_in,
                   const bool transpose_out, const int leading_dimension_in,
                   const int leading_dimension_out);

/*******************************************************************************
 * \brief   Performs a (double precision complex) 1D-FFT on the GPU.
 * \author  Frederick Stein
 ******************************************************************************/
void fft_gpu_ff(const double *zin, double *zout, const bool dir, const int n[2],
                const int m, const bool transpose_in, const bool transpose_out);

/*******************************************************************************
 * \brief   Performs a (double precision complex) R2C 1D-FFT on the GPU.
 * \author  Frederick Stein
 ******************************************************************************/
void fft_r2c_gpu_ff(const double *zin, double *zout, const bool dir,
                    const int n[2], const int m, const bool transpose_in,
                    const bool transpose_out);

/*******************************************************************************
 * \brief   Performs a (double precision complex) 3D-FFT on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void fft_gpu_fff(const double *zin, double *zout, const bool dir,
                 const int *npts);

/*******************************************************************************
 * \brief   Performs a 3D-R2C/C2R-FFT, on the GPU.
 * \author  Andreas Gloess, Ole Schuett
 ******************************************************************************/
void fft_r2c_gpu_fff(const double *zin, double *zout, const bool dir,
                     const int *npts);

// EOF
