/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include "fft_lib.h"
#include "../offload/offload_runtime.h"
#include "fft_lib_fftw.h"
#include "fft_timer.h"
#include "fft_utils.h"
#include "gpu/fft_gpu.h"

#include "../base/cp_assert.h"
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
const int fft_lib_default = FFT_LIB_GPU;
#elif defined(__FFTW3)
const int fft_lib_default = FFT_LIB_FFTW;
#else
#error "The FFT backend needs at least the FFTW3 backend."
const int fft_lib_default = -1; // No backend available
#endif
fft_lib fft_lib_choice = fft_lib_default;
bool fft_lib_initialized = false;

/*******************************************************************************
 * \brief Start an FFT timer from a printf-style label.
 ******************************************************************************/
static int fft_timer_start(const char *format, ...) {
  char routine_name[FFT_MAX_STRING_LENGTH + 1];
  va_list arguments;
  va_start(arguments, format);
  vsnprintf(routine_name, sizeof(routine_name), format, arguments);
  va_end(arguments);
  return fft_start_timer(routine_name);
}

/*******************************************************************************
 * \brief Forward a void FFT call to the selected backend.
 ******************************************************************************/
#define FFT_DISPATCH_VOID(gpu_call, fftw_call)                                 \
  do {                                                                         \
    switch (fft_lib_choice) {                                                  \
    case FFT_LIB_GPU:                                                          \
      gpu_call;                                                                \
      break;                                                                   \
    case FFT_LIB_FFTW:                                                         \
      fftw_call;                                                               \
      break;                                                                   \
    default:                                                                   \
      cp_assert(0 && "Unknown FFT library.");                                  \
    }                                                                          \
  } while (0)

/*******************************************************************************
 * \brief Initialize the accelerated FFT library.
 * \author Frederick Stein
 ******************************************************************************/
int fft_init_acc_lib() {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
  fft_gpu_init();
  return 0;
#else
  return 1;
#endif
}

/*******************************************************************************
 * \brief Initialize the FFT library (if not done externally).
 * \param lib FFT backend to use.
 * \param fftw_planning_flag FFTW planning mode.
 * \param use_guru_interface Whether the FFTW Guru interface may be used.
 * \param wisdom_file Path to the FFTW wisdom file.
 * \author Frederick Stein
 ******************************************************************************/
int fft_init_lib(const fft_lib lib, const int fftw_planning_flag,
                 const bool use_guru_interface, const char *wisdom_file) {
  if (fft_lib_initialized) {
    return 0;
  }
  fft_lib_initialized = true;
  if (lib != FFT_LIB_DEFAULT)
    fft_lib_choice = lib;
  int status = 1;
  switch (fft_lib_choice) {
  case FFT_LIB_FFTW:
    status =
        fft_fftw_init_lib(fftw_planning_flag, use_guru_interface, wisdom_file);
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
    // Distributed FFTs always use the GPU backend, so it must be initialized
    // even when FFTW3 is selected for the local FFTs.
    if (status == 0) {
      status = fft_init_acc_lib();
    }
#endif
    break;
  case FFT_LIB_GPU:
    status = fft_init_acc_lib();
    break;
  default:
    cp_assert(0 && "Unknown FFT library.");
  }
  if (status != 0) {
    fft_finalize_lib(NULL);
    return status;
  }
  fft_init_timer();
  return status;
}

/*******************************************************************************
 * \brief Finalize the accelerated FFT library.
 * \author Frederick Stein
 ******************************************************************************/
void fft_finalize_acc_lib() {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
  fft_gpu_finalize();
#endif
}

/*******************************************************************************
 * \brief Finalize the FFT library (if not done externally).
 * \param wisdom_file Path to the FFTW wisdom file.
 * \author Frederick Stein
 ******************************************************************************/
void fft_finalize_lib(const char *wisdom_file) {
  fft_fftw_finalize_lib(wisdom_file);
  fft_finalize_acc_lib();
  fft_lib_initialized = false;
}

/*******************************************************************************
 * \brief Get the default library (GPU if offloading was enabled, else FFTW3).
 * \author Frederick Stein
 ******************************************************************************/
int fft_lib_default_library() { return fft_lib_default; }

/*******************************************************************************
 * \brief Whether a compound MPI implementation is available.
 * \author Frederick Stein
 ******************************************************************************/
int fft_lib_backend_in_use() { return fft_lib_choice; }

/*******************************************************************************
 * \brief Whether a compound MPI implementation is available.
 * \author Frederick Stein
 ******************************************************************************/
bool fft_lib_has_guru_interface() {
  switch (fft_lib_choice) {
  case FFT_LIB_FFTW:
    return fft_fftw_lib_has_guru_interface();
  case FFT_LIB_GPU:
    return false;
  default:
    cp_assert(0 && "Unknown FFT library.");
    return false;
  }
}

/*******************************************************************************
 * \brief Allocate buffer of type double.
 * \param length Number of elements.
 * \param buffer Buffer pointer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_allocate_double(const int length, double **buffer) {
  FFT_DISPATCH_VOID((fft_gpu_allocate_double(length, buffer)),
                    (fft_fftw_allocate_double(length, buffer)));
}

/*******************************************************************************
 * \brief Allocate buffer of type double complex.
 * \param length Number of elements.
 * \param buffer Buffer pointer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_allocate_complex(const int length, double complex **buffer) {
  FFT_DISPATCH_VOID((fft_gpu_allocate_complex(imax(1, length), buffer)),
                    (fft_fftw_allocate_complex(imax(1, length), buffer)));
}

/*******************************************************************************
 * \brief Allocate buffer of type double.
 * \param buffer Buffer pointer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_free_double(double *buffer) {
  if (fft_lib_choice == FFT_LIB_FFTW) {
    fft_fftw_free_double(buffer);
  } else if (fft_lib_choice == FFT_LIB_GPU) {
    fft_gpu_free_double(buffer);
  } else {
    cp_assert(0 && "Unknown FFT library.");
  }
}

/*******************************************************************************
 * \brief Allocate buffer of type double complex.
 * \param buffer Buffer pointer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_free_complex(double complex *buffer) {
  if (fft_lib_choice == FFT_LIB_FFTW) {
    fft_fftw_free_complex(buffer);
  } else if (fft_lib_choice == FFT_LIB_GPU) {
    fft_gpu_free_complex(buffer);
  } else {
    cp_assert(0 && "Unknown FFT library.");
  }
}

/*******************************************************************************
 * \brief Register a local 1D C2C FFT.
 * \param dir FFT direction (forward or backward).
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param leading_dimension_in Leading dimension of the input array.
 * \param leading_dimension_out Leading dimension of the output array.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_1d_local(const bool dir, const int fft_size,
                           const int number_of_ffts, const bool transpose_in,
                           const bool transpose_out,
                           const int leading_dimension_in,
                           const int leading_dimension_out,
                           double complex *grid_in, double complex *grid_out) {
  const int handle = fft_timer_start("register_1d_local");
  const int handle2 = fft_timer_start(
      "register_1d_l_%s_%s_%i_%i", dir ? "fw" : "bw",
      grid_in == grid_out ? "i" : "o", fft_size, number_of_ffts);
  FFT_DISPATCH_VOID(
      (fft_gpu_register_1d_local(dir, fft_size, number_of_ffts, transpose_in,
                                 transpose_out, leading_dimension_in,
                                 leading_dimension_out, grid_in == grid_out)),
      (fft_fftw_register_1d_local(dir, fft_size, number_of_ffts, transpose_in,
                                  transpose_out, leading_dimension_in,
                                  leading_dimension_out, grid_in, grid_out)));
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Register a local 1D R2C FFT.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param leading_dimension_in Leading dimension of the input array.
 * \param leading_dimension_out Leading dimension of the output array.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_1d_r2c_local(const int fft_size, const int number_of_ffts,
                               const bool transpose_in,
                               const bool transpose_out,
                               const int leading_dimension_in,
                               const int leading_dimension_out, double *grid_in,
                               double complex *grid_out) {
  const int handle = fft_timer_start("register_1d_r2c_local");
  const int handle2 = fft_timer_start("register_1d_r2c_l_%s_%i_%i",
                                      grid_in == (double *)grid_out ? "i" : "o",
                                      fft_size, number_of_ffts);
  FFT_DISPATCH_VOID(
      (fft_gpu_register_1d_r2c_local(
          true, fft_size, number_of_ffts, transpose_in, transpose_out,
          leading_dimension_in, leading_dimension_out,
          grid_in == (double *)grid_out)),
      (fft_fftw_register_1d_r2c_local(
          fft_size, number_of_ffts, transpose_in, transpose_out,
          leading_dimension_in, leading_dimension_out, grid_in, grid_out)));
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Register a local 1D C2R FFT.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param leading_dimension_in Leading dimension of the input array.
 * \param leading_dimension_out Leading dimension of the output array.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_1d_c2r_local(const int fft_size, const int number_of_ffts,
                               const bool transpose_in,
                               const bool transpose_out,
                               const int leading_dimension_in,
                               const int leading_dimension_out,
                               double complex *grid_in, double *grid_out) {
  const int handle = fft_timer_start("register_1d_c2r_local");
  const int handle2 = fft_timer_start("register_1d_c2r_l_%s_%i_%i",
                                      (double *)grid_in == grid_out ? "i" : "o",
                                      fft_size, number_of_ffts);
  FFT_DISPATCH_VOID(
      (fft_gpu_register_1d_r2c_local(
          false, fft_size, number_of_ffts, transpose_in, transpose_out,
          leading_dimension_in, leading_dimension_out,
          grid_in == (double complex *)grid_out)),
      (fft_fftw_register_1d_c2r_local(
          fft_size, number_of_ffts, transpose_in, transpose_out,
          leading_dimension_in, leading_dimension_out, grid_in, grid_out)));
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Register a local 2D FFT.
 * \param dir FFT direction (forward or backward).
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_2d_local(const bool dir, const int fft_size[2],
                           const int number_of_ffts, const bool transpose_in,
                           const bool transpose_out, double complex *grid_in,
                           double complex *grid_out) {
  const int handle = fft_timer_start("register_2d_local");
  const int handle2 =
      fft_timer_start("register_2d_l_%s_%s_%i_%i_%i", dir ? "fw" : "bw",
                      grid_in == grid_out ? "i" : "o", fft_size[0], fft_size[1],
                      number_of_ffts);
  FFT_DISPATCH_VOID(
      (fft_gpu_register_2d_local(dir, fft_size, number_of_ffts, transpose_in,
                                 transpose_out, grid_in == grid_out)),
      (fft_fftw_register_2d_local(dir, fft_size, number_of_ffts, transpose_in,
                                  transpose_out, grid_in, grid_out)));
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Register a local 2D R2C FFT.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_2d_r2c_local(const int fft_size[2], const int number_of_ffts,
                               const bool transpose_in,
                               const bool transpose_out, double *grid_in,
                               double complex *grid_out) {
  const int handle = fft_timer_start("register_2d_r2c_local");
  const int handle2 = fft_timer_start("register_2d_r2c_l_%s_%i_%i_%i",
                                      grid_in == (double *)grid_out ? "i" : "o",
                                      fft_size[0], fft_size[1], number_of_ffts);
  FFT_DISPATCH_VOID(
      (fft_gpu_register_2d_r2c_local(true, fft_size, number_of_ffts,
                                     transpose_in, transpose_out,
                                     grid_in == (double *)grid_out)),
      (fft_fftw_register_2d_r2c_local(fft_size, number_of_ffts, transpose_in,
                                      transpose_out, grid_in, grid_out)));
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Register a local 2D C2R FFT.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_2d_c2r_local(const int fft_size[2], const int number_of_ffts,
                               const bool transpose_in,
                               const bool transpose_out,
                               double complex *grid_in, double *grid_out) {
  const int handle = fft_timer_start("register_2d_c2r_local");
  const int handle2 = fft_timer_start("register_2d_c2r_l_%s_%i_%i_%i",
                                      (double *)grid_in == grid_out ? "i" : "o",
                                      fft_size[0], fft_size[1], number_of_ffts);
  FFT_DISPATCH_VOID(
      (fft_gpu_register_2d_r2c_local(false, fft_size, number_of_ffts,
                                     transpose_in, transpose_out,
                                     grid_in == (double complex *)grid_out)),
      (fft_fftw_register_2d_c2r_local(fft_size, number_of_ffts, transpose_in,
                                      transpose_out, grid_in, grid_out)));
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Register a local 3D FFT.
 * \note fft_3d_bw_local(grid_gs, grid_rs, n) is the reverse to
 * fft_3d_rw_local(grid_rs, grid_gs, n) (ignoring normalization).
 * \param dir FFT direction (forward or backward).
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_3d_local(const bool dir, const int fft_size[3],
                           const int number_of_ffts, const bool transpose_in,
                           const bool transpose_out, double complex *grid_in,
                           double complex *grid_out) {
  const int handle = fft_timer_start("register_3d_local");
  const int handle2 = fft_timer_start(
      "register_3d_l_%s_%s_%i_%i_%i", dir ? "fw" : "bw",
      grid_in == grid_out ? "i" : "o", fft_size[0], fft_size[1], fft_size[2]);
  FFT_DISPATCH_VOID(
      (fft_gpu_register_3d_local(dir, fft_size, grid_in == grid_out,
                                 number_of_ffts, transpose_in, transpose_out)),
      (fft_fftw_register_3d_local(dir, fft_size, number_of_ffts, transpose_in,
                                  transpose_out, grid_in, grid_out)));
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Register a local 3D R2C FFT.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_3d_r2c_local(const int fft_size[3], const int number_of_ffts,
                               const bool transpose_in,
                               const bool transpose_out, double *grid_in,
                               double complex *grid_out) {
  const int handle = fft_timer_start("register_3d_r2c_local");
  const int handle2 = fft_timer_start("register_3d_r2c_l_%s_%i_%i_%i",
                                      grid_in == (double *)grid_out ? "i" : "o",
                                      fft_size[0], fft_size[1], fft_size[2]);
  FFT_DISPATCH_VOID(
      (fft_gpu_register_3d_r2c_local(
          true, fft_size, grid_in == (double *)grid_out, number_of_ffts,
          transpose_in, transpose_out)),
      (fft_fftw_register_3d_r2c_local(fft_size, number_of_ffts, transpose_in,
                                      transpose_out, grid_in, grid_out)));
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Register a local 3D C2R FFT.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_3d_c2r_local(const int fft_size[3], const int number_of_ffts,
                               const bool transpose_in,
                               const bool transpose_out,
                               double complex *grid_in, double *grid_out) {
  const int handle = fft_timer_start("register_3d_c2r_local");
  const int handle2 = fft_timer_start("register_3d_c2r_l_%s_%i_%i_%i",
                                      (double *)grid_in == grid_out ? "i" : "o",
                                      fft_size[0], fft_size[1], fft_size[2]);
  FFT_DISPATCH_VOID(
      (fft_gpu_register_3d_r2c_local(
          false, fft_size, grid_in == (double complex *)grid_out,
          number_of_ffts, transpose_in, transpose_out)),
      (fft_fftw_register_3d_c2r_local(fft_size, number_of_ffts, transpose_in,
                                      transpose_out, grid_in, grid_out)));
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Perform a local 1D FFT.
 * \param dir FFT direction (forward or backward).
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param leading_dimension_in Leading dimension of the input array.
 * \param leading_dimension_out Leading dimension of the output array.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_1d_local(const bool dir, const int fft_size, const int number_of_ffts,
                  const bool transpose_in, const bool transpose_out,
                  const int leading_dimension_in,
                  const int leading_dimension_out, double complex *grid_in,
                  double complex *grid_out) {
  const int handle = fft_timer_start("1d_local");
  const int handle2 = fft_timer_start("1d_l_%s_%s_%i_%i", dir ? "fw" : "bw",
                                      grid_in == grid_out ? "i" : "o", fft_size,
                                      number_of_ffts);
  FFT_DISPATCH_VOID(
      (fft_gpu_1d_local((const double *)grid_in, (double *)grid_out, dir,
                        fft_size, number_of_ffts, transpose_in, transpose_out,
                        leading_dimension_in, leading_dimension_out)),
      (fft_fftw_1d_local(dir, fft_size, number_of_ffts, transpose_in,
                         transpose_out, leading_dimension_in,
                         leading_dimension_out, grid_in, grid_out)));
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Perform a local 1D R2C FFT.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param leading_dimension_in Leading dimension of the input array.
 * \param leading_dimension_out Leading dimension of the output array.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_1d_r2c_local(const int fft_size, const int number_of_ffts,
                      const bool transpose_in, const bool transpose_out,
                      const int leading_dimension_in,
                      const int leading_dimension_out, double *grid_in,
                      double complex *grid_out) {
  const int handle = fft_timer_start("1d_r2c_local");
  const int handle2 = fft_timer_start("1d_r2c_l_%s_%i_%i",
                                      grid_in == (double *)grid_out ? "i" : "o",
                                      fft_size, number_of_ffts);
  FFT_DISPATCH_VOID(
      (fft_gpu_1d_r2c_local((const double *)grid_in, (double *)grid_out, true,
                            fft_size, number_of_ffts, transpose_in,
                            transpose_out, leading_dimension_in,
                            leading_dimension_out)),
      (fft_fftw_1d_r2c_local(fft_size, number_of_ffts, transpose_in,
                             transpose_out, leading_dimension_in,
                             leading_dimension_out, grid_in, grid_out)));
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Perform a local 1D C2R FFT.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param leading_dimension_in Leading dimension of the input array.
 * \param leading_dimension_out Leading dimension of the output array.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_1d_c2r_local(const int fft_size, const int number_of_ffts,
                      const bool transpose_in, const bool transpose_out,
                      const int leading_dimension_in,
                      const int leading_dimension_out, double complex *grid_in,
                      double *grid_out) {
  const int handle = fft_timer_start("1d_c2r_local");
  const int handle2 = fft_timer_start("1d_c2r_l_%s_%i_%i",
                                      (double *)grid_in == grid_out ? "i" : "o",
                                      fft_size, number_of_ffts);
  FFT_DISPATCH_VOID(
      (fft_gpu_1d_r2c_local((const double *)grid_in, (double *)grid_out, false,
                            fft_size, number_of_ffts, transpose_in,
                            transpose_out, leading_dimension_in,
                            leading_dimension_out)),
      (fft_fftw_1d_c2r_local(fft_size, number_of_ffts, transpose_in,
                             transpose_out, leading_dimension_in,
                             leading_dimension_out, grid_in, grid_out)));
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Perform a local 2D FFT.
 * \param dir FFT direction (forward or backward).
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_2d_local(const bool dir, const int fft_size[2],
                  const int number_of_ffts, const bool transpose_in,
                  const bool transpose_out, double complex *grid_in,
                  double complex *grid_out) {
  const int handle = fft_timer_start("2d_local");
  const int handle2 = fft_timer_start("2d_l_%s_%s_%i_%i_%i", dir ? "fw" : "bw",
                                      grid_in == grid_out ? "i" : "o",
                                      fft_size[0], fft_size[1], number_of_ffts);
  FFT_DISPATCH_VOID(
      (fft_gpu_2d_local((const double *)grid_in, (double *)grid_out, dir,
                        fft_size, number_of_ffts, transpose_in, transpose_out)),
      (fft_fftw_2d_local(dir, fft_size, number_of_ffts, transpose_in,
                         transpose_out, grid_in, grid_out)));
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Perform a local 2D R2C FFT.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_2d_r2c_local(const int fft_size[2], const int number_of_ffts,
                      const bool transpose_in, const bool transpose_out,
                      double *grid_in, double complex *grid_out) {
  const int handle = fft_timer_start("2d_r2c_local");
  const int handle2 = fft_timer_start("2d_r2c_l_%s_%i_%i_%i",
                                      grid_in == (double *)grid_out ? "i" : "o",
                                      fft_size[0], fft_size[1], number_of_ffts);
  FFT_DISPATCH_VOID(
      (fft_gpu_2d_r2c_local((const double *)grid_in, (double *)grid_out, true,
                            fft_size, number_of_ffts, transpose_in,
                            transpose_out)),
      (fft_fftw_2d_r2c_local(fft_size, number_of_ffts, transpose_in,
                             transpose_out, grid_in, grid_out)));
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Perform a local 2D C2R FFT.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_2d_c2r_local(const int fft_size[2], const int number_of_ffts,
                      const bool transpose_in, const bool transpose_out,
                      double complex *grid_in, double *grid_out) {
  const int handle = fft_timer_start("2d_c2r_local");
  const int handle2 = fft_timer_start("2d_c2r_l_%s_%i_%i_%i",
                                      (double *)grid_in == grid_out ? "i" : "o",
                                      fft_size[0], fft_size[1], number_of_ffts);
  FFT_DISPATCH_VOID(
      (fft_gpu_2d_r2c_local((const double *)grid_in, (double *)grid_out, false,
                            fft_size, number_of_ffts, transpose_in,
                            transpose_out)),
      (fft_fftw_2d_c2r_local(fft_size, number_of_ffts, transpose_in,
                             transpose_out, grid_in, grid_out)));
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Perform a local 3D FFT.
 * \param dir FFT direction (forward or backward).
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_3d_local(const bool dir, const int fft_size[3],
                  const int number_of_ffts, const bool transpose_in,
                  const bool transpose_out, double complex *grid_in,
                  double complex *grid_out) {
  const int handle = fft_timer_start("3d_local");
  const int handle2 = fft_timer_start("3d_l_%s_%s_%i_%i_%i", dir ? "fw" : "bw",
                                      grid_in == grid_out ? "i" : "o",
                                      fft_size[0], fft_size[1], fft_size[2]);
  FFT_DISPATCH_VOID(
      (fft_gpu_3d_local((const double *)grid_in, (double *)grid_out, dir,
                        fft_size, number_of_ffts, transpose_in, transpose_out)),
      (fft_fftw_3d_local(dir, fft_size, number_of_ffts, transpose_in,
                         transpose_out, grid_in, grid_out)));
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Perform a local 3D R2C FFT.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_3d_r2c_local(const int fft_size[3], const int number_of_ffts,
                      const bool transpose_in, const bool transpose_out,
                      double *grid_in, double complex *grid_out) {
  const int handle = fft_timer_start("3d_r2c_local");
  const int handle2 = fft_timer_start("3d_r2c_l_%s_%i_%i_%i",
                                      grid_in == (double *)grid_out ? "i" : "o",
                                      fft_size[0], fft_size[1], fft_size[2]);
  FFT_DISPATCH_VOID(
      (fft_gpu_3d_r2c_local((const double *)grid_in, (double *)grid_out, true,
                            fft_size, number_of_ffts, transpose_in,
                            transpose_out)),
      (fft_fftw_3d_r2c_local(fft_size, number_of_ffts, transpose_in,
                             transpose_out, grid_in, grid_out)));
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Perform a local 3D C2R FFT.
 * \param fft_size Size of the transform in each dimension.
 * \param number_of_ffts Number of independent (batched) transforms.
 * \param transpose_in Whether the input array uses the transposed layout.
 * \param transpose_out Whether the output array uses the transposed layout.
 * \param grid_in Input grid buffer.
 * \param grid_out Output grid buffer.
 * \author Frederick Stein
 ******************************************************************************/
void fft_3d_c2r_local(const int fft_size[3], const int number_of_ffts,
                      const bool transpose_in, const bool transpose_out,
                      double complex *grid_in, double *grid_out) {
  const int handle = fft_timer_start("3d_c2r_local");
  const int handle2 = fft_timer_start("3d_c2r_l_%s_%i_%i_%i",
                                      (double *)grid_in == grid_out ? "i" : "o",
                                      fft_size[0], fft_size[1], fft_size[2]);
  FFT_DISPATCH_VOID(
      (fft_gpu_3d_r2c_local((const double *)grid_in, (double *)grid_out, false,
                            fft_size, number_of_ffts, transpose_in,
                            transpose_out)),
      (fft_fftw_3d_c2r_local(fft_size, number_of_ffts, transpose_in,
                             transpose_out, grid_in, grid_out)));
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

// EOF
