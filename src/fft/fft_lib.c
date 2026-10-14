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
#include "fpga/fft_fpga.h"
#include "gpu/fft_gpu.h"

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
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
 * \brief Initialize the accelerated FFT library (FPGA+GPU).
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
#if defined(__FFT_FPGA)
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
#error                                                                         \
    "OFFLOAD and FPGA cannot be configured concurrently! Recompile with -D__NO_OFFLOAD_FFT."
  CPABORT("OFFLOAD and FPGA cannot be configured concurrently! Recompile with "
          "-D__NO_OFFLOAD_FFT.")
#endif
  const int stat = fft_fpga_initialize()
      assert(stat == 0 && "Initialization of FPGA failed!");
#endif
  int status = 1;
  switch (fft_lib_choice) {
  case FFT_LIB_FFTW:
    status =
        fft_fftw_init_lib(fftw_planning_flag, use_guru_interface, wisdom_file);
    printf("Using FFTW library.\n");
    break;
  case FFT_LIB_GPU:
    status = fft_init_acc_lib();
    printf("Using GPU-FFT library.\n");
    break;
  default:
    assert(0 && "Unknown FFT library.");
  }
  if (status != 0) {
    fft_finalize_lib(NULL);
  }
  return status;
}

/*******************************************************************************
 * \brief Finalize the accelerated FFT library.
 * \author Frederick Stein
 ******************************************************************************/
void fft_finalize_acc_lib() {
#if defined(__FFT_FPGA)
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
#error                                                                         \
    "OFFLOAD and FPGA cannot be configured concurrently! Recompile with -D__NO_OFFLOAD_FFT."
  CPABORT("OFFLOAD and FPGA cannot be configured concurrently! Recompile with "
          "-D__NO_OFFLOAD_FFT.")
#endif
  fft_fpga_final_();
#endif
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_FFT)
  fft_gpu_finalize();
#endif
}

/*******************************************************************************
 * \brief Finalize the FFT library (if not done externally).
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
    assert(0 && "Unknown FFT library.");
    return false;
  }
}

/*******************************************************************************
 * \brief Allocate buffer of type double.
 * \author Frederick Stein
 ******************************************************************************/
void fft_allocate_double(const int length, double **buffer) {
  switch (fft_lib_choice) {
  case FFT_LIB_FFTW:
    fft_fftw_allocate_double(length, buffer);
    break;
  case FFT_LIB_GPU:
    fft_gpu_allocate_double(length, buffer);
    break;
  default:
    assert(0 && "Unknown FFT library.");
  }
}

/*******************************************************************************
 * \brief Allocate buffer of type double complex.
 * \author Frederick Stein
 ******************************************************************************/
void fft_allocate_complex(const int length, double complex **buffer) {
  switch (fft_lib_choice) {
  case FFT_LIB_FFTW:
    fft_fftw_allocate_complex(imax(1, length), buffer);
    break;
  case FFT_LIB_GPU:
    fft_gpu_allocate_complex(imax(1, length), buffer);
    break;
  default:
    assert(0 && "Unknown FFT library.");
  }
}

/*******************************************************************************
 * \brief Allocate buffer of type double.
 * \author Frederick Stein
 ******************************************************************************/
void fft_free_double(double *buffer) {
  if (fft_lib_choice == FFT_LIB_FFTW) {
    fft_fftw_free_double(buffer);
  } else if (fft_lib_choice == FFT_LIB_GPU) {
    fft_gpu_free_double(buffer);
  } else {
    assert(0 && "Unknown FFT library.");
  }
}

/*******************************************************************************
 * \brief Allocate buffer of type double complex.
 * \author Frederick Stein
 ******************************************************************************/
void fft_free_complex(double complex *buffer) {
  if (fft_lib_choice == FFT_LIB_FFTW) {
    fft_fftw_free_complex(buffer);
  } else if (fft_lib_choice == FFT_LIB_GPU) {
    fft_gpu_free_complex(buffer);
  } else {
    assert(0 && "Unknown FFT library.");
  }
}

/*******************************************************************************
 * \brief Register a local 1D C2C FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_1d_local(const bool dir, const int fft_size,
                           const int number_of_ffts, const bool transpose_in,
                           const bool transpose_out,
                           const int leading_dimension_in,
                           const int leading_dimension_out,
                           double complex *grid_in, double complex *grid_out) {
  char routine_name[FFT_MAX_STRING_LENGTH + 1];
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH, "fft_register_1d_local");
  const int handle = fft_start_timer(routine_name);
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH,
           "fft_register_1d_local_%s_%s_%i_%i", dir ? "fw" : "bw",
           grid_in == grid_out ? "i" : "o", fft_size, number_of_ffts);
  const int handle2 = fft_start_timer(routine_name);
  switch (fft_lib_choice) {
  case FFT_LIB_GPU:
    fft_register_gpu_f(dir, fft_size, number_of_ffts, transpose_in,
                       transpose_out, leading_dimension_in,
                       leading_dimension_out, grid_in == grid_out);
    break;
  case FFT_LIB_FFTW:
    fft_fftw_register_1d_local(dir, fft_size, number_of_ffts, transpose_in,
                               transpose_out, leading_dimension_in,
                               leading_dimension_out, grid_in, grid_out);
    break;
  default:
    assert(0 && "Unknown FFT library.");
  }
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Register a local 1D R2C FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_1d_r2c_local(const int fft_size, const int number_of_ffts,
                               const bool transpose_in,
                               const bool transpose_out,
                               const int leading_dimension_in,
                               const int leading_dimension_out, double *grid_in,
                               double complex *grid_out) {
  char routine_name[FFT_MAX_STRING_LENGTH + 1];
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH, "fft_register_1d_r2c_local");
  const int handle = fft_start_timer(routine_name);
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH,
           "fft_register_1d_r2c_local_%s_%i_%i",
           grid_in == (double *)grid_out ? "i" : "o", fft_size, number_of_ffts);
  const int handle2 = fft_start_timer(routine_name);
  switch (fft_lib_choice) {
  case FFT_LIB_GPU:
    fft_register_r2c_gpu_f(true, fft_size, number_of_ffts, transpose_in,
                           transpose_out, leading_dimension_in,
                           leading_dimension_out,
                           grid_in == (double *)grid_out);
    break;
  case FFT_LIB_FFTW:
    fft_fftw_register_1d_r2c_local(fft_size, number_of_ffts, transpose_in,
                                   transpose_out, leading_dimension_in,
                                   leading_dimension_out, grid_in, grid_out);
    break;
  default:
    assert(0 && "Unknown FFT library.");
  }
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Register a local 1D C2R FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_1d_c2r_local(const int fft_size, const int number_of_ffts,
                               const bool transpose_in,
                               const bool transpose_out,
                               const int leading_dimension_in,
                               const int leading_dimension_out,
                               double complex *grid_in, double *grid_out) {
  char routine_name[FFT_MAX_STRING_LENGTH + 1];
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH, "fft_register_1d_c2r_local");
  const int handle = fft_start_timer(routine_name);
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH,
           "fft_register_1d_c2r_local_%s_%i_%i",
           (double *)grid_in == grid_out ? "i" : "o", fft_size, number_of_ffts);
  const int handle2 = fft_start_timer(routine_name);
  switch (fft_lib_choice) {
  case FFT_LIB_GPU:
    fft_register_r2c_gpu_f(false, fft_size, number_of_ffts, transpose_in,
                           transpose_out, leading_dimension_in,
                           leading_dimension_out,
                           grid_in == (double complex *)grid_out);
    break;
  case FFT_LIB_FFTW:
    fft_fftw_register_1d_c2r_local(fft_size, number_of_ffts, transpose_in,
                                   transpose_out, leading_dimension_in,
                                   leading_dimension_out, grid_in, grid_out);
    break;
  default:
    assert(0 && "Unknown FFT library.");
  }
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Register a local 2D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_2d_local(const bool dir, const int fft_size[2],
                           const int number_of_ffts, const bool transpose_in,
                           const bool transpose_out, double complex *grid_in,
                           double complex *grid_out) {
  char routine_name[FFT_MAX_STRING_LENGTH + 1];
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH, "fft_register_2d_local");
  const int handle = fft_start_timer(routine_name);
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH,
           "fft_register_2d_local_%s_%s_%i_%i_%i", dir ? "fw" : "bw",
           grid_in == grid_out ? "i" : "o", fft_size[0], fft_size[1],
           number_of_ffts);
  const int handle2 = fft_start_timer(routine_name);
  switch (fft_lib_choice) {
  case FFT_LIB_GPU:
    fft_register_gpu_ff(dir, fft_size, number_of_ffts, transpose_in,
                        transpose_out, grid_in == grid_out);
    break;
  case FFT_LIB_FFTW:
    fft_fftw_register_2d_local(dir, fft_size, number_of_ffts, transpose_in,
                               transpose_out, grid_in, grid_out);
    break;
  default:
    assert(0 && "Unknown FFT library.");
  }
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Register a local 2D R2C FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_2d_r2c_local(const int fft_size[2], const int number_of_ffts,
                               const bool transpose_in,
                               const bool transpose_out, double *grid_in,
                               double complex *grid_out) {
  char routine_name[FFT_MAX_STRING_LENGTH + 1];
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH, "fft_register_2d_r2c_local");
  const int handle = fft_start_timer(routine_name);
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH,
           "fft_register_2d_r2c_local_%s_%i_%i_%i",
           grid_in == (double *)grid_out ? "i" : "o", fft_size[0], fft_size[1],
           number_of_ffts);
  const int handle2 = fft_start_timer(routine_name);
  switch (fft_lib_choice) {
  case FFT_LIB_GPU:
    fft_register_r2c_gpu_ff(true, fft_size, number_of_ffts, transpose_in,
                            transpose_out, grid_in == (double *)grid_out);
    break;
  case FFT_LIB_FFTW:
    fft_fftw_register_2d_r2c_local(fft_size, number_of_ffts, transpose_in,
                                   transpose_out, grid_in, grid_out);
    break;
  default:
    assert(0 && "Unknown FFT library.");
  }
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Register a local 2D C2R FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_2d_c2r_local(const int fft_size[2], const int number_of_ffts,
                               const bool transpose_in,
                               const bool transpose_out,
                               double complex *grid_in, double *grid_out) {
  char routine_name[FFT_MAX_STRING_LENGTH + 1];
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH, "fft_register_2d_c2r_local");
  const int handle = fft_start_timer(routine_name);
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH,
           "fft_register_2d_c2r_local_%s_%i_%i_%i",
           (double *)grid_in == grid_out ? "i" : "o", fft_size[0], fft_size[1],
           number_of_ffts);
  const int handle2 = fft_start_timer(routine_name);
  switch (fft_lib_choice) {
  case FFT_LIB_GPU:
    fft_register_r2c_gpu_ff(false, fft_size, number_of_ffts, transpose_in,
                            transpose_out,
                            grid_in == (double complex *)grid_out);
    break;
  case FFT_LIB_FFTW:
    fft_fftw_register_2d_c2r_local(fft_size, number_of_ffts, transpose_in,
                                   transpose_out, grid_in, grid_out);
    break;
  default:
    assert(0 && "Unknown FFT library.");
  }
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Register a local 3D FFT.
 * \note fft_3d_bw_local(grid_gs, grid_rs, n) is the reverse to
 * fft_3d_rw_local(grid_rs, grid_gs, n) (ignoring normalization).
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_3d_local(const bool dir, const int fft_size[3],
                           double complex *grid_in, double complex *grid_out) {
  char routine_name[FFT_MAX_STRING_LENGTH + 1];
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH, "fft_register_3d_local");
  const int handle = fft_start_timer(routine_name);
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH,
           "fft_register_3d_local_%s_%s_%i_%i_%i", dir ? "fw" : "bw",
           grid_in == grid_out ? "i" : "o", fft_size[0], fft_size[1],
           fft_size[2]);
  const int handle2 = fft_start_timer(routine_name);
  switch (fft_lib_choice) {
  case FFT_LIB_GPU:
    fft_register_gpu_fff(dir, fft_size, grid_in == grid_out);
    break;
  case FFT_LIB_FFTW:
    fft_fftw_register_3d_local(dir, fft_size, grid_in, grid_out);
    break;
  default:
    assert(0 && "Unknown FFT library.");
  }
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Register a local 3D R2C FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_3d_r2c_local(const int fft_size[3], double *grid_in,
                               double complex *grid_out) {
  char routine_name[FFT_MAX_STRING_LENGTH + 1];
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH, "fft_register_3d_r2c_local");
  const int handle = fft_start_timer(routine_name);
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH,
           "fft_register_3d_r2c_local_%s_%i_%i_%i",
           grid_in == (double *)grid_out ? "i" : "o", fft_size[0], fft_size[1],
           fft_size[2]);
  const int handle2 = fft_start_timer(routine_name);
  switch (fft_lib_choice) {
  case FFT_LIB_GPU:
    fft_register_r2c_gpu_fff(true, fft_size, grid_in == (double *)grid_out);
    break;
  case FFT_LIB_FFTW:
    fft_fftw_register_3d_r2c_local(fft_size, grid_in, grid_out);
    break;
  default:
    assert(0 && "Unknown FFT library.");
  }
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Register a local 3D C2R FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_register_3d_c2r_local(const int fft_size[3], double complex *grid_in,
                               double *grid_out) {
  char routine_name[FFT_MAX_STRING_LENGTH + 1];
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH, "fft_register_3d_c2r_local");
  const int handle = fft_start_timer(routine_name);
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH,
           "fft_register_3d_c2r_local_%s_%i_%i_%i",
           (double *)grid_in == grid_out ? "i" : "o", fft_size[0], fft_size[1],
           fft_size[2]);
  const int handle2 = fft_start_timer(routine_name);
  switch (fft_lib_choice) {
  case FFT_LIB_GPU:
    fft_register_r2c_gpu_fff(false, fft_size,
                             grid_in == (double complex *)grid_out);
    break;
  case FFT_LIB_FFTW:
    fft_fftw_register_3d_c2r_local(fft_size, grid_in, grid_out);
    break;
  default:
    assert(0 && "Unknown FFT library.");
  }
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Perform a local 1D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_1d_local(const bool dir, const int fft_size, const int number_of_ffts,
                  const bool transpose_in, const bool transpose_out,
                  const int leading_dimension_in,
                  const int leading_dimension_out, double complex *grid_in,
                  double complex *grid_out) {
  char routine_name[FFT_MAX_STRING_LENGTH + 1];
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH, "fft_1d_local");
  const int handle = fft_start_timer(routine_name);
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH, "fft_1d_local_%s_%s_%i_%i",
           dir ? "fw" : "bw", grid_in == grid_out ? "i" : "o", fft_size,
           number_of_ffts);
  const int handle2 = fft_start_timer(routine_name);
  switch (fft_lib_choice) {
  case FFT_LIB_GPU:
    fft_gpu_f((const double *)grid_in, (double *)grid_out, dir, fft_size,
              number_of_ffts, transpose_in, transpose_out, leading_dimension_in,
              leading_dimension_out);
    break;
  case FFT_LIB_FFTW:
    fft_fftw_1d_local(dir, fft_size, number_of_ffts, transpose_in,
                      transpose_out, leading_dimension_in,
                      leading_dimension_out, grid_in, grid_out);
    break;
  default:
    assert(0 && "Unknown FFT library.");
  }
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Perform a local 1D R2C FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_1d_r2c_local(const int fft_size, const int number_of_ffts,
                      const bool transpose_in, const bool transpose_out,
                      const int leading_dimension_in,
                      const int leading_dimension_out, double *grid_in,
                      double complex *grid_out) {
  char routine_name[FFT_MAX_STRING_LENGTH + 1];
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH, "fft_1d_r2c_local");
  const int handle = fft_start_timer(routine_name);
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH, "fft_1d_r2c_local_%s_%i_%i",
           grid_in == (double *)grid_out ? "i" : "o", fft_size, number_of_ffts);
  const int handle2 = fft_start_timer(routine_name);
  switch (fft_lib_choice) {
  case FFT_LIB_GPU:
    fft_r2c_gpu_f((const double *)grid_in, (double *)grid_out, true, fft_size,
                  number_of_ffts, transpose_in, transpose_out,
                  leading_dimension_in, leading_dimension_out);
    break;
  case FFT_LIB_FFTW:
    fft_fftw_1d_r2c_local(fft_size, number_of_ffts, transpose_in, transpose_out,
                          leading_dimension_in, leading_dimension_out, grid_in,
                          grid_out);
    break;
  default:
    assert(0 && "Unknown FFT library.");
  }
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Perform a local 1D C2R FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_1d_c2r_local(const int fft_size, const int number_of_ffts,
                      const bool transpose_in, const bool transpose_out,
                      const int leading_dimension_in,
                      const int leading_dimension_out, double complex *grid_in,
                      double *grid_out) {
  char routine_name[FFT_MAX_STRING_LENGTH + 1];
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH, "fft_1d_c2r_local");
  const int handle = fft_start_timer(routine_name);
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH, "fft_1d_c2r_local_%s_%i_%i",
           (double *)grid_in == grid_out ? "i" : "o", fft_size, number_of_ffts);
  const int handle2 = fft_start_timer(routine_name);
  switch (fft_lib_choice) {
  case FFT_LIB_GPU:
    fft_r2c_gpu_f((const double *)grid_in, (double *)grid_out, false, fft_size,
                  number_of_ffts, transpose_in, transpose_out,
                  leading_dimension_in, leading_dimension_out);
    break;
  case FFT_LIB_FFTW:
    fft_fftw_1d_c2r_local(fft_size, number_of_ffts, transpose_in, transpose_out,
                          leading_dimension_in, leading_dimension_out, grid_in,
                          grid_out);
    break;
  default:
    assert(0 && "Unknown FFT library.");
  }
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Perform a local 2D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_2d_local(const bool dir, const int fft_size[2],
                  const int number_of_ffts, const bool transpose_in,
                  const bool transpose_out, double complex *grid_in,
                  double complex *grid_out) {
  char routine_name[FFT_MAX_STRING_LENGTH + 1];
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH, "fft_2d_local");
  const int handle = fft_start_timer(routine_name);
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH, "fft_2d_local_%s_%s_%i_%i_%i",
           dir ? "fw" : "bw", grid_in == grid_out ? "i" : "o", fft_size[0],
           fft_size[1], number_of_ffts);
  const int handle2 = fft_start_timer(routine_name);
  switch (fft_lib_choice) {
  case FFT_LIB_GPU:
    fft_gpu_ff((const double *)grid_in, (double *)grid_out, dir, fft_size,
               number_of_ffts, transpose_in, transpose_out);
    break;
  case FFT_LIB_FFTW:
    fft_fftw_2d_local(dir, fft_size, number_of_ffts, transpose_in,
                      transpose_out, grid_in, grid_out);
    break;
  default:
    assert(0 && "Unknown FFT library.");
  }
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Perform a local 2D R2C FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_2d_r2c_local(const int fft_size[2], const int number_of_ffts,
                      const bool transpose_in, const bool transpose_out,
                      double *grid_in, double complex *grid_out) {
  char routine_name[FFT_MAX_STRING_LENGTH + 1];
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH, "fft_2d_r2c_local");
  const int handle = fft_start_timer(routine_name);
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH, "fft_2d_r2c_local_%s_%i_%i_%i",
           grid_in == (double *)grid_out ? "i" : "o", fft_size[0], fft_size[1],
           number_of_ffts);
  const int handle2 = fft_start_timer(routine_name);
  switch (fft_lib_choice) {
  case FFT_LIB_GPU:
    fft_r2c_gpu_ff((const double *)grid_in, (double *)grid_out, true, fft_size,
                   number_of_ffts, transpose_in, transpose_out);
    break;
  case FFT_LIB_FFTW:
    fft_fftw_2d_r2c_local(fft_size, number_of_ffts, transpose_in, transpose_out,
                          grid_in, grid_out);
    break;
  default:
    assert(0 && "Unknown FFT library.");
  }
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Perform a local 2D C2R FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_2d_c2r_local(const int fft_size[2], const int number_of_ffts,
                      const bool transpose_in, const bool transpose_out,
                      double complex *grid_in, double *grid_out) {
  char routine_name[FFT_MAX_STRING_LENGTH + 1];
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH, "fft_2d_c2r_local");
  const int handle = fft_start_timer(routine_name);
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH, "fft_2d_c2r_local_%s_%i_%i_%i",
           (double *)grid_in == grid_out ? "i" : "o", fft_size[0], fft_size[1],
           number_of_ffts);
  const int handle2 = fft_start_timer(routine_name);
  switch (fft_lib_choice) {
  case FFT_LIB_GPU:
    fft_r2c_gpu_ff((const double *)grid_in, (double *)grid_out, false, fft_size,
                   number_of_ffts, transpose_in, transpose_out);
    break;
  case FFT_LIB_FFTW:
    fft_fftw_2d_c2r_local(fft_size, number_of_ffts, transpose_in, transpose_out,
                          grid_in, grid_out);
    break;
  default:
    assert(0 && "Unknown FFT library.");
  }
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Perform a local 3D FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_3d_local(const bool dir, const int fft_size[3],
                  double complex *grid_in, double complex *grid_out) {
  char routine_name[FFT_MAX_STRING_LENGTH + 1];
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH, "fft_3d_local");
  const int handle = fft_start_timer(routine_name);
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH, "fft_3d_local_%s_%s_%i_%i_%i",
           dir ? "fw" : "bw", grid_in == grid_out ? "i" : "o", fft_size[0],
           fft_size[1], fft_size[2]);
  const int handle2 = fft_start_timer(routine_name);
#if defined(__FFT_FPGA)
  if (fft_fpga_check_bitstream_(get_data_dir(), fft_size)) {
    const int number_of_elements = product3(fft_size);
#if (__FFT_FPGA_SP && __FFT_FPGA)
    float complex *grid_sp = calloc(number_of_elements, sizeof(float complex));
    for (int i = 0; i < number_of_elements; i++)
      grid_sp[i] = (float complex)grid_in[i];
    fft_fpga_fft3d_sp_(dir ? 1 : -1, fft_size, grid_sp);
    for (int i = 0; i < number_of_elements; i++)
      grid_out[i] = (double complex)grid_sp[i];
#else
    memcpy(grid_out, grid_in, number_of_elements * sizeof(double complex));
    fft_fpga_fft3d_dp_(dir ? 1 : -1, fft_size, grid_out);
#endif
  } else {
#endif
    switch (fft_lib_choice) {
    case FFT_LIB_GPU:
      fft_gpu_fff((const double *)grid_in, (double *)grid_out, dir, fft_size);
      break;
    case FFT_LIB_FFTW:
      fft_fftw_3d_local(dir, fft_size, grid_in, grid_out);
      break;
    default:
      assert(0 && "Unknown FFT library.");
    }
#if defined(__FFT_FPGA)
  }
#endif
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Perform a local 3D R2C FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_3d_r2c_local(const int fft_size[3], double *grid_in,
                      double complex *grid_out) {
  char routine_name[FFT_MAX_STRING_LENGTH + 1];
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH, "fft_3d_r2c_local");
  const int handle = fft_start_timer(routine_name);
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH, "fft_3d_r2c_local_%s_%i_%i_%i",
           grid_in == (double *)grid_out ? "i" : "o", fft_size[0], fft_size[1],
           fft_size[2]);
  const int handle2 = fft_start_timer(routine_name);
#if defined(__FFT_FPGA)
  if (fft_fpga_check_bitstream_(get_data_dir(), fft_size)) {
    const int number_of_elements = product3(fft_size);
#if (__FFT_FPGA_SP && __FFT_FPGA)
    float complex *grid_sp = calloc(number_of_elements, sizeof(float complex));
    for (int i = 0; i < number_of_elements; i++)
      grid_sp[i] = CMPLXF((float)grid_in[i], 0.0f);
    fft_fpga_fft3d_sp_(1, fft_size, grid_sp);
    for (int i = 0; i < number_of_elements; i++)
      grid_out[i] = (double complex)grid_sp[i];
#else
    for (int i = 0; i < number_of_elements; i++)
      grid_out[i] = CMPLX(grid_in[i], 0.0);
    fft_fpga_fft3d_dp_(1, fft_size, grid_out);
#endif
  } else {
#endif
    switch (fft_lib_choice) {
    case FFT_LIB_GPU:
      fft_r2c_gpu_fff((const double *)grid_in, (double *)grid_out, true,
                      fft_size);
      break;
    case FFT_LIB_FFTW:
      fft_fftw_3d_r2c_local(fft_size, grid_in, grid_out);
      break;
    default:
      assert(0 && "Unknown FFT library.");
    }
#if defined(__FFT_FPGA)
  }
#endif
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

/*******************************************************************************
 * \brief Perform a local 3D C2R FFT.
 * \author Frederick Stein
 ******************************************************************************/
void fft_3d_c2r_local(const int fft_size[3], double complex *grid_in,
                      double *grid_out) {
  char routine_name[FFT_MAX_STRING_LENGTH + 1];
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH, "fft_3d_c2r_local");
  const int handle = fft_start_timer(routine_name);
  memset(routine_name, '\0', FFT_MAX_STRING_LENGTH + 1);
  snprintf(routine_name, FFT_MAX_STRING_LENGTH, "fft_3d_c2r_local_%s_%i_%i_%i",
           (double *)grid_in == grid_out ? "i" : "o", fft_size[0], fft_size[1],
           fft_size[2]);
  const int handle2 = fft_start_timer(routine_name);
#if defined(__FFT_FPGA)
  if (fft_fpga_check_bitstream_(get_data_dir(), fft_size)) {
    const int number_of_elements =
        (fft_size[0] / 2 + 1) * fft_size[1] * fft_size[2];
#if (__FFT_FPGA_SP && __FFT_FPGA)
    float complex *grid_sp = calloc(number_of_elements, sizeof(float complex));
    for (int i = 0; i < number_of_elements; i++)
      grid_sp[i] = (float complex)grid_in[i];
    fft_fpga_fft3d_sp_(-1, fft_size, grid_sp);
    for (int i = 0; i < number_of_elements; i++)
      grid_out[i] = (double)crealf(grid_sp[i]);
#else
    memcpy((double complex *)grid_out, grid_in,
           number_of_elements * sizeof(double));
    fft_fpga_fft3d_dp_(-1, fft_size, grid_out);
#endif
  } else {
#endif
    switch (fft_lib_choice) {
    case FFT_LIB_GPU:
      fft_r2c_gpu_fff((const double *)grid_in, (double *)grid_out, false,
                      fft_size);
      break;
    case FFT_LIB_FFTW:
      fft_fftw_3d_c2r_local(fft_size, grid_in, grid_out);
      break;
    default:
      assert(0 && "Unknown FFT library.");
    }
#if defined(__FFT_FPGA)
  }
#endif
  fft_stop_timer(handle2);
  fft_stop_timer(handle);
}

// EOF
