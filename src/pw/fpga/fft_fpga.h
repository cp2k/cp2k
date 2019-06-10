/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2019  CP2K developers group                         *
 *****************************************************************************/

#ifndef FFT_FPGA_H
#define FFT_FPGA_H
/******************************************************************************
 *  Authors: Arjun Ramaswami
 *****************************************************************************/
#if defined ( __PW_FPGA )

typedef struct {
  double x;
  double y;
} double2;

typedef struct {
  float x;
  float y;
} float2;

#ifdef __PW_FPGA_SP
    typedef float2 cmplx;
#else
    typedef double2 cmplx;
#endif

// Initialize FPGA
int pw_fpga_initialize_();

// Finalize FPGA
void pw_fpga_final_();

// Single precision FFT3d procedure
void pw_fpga_fft3d_sp_(int data_path_len, char *data_path, int direction, int N[3], cmplx *din);

// Double precision FFT3d procedure
void pw_fpga_fft3d_dp_(int data_path_len, char *data_path, int direction, int N[3], cmplx *din);

// Check fpga bitstream present in directory
bool pw_fpga_check_bitstream_(int N[3]);

#endif
#endif
