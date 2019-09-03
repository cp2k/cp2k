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
void pw_fpga_fft3d_sp_(int direction, int N[3], cmplx *din);

// Double precision FFT3d procedure
void pw_fpga_fft3d_dp_(int direction, int N[3], cmplx *din);

// Check fpga bitstream present in directory
int pw_fpga_check_bitstream_(char *data_path, int N[3]);

// host variables
static cl_platform_id platform = NULL;
static cl_device_id *devices;
static cl_device_id device = NULL;
static cl_context context = NULL;
static cl_program program = NULL;

static cl_command_queue queue1 = NULL, queue2 = NULL, queue3 = NULL;
static cl_command_queue queue4 = NULL, queue5 = NULL, queue6 = NULL;

#endif
#endif
