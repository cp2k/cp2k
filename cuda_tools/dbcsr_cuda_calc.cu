/******************************************************************************
*  CP2K: A general program to perform molecular dynamics simulations
*  Copyright (C) 2000 - 2013  Urban Borstnik and the CP2K developers group
*****************************************************************************/
#include <cuda_runtime.h>
#include <stdio.h>
#include <sm_11_atomic_functions.h>
#include "dbcsr_cuda.h"
#include "dbcsr_kernels/stack_mm_r.h"
#include "dbcsr_kernels/stack_mm_d.h"
#include "dbcsr_kernels/stack_mm_c.h"
#include "dbcsr_kernels/stack_mm_z.h"
#include "dbcsr_kernels/stack_mm_mnk_d.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_5_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_5_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_5_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_5_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_5_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_5_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_8_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_8_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_8_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_8_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_8_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_8_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_13_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_13_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_13_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_13_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_13_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_13_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_16_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_16_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_16_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_16_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_16_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_16_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_23_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_23_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_23_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_23_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_23_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_23_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_26_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_26_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_26_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_26_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_26_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_5_26_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_5_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_5_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_5_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_5_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_5_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_5_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_8_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_8_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_8_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_8_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_8_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_8_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_13_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_13_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_13_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_13_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_13_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_13_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_16_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_16_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_16_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_16_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_16_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_16_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_23_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_23_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_23_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_23_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_23_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_23_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_26_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_26_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_26_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_26_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_26_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_8_26_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_5_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_5_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_5_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_5_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_5_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_5_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_8_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_8_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_8_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_8_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_8_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_8_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_13_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_13_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_13_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_13_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_13_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_13_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_16_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_16_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_16_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_16_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_16_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_16_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_23_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_23_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_23_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_23_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_23_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_23_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_26_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_26_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_26_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_26_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_26_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_13_26_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_5_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_5_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_5_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_5_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_5_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_5_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_8_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_8_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_8_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_8_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_8_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_8_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_13_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_13_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_13_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_13_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_13_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_13_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_16_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_16_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_16_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_16_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_16_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_16_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_23_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_23_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_23_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_23_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_23_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_23_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_26_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_26_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_26_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_26_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_26_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_16_26_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_5_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_5_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_5_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_5_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_5_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_5_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_8_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_8_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_8_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_8_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_8_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_8_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_13_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_13_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_13_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_13_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_13_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_13_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_16_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_16_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_16_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_16_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_16_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_16_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_23_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_23_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_23_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_23_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_23_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_23_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_26_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_26_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_26_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_26_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_26_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_23_26_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_5_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_5_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_5_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_5_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_5_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_5_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_8_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_8_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_8_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_8_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_8_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_8_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_13_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_13_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_13_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_13_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_13_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_13_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_16_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_16_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_16_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_16_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_16_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_16_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_23_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_23_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_23_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_23_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_23_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_23_26_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_26_5_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_26_8_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_26_13_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_26_16_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_26_23_2.h"
#include "dbcsr_kernels/stack_mm_mnk_kepler_NxNd_26_26_26_2.h"


#define dbcsr_type_real_4     1
#define dbcsr_type_real_8     3
#define dbcsr_type_complex_4  5
#define dbcsr_type_complex_8  7

//    static const int verbose_print = 0;


/**
 * \brief Bridge routine to call appropriate CUDA kernel.
 */
extern "C" int
dc_do_stack_cu (int *param_stack, int stack_size, int nparams, int datatype,
		void *a_data, void *b_data, void *c_data,
		int m_max, int n_max, int k_max, int def_mnk,
		cudaStream_t * stream)
{
  int stat;

/*  if (verbose_print)
    {
      printf ("A data %p.\n", a_data);
      printf ("B data %p.\n", b_data);
      printf ("C data %p.\n", c_data);
      printf ("params %p.\n", param_stack);
    }
*/
  // printf("Got m,n,k: %d %d %d; %d.\n", m_max, n_max, k_max, stack_size);
  switch (datatype)
    {
    case dbcsr_type_real_4:
      stat =
	launch_stack_mm_r (param_stack, stack_size, *stream, m_max, n_max,
			   k_max, (float *) a_data, (float *) b_data,
			   (float *) c_data);
      break;
    case dbcsr_type_real_8:
      if (def_mnk)
	{
	  if (m_max == 5 && n_max == 5 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_5_5_2 (param_stack,
							 stack_size, *stream,
							 m_max, n_max, k_max,
							 (double *) a_data,
							 (double *) b_data,
							 (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 5 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_5_8_2 (param_stack,
							 stack_size, *stream,
							 m_max, n_max, k_max,
							 (double *) a_data,
							 (double *) b_data,
							 (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 5 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_5_13_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 5 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_5_16_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 5 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_5_23_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 5 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_5_26_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 8 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_8_5_2 (param_stack,
							 stack_size, *stream,
							 m_max, n_max, k_max,
							 (double *) a_data,
							 (double *) b_data,
							 (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 8 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_8_8_2 (param_stack,
							 stack_size, *stream,
							 m_max, n_max, k_max,
							 (double *) a_data,
							 (double *) b_data,
							 (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 8 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_8_13_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 8 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_8_16_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 8 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_8_23_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 8 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_8_26_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 13 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_13_5_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 13 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_13_8_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 13 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_13_13_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 13 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_13_16_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 13 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_13_23_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 13 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_13_26_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 16 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_16_5_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 16 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_16_8_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 16 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_16_13_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 16 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_16_16_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 16 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_16_23_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 16 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_16_26_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 23 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_23_5_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 23 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_23_8_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 23 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_23_13_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 23 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_23_16_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 23 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_23_23_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 23 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_23_26_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 26 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_26_5_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 26 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_26_8_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 26 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_26_13_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 26 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_26_16_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 26 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_26_23_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 5 && n_max == 26 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_5_26_26_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 5 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_5_5_2 (param_stack,
							 stack_size, *stream,
							 m_max, n_max, k_max,
							 (double *) a_data,
							 (double *) b_data,
							 (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 5 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_5_8_2 (param_stack,
							 stack_size, *stream,
							 m_max, n_max, k_max,
							 (double *) a_data,
							 (double *) b_data,
							 (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 5 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_5_13_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 5 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_5_16_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 5 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_5_23_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 5 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_5_26_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 8 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_8_5_2 (param_stack,
							 stack_size, *stream,
							 m_max, n_max, k_max,
							 (double *) a_data,
							 (double *) b_data,
							 (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 8 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_8_8_2 (param_stack,
							 stack_size, *stream,
							 m_max, n_max, k_max,
							 (double *) a_data,
							 (double *) b_data,
							 (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 8 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_8_13_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 8 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_8_16_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 8 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_8_23_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 8 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_8_26_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 13 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_13_5_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 13 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_13_8_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 13 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_13_13_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 13 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_13_16_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 13 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_13_23_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 13 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_13_26_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 16 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_16_5_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 16 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_16_8_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 16 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_16_13_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 16 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_16_16_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 16 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_16_23_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 16 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_16_26_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 23 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_23_5_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 23 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_23_8_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 23 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_23_13_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 23 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_23_16_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 23 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_23_23_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 23 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_23_26_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 26 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_26_5_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 26 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_26_8_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 26 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_26_13_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 26 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_26_16_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 26 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_26_23_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 8 && n_max == 26 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_8_26_26_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 5 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_5_5_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 5 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_5_8_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 5 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_5_13_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 5 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_5_16_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 5 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_5_23_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 5 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_5_26_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 8 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_8_5_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 8 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_8_8_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 8 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_8_13_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 8 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_8_16_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 8 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_8_23_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 8 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_8_26_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 13 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_13_5_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 13 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_13_8_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 13 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_13_13_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 13 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_13_16_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 13 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_13_23_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 13 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_13_26_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 16 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_16_5_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 16 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_16_8_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 16 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_16_13_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 16 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_16_16_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 16 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_16_23_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 16 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_16_26_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 23 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_23_5_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 23 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_23_8_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 23 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_23_13_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 23 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_23_16_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 23 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_23_23_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 23 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_23_26_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 26 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_26_5_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 26 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_26_8_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 26 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_26_13_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 26 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_26_16_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 26 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_26_23_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 13 && n_max == 26 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_13_26_26_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 5 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_5_5_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 5 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_5_8_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 5 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_5_13_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 5 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_5_16_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 5 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_5_23_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 5 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_5_26_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 8 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_8_5_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 8 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_8_8_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 8 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_8_13_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 8 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_8_16_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 8 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_8_23_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 8 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_8_26_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 13 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_13_5_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 13 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_13_8_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 13 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_13_13_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 13 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_13_16_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 13 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_13_23_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 13 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_13_26_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 16 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_16_5_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 16 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_16_8_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 16 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_16_13_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 16 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_16_16_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 16 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_16_23_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 16 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_16_26_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 23 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_23_5_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 23 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_23_8_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 23 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_23_13_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 23 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_23_16_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 23 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_23_23_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 23 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_23_26_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 26 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_26_5_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 26 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_26_8_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 26 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_26_13_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 26 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_26_16_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 26 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_26_23_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 16 && n_max == 26 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_16_26_26_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 5 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_5_5_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 5 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_5_8_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 5 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_5_13_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 5 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_5_16_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 5 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_5_23_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 5 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_5_26_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 8 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_8_5_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 8 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_8_8_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 8 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_8_13_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 8 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_8_16_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 8 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_8_23_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 8 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_8_26_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 13 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_13_5_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 13 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_13_8_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 13 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_13_13_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 13 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_13_16_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 13 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_13_23_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 13 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_13_26_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 16 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_16_5_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 16 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_16_8_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 16 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_16_13_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 16 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_16_16_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 16 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_16_23_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 16 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_16_26_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 23 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_23_5_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 23 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_23_8_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 23 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_23_13_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 23 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_23_16_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 23 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_23_23_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 23 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_23_26_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 26 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_26_5_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 26 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_26_8_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 26 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_26_13_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 26 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_26_16_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 26 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_26_23_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 23 && n_max == 26 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_23_26_26_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 5 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_5_5_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 5 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_5_8_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 5 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_5_13_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 5 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_5_16_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 5 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_5_23_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 5 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_5_26_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 8 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_8_5_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 8 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_8_8_2 (param_stack,
							  stack_size, *stream,
							  m_max, n_max, k_max,
							  (double *) a_data,
							  (double *) b_data,
							  (double *) c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 8 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_8_13_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 8 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_8_16_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 8 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_8_23_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 8 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_8_26_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 13 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_13_5_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 13 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_13_8_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 13 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_13_13_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 13 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_13_16_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 13 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_13_23_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 13 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_13_26_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 16 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_16_5_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 16 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_16_8_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 16 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_16_13_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 16 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_16_16_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 16 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_16_23_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 16 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_16_26_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 23 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_23_5_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 23 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_23_8_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 23 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_23_13_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 23 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_23_16_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 23 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_23_23_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 23 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_23_26_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 26 && k_max == 5)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_26_5_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 26 && k_max == 8)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_26_8_2 (param_stack,
							   stack_size,
							   *stream, m_max,
							   n_max, k_max,
							   (double *) a_data,
							   (double *) b_data,
							   (double *) c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 26 && k_max == 13)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_26_13_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 26 && k_max == 16)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_26_16_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 26 && k_max == 23)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_26_23_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  if (m_max == 26 && n_max == 26 && k_max == 26)
	    {
	      stat =
		launch_stack_mm_mnk_kepler_NxNd_26_26_26_2 (param_stack,
							    stack_size,
							    *stream, m_max,
							    n_max, k_max,
							    (double *) a_data,
							    (double *) b_data,
							    (double *)
							    c_data);

	      break;
	    }
	  stat =
	    launch_stack_mm_mnk_d (param_stack, stack_size, *stream, m_max,
				   n_max, k_max, (double *) a_data,
				   (double *) b_data, (double *) c_data);

	  break;
	}
      stat =
	launch_stack_mm_d (param_stack, stack_size, *stream, m_max, n_max,
			   k_max, (double *) a_data, (double *) b_data,
			   (double *) c_data);
      break;
    case dbcsr_type_complex_4:
      stat =
	launch_stack_mm_c (param_stack, stack_size, *stream, m_max, n_max,
			   k_max, (float *) a_data, (float *) b_data,
			   (float *) c_data);
      break;
    case dbcsr_type_complex_8:
      stat =
	launch_stack_mm_z (param_stack, stack_size, *stream, m_max, n_max,
			   k_max, (double *) a_data, (double *) b_data,
			   (double *) c_data);
      break;
    }
  if (stat != 0)
    return (stat);
  if (cuda_error_check (cudaGetLastError ()))
    return 1;
  return 0;
};

//EOF
