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
#include "libcusmm/include/libcusmm.h"

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
        if (def_mnk){
          stat = libcusmm_process_d (param_stack, stack_size, *stream, m_max,
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
    return 100;
  return 0;
};

//EOF
