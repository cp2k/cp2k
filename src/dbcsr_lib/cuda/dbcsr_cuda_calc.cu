/******************************************************************************
*  CP2K: A general program to perform molecular dynamics simulations
*  Copyright (C) 2000 - 2013  Urban Borstnik and the CP2K developers group
*****************************************************************************/
#include <cuda_runtime.h>
#include <stdio.h>
#include "dbcsr_cuda.h"
#include "libcusmm/libcusmm.h"

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
                cudaStream_t * stream){

  if(datatype==dbcsr_type_real_8 && def_mnk==1){
      return(libcusmm_process_d (param_stack, stack_size, *stream, m_max, n_max, k_max,
                                 (double *) a_data, (double *) b_data, (double *) c_data));
  }

  // not appropriate kernel was found
  return(-1);
};

//EOF
