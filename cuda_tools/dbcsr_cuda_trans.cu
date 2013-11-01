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


/**
 * \brief Bridge routine to call appropriate CUDA transpose kernel.
 */
extern "C" int
dc_do_transpose_stack_cu (int *trs_stack, int offset, int nblks, void *buffer,
                         int datatype, int m, int n, cudaStream_t * stream) {


   if(datatype != dbcsr_type_real_8) return 0;

   return libcusmm_transpose_d(trs_stack, offset, nblks, (double*) buffer, m, n, stream);
};

//EOF
