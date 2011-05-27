/******************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations
 *  Copyright (C) 2000 - 2011  Urban Borstnik and the CP2K developers group
 *****************************************************************************/

#include <cuda_runtime.h>
#include <stdio.h>
#include <sm_11_atomic_functions.h>

#include "dbcsr_cuda.h"


int cuda_error_check (cudaError_t cudaError) {
  if (cudaError != cudaSuccess) {
    printf("CUDA Error: %s\n", cudaGetErrorString(cudaError));
    return 1;
  }
  return 0;
};


/**
 * \var cache  Per-threadblock cache for A and B data.
 */
extern __shared__ double cache[];


/* These file are included here to avoid linking issues. */

__global__ void stack_mm_r
                   (const int *__restrict__ param_stack,
		    int stack_size, int nparams,
		    const float *__restrict__ a_data,
		    const float *__restrict__ b_data,
		    float *__restrict__ c_data,
		    int *__restrict__ c_locks);
__global__ void stack_mm_d
                   (const int *__restrict__ param_stack,
		    int stack_size, int nparams,
		    const double *__restrict__ a_data,
		    const double *__restrict__ b_data,
		    double *__restrict__ c_data,
		    int *__restrict__ c_locks);
__global__ void stack_mm_c
                   (const int *__restrict__ param_stack,
		    int stack_size, int nparams,
		    const float *__restrict__ a_data,
		    const float *__restrict__ b_data,
		    float *__restrict__ c_data,
		    int *__restrict__ c_locks);
__global__ void stack_mm_z
                   (const int *__restrict__ param_stack,
		    int stack_size, int nparams,
		    const double *__restrict__ a_data,
		    const double *__restrict__ b_data,
		    double *__restrict__ c_data,
		    int *__restrict__ c_locks);


/**
 * \brief Bridge routine to call appropriate CUDA kernel.
 */
extern "C" int dc_do_stack_cu(int *param_stack, int stack_size, int nparams,
			     int which_data,
			     void *a_data, void *b_data, void *c_data,
			     int *c_locks,
			     int m_max, int n_max, int k_max) {
  int maxt;
  cudaError_t cErr;
  int myDevice;
  size_t shared_size;
  struct cudaDeviceProp devProperties;

  maxt = m_max * n_max;

  cErr = cudaGetDevice(&myDevice);
  if (cuda_error_check (cErr)) return 1;

  cErr = cudaGetDeviceProperties(&devProperties, myDevice);
  if (cuda_error_check (cErr)) return 1;

  if (maxt > devProperties.maxThreadsPerBlock)
    return 3;

  switch (which_data) {
    /* The data type identifier numbers correspond to the values
       defined in dbcsr_types.F. */
  case 1:
    /* Real, single precision */
    shared_size = (m_max*k_max + k_max*n_max)*sizeof(float);
    if (shared_size > devProperties.sharedMemPerBlock) return 4;
    stack_mm_r <<< stack_size, maxt, shared_size >>>
      (param_stack, stack_size, nparams,
       (float *) a_data, (float *) b_data, (float *) c_data,
       c_locks);
    break;
  case 3:
    /* Real, double precision */
    shared_size = (m_max*k_max + k_max*n_max)*sizeof(double);
    if (shared_size > devProperties.sharedMemPerBlock) return 4;
    stack_mm_d <<< stack_size, maxt, shared_size >>>
      (param_stack, stack_size, nparams,
       (double *) a_data, (double *) b_data, (double *) c_data,
       c_locks);
    break;
  case 5:
    /* Complex, single precision */
    shared_size = (m_max*k_max + k_max*n_max)*sizeof(float)*2;
    if (shared_size > devProperties.sharedMemPerBlock) return 4;
    stack_mm_c <<< stack_size, maxt, shared_size >>>
      (param_stack, stack_size, nparams,
       (float *) a_data, (float *) b_data, (float *) c_data,
       c_locks);
    break;
  case 7:
    /* Complex, single precision */
    shared_size = (m_max*k_max + k_max*n_max)*sizeof(double)*2;
    if (shared_size > devProperties.sharedMemPerBlock) return 4;
    stack_mm_z <<< stack_size, maxt, shared_size >>>
      (param_stack, stack_size, nparams,
       (double *) a_data, (double *) b_data, (double *) c_data,
       c_locks);
    break;
  default:
    return 2;
  }
  if (cuda_error_check (cudaGetLastError())) return 1;

  return 0;
};


extern "C" int dc_thread_sync_cu() {
  cudaError_t cErr;

  cErr = cudaThreadSynchronize ();
  if (cuda_error_check (cErr)) return 1;
  return 0;
}    
