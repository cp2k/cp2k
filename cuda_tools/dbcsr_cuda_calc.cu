/******************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations
 *  Copyright (C) 2000 - 2011  Urban Borstnik and the CP2K developers group
 *****************************************************************************/

#include <cuda_runtime.h>
#include <stdio.h>
#include <sm_11_atomic_functions.h>

#include "dbcsr_cuda.h"

#define MAX(a,b) ((a > b) ? (a) : (b))


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

__global__ void stack_mm_d
                   (const int* __restrict__ param_stack,
		    int stack_size, int nparams,
		    const double* __restrict__ a_data,
		    const double* __restrict__ b_data,
		    double* __restrict__ c_data,
		    int* __restrict__ c_locks) {

  /**
   *  \var sp        which stack member this thread block is processing
   *  \var our_params  cache for this thread block's multiplication parameters
   *  \var m, n, k   dimensions of the blocks (C is m*n, A is m*k, B is k*n)
   *  \var mn, mk, kn  product of the block dimensions
   *  \var l         multiplication loop index
   *  \var cr        C matrix element (fortran-style) of this thread
   *  \var c, r      C matrix row, column of this thread
   *  \var myc       C matrix accumulator
   *  \var buff      cache for A and B data
   */ 

  int sp, lock_owner, c_id, sp_one;
  int r, c, l;
  int m, n, k;
  int mn, cr, mk, kn;
  double myc;
  __shared__ int our_params[7];
  double *buff;


  /* Setup shared memory. */
  buff = (double *) cache;

  /* Determine who I am. */
  sp = blockIdx.x;
  cr = threadIdx.x;

  /* Load in the parameters. */
  if (blockDim.x >= 7) {
    if (cr < 7) {
      our_params[cr] = param_stack[cr+7*sp];
    }
  } else if (cr == 0) {
    for (l = 0; l < 7; l++) {
      our_params[l] = param_stack[l+7*sp];
    }
  }
  syncthreads();
  m = our_params[0];
  n = our_params[1];
  k = our_params[2];

  /* Load in the buffers. */
  mk = m*k;
  kn = k*n;
  if (cr < mk)
    buff[cr] = a_data[our_params[3]-1+cr];
  if (cr >= blockDim.x - kn)
    buff[mk+cr-(blockDim.x-kn)] = b_data[our_params[4]-1+(cr-(blockDim.x-kn))];


  /* Calculate who I am. */
  syncthreads();

  mn = m*n;

  /* Do actual multiplication. */
  if (cr < mn) {
    r = cr % m;
    c = cr / m;
    myc = 0.0;

    for (l = 0; l < k; l++) {
      myc = myc +
	buff[   l*m+r] *
	buff[mk+c*k+l];
    }
  }

  /* Lock the C block. */
  syncthreads();
  if (cr == 0) {
    sp_one = sp + 1;
    c_id = our_params[6]-1;
    lock_owner = 0;
    while ((lock_owner != sp_one))
      lock_owner = atomicCAS (&(c_locks[c_id]), 0, sp_one);
  }

  /* Add our results to the C block. */
  syncthreads();
  if (cr < mn) {
    c_data[our_params[5]-1+cr] += myc;
  }

  /* Release the lock on the C block. */
  syncthreads();
  if (cr == 0)
    c_locks[c_id] = 0;

  if (mn == 16) return;

};

/**
 * \var cache  Per-threadblock cache for A and B data.
 */
extern __shared__ double cache[];

__global__ void stack_mm_r
                   (const int* __restrict__ param_stack,
		    int stack_size, int nparams,
		    const float* __restrict__ a_data,
		    const float* __restrict__ b_data,
		    float* __restrict__ c_data,
		    int* __restrict__ c_locks) {

  /**
   *  \var sp        which stack member this thread block is processing
   *  \var our_params  cache for this thread block's multiplication parameters
   *  \var m, n, k   dimensions of the blocks (C is m*n, A is m*k, B is k*n)
   *  \var mn, mk, kn  product of the block dimensions
   *  \var l         multiplication loop index
   *  \var cr        C matrix element (fortran-style) of this thread
   *  \var c, r      C matrix row, column of this thread
   *  \var myc       C matrix accumulator
   *  \var buff      cache for A and B data
   */ 

  int sp, lock_owner, c_id, sp_one;
  int r, c, l;
  int m, n, k;
  int mn, cr, mk, kn;
  float myc;
  __shared__ int our_params[7];
  float *buff;


  /* Setup shared memory. */
  buff = (float *) cache;

  /* Determine who I am. */
  sp = blockIdx.x;
  cr = threadIdx.x;

  /* Load in the parameters. */
  if (blockDim.x >= 7) {
    if (cr < 7) {
      our_params[cr] = param_stack[cr+7*sp];
    }
  } else if (cr == 0) {
    for (l = 0; l < 7; l++) {
      our_params[l] = param_stack[l+7*sp];
    }
  }
  syncthreads();
  m = our_params[0];
  n = our_params[1];
  k = our_params[2];

  /* Load in the buffers. */
  mk = m*k;
  kn = k*n;
  if (cr < mk)
    buff[cr] = a_data[our_params[3]-1+cr];
  if (cr >= blockDim.x - kn)
    buff[mk+cr-(blockDim.x-kn)] = b_data[our_params[4]-1+(cr-(blockDim.x-kn))];


  /* Calculate who I am. */
  syncthreads();

  mn = m*n;

  /* Do actual multiplication. */
  if (cr < mn) {
    r = cr % m;
    c = cr / m;
    myc = 0.0;

    for (l = 0; l < k; l++) {
      myc = myc +
	buff[   l*m+r] *
	buff[mk+c*k+l];
    }
  }

  /* Lock the C block. */
  syncthreads();
  if (cr == 0) {
    sp_one = sp + 1;
    c_id = our_params[6]-1;
    lock_owner = 0;
    while ((lock_owner != sp_one))
      lock_owner = atomicCAS (&(c_locks[c_id]), 0, sp_one);
  }

  /* Add our results to the C block. */
  syncthreads();
  if (cr < mn) {
    c_data[our_params[5]-1+cr] += myc;
  }

  /* Release the lock on the C block. */
  syncthreads();
  if (cr == 0)
    c_locks[c_id] = 0;

  if (mn == 16) return;

};


__global__ void stack_mm_z
                   (const int* __restrict__ param_stack,
		    int stack_size, int nparams,
		    const double* __restrict__ a_data,
		    const double* __restrict__ b_data,
		    double* __restrict__ c_data,
		    int* __restrict__ c_locks) {
  /**
   *  \var sp        which stack member this thread block is processing
   *  \var our_params  cache for this thread block's multiplication parameters
   *  \var m, n, k   dimensions of the blocks (C is m*n, A is m*k, B is k*n)
   *  \var mn, mk, kn  product of the block dimensions
   *  \var l         multiplication loop index
   *  \var cr        C matrix element (fortran-style) of this thread
   *  \var c, r      C matrix row, column of this thread
   *  \var myc       C matrix accumulator
   *  \var buff      cache for A and B data
   */ 
  int sp, lock_owner, c_id, sp_one;
  int r, c, l;
  int m, n, k;
  int mn, cr, mk, kn;
  double myc_r, myc_i;
  __shared__ int our_params[7];
  double *buff;


  /* Setup shared memory. */
  buff = (double *) cache;

  /* Determine who I am. */
  sp = blockIdx.x;
  cr = threadIdx.x;

  /* Load in the parameters. */
  if (blockDim.x >= 7) {
    if (cr < 7) {
      our_params[cr] = param_stack[cr+7*sp];
    }
  } else if (cr == 0) {
    for (l = 0; l < 7; l++) {
      our_params[l] = param_stack[l+7*sp];
    }
  }
  syncthreads();
  m = our_params[0];
  n = our_params[1];
  k = our_params[2];

  /* Load in the buffers. */
  mk = m*k;
  kn = k*n;
  if (cr < mk) {
    buff[2*cr] = a_data[2*(our_params[3]-1+cr)];
    buff[2*(mk+cr)+1] = a_data[2*(our_params[3]-1+cr)+1];
  }
  if (cr >= blockDim.x - kn)
    buff[2*(mk+(cr-(blockDim.x-kn)))] = b_data[2*(our_params[4]-1+(cr-(blockDim.x-kn)))];
    buff[2*(mk+(cr-(blockDim.x-kn)))+1] = b_data[2*(our_params[4]-1+(cr-(blockDim.x-kn)))+1];

  /* Calculate who I am. */
  syncthreads();
  mn = m*n;

  /* Do actual multiplication. */
  if (cr < mn) {
    r = cr % m;
    c = cr / m;
    myc_r = 0.0;
    myc_i = 0.0;

    for (l = 0; l < k; l++) {
      myc_r = myc_r +
	buff[2*(   l*m+r)] *
	buff[2*(mk+c*k+l)] -
	buff[2*(   l*m+r)+1] *
	buff[2*(mk+c*k+l)+1];
      myc_i = myc_i +
	buff[2*(   l*m+r)] *
	buff[2*(mk+c*k+l)+1] +
	buff[2*(   l*m+r)+1] *
	buff[2*(mk+c*k+l)];
    }
  }

  /* Lock the C block. */
  syncthreads();
  if (cr == 0) {
    sp_one = sp + 1;
    c_id = our_params[6]-1;
    lock_owner = 0;
    while ((lock_owner != sp_one))
      lock_owner = atomicCAS (&(c_locks[c_id]), 0, sp_one);
  }

  /* Add our results to the C block. */
  syncthreads();
  if (cr < mn) {
    c_data[2*(our_params[5]-1+cr)] += myc_r;
    c_data[2*(our_params[5]-1+cr)+1] += myc_i;
  }

  /* Release the lock on the C block. */
  syncthreads();
  if (cr == 0)
    c_locks[c_id] = 0;

};

__global__ void stack_mm_c
                   (const int* __restrict__ param_stack,
		    int stack_size, int nparams,
		    const float* __restrict__ a_data,
		    const float* __restrict__ b_data,
		    float* __restrict__ c_data,
		    int* __restrict__ c_locks) {
  /**
   *  \var sp        which stack member this thread block is processing
   *  \var our_params  cache for this thread block's multiplication parameters
   *  \var m, n, k   dimensions of the blocks (C is m*n, A is m*k, B is k*n)
   *  \var mn, mk, kn  product of the block dimensions
   *  \var l         multiplication loop index
   *  \var cr        C matrix element (fortran-style) of this thread
   *  \var c, r      C matrix row, column of this thread
   *  \var myc       C matrix accumulator
   *  \var buff      cache for A and B data
   */ 
  int sp, lock_owner, c_id, sp_one;
  int r, c, l;
  int m, n, k;
  int mn, cr, mk, kn;
  float myc_r, myc_i;
  __shared__ int our_params[7];
  float *buff;


  /* Setup shared memory. */
  buff = (float *) cache;

  /* Determine who I am. */
  sp = blockIdx.x;
  cr = threadIdx.x;

  /* Load in the parameters. */
  if (blockDim.x >= 7) {
    if (cr < 7) {
      our_params[cr] = param_stack[cr+7*sp];
    }
  } else if (cr == 0) {
    for (l = 0; l < 7; l++) {
      our_params[l] = param_stack[l+7*sp];
    }
  }
  syncthreads();
  m = our_params[0];
  n = our_params[1];
  k = our_params[2];

  /* Load in the buffers. */
  mk = m*k;
  kn = k*n;
  if (cr < mk) {
    buff[2*cr] = a_data[2*(our_params[3]-1+cr)];
    buff[2*(mk+cr)+1] = a_data[2*(our_params[3]-1+cr)+1];
  }
  if (cr >= blockDim.x - kn)
    buff[2*(mk+(cr-(blockDim.x-kn)))] = b_data[2*(our_params[4]-1+(cr-(blockDim.x-kn)))];
    buff[2*(mk+(cr-(blockDim.x-kn)))+1] = b_data[2*(our_params[4]-1+(cr-(blockDim.x-kn)))+1];

  /* Calculate who I am. */
  syncthreads();
  mn = m*n;

  /* Do actual multiplication. */
  if (cr < mn) {
    r = cr % m;
    c = cr / m;
    myc_r = 0.0;
    myc_i = 0.0;

    for (l = 0; l < k; l++) {
      myc_r = myc_r +
	buff[2*(   l*m+r)] *
	buff[2*(mk+c*k+l)] -
	buff[2*(   l*m+r)+1] *
	buff[2*(mk+c*k+l)+1];
      myc_i = myc_i +
	buff[2*(   l*m+r)] *
	buff[2*(mk+c*k+l)+1] +
	buff[2*(   l*m+r)+1] *
	buff[2*(mk+c*k+l)];
    }
  }

  /* Lock the C block. */
  syncthreads();
  if (cr == 0) {
    sp_one = sp + 1;
    c_id = our_params[6]-1;
    lock_owner = 0;
    while ((lock_owner != sp_one))
      lock_owner = atomicCAS (&(c_locks[c_id]), 0, sp_one);
  }

  /* Add our results to the C block. */
  syncthreads();
  if (cr < mn) {
    c_data[2*(our_params[5]-1+cr)] += myc_r;
    c_data[2*(our_params[5]-1+cr)+1] += myc_i;
  }

  /* Release the lock on the C block. */
  syncthreads();
  if (cr == 0)
    c_locks[c_id] = 0;

};



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
  int test[1];

  maxt = MAX(MAX(m_max*n_max, m_max*k_max), k_max*n_max);

  cErr = cudaGetDevice(&myDevice);
  if (cuda_error_check (cErr)) return 1;

  cErr = cudaGetDeviceProperties(&devProperties, myDevice);
  if (cuda_error_check (cErr)) return 1;

  if (maxt > devProperties.maxThreadsPerBlock)
    return 3;

  //cErr = cudaMemset ((void *) c_locks, (int) 0, sizeof(int));
  //if (cuda_error_check (cErr)) return 1;

//  printf("nparams %d, %d, m, n, k max %d %d %d\n", nparams, stack_size, m_max, n_max, k_max);
//
//  cErr = cudaMemcpy((void *) test, (void *) c_locks, sizeof(int), cudaMemcpyDeviceToHost);
//  if (cuda_error_check (cErr)) return 1;
//  printf("c_locks[0]: %d\n", c_locks);
//
  /*
  printf("param.stack %p\n", param_stack);
  printf("a_data %p\n", a_data);
  printf("b_data %p\n", b_data);
  printf("c_data %p\n", c_data);
  printf("c_lock %p\n", c_locks);
  */
  //  test = (int *) malloc ((size_t) (stack_size*7*sizeof(int)));
//  cErr = cudaMemcpy((void *) test, (void *) param_stack, 7*stack_size*sizeof(int), cudaMemcpyDeviceToHost);
//  if (cuda_error_check (cErr)) return 1;
//  max_c = 0;
//  for (sp = 0; sp < stack_size; sp++)
//    max_c = MAX(max_c, test[sp*7+6]);
//  printf("max_c %d %d\n", max_c, stack_size);
//  test = (int *) realloc(test, (size_t) (max_c*sizeof(int)));
//  cErr = cudaMemcpy((void *) test, (void *) c_locks, (size_t) (max_c*sizeof(int)), cudaMemcpyDeviceToHost);
//  if (cuda_error_check (cErr)) return 1;
//  n = 0;
//  for (sp = 0; sp < max_c; sp++)
//    if (test[sp] != 0) n++;
//  printf("nz %d\n", n);
// // free(test);

//  cErr = cudaMemset ((void *) c_locks, (int) 0, (size_t) (sizeof(int)*max_c));
//  if (cuda_error_check (cErr)) return 1;  

  switch (which_data) {
  case 1:
    shared_size = (m_max*k_max + k_max*n_max)*sizeof(float);
    if (shared_size > devProperties.sharedMemPerBlock) return 4;
    stack_mm_r <<< stack_size, maxt >>>
      (param_stack, stack_size, nparams,
       (float *) a_data, (float *) b_data, (float *) c_data,
       c_locks);
    break;
  case 3:
    shared_size = (m_max*k_max + k_max*n_max)*sizeof(double);
    if (shared_size > devProperties.sharedMemPerBlock) return 4;
    stack_mm_d <<< stack_size, maxt, shared_size >>>
      (param_stack, stack_size, nparams,
       (double *) a_data, (double *) b_data, (double *) c_data,
       c_locks);
    break;
  case 5:
    shared_size = (m_max*k_max + k_max*n_max)*sizeof(float)*2;
    if (shared_size > devProperties.sharedMemPerBlock) return 4;
    stack_mm_c <<< stack_size, maxt >>>
      (param_stack, stack_size, nparams,
       (float *) a_data, (float *) b_data, (float *) c_data,
       c_locks);
    break;
  case 7:
    shared_size = (m_max*k_max + k_max*n_max)*sizeof(double)*2;
    if (shared_size > devProperties.sharedMemPerBlock) return 4;
    stack_mm_z <<< stack_size, maxt >>>
      (param_stack, stack_size, nparams,
       (double *) a_data, (double *) b_data, (double *) c_data,
       c_locks);
    break;
  default:
    return 2;
  }
  if (cuda_error_check (cudaGetLastError())) return 1;

  
//  printf("syncing...\n");
//  cErr = cudaThreadSynchronize ();
//  if (cuda_error_check (cErr)) return 1;


  /*
  printf("dl...\n");

  test = (int *) malloc ((size_t) (stack_size*7*sizeof(int)));
  cErr = cudaMemcpy((void *) test, (void *) param_stack,
		    (size_t) (7*stack_size*sizeof(int)),
		    cudaMemcpyDeviceToHost);
  if (cuda_error_check (cErr)) return 1;
  for (sp = 0; sp < 6; sp++)
    printf ("%d ", test[sp]);
  printf("%d\n", test[6]);
*/
  //  n = 0;
//  for (sp = 0; sp < max_c; sp++)
//    if (test[sp] != 0) n++;
//  printf("nz %d\n", n);
//
//  free(test);

  return 0;
};


extern "C" int dc_thread_sync_cu() {
  cudaError_t cErr;

  cErr = cudaThreadSynchronize ();
  if (cuda_error_check (cErr)) return 1;
  return 0;
}    
