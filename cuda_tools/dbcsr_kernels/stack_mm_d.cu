/******************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations
 *  Copyright (C) 2000 - 2013 the CP2K developers group
 *****************************************************************************/

// DBCSR_KERNEL datatype=dbcsr_type_real_8, homogeneous_only=False

#include "dbcsr_kernel.h"
#include "stack_mm_d.h"
#include <stdio.h>

//==============================================================================
__global__ void stack_mm_d
  (const int *__restrict__ param_stack,
   int stack_size,
   const double *__restrict__ a_data,
   const double *__restrict__ b_data,
   double *__restrict__ c_data)
{

  /**
   *  \var sp        which stack member this thread block is processing
                     (= CUDA thread block)
   *  \var sp_one    translated stack (=sp+1)
   *  \var tn        thread number (of CUDA thread block)
   *  \var nt        number of threads (size of CUDA thread block)
   *  \var m, n, k   dimensions of the blocks (C is m*n, A is m*k, B is k*n)
   *  \var mn, mk, kn  product of the block dimensions
   *  \var l         multiplication loop index
   *  \var c, r      C matrix row, column of this thread
   *  \var myc       C matrix accumulator
   *  \var buff      cache for A and B data
   */

  int sp;
  int tn;
  int r, c, l;
  int m, n, k;
  int mn;
  double myc;
  const double *buff_l, *buff_r;

  int psp, c_loc;

  /* Setup shared memory. */
  //buff = (double *) cache;

  /* Determine who I am. */
  sp = blockIdx.x;
  tn = threadIdx.x;

  psp = 7 * sp;
  m = param_stack[psp];
  n = param_stack[psp + 1];
  k = param_stack[psp + 2];

  buff_l = &(a_data[param_stack[psp + 3] - 1]);
  buff_r = &(b_data[param_stack[psp + 4] - 1]);

  /* Calculate who I am. */

  mn = m * n;

  /* Do actual multiplication. */
  if (tn < mn)
    {
      r = tn % m;
      c = tn / m;
      myc = 0.0l;

      for (l = 0; l < k; l++)
	{
	  myc = myc + buff_l[l * m + r] * buff_r[c * k + l];
	}
    }

  /* Add results to global C block. */
  c_loc = param_stack[psp + 5] - 1;
  if (tn < mn)
    atomicAdd (&c_data[c_loc + tn], myc);

}


//==============================================================================
int launch_stack_mm_d(int *param_stack, int stack_size, cudaStream_t stream,
    int m_max, int n_max, int k_max,
    double *a_data, double *b_data, double *c_data){

    //printf("launch_stack_mm_d: m,n,k: %d %d %d; %d.\n", m_max, n_max, k_max, stack_size);
     int shared_size = (m_max * k_max + k_max * n_max) * sizeof (double);
     if (shared_size > devProperties.sharedMemPerBlock)
          return 4;

     int maxt = m_max * n_max;
     if (maxt > devProperties.maxThreadsPerBlock)
         return 3;

     stack_mm_d <<< stack_size, maxt, shared_size, stream >>>
        (param_stack, stack_size, a_data, b_data, c_data);
     return(0);
}
