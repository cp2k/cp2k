/******************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations
 *  Copyright (C) 2000 - 2013 the CP2K developers group
 *****************************************************************************/

#include "dbcsr_kernel.h"
#include "stack_mm_mnk_sq23_d_slow.h"

extern __shared__ double cache[];

/* The following are defined in "dbcsr_cuda.h" */
/* SQUARESIZE, BLOCK, TDIM, NUMTHREADS23SQ, BLOCKSPLIT23SQ, TOTALTHREADS23SQ */

/******************************************************************************/
/* The stack_mm_mnk_sq23_d routine was added by Neil Stringfellow */
/* This is a specialised routine for square matrices of size 23 (m==n==k==23) */
/******************************************************************************/

/* The following is only valid for square matrices of size 23 (m==n==k==23) */
/* Should be called with 64 threads only (blockDim.x==64) */
/* N.B. This is checked for in the calling routine, not down here */

/* Why 64 threads I here you ask ?  What about occupancy !!!!! */
/* !!! Update - using BLOCKSPLIT23SQ we can have 128 threads !!! */
/* I'm glad you asked ... */
/* With 64 threads we can provide the device with lots of work, and with a small
   number of threads we have a better chance of having multiple blocks being
   executed concurrently so that they will overlap core floating point computation
   and global memory traffic, hiding the memory bandwidth pressure that we saw in
   the original implemenation.
   The number of thread-blocks that can be executed concurrently is probably more
   driven by the shared memory requirements per thread-block which is a little
   over 8Kbytes, and therefore we should be able to get up to 5 blocks working at
   the same time per SM.
   If indeed we can have a maximum of 5 thread-blocks each of 64 threads working
   concurrently then we have at most 320 threads, and these threads then can use
   the full complement of 63 32-bit registers available as a maximum per thread.
   This then means that we can hold up to 31 double precision numbers in registers
   and allows the inner kernel (see the code) to operate C=Matmul(A,B) on 3x3
   square matrices entirely out of registers (3 arrays each of 3x3 = 27 registers
   required for the inner kernel */

/* Needs 24x24*8*2 bytes of shared memory when being called */
__global__ void
stack_mm_mnk_sq23_d_slow (const int *__restrict__ param_stack,
		     const int careful, const int nruns,
		     const int m, const int n, const int k,
		     const int liter,
		     const double *__restrict__ a_data,
		     const double *__restrict__ b_data,
		     double *__restrict__ c_data, int *__restrict__ c_locks,
		     int lock_offset)
{

  const int mn = m * n;
  const int mk = m * k;
  const int kn = n * k;
  const int r = threadIdx.x % m;
  const int c = threadIdx.x / m;
  double myc;
  const double *__restrict__ buff_l, *__restrict__ buff_r;

  int psp, c_loc;

  double *buff;

  buff = (double *) cache;
  buff_l = buff;
  buff_r = &(buff[mk]);

  int nrun = GROUPING;

  if (blockIdx.x == gridDim.x-1)
    nrun = nruns;

  /* Set the partial sum to zero */
  myc = 0.0;
  for (int run = 0; run < nrun; run++)
    {
      psp = 7 * (blockIdx.x * GROUPING + run);

      // load a and b matrix into smem
      for(int i = threadIdx.x; i < mk; i += blockDim.x){
          buff[i] = a_data[param_stack[psp + 3] - 1 + i];
      }

      for(int i = threadIdx.x; i < kn; i += blockDim.x){
          buff[mk + i] = b_data[param_stack[psp + 4] - 1 + i];
      }
      syncthreads ();

      /* Do actual multiplication. */
      if (threadIdx.x < mn)
        {

          for (int l = 0; l < k; l++)
            {
              myc = myc + buff_l[l * m + r] * buff_r[c * k + l];
              //printf("GPU : a[%d, %d] = %g, b[%d, %d] = %g\n", r, l, buff_l[l * m + r] , l, c, buff_r[c * k + l]);
            }

        }
      if (run == nrun - 1
          || param_stack[psp + 6] - 1 != param_stack[psp + 6 + 7] - 1)
        {
          /* Add results to global C block. */
          c_loc = param_stack[psp + 5] - 1;
          if (threadIdx.x < mn)
            atomicAdd (&c_data[c_loc + threadIdx.x], myc);

          myc = 0.0l;
        }

      syncthreads ();

    }

}

