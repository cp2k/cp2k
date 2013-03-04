/******************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations
 *  Copyright (C) 2000 - 2013  Urban Borstnik and the CP2K developers group
 *****************************************************************************/

#include "dbcsr_cuda.h"

extern __shared__ double cache[];

__global__ void stack_mm_r
  (const int *__restrict__ param_stack,
   int stack_size, int nparams,
   const float *__restrict__ a_data,
   const float *__restrict__ b_data,
   float *__restrict__ c_data, int *__restrict__ c_locks, int lock_offset)
{

  /**
   *  \var sp        which stack member this thread block is processing
                     (= CUDA thread block)
   *  \var sp_one    translated stack (=sp+1)
   *  \var tn        thread number (of CUDA thread block)
   *  \var nt        number of threads (size of CUDA thread block)
   *  \var our_params  cache for this thread block's multiplication parameters
   *  \var m, n, k   dimensions of the blocks (C is m*n, A is m*k, B is k*n)
   *  \var mn, mk, kn  product of the block dimensions
   *  \var l         multiplication loop index
   *  \var c, r      C matrix row, column of this thread
   *  \var myc       C matrix accumulator
   *  \var buff      cache for A and B data
   *  \var c_id      translated C block number (used in locking)
   *  \var lock_owner  current C block owner (used in locking)
   */

  int sp;			//, lock_owner, c_id, sp_one;
  int tn, nt;
  int r, c, l;
  int m, n, k;
  int mn, mk, kn;
  float myc;
  __shared__ int our_params[7];
  float *buff;


  /* Setup shared memory. */
  buff = (float *) cache;

  /* Determine who I am. */
  sp = blockIdx.x;
  tn = threadIdx.x;
  nt = blockDim.x;

  /* Load in the parameters. */
  for (l = 0; l <= 6 / nt; l++)
    {
      r = tn + nt * l;
      if (r < 7)
	our_params[r] = param_stack[7 * sp + r];
    }

  syncthreads ();
  m = our_params[0];
  n = our_params[1];
  k = our_params[2];

  /* Load in the buffers. */
  mk = m * k;
  kn = k * n;
  for (l = 0; l <= (mk - 1) / nt; l++)
    {
      r = tn + nt * l;
      if (r < mk)
	buff[r] = a_data[our_params[3] - 1 + r];
    }
  for (l = 0; l <= (kn - 1) / nt; l++)
    {
      r = tn + nt * l;
      if (r < kn)
	buff[mk + r] = b_data[our_params[4] - 1 + r];
    }

  /* Calculate who I am. */
  syncthreads ();

  mn = m * n;

  /* Do actual multiplication. */
  if (tn < mn)
    {
      r = tn % m;
      c = tn / m;
      myc = 0.0f;

      for (l = 0; l < k; l++)
	{
	  myc = myc + buff[l * m + r] * buff[mk + c * k + l];
	}
    }

  /* Add results to global C block. */
  if (tn < mn)
    atomicAdd (&c_data[our_params[5] - 1 + tn], myc);

  // /* Lock the C block. */
  // syncthreads();
  // if (tn == 0) {
  //   sp_one = lock_offset + sp + 1;
  //   c_id = our_params[6]-1;
  //   lock_owner = 0;
  //   while ((lock_owner != sp_one))
  //     lock_owner = atomicCAS (&(c_locks[c_id]), 0, sp_one);
  // }
  //
  // /* Add our results to the C block. */
  // syncthreads();
  // if (tn < mn) {
  //   c_data[our_params[5]-1+tn] += myc;
  // }
  //
  // /* Release the lock on the C block. */
  // syncthreads();
  // if (tn == 0) {
  //   c_locks[c_id] = 0;
  //   //threadfence();
  // }

};
