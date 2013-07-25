/******************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations
 *  Copyright (C) 2000 - 2013 the CP2K developers group
 *****************************************************************************/


#include "dbcsr_kernel.h"
#include "stack_mm_mnk_sq5_d.h"


__global__ void
//__launch_bounds__(32, 16)
stack_mm_mnk_sq5_d (const int *__restrict__ param_stack,
		    const int careful, const int nruns,
		    const int m, const int n, const int k,
		    const int liter,
		    const double *__restrict__ a_data,
		    const double *__restrict__ b_data,
		    double *__restrict__ c_data, int *__restrict__ c_locks,
		    int lock_offset)
{

	/**
	 *  \var sp        which stack member this thread block is processing
	 (= CUDA thread block)
	 *  \var psp       pointer to first element of parameters
	 *  \var c_loc     pointer to C data
	 *  \var run       run number
         *  \var nrun      number of runs
	 *  \var my_id     my ID for locking
	 *  \var tn        thread number (of CUDA thread block)
	 *  \var mn        product of the block dimensions
	 *  \var l         multiplication loop index
	 *  \var c, r      C matrix row, column of this thread
	 *  \var myc       C matrix accumulator
	 *  \var buff_l    cache for A data
	 *  \var buff_r    cache for B data
	 *  \var c_id      translated C block number (used in locking)
	 *  \var lock_owner  current C block owner (used in locking)
	 */

  __shared__ double c_s[128];

  //int  c_id, my_id;
  double myc;

  int psp;

  int run, nrun;


  nrun = GROUPING;
  if (blockIdx.x == careful)
    nrun = nruns;

  /* Set the partial sum to zero (this used to be done in the inner 
     loop, but now we might carry it over loops) */
  unsigned int tid = threadIdx.x;
  unsigned int wid = (threadIdx.x / 32);
  unsigned int quad = wid * 7;
  unsigned int lid = tid % 32;

  int j = lid / 5;		// column index
  int i = lid - j * 5;		// row index

  myc = 0.0;

  /* load the first four parameter sets */
  psp = 7 * (blockIdx.x * GROUPING);

  param_stack += psp + lid;

  int param_r = 0;
  if (wid < nrun + wid)
    {
      param_r = __ldg (param_stack);
    }

  //for (run = 0; run < nrun; run+=4) {
  for (run = wid; run < nrun + wid; run += 4)
    {
      //int param_r = 0;

      //param_r = __ldg(param_stack+psp+lid);
      //param_r = __ldg(param_stack);

      // load matrix elements
      int srcA = __shfl (param_r, quad + 3) - 1;
      int srcB = __shfl (param_r, quad + 4) - 1;

      double mya = 0.0;
      double myb = 0.0;
      //if( (run + wid < nrun) && (lid < 25)) {
      if ((run < nrun) && (lid < 25))
	{
	  mya = a_data[srcA + lid];
	  myb = b_data[srcB + lid];
	}
      // initialization and first product
      double tmpA = dbl_shfl (mya, i);
      double tmpB = dbl_shfl (myb, j * 5);
      myc += tmpA * tmpB;

      tmpA = dbl_shfl (mya, i + 5);
      tmpB = dbl_shfl (myb, j * 5 + 1);
      myc += tmpA * tmpB;

      tmpA = dbl_shfl (mya, i + 2 * 5);
      tmpB = dbl_shfl (myb, j * 5 + 2);
      myc += tmpA * tmpB;

      tmpA = dbl_shfl (mya, i + 3 * 5);
      tmpB = dbl_shfl (myb, j * 5 + 3);
      myc += tmpA * tmpB;

      tmpA = dbl_shfl (mya, i + 4 * 5);
      tmpB = dbl_shfl (myb, j * 5 + 4);
      myc += tmpA * tmpB;


      bool flush_c = true;


      // if we are in the last quadrant, so load the new parameter
      // set now
      int c_loc = __shfl (param_r, quad + 5) - 1;
      //c_id = __shfl(param_r, quad+6) - 1;

      //if(run + wid >= nrun) flush_c = false;
      if (run >= nrun)
	flush_c = false;


      if (wid == 1 || wid == 3)
	{
//          if( (wid && 0x1) == 1 ){
	  if (c_loc == __shfl (param_r, quad + 5 - 7) - 1)
	    {
	      c_s[tid] = myc;
	      flush_c = false;
	      myc = 0.0;
//             if(lid == 0) printf("reducing result in smem\n");
	    }
	}

      syncthreads ();

      if (wid == 0 || wid == 2)
	{
//          if((wid && 0x1) == 0){
	  if (c_loc == __shfl (param_r, quad + 5 + 7) - 1)
	    myc += c_s[tid + 32];
	  if (wid == 2)
	    {
	      if (c_loc == __shfl (param_r, quad + 5 - 14) - 1)
		{
		  c_s[tid] = myc;
		  flush_c = false;
		  myc = 0.0;
		}
	    }
	}

      syncthreads ();

      if (wid == 0)
	{
	  if (c_loc == __shfl (param_r, quad + 5 + 14) - 1)
	    myc += c_s[tid + 64];
	}

      /* Only update c_date if we are in the last iteration, or if the next C-block will be
         different to this C-block */
      /* param_stack[psp+6] is the current C-block ID, so adding 7 means that param_stack[psp+6+7]
         should be the next C-block ID */

      param_stack += 28;
      param_r = __ldg (param_stack);

      if (wid == 0 && run + 4 < nrun + wid)
	{
	  if (c_loc == __shfl (param_r, quad + 5) - 1)
	    flush_c = false;
	}

      if (flush_c)
	{


/*
             if (lid == 0) {
                // TODO: make shure 2^16 is enough as lock_offset stepsize
                 my_id = lock_offset + 4*blockIdx.x  + wid + 1;
                 int lock_owner = 0;
                    while ((lock_owner != my_id))
                         lock_owner = atomicCAS (&(c_locks[c_id]), 0, my_id);
             }
*/

	  if (lid < 25)
	    {
	      //c_data[c_loc+lid] += myc;
	      atomicAdd (&c_data[c_loc + lid], myc);
	    }

	  /* Release the lock on the C block. */
//             if (lid == 0) {
//                  c_locks[c_id] = 0;
//             }
	  /* If we have another C-block then we need to reset 
	     our partial sum to zero for the new C-block */
	  myc = 0.0;

	}

      //psp += 28;
      //param_stack += 28;
    }

};


