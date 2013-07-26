/******************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations
 *  Copyright (C) 2000 - 2013 the CP2K developers group
 *****************************************************************************/

// DBCSR_KERNEL datatype=dbcsr_type_real_8, homogeneous_only=True

#include "dbcsr_kernel.h"
#include "stack_mm_mnk_d.h"

extern __shared__ double cache[];

__global__ void
stack_mm_mnk_d (const int *__restrict__ param_stack,
		const int careful, const int nruns,
		const int m, const int n, const int k,
		const int liter,
		const double *__restrict__ a_data,
		const double *__restrict__ b_data,
		double *__restrict__ c_data)
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

  //int lock_owner, c_id, my_id;
  const int mn = m * n;
  const int mk = m * k;
  const int kn = n * k;
  const int r = threadIdx.x % m;
  const int c = threadIdx.x / m;
  int l, i;
  double myc;
  const double *__restrict__ buff_l, *__restrict__ buff_r;

  int psp, c_loc;

  int run, nrun;

  double *buff;

  buff = (double *) cache;
  buff_l = buff;
  buff_r = &(buff[mk]);

  nrun = GROUPING;
  if (blockIdx.x == careful)
    nrun = nruns;

  /* Set the partial sum to zero (this used to be done in the inner loop, but now we might carry it over loops) */
  myc = 0.0l;

  for (run = 0; run < nrun; run++)
    {
      psp = 7 * (blockIdx.x * GROUPING + run);

/*
      for(int i = threadIdx.x; i < mk; i += blockDim.x){
          buff[i] = a_data[param_stack[psp + 3] - 1 + i];
      }

      for(int i = threadIdx.x; i < kn; i += blockDim.x){
          buff[mk + i] = b_data[param_stack[psp + 4] - 1 + i];
      }
*/

      for (l = 0; l <= liter; l++)
	{
	  i = threadIdx.x + blockDim.x * l;
	  if (i < mk)
	    buff[i] = a_data[param_stack[psp + 3] - 1 + i];
	  if (i < kn)
	    buff[mk + i] = b_data[param_stack[psp + 4] - 1 + i];
	}

      syncthreads ();

      /* Do actual multiplication. */
      if (threadIdx.x < mn)
	{

	  for (l = 0; l < k; l++)
	    {
	      myc = myc + buff_l[l * m + r] * buff_r[c * k + l];
	    }

	}

      /* Only update c_date if we are in the last iteration, or if the next C-block will be
         different to this C-block */
      /* param_stack[psp+6] is the current C-block ID, so adding 7 means that param_stack[psp+6+7]
         should be the next C-block ID */
      if (run == nrun - 1
	  || param_stack[psp + 6] - 1 != param_stack[psp + 6 + 7] - 1)
	{
	  /* Add results to global C block. */
	  c_loc = param_stack[psp + 5] - 1;
	  if (threadIdx.x < mn)
	    atomicAdd (&c_data[c_loc + threadIdx.x], myc);


	  /* If we have another C-block then we need to reset our partial sum to zero for the new C-block */
	  myc = 0.0l;
	}

      syncthreads ();

    }


};

//==============================================================================
int launch_stack_mm_mnk_d(int *param_stack, int stack_size, cudaStream_t stream,
    int m_max, int n_max, int k_max,
    double *a_data, double *b_data, double *c_data){

     int shared_size = (m_max * k_max + k_max * n_max) * sizeof (double);
     if (shared_size > devProperties.sharedMemPerBlock)
          return 4;

     int maxt = m_max * n_max;
     if (maxt > devProperties.maxThreadsPerBlock)
       return 3;

     //int mn = m_max * n_max;
     int mk = m_max * k_max;
     int nk = n_max * k_max;
     int careful = (stack_size / GROUPING);
     int nruns = stack_size - careful * GROUPING;
     int maxb = MAX (mk, nk);
     int liter = (maxb - 1) / maxt;

     stack_mm_mnk_d <<< (stack_size + GROUPING - 1) / GROUPING, maxt,
                shared_size, stream >>> (param_stack, careful, nruns, m_max,
                    n_max, k_max, liter, a_data, b_data, c_data);

     return(0);
}
