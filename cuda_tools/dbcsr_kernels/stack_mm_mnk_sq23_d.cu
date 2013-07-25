/******************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations
 *  Copyright (C) 2000 - 2013 the CP2K developers group
 *****************************************************************************/


#include "dbcsr_kernel.h"
#include "stack_mm_mnk_sq23_d.h"

__global__ void
stack_mm_mnk_sq23_d(const int *__restrict__ param_stack,
		     const int careful, const int nruns,
		     const int p_m, const int p_n, const int p_k,
		     const int liter,
		     const double *__restrict__ a_data,
		     const double *__restrict__ b_data,
		     double *__restrict__ c_data, int *__restrict__ c_locks,
		     int lock_offset)
{

  const int mn = 23 * 23;
  const int mk = 23 * 23;
  const int kn = 23 * 23;

  const int cmax = (23 % 2 == 0)?  23/2: 23/2+1;
  const int rmax = (23 % 2 == 0)?  23/2: 23/2+1;

  const int r = threadIdx.x % rmax;
  const int c = threadIdx.x / cmax;

  double myc[2*2];

  const double *__restrict__ buff_l, *__restrict__ buff_r;

  int psp, c_loc;

  __shared__ double buff[mk + kn];
  __shared__ int param_stack_s[4*GROUPING];

  buff_l = buff;
  buff_r = &(buff[mk]);

  int nrun = GROUPING;

  psp = 7 * (blockIdx.x * GROUPING);
  if (blockIdx.x == careful)
    nrun = nruns;


  /* Set the partial sum to zero */
  for (int i = 0; i < 2*2; i++)
    myc[i] = 0.0;

  // load and pack task data into smem
  for(int i = threadIdx.x; i < 7*GROUPING; i += blockDim.x){
    int p_tmp = __ldg(param_stack + psp + i);
    if (i % 7 > 2)
      param_stack_s[(i / 7)*4 + i % 7 - 3] = p_tmp-1;
  }

  for (int run = 0; run < nrun; run++)
  {
    psp = run*4;
    syncthreads ();

    int srcA = param_stack_s[psp];     //__shfl (param_r, 3) - 1;
    int srcB = param_stack_s[psp + 1]; //__shfl (param_r, 4) - 1;

    // load a and b matrix into smem
    for(int i = threadIdx.x; i < mk; i += blockDim.x){
        buff[i] = __ldg(a_data + srcA + i);
    }

    for(int i = threadIdx.x; i < kn; i += blockDim.x){
        buff[mk + i] = __ldg(b_data +srcB + i);
    }
    syncthreads ();

    /* Do actual multiplication. */
    if (c < cmax  && r < rmax)
    {
      for (int l = 0; l < 23; l++) {
        for (int i = 0; i < 2; i++)
          for (int j = 0; j < 2; j++)
              myc[2*i+j] = myc[2*i+j] + buff_l[l * 23 + 2*r+i] * buff_r[(2*c+j) * 23 + l];
      }
    }

    if (run == nrun - 1
        || param_stack_s[psp + 3] != param_stack_s[psp + 3 + 4] )
    {

      // decompress results to a:
      syncthreads();
      if (c < cmax  && r < rmax)
      {
        for (int i = 0; i < 2; i++)
          for (int j = 0; j < 2; j++)
            if (2*r+i < 23 && 2*c+j < 23)
              buff[(2*c+j) * 23 + 2*r+i] = myc[2*i+j];
      }
      syncthreads();
      /* Add results to global C block. */
      //c_loc = param_stack_s[psp + 3]; //__shfl (param_r, 5) - 1;
      c_loc = param_stack_s[psp + 2]; //__shfl (param_r, 5) - 1;
      for(int i = threadIdx.x; i < mn; i += blockDim.x)
        atomicAdd (&c_data[c_loc + i], buff[i]);

      for (int i = 0; i < 2*2; i++)
        myc[i] = 0.0;
    }
    //syncthreads ();
  }

};

