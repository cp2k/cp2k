/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2013 the CP2K developers group                      *
 *  Authors: Peter Messmer <pmessmer@nvidia.com>,                            *
 *           Nikolay Markovskiy <nmarkovskiy@nvidia.com>                     *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include "cusmm_common.h"

#define GROUPING 16

template < int m,  int n,  int k, int N>
__global__ void
cusmm_kernel_default(const int *__restrict__ param_stack,
                     const int careful, const int nruns,
                     const int p_m, const int p_n, const int p_k,
                     const double *__restrict__ a_data,
                     const double *__restrict__ b_data,
                     double *__restrict__ c_data)
{
  const int mn = m * n;
  const int mk = m * k;
  const int kn = n * k;


  const int buf_sz = (mk + kn < mn)? mn : mk + kn;

  const int cmax = (n % N == 0)?  n/N: n/N+1;
  const int rmax = (m % N == 0)?  m/N: m/N+1;

  const int r = threadIdx.x % rmax;
  const int c = threadIdx.x / rmax;


  double myc[N*N];

  const double *__restrict__ buff_l, *__restrict__ buff_r;

  int psp, c_loc;

  __shared__ double buff[buf_sz];
  __shared__ int param_stack_s[4*GROUPING];

  buff_l = buff;
  buff_r = &(buff[mk]);

  int nrun = GROUPING;

  psp = 7 * (blockIdx.x * GROUPING);
  if (blockIdx.x == careful)
    nrun = nruns;


  /* Set the partial sum to zero */
  for (int i = 0; i < N*N; i++)
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

    int srcA = param_stack_s[psp];     
    int srcB = param_stack_s[psp + 1]; 

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
      for (int l = 0; l < k; l++) {
        for (int i = 0; i < N; i++)
          for (int j = 0; j < N; j++)
              myc[N*i+j] = myc[N*i+j] + buff_l[l * m + N*r+i] * buff_r[(N*c+j) * k + l];
      }
    }

    if (run == nrun - 1
        || param_stack_s[psp + 3] != param_stack_s[psp + 3 + 4] )
    {

      // decompress results to a:
      syncthreads();
      if (c < cmax  && r < rmax)
      {
        for (int i = 0; i < N; i++)
          for (int j = 0; j < N; j++)
            if (N*r+i < m && N*c+j < n)
              buff[(N*c+j) * m + N*r+i] = myc[N*i+j];
      }
      syncthreads();
      /* Add results to global C block. */
      c_loc = param_stack_s[psp + 2];
      for(int i = threadIdx.x; i < mn; i += blockDim.x)
        atomicAdd (&c_data[c_loc + i], buff[i]);

      for (int i = 0; i < N*N; i++)
        myc[i] = 0.0;
    }
  }

}

