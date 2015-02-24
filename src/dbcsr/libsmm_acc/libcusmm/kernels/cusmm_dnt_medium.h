/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015  CP2K developers group                         *
 *****************************************************************************/

/*****************************************************************************
 *  Authors: Peter Messmer <pmessmer@nvidia.com>,                            *
 *           Nikolay Markovskiy <nmarkovskiy@nvidia.com>                     *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include "cusmm_common.h"

#define PARAMS_IN_SMEM_OPT

template < int m,  int n,  int k, int M, int N, int bsize, int grouping, int minblocks>
__global__ void
__launch_bounds__(bsize, minblocks)
cusmm_dnt_medium(const int* __restrict__ param_stack, int  careful,  int nruns,
     double* a_data, double* b_data, double* c_data){

  const int mn = m * n;
  const int mk = m * k;
  const int kn = n * k;

  const int buf_sz = (mk + kn < mn)? mn: mk + kn;

  const int cmax = (n % N == 0)?  n/N: n/N+1;
  const int rmax = (m % M == 0)?  m/M: m/M+1;

  const int c = threadIdx.x / rmax;
  const int r = threadIdx.x - rmax*c;

  double myc[N*M];

  const int load_unroll_factor_1 = mk / bsize + 1; //6;
  const int load_unroll_factor_2 = kn / bsize + 1; //3;

  const int n_mkloads = mk / (load_unroll_factor_1*bsize);
  const int n_knloads = kn / (load_unroll_factor_2*bsize);

  const int mkloads_remained = (mk - n_mkloads*load_unroll_factor_1*bsize) / bsize;
  const int knloads_remained = (kn - n_knloads*load_unroll_factor_2*bsize) / bsize;

  const int mkloads_tail = ((mkloads_remained + n_mkloads*load_unroll_factor_1)*bsize == mk)? 0:1;
  const int knloads_tail = ((knloads_remained + n_knloads*load_unroll_factor_2)*bsize == kn)? 0:1;

  const int  m_loads_to_finish = n_mkloads*(load_unroll_factor_1*bsize) + mkloads_remained*bsize;
  const int  n_loads_to_finish = n_knloads*(load_unroll_factor_2*bsize) + knloads_remained*bsize;
  const int left_to_finish_1 =  m_loads_to_finish + threadIdx.x;
  const int left_to_finish_2 =  n_loads_to_finish + threadIdx.x;

  double lba_r[load_unroll_factor_1 + mkloads_remained + mkloads_tail];
  double lbb_r[load_unroll_factor_2 + knloads_remained + knloads_tail];

  const double *__restrict__ buff_l, *__restrict__ buff_r;

  int psp;

  bool is_loaded = false;

  __shared__ double buff[buf_sz];

#ifdef PARAMS_IN_SMEM_OPT
  const int paramNum = 4;
//  const int paramNum = 3;
  __shared__ int param_stack_s[paramNum*grouping];
#endif

  buff_l = buff;
  buff_r = &(buff[mk]);

  int nrun = grouping;

  psp = 7 * (blockIdx.x * grouping);
  //if (blockIdx.x == gridDim.x-1)
  if (blockIdx.x == careful)
    nrun = nruns;


  /* Set the partial sum to zero */
  for (int i = 0; i < M*N; i++)
    myc[i] = 0.0;

#if defined(PARAMS_IN_SMEM_OPT) || defined(FLUSH_OPT)
  psp = 7 * (blockIdx.x * grouping);
#endif

#ifdef PARAMS_IN_SMEM_OPT
  // load and pack task data into smem
  for(int i = threadIdx.x; i < 7*nrun; i += blockDim.x){
    int p_tmp = __ldg(param_stack + psp + i);
    //if ((i % 7 > 2) && (i % 7 < 6))
    if ((i % 7 > 2))
      param_stack_s[(i / 7)*paramNum + i % 7 - 3] = p_tmp-1;
  }
  syncthreads ();
#endif

#ifdef FLUSH_OPT
  unsigned int flushMap = 0;
  for(int i = 0; i < nrun-1; i ++){
    if (param_stack_s[paramNum*i + 2] != param_stack_s[paramNum*i + 2 + paramNum])
      flushMap |= (1<<i); // Storing positions when to flush
  }

  flushMap |= (1<< (nrun-1));
#endif

  for (int run = 0; run < nrun; run++)
  {
#ifdef PARAMS_IN_SMEM_OPT
  psp = run*paramNum;
#else
  psp = 7 * (bsize * grouping + run);
#endif

    int srcA, srcB;

    if (!is_loaded)
    {
#ifdef PARAMS_IN_SMEM_OPT
      srcA = param_stack_s[psp];
      srcB = param_stack_s[psp+1];
#else
      srcA = param_stack[psp+3]-1;
      srcB = param_stack[psp+4]-1;
#endif
#pragma unroll 3
    for(int i = threadIdx.x; i < n_mkloads*(load_unroll_factor_1*bsize); i += load_unroll_factor_1*bsize)
    {
#pragma unroll
     for (int l = 0; l < load_unroll_factor_1; l++)
      {
        lba_r[l] = __ldg(a_data + srcA + i + l*bsize);
        if (m == n)
          lbb_r[l] = __ldg(b_data +srcB + i + l*bsize);
      }
    }

#pragma unroll 3
    for(int i = n_mkloads*(load_unroll_factor_1*bsize) + threadIdx.x; i < m_loads_to_finish; i += mkloads_remained*bsize)
    {
#pragma unroll
      for (int l = 0; l < mkloads_remained; l++)
      {
        lba_r[l] = __ldg(a_data + srcA + i + l*bsize);
        if (m == n)
          lbb_r[l] = __ldg(b_data +srcB + i +l*bsize );
      }
    }

    if (left_to_finish_1 < mk)
    {
      lba_r[ load_unroll_factor_1 + mkloads_remained] = __ldg(a_data + srcA + left_to_finish_1);
      if (m == n)
        lbb_r[load_unroll_factor_2 + knloads_remained] = __ldg(b_data +srcB + left_to_finish_2);
    }

    if (m != n)
    {
#pragma unroll 3
      for(int i = threadIdx.x; i < n_knloads*(load_unroll_factor_2*bsize); i += load_unroll_factor_2*bsize)
#pragma unroll
        for (int l = 0; l < load_unroll_factor_2; l++)
          lbb_r[l] = __ldg(b_data +srcB + i + l*bsize);

#pragma unroll 3
      for(int i = n_knloads*(load_unroll_factor_2*bsize) + threadIdx.x; i < n_loads_to_finish; i += knloads_remained*bsize)
#pragma unroll
        for (int l = 0; l < knloads_remained; l++)
          lbb_r[l] = __ldg(b_data +srcB + i +l*bsize );

      if (left_to_finish_2 < kn)
        lbb_r[load_unroll_factor_2 + knloads_remained] = __ldg(b_data +srcB + left_to_finish_2);
    }
    }
    else
     {
       is_loaded = false;
     }

// Copy to smem:
#pragma unroll 3
    for(int i = n_mkloads*(load_unroll_factor_1*bsize) + threadIdx.x; i < m_loads_to_finish; i += mkloads_remained*bsize)
    {
#pragma unroll
      for (int l = 0; l < mkloads_remained; l++)
      {
        buff[i + l*bsize] = lba_r[l];
        if (m == n)
          buff[i+mk + l*bsize] = lbb_r[l];
      }
    }

#pragma unroll 3
    for(int i = threadIdx.x; i < n_mkloads*(load_unroll_factor_1*bsize); i += load_unroll_factor_1*bsize)
    {
#pragma unroll
      for (int l = 0; l < load_unroll_factor_1; l++)
      {
        buff[i + l*bsize] = lba_r[l];
        if (m==n)
          buff[i + mk + l*bsize] = lbb_r[l];
      }
    }

    if (left_to_finish_1 < mk)
    {
      buff[left_to_finish_1] = lba_r[load_unroll_factor_1 + mkloads_remained];
      if (m==n)
        buff[mk + left_to_finish_2] = lbb_r[load_unroll_factor_2 + knloads_remained];
    }

    if (m != n)
    {
#pragma unroll 3
      for(int i = n_knloads*(load_unroll_factor_2*bsize) + threadIdx.x; i < n_loads_to_finish; i += knloads_remained*bsize)
#pragma unroll
        for (int l = 0; l < knloads_remained; l++)
          buff[i+mk + l*bsize] = lbb_r[l];
#pragma unroll 3
      for(int i = threadIdx.x; i < n_knloads*(load_unroll_factor_2*bsize); i += load_unroll_factor_2*bsize)
#pragma unroll
       for (int l = 0; l < load_unroll_factor_2; l++)
          buff[i + mk + l*bsize] = lbb_r[l];
      if (left_to_finish_2 < kn)
        buff[mk + left_to_finish_2] = lbb_r[load_unroll_factor_2 + knloads_remained];
    }

    syncthreads ();

    int next_run = run+1;
    if (next_run >= nrun)
     {
        is_loaded=true;
     }

    if (!is_loaded || (run == 0 && nrun>1))
//     if (!is_loaded || run == 0 )
      {
#ifdef PARAMS_IN_SMEM_OPT
        int next_psp = next_run*paramNum;
#else
        next_psp = 7 * (bsize * grouping + next_run);
#endif

#ifdef PARAMS_IN_SMEM_OPT
        srcA = param_stack_s[next_psp];
        srcB = param_stack_s[next_psp+1];
#else
        srcA = param_stack[next_psp+3]-1;
        srcB = param_stack[next_psp+4]-1;
#endif

#pragma unroll 3
    for(int i = threadIdx.x; i < n_mkloads*(load_unroll_factor_1*bsize); i += load_unroll_factor_1*bsize)
    {
#pragma unroll
      for (int l = 0; l < load_unroll_factor_1; l++)
      {
        lba_r[l] = __ldg(a_data + srcA + i + l*bsize);
        if (m == n)
          lbb_r[l] = __ldg(b_data +srcB + i + l*bsize);
      }
    }

#pragma unroll 3
    for(int i = n_mkloads*(load_unroll_factor_1*bsize) + threadIdx.x; i < m_loads_to_finish; i += mkloads_remained*bsize)
    {
#pragma unroll
      for (int l = 0; l < mkloads_remained; l++)
      {
        lba_r[l] = __ldg(a_data + srcA + i + l*bsize);
        if (m == n)
          lbb_r[l] = __ldg(b_data +srcB + i +l*bsize );
      }
    }

    if (left_to_finish_1 < mk)
    {
      lba_r[ load_unroll_factor_1 + mkloads_remained] = __ldg(a_data + srcA + left_to_finish_1);
      if (m == n)
        lbb_r[load_unroll_factor_2 + knloads_remained] = __ldg(b_data +srcB + left_to_finish_2);
    }

    if (m != n)
    {
#pragma unroll 3
      for(int i = threadIdx.x; i < n_knloads*(load_unroll_factor_2*bsize); i += load_unroll_factor_2*bsize)
#pragma unroll
        for (int l = 0; l < load_unroll_factor_2; l++)
          lbb_r[l] = __ldg(b_data +srcB + i + l*bsize);

#pragma unroll 3
      for(int i = n_knloads*(load_unroll_factor_2*bsize) + threadIdx.x; i < n_loads_to_finish; i += knloads_remained*bsize)
#pragma unroll
        for (int l = 0; l < knloads_remained; l++)
          lbb_r[l] = __ldg(b_data +srcB + i +l*bsize );

      if (left_to_finish_2 < kn)
        lbb_r[load_unroll_factor_2 + knloads_remained] = __ldg(b_data +srcB + left_to_finish_2);
    }
    is_loaded = true;
      }

    /* Do actual multiplication. */
    if (c < cmax  && r < rmax)
    {
      for (int l = 0; l < k; l++) {
        for (int i = 0; i < M; i++)
          for (int j = 0; j < N; j++)
              myc[N*i+j] = myc[N*i+j] + buff_l[l * m + M*r+i] * buff_r[l*n+N*c+j];
      //        myc[N*i+j] = myc[N*i+j] + buff_l[l * m + M*r+i] * buff_r[(N*c+j) * k + l];
      }
    }

    if (
#if defined(FLUSH_OPT)
      flushMap & (1 << run)
#else
#if defined(PARAMS_IN_SMEM_OPT)
//      run == nrun-1 || param_stack_s[psp + 2] != param_stack_s[psp + 2 + paramNum]
    run == nrun-1 || param_stack_s[psp + 3] != param_stack_s[psp + 3 + paramNum]
#else
      param_stack[psp + 6] -1 != param_stack[psp + 6 + 7] -1
#endif
#endif
    )
    {

      // decompress results to a:
      /* Add results to global C block. */
#ifdef PARAMS_IN_SMEM_OPT
     int  c_loc = param_stack_s[psp+2];
#else
     int  c_loc = param_stack[psp + 5] - 1;
#endif

     if (M > 1 || N > 1)
     {
       syncthreads();
       if (c < cmax  && r < rmax)
       {
         for (int i = 0; i < M; i++)
           for (int j = 0; j < N; j++)
             if (M*r+i < m && N*c+j < n)
             {
               buff[(N*c+j) * m + M*r+i] = myc[N*i+j];
             }
       }
       syncthreads();
       for(int i = threadIdx.x; i < mn; i += bsize) {
         atomicAdd (&c_data[c_loc + i], buff[i]); //c_data[c_loc + i] += buff[i];
       }
     }
     else
     {
       for(int i = threadIdx.x; i < m*n; i += bsize)
         atomicAdd (&c_data[c_loc + i], myc[0]);
     }
     for (int i = 0; i < N*M; i++)
       myc[i] = 0.0;
    }
    syncthreads ();
  }

}
