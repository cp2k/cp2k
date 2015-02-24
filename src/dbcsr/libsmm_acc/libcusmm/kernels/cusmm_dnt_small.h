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

// optimized for 13x13x13 products
template < int m,  int n,  int k, int M, int N, int blockdim, int grouping, int minblocks>
__global__ void
__launch_bounds__(blockdim, minblocks)
cusmm_dnt_small(const int* __restrict__ param_stack, int  careful,  int nruns,
     double* a_data, double* b_data, double* c_data){

  const int mn = m * n;
  const int mk = m * k;
  const int kn = n * k;

  const int buf_sz = (mk + kn < mn)? mn: mk + kn;

  const int cmax = (n % N == 0)?  n/N: n/N+1;
  const int rmax = (m % M == 0)?  m/M: m/M+1;

  const int c = threadIdx.x / rmax;
  const int r = threadIdx.x - c  * rmax;

  double myc[M*N];


  int psp, c_loc;

  const int paramNum = 3;

  __shared__ double buff[buf_sz];
  __shared__ int param_stack_s[paramNum*grouping];

  const double *buff_l = buff;
  double *buff_r = &(buff[mk]);

  int nrun = grouping;

  psp = 7 * (blockIdx.x * grouping);
  if (blockIdx.x == careful)
    nrun = nruns;

  /* Set the partial sum to zero */
  for (int i = 0; i < M*N; i++)
    myc[i] = 0.0;

  // load and pack task data into smem
  for(int i = threadIdx.x; i < 7*nrun; i += blockDim.x){
    int p_tmp = __ldg(param_stack + psp + i);
    if ((i % 7 > 2) && (i%7 < 6)) 
      param_stack_s[(i / 7)*paramNum + i % 7 - 3] = p_tmp-1;
  }

  
/*
  unsigned int flushMap = 1;
//  for(int i=grouping-2; i>=0; i--){
  for(int i=nrun-2; i>=0; i--){
    unsigned int tmp1 = (flushMap << 1);
    unsigned int tmp2 = (param_stack_s[paramNum*i + 2] != param_stack_s[paramNum*i + 2 + paramNum]);
    flushMap = tmp1 | tmp2;
  }
*/

 for (int run = 0; run < nrun; run++)
  {
    psp = run*paramNum;
    syncthreads ();

    int srcA = param_stack_s[psp];
    int srcB = param_stack_s[psp + 1];

    // load a and b matrix into smem

    if( m == n) {
       for(int i = threadIdx.x; i < mk; i += blockDim.x){
          buff[i] = __ldg(a_data +srcA + i);
          //buff[i] = a_data[srcA + i];
          buff_r[i] = __ldg(b_data +srcB + i);
          //buff_r[i] = b_data[srcB + i];
       }
    }  else  { 

       for(int i = threadIdx.x; i < mk; i += blockDim.x){
          buff[i] = __ldg(a_data +srcA + i);
          //buff[i] = a_data[srcA + i];
       }

       for(int i = threadIdx.x; i < kn; i += blockDim.x){
          //buff_r[i] = __ldg(b_data +srcB + i);
          buff_r[i] = b_data[srcB + i];
       }
    }

    syncthreads ();

    /* Do actual multiplication. */
    if (c < cmax )
    {
      for (int l = 0; l < k; l++) {
        for (int i = 0; i < N; i++)
          for (int j = 0; j < M; j++){
              //myc[M*i+j] += buff_l[l * m + M*r+j] * buff_r[(N*c+i) * k + l];
              myc[M*i+j] += buff_l[l * m + M*r+j] * buff_r[l*n + N*c + i];
          } 
      }
    
    }


//    if ( flushMap & (1 << run))
    if (run == nrun-1 || (param_stack_s[psp + 2] != param_stack_s[ psp + 2 + paramNum]))
    {

      if( M > 1 || N > 1) {
         // decompress results to a:
        syncthreads();
        if (c < cmax  && r < rmax) {
          for (int i = 0; i < N; i++)
            for (int j = 0; j < M; j++){
              if (M*r+j < m && N*c+i < n){
                buff[(N*c+i) * m + M*r+j] = myc[M*i+j];
                myc[M*i+j] = 0.0;
              }
            }
        }
        syncthreads();
        /* Add results to global C block. */
        c_loc = param_stack_s[psp + 2]; 
        for(int i = threadIdx.x; i < mn; i += blockDim.x)
          atomicAdd (&c_data[c_loc + i], buff[i]);
      } else {
        c_loc = param_stack_s[psp + 2]; 
        for(int i = threadIdx.x; i < mn; i += blockDim.x)
          atomicAdd (&c_data[c_loc + i], myc[0]);
        myc[0]= 0.0;
      }
    }
  }

}

