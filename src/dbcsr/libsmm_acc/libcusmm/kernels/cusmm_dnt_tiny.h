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

// optimized for  that fit into shared memory
template < int m,  int n,  int k, int splitThread, int blockdim, int grouping, int minblocks>
__global__ void
__launch_bounds__(blockdim, minblocks)
cusmm_dnt_tiny(const int* __restrict__ param_stack, int careful, int nruns,
     double* a_data, double* b_data, double* c_data){

  const int mn = m * n;
  const int mk = m * k;
  const int kn = n * k;

  const int buf_sz = 2*(mk + kn);

  const int rmax = m;

  //const int r = threadIdx.x % rmax;
  const int c = threadIdx.x / rmax;
  const int r = threadIdx.x - c  * rmax;

  double myc;

  int psp, c_loc;

  const int paramNum = 3;

  __shared__ double buff[buf_sz];
  __shared__ int param_stack_s[paramNum*grouping];

  double *buff_l1 =   buff;
  double *buff_r1 = &(buff[mk]);
  double *buff_l2 = &(buff[mk+kn]);
  double *buff_r2 = &(buff[2*mk+kn]);

  int nrun = grouping;

  psp = 7 * (blockIdx.x * grouping);
  if (blockIdx.x == careful)
    nrun = nruns;

  /* Set the partial sum to zero */
  myc = 0.0;

  // load and pack task data into smem
  for(int i = threadIdx.x; i < 7*nrun; i += blockDim.x){
    int p_tmp = __ldg(param_stack + psp + i);
    if ((i % 7 > 2) && (i%7 < 6)) 
      param_stack_s[(i / 7)*paramNum + i % 7 - 3] = p_tmp-1;
  }

  syncthreads();
  
  unsigned int flushMap = 1;
  for(int i=nrun-2; i>=0; i--){
    unsigned int tmp1 = (flushMap << 1);
    unsigned int tmp2 = (param_stack_s[paramNum*i + 2] != param_stack_s[paramNum*i + 2 + paramNum]);
    flushMap = tmp1 | tmp2;
  }

//  printf("flushMap = %d\n", flushMap);
 for (int run = 0; run < nrun; run++)
  {
    psp = run*paramNum;

    if(threadIdx.x >= splitThread){
      int srcA = param_stack_s[psp];
      int srcB = param_stack_s[psp + 1];

      // load a and b matrix into smem

      if( m == n) {
         for(int i = threadIdx.x-splitThread; i < mk; i += (blockDim.x-splitThread)){
            buff_l1[i] = __ldg(a_data +srcA + i);
            buff_r1[i] = __ldg(b_data +srcB + i);
         }
      }  else  { 
  
         for(int i = threadIdx.x-splitThread; i < mk; i += (blockDim.x-splitThread)){
            buff_l1[i] = __ldg(a_data +srcA + i);
         }

         for(int i = threadIdx.x-splitThread; i < kn; i += (blockDim.x-splitThread)){
            buff_r1[i] = __ldg(b_data +srcB + i);
         }
      }
    }
 //   printf("%d arrived \n", threadIdx.x);
    syncthreads ();

    if(threadIdx.x < mn){

      /* Do actual multiplication. */
      for (int l = 0; l < k; l++) {
          //myc += buff_l1[l * m + r] * buff_r1[c * k + l];
          myc += buff_l1[l * m + r] * buff_r1[l*n + c];
      } 

      //if ( run== nrun-1  || flushMap & (1 << run)) {
      if ( flushMap & (1 << run)) {
        c_loc = param_stack_s[psp + 2];
//        if(threadIdx.x < mn) {
         atomicAdd (&c_data[c_loc + threadIdx.x], myc);
//        }
        myc  =0.0;
      }
    }

    run++;
    psp = run*paramNum;
    
    if(threadIdx.x >= splitThread && run < nrun){
        int srcA = param_stack_s[psp];
        int srcB = param_stack_s[psp + 1];

        if( m == n) {
           for(int i = threadIdx.x-splitThread; i < mk; i += (blockDim.x-splitThread)){
              buff_l2[i] = __ldg(a_data +srcA + i);
              buff_r2[i] = __ldg(b_data +srcB + i);
           }
        }  else  {

           for(int i = threadIdx.x-splitThread; i < mk; i += (blockDim.x-splitThread)){
              buff_l2[i] = __ldg(a_data +srcA + i);
           }

           for(int i = threadIdx.x-splitThread; i < kn; i += (blockDim.x-splitThread)){
              buff_r2[i] = __ldg(b_data +srcB + i);
           }
        }
    }
     
    syncthreads();

    if(threadIdx.x < mn && run < nrun){

      /* Do actual multiplication. */
      for (int l = 0; l < k; l++) {
          //myc += buff_l2[l * m + r] * buff_r2[c * k + l];
          myc += buff_l2[l * m + r] * buff_r2[l*n+c];
      }

      if ( flushMap & (1 << run)) {
        c_loc = param_stack_s[psp + 2];
        atomicAdd(&c_data[c_loc + threadIdx.x], myc);
        myc  =0.0;
      }
    }

  }

}

