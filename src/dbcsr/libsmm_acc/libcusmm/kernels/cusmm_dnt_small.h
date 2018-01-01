/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2018  CP2K developers group                         *
 *****************************************************************************/

/*****************************************************************************
 *  Authors: Peter Messmer <pmessmer@nvidia.com>,                            *
 *           Nikolay Markovskiy <nmarkovskiy@nvidia.com>                     *
 *  Simplified, cleanup & refacture: Based on 'cusmm_dnt_small.h'            *
 *           Andreas Gloess <andreas.gloess@chem.uzh.ch>                     *
 *                                                                           *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include "cusmm_common.h"

/*
gridDim.x = (stack_size + (grouping-1))/grouping
blockIdx.x = {0, ..., gridDim.x-1}
blockDim.x = threads
threadIdx.x = {0, ..., threads-1}
*/

// optimized for 13x13x13 products
template <int m, int n, int k, int M, int N, int threads, int grouping, int minblocks>
__global__ void
__launch_bounds__(threads, minblocks)
cusmm_dnt_small(const int* __restrict__ param_stack, int stack_size,
     const double* __restrict__ a_data, const double* __restrict__ b_data, double* c_data){

  const int mn = m * n;
  const int mk = m * k;
  const int kn = n * k;

  const int bidx = blockIdx.x;
  const int tidx = threadIdx.x;

  const int cmax = (n % N == 0)?  n / N: n / N + 1;
  const int rmax = (m % M == 0)?  m / M: m / M + 1;

  // buff_l and buff_r can overlap in the multiplication step
  // (still better than 'if' in inner loop, see ^ref1) 
  const int buf_tmp = (mk + k * N * cmax < M * rmax * k + 1)? M * rmax * k + 1: mk + k * N * cmax;
  const int buf_sz = (buf_tmp < mn)? mn: buf_tmp;

  const int c = tidx / rmax;
  const int r = tidx - c * rmax;

  const int  npar = 3;
  const int  warp_size = 32;
  const bool need_sync = (mn > warp_size || mk > warp_size || kn > warp_size || threads > warp_size);

  int    nrun, psp;
  double myc[M * N];

  __shared__ int    param_stack_s[npar * grouping];
  __shared__ double buff[buf_sz];

  double *buff_l = buff;
  double *buff_r = &buff[mk];

  nrun = grouping;
  if (((bidx + 1) * grouping) > stack_size) nrun = stack_size % grouping;

  /* Set the partial sum to zero */
#pragma unroll
  for (int i = 0; i < M * N; i++)
    myc[i] = 0.0;

  /* load and pack stack data for current block into smem */
  psp = bidx * npar * grouping;
#pragma unroll
  for (int i = tidx; i < nrun; i += threads){
    param_stack_s[i * npar    ] = __ldg(&param_stack[psp + i * npar    ]) - 1;
    param_stack_s[i * npar + 1] = __ldg(&param_stack[psp + i * npar + 1]) - 1;
    param_stack_s[i * npar + 2] = __ldg(&param_stack[psp + i * npar + 2]) - 1;
  }

  if (need_sync) syncthreads ();

  for (int run = 0; run < nrun; run++){
    /* pointer offsets */
    psp = run * npar;
    int srcA = param_stack_s[psp];
    int srcB = param_stack_s[psp + 1];
    int srcC = param_stack_s[psp + 2]; 

    /* load a and b matrix into smem */
    if ( m == n){
#pragma unroll
      for (int i = tidx; i < mk; i += threads){
         buff_l[i] = __ldg(&a_data[srcA + i]);
         buff_r[i] = __ldg(&b_data[srcB + i]);
      }
    } else { 
#pragma unroll
      for (int i = tidx; i < mk; i += threads){
         buff_l[i] = __ldg(&a_data[srcA + i]);
      }
#pragma unroll
      for (int i = tidx; i < kn; i += threads){
         buff_r[i] = __ldg(&b_data[srcB + i]);
      }
    }

    if (need_sync) syncthreads ();

    /* Do actual multiplication. */
    if (c < cmax && r < rmax){
      for (int l = 0; l < k; l++){
#pragma unroll
        for (int i = 0; i < N; i++){
#pragma unroll
          for (int j = 0; j < M; j++){
//            if ((M * r + j < m) && (N * c + i < n))          // ^ref1
               myc[M * i + j] += buff_l[l * m + M * r + j] * buff_r[l * n + N * c + i];
          }
        }
      }
    }

    /* last loop or C_idx for next stack entry is different */
    if ((run == (nrun - 1)) || (srcC != param_stack_s[psp + 2 + npar])){
      if ((M > 1) || (N > 1)){
        if (need_sync) syncthreads ();

        /* decompress results to buffer */
        if (c < cmax && r < rmax){
#pragma unroll
          for (int i = 0; i < N; i++){
#pragma unroll
            for (int j = 0; j < M; j++){
              if (M * r + j < m && N * c + i < n){
                 buff[(N * c + i) * m + M * r + j] = myc[M * i + j];
                 myc[M * i + j] = 0.0;
              }
            }
          }
        }

        if (need_sync) syncthreads ();

        /* Add results to global C block. */
#pragma unroll
        for (int i = tidx; i < mn; i += threads)
          atomicAdd (&c_data[srcC + i], buff[i]);
      } else {
#pragma unroll
        for (int i = tidx; i < mn; i += threads)
          atomicAdd (&c_data[srcC + i], myc[0]);
        myc[0]= 0.0;
      }
    }

    if (need_sync) syncthreads ();
  }
}

