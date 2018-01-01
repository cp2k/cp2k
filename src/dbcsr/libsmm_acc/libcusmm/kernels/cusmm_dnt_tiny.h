/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2018  CP2K developers group                         *
 *****************************************************************************/

/*****************************************************************************
 *  Authors: Peter Messmer <pmessmer@nvidia.com>,                            *
 *           Nikolay Markovskiy <nmarkovskiy@nvidia.com>                     *
 *  Simplified, cleanup & refacture: Based on 'cusmm_dnt_tiny.h'             *
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

// optimized for matrices that fit into shared memory
template <int m, int n, int k, int threads, int grouping, int minblocks>
__global__ void
__launch_bounds__(threads, minblocks)
cusmm_dnt_tiny(const int* __restrict__ param_stack, int stack_size,
     const double* __restrict__ a_data, const double* __restrict__ b_data, double* c_data){

  const int mn = m * n;
  const int mk = m * k;
  const int kn = n * k;

  const int bidx = blockIdx.x;
  const int tidx = threadIdx.x;

  const int c = tidx / m;
  const int r = tidx - c * m;

  const int  npar      = 3;
  const int  warp_size = 32;
  const bool need_sync = (mn > warp_size || mk > warp_size || kn > warp_size || threads > warp_size);

  int    psp, nrun;
  double myc;

  __shared__ int    param_stack_s[npar * grouping];
  __shared__ double buff_a[mk];
  __shared__ double buff_b[kn];

  nrun = grouping;
  if (((bidx + 1) * grouping) > stack_size) nrun = stack_size % grouping;

  /* Set the partial sum to zero */
  myc = 0.0;

  /* load and pack stack data for current block into smem */
  psp = bidx * npar * grouping;
#pragma unroll
  for (int i = tidx; i < nrun; i += threads){
    param_stack_s[i * npar    ] = __ldg(&param_stack[psp + i * npar    ]) - 1;
    param_stack_s[i * npar + 1] = __ldg(&param_stack[psp + i * npar + 1]) - 1;
    param_stack_s[i * npar + 2] = __ldg(&param_stack[psp + i * npar + 2]) - 1;
  }

  if (need_sync) syncthreads ();
  
  for (int run = 0; run < nrun; run++)
  {
    /* pointer offsets */
    psp = npar * run;
    int srcA = param_stack_s[psp    ];
    int srcB = param_stack_s[psp + 1];
    int srcC = param_stack_s[psp + 2];
  
    /* load a and b matrix for current block and stack into smem (no comp/load overlap!) */
    if (m == n) {
#pragma unroll
      for (int i = tidx; i < mk; i += threads){
        buff_a[i] = __ldg(&a_data[srcA + i]);
        buff_b[i] = __ldg(&b_data[srcB + i]);
      }
    } else {
#pragma unroll
      for (int i = tidx; i < mk; i += threads){
        buff_a[i] = __ldg(&a_data[srcA + i]);
      }
#pragma unroll
      for (int i = tidx; i < kn; i += threads){
        buff_b[i] = __ldg(&b_data[srcB + i]);
      }
    }

    if (need_sync) syncthreads ();

    if (tidx < mn) {
      /* do actual multiplication in shared memory */
#pragma unroll
      for (int l = 0; l < k; l++) {
        myc += buff_a[l * m + r] * buff_b[l * n + c];
      } 
      /* store result in global memory */
      atomicAdd (&c_data[srcC + tidx], myc);
      myc = 0.0;
    }

    if (need_sync) syncthreads ();
  }

}
