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

template < int m,  int n,  int k, int M, int N, int threads, int grouping, int minblocks>
__global__ void
__launch_bounds__(threads, minblocks)
cusmm_dnt_medium(const int* __restrict__ param_stack, int stack_size,
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
  const int r = tidx - rmax * c;

  const int load_unroll_factor_1 = mk / threads + 1;
  const int load_unroll_factor_2 = kn / threads + 1;

  const int n_mkloads = mk / (load_unroll_factor_1 * threads);
  const int n_knloads = kn / (load_unroll_factor_2 * threads);

  const int mkloads_remained = (mk - n_mkloads * load_unroll_factor_1 * threads) / threads;
  const int knloads_remained = (kn - n_knloads * load_unroll_factor_2 * threads) / threads;

  const int mkloads_tail = ((mkloads_remained + n_mkloads * load_unroll_factor_1) * threads == mk)? 0: 1;
  const int knloads_tail = ((knloads_remained + n_knloads * load_unroll_factor_2) * threads == kn)? 0: 1;

  const int m_loads_to_finish = n_mkloads * (load_unroll_factor_1 * threads) + mkloads_remained * threads;
  const int n_loads_to_finish = n_knloads * (load_unroll_factor_2 * threads) + knloads_remained * threads;
  const int left_to_finish_1 = m_loads_to_finish + tidx;
  const int left_to_finish_2 = n_loads_to_finish + tidx;

  const int  npar = 3;
  const int  warp_size = 32;
  const bool need_sync = (mn > warp_size || mk > warp_size || kn > warp_size || threads > warp_size);

  double myc[N * M];
  double lba_r[load_unroll_factor_1 + mkloads_remained + mkloads_tail];
  double lbb_r[load_unroll_factor_2 + knloads_remained + knloads_tail];

  const double* __restrict__ buff_l;
  const double* __restrict__ buff_r;

  int  psp, nrun;
  bool is_loaded = false;

  __shared__ int    param_stack_s[npar * grouping];
  __shared__ double buff[buf_sz];

  buff_l = buff;
  buff_r = &(buff[mk]);

  nrun = grouping;
  if (((bidx + 1) * grouping) > stack_size) nrun = stack_size % grouping;

  /* set the partial sum to zero */
  for (int i = 0; i < M * N; i++)
    myc[i] = 0.0;


  /* load and pack stack data for current block into smem */
  psp = bidx * npar * grouping;
#pragma unroll 3
  for (int i = tidx; i < nrun; i += threads){
    param_stack_s[i * npar    ] = __ldg(&param_stack[psp + i * npar    ]) - 1;
    param_stack_s[i * npar + 1] = __ldg(&param_stack[psp + i * npar + 1]) - 1;
    param_stack_s[i * npar + 2] = __ldg(&param_stack[psp + i * npar + 2]) - 1;
  }

  if (need_sync) syncthreads ();

  for (int run = 0; run < nrun; run++){
    psp = run * npar;

    int srcA, srcB;

    if (!is_loaded){

      srcA = param_stack_s[psp    ];
      srcB = param_stack_s[psp + 1];

    if (m == n){
#pragma unroll 3
      for (int i = tidx; i < n_mkloads * (load_unroll_factor_1 * threads); i += load_unroll_factor_1 * threads){
#pragma unroll
        for (int l = 0; l < load_unroll_factor_1; l++){
          lba_r[l] = __ldg(&a_data[srcA + i + l * threads]);
          lbb_r[l] = __ldg(&b_data[srcB + i + l * threads]);
        }
      }

#pragma unroll 3
      for (int i = n_mkloads * (load_unroll_factor_1 * threads) + tidx; i < m_loads_to_finish; i += mkloads_remained * threads){
#pragma unroll
        for (int l = 0; l < mkloads_remained; l++){
          lba_r[l] = __ldg(&a_data[srcA + i + l * threads]);
          lbb_r[l] = __ldg(&b_data[srcB + i + l * threads]);
        }
      }

      if (left_to_finish_1 < mk){
        lba_r[load_unroll_factor_1 + mkloads_remained] = __ldg(&a_data[srcA + left_to_finish_1]);
        lbb_r[load_unroll_factor_2 + knloads_remained] = __ldg(&b_data[srcB + left_to_finish_2]);
      }
    } else {
#pragma unroll 3
      for (int i = tidx; i < n_mkloads * (load_unroll_factor_1 * threads); i += load_unroll_factor_1 * threads){
#pragma unroll
        for (int l = 0; l < load_unroll_factor_1; l++){
          lba_r[l] = __ldg(&a_data[srcA + i + l * threads]);
        }
      }

#pragma unroll 3
      for (int i = tidx; i < n_knloads * (load_unroll_factor_2 * threads); i += load_unroll_factor_2 * threads){
#pragma unroll
        for (int l = 0; l < load_unroll_factor_2; l++){
          lbb_r[l] = __ldg(&b_data[srcB + i + l * threads]);
        }
      }

#pragma unroll 3
      for (int i = n_mkloads * (load_unroll_factor_1 * threads) + tidx; i < m_loads_to_finish; i += mkloads_remained * threads){
#pragma unroll
        for (int l = 0; l < mkloads_remained; l++){
          lba_r[l] = __ldg(&a_data[srcA + i + l * threads]);
        }
      }

#pragma unroll 3
      for (int i = n_knloads * (load_unroll_factor_2 * threads) + tidx; i < n_loads_to_finish; i += knloads_remained * threads){
#pragma unroll
        for (int l = 0; l < knloads_remained; l++){
          lbb_r[l] = __ldg(&b_data[srcB + i + l * threads]);
        }
      }

      if (left_to_finish_1 < mk) lba_r[load_unroll_factor_1 + mkloads_remained] = __ldg(&a_data[srcA + left_to_finish_1]);
      if (left_to_finish_2 < kn) lbb_r[load_unroll_factor_2 + knloads_remained] = __ldg(&b_data[srcB + left_to_finish_2]);
    }

      /* hard sync necessary */
      syncthreads ();
    } else {
       is_loaded = false;
    }


    /* copy to smem */
    if (m == n){
#pragma unroll 3
      for (int i = n_mkloads * (load_unroll_factor_1 * threads) + tidx; i < m_loads_to_finish; i += mkloads_remained * threads){
#pragma unroll
        for (int l = 0; l < mkloads_remained; l++){
          buff[i + l * threads] = lba_r[l];
          buff[i + mk + l * threads] = lbb_r[l];
        }
      }

#pragma unroll 3
      for (int i = tidx; i < n_mkloads * (load_unroll_factor_1 * threads); i += load_unroll_factor_1 * threads){
#pragma unroll
        for (int l = 0; l < load_unroll_factor_1; l++){
          buff[i + l * threads] = lba_r[l];
          buff[i + mk + l * threads] = lbb_r[l];
        }
      }

      if (left_to_finish_1 < mk){
        buff[left_to_finish_1] = lba_r[load_unroll_factor_1 + mkloads_remained];
        buff[mk + left_to_finish_2] = lbb_r[load_unroll_factor_2 + knloads_remained];
      }
    } else {
#pragma unroll 3
      for (int i = n_mkloads * (load_unroll_factor_1 * threads) + tidx; i < m_loads_to_finish; i += mkloads_remained * threads){
#pragma unroll
        for (int l = 0; l < mkloads_remained; l++){
          buff[i + l * threads] = lba_r[l];
        }
      }

#pragma unroll 3
      for (int i = n_knloads * (load_unroll_factor_2 * threads) + tidx; i < n_loads_to_finish; i += knloads_remained * threads){
#pragma unroll
        for (int l = 0; l < knloads_remained; l++){
          buff[i + mk + l * threads] = lbb_r[l];
        }
      }

#pragma unroll 3
      for (int i = tidx; i < n_mkloads * (load_unroll_factor_1 * threads); i += load_unroll_factor_1 * threads){
#pragma unroll
        for (int l = 0; l < load_unroll_factor_1; l++){
          buff[i + l * threads] = lba_r[l];
        }
      }

#pragma unroll 3
      for (int i = tidx; i < n_knloads * (load_unroll_factor_2 * threads); i += load_unroll_factor_2 * threads){
#pragma unroll
        for (int l = 0; l < load_unroll_factor_2; l++){
          buff[i + mk + l * threads] = lbb_r[l];
        }
      }

      if (left_to_finish_1 < mk) buff[left_to_finish_1] = lba_r[load_unroll_factor_1 + mkloads_remained];
      if (left_to_finish_2 < kn) buff[mk + left_to_finish_2] = lbb_r[load_unroll_factor_2 + knloads_remained];
    }

    if (need_sync) syncthreads ();


    int next_run = run + 1;
    if (next_run >= nrun) is_loaded = true;

    if (!is_loaded || (run == 0 && nrun > 1)){
      int next_psp = next_run * npar;
      srcA = param_stack_s[next_psp    ];
      srcB = param_stack_s[next_psp + 1];


      if (m == n){
#pragma unroll 3
        for (int i = tidx; i < n_mkloads * (load_unroll_factor_1 * threads); i += load_unroll_factor_1 * threads){
#pragma unroll
          for (int l = 0; l < load_unroll_factor_1; l++){
            lba_r[l] = __ldg(&a_data[srcA + i + l * threads]);
            lbb_r[l] = __ldg(&b_data[srcB + i + l * threads]);
          }
        }

#pragma unroll 3
        for (int i = n_mkloads * (load_unroll_factor_1 * threads) + tidx; i < m_loads_to_finish; i += mkloads_remained * threads){
#pragma unroll
          for (int l = 0; l < mkloads_remained; l++){
            lba_r[l] = __ldg(&a_data[srcA + i + l * threads]);
            lbb_r[l] = __ldg(&b_data[srcB + i +l * threads]);
          }
        }

        if (left_to_finish_1 < mk){
          lba_r[load_unroll_factor_1 + mkloads_remained] = __ldg(&a_data[srcA + left_to_finish_1]);
          lbb_r[load_unroll_factor_2 + knloads_remained] = __ldg(&b_data[srcB + left_to_finish_2]);
        }
      } else {
#pragma unroll 3
        for (int i = tidx; i < n_mkloads * (load_unroll_factor_1 * threads); i += load_unroll_factor_1 * threads){
#pragma unroll
          for (int l = 0; l < load_unroll_factor_1; l++){
            lba_r[l] = __ldg(&a_data[srcA + i + l * threads]);
          }
        }

#pragma unroll 3
        for (int i = tidx; i < n_knloads * (load_unroll_factor_2 * threads); i += load_unroll_factor_2 * threads){
#pragma unroll
          for (int l = 0; l < load_unroll_factor_2; l++){
            lbb_r[l] = __ldg(&b_data[srcB + i + l * threads]);
          }
        }

#pragma unroll 3
        for (int i = n_mkloads * (load_unroll_factor_1 * threads) + tidx; i < m_loads_to_finish; i += mkloads_remained * threads){
#pragma unroll
          for (int l = 0; l < mkloads_remained; l++){
            lba_r[l] = __ldg(&a_data[srcA + i + l * threads]);
          }
        }

#pragma unroll 3
        for (int i = n_knloads * (load_unroll_factor_2 * threads) + tidx; i < n_loads_to_finish; i += knloads_remained * threads){
#pragma unroll
          for (int l = 0; l < knloads_remained; l++){
            lbb_r[l] = __ldg(&b_data[srcB + i + l * threads]);
          }
        }

        if (left_to_finish_1 < mk) lba_r[load_unroll_factor_1 + mkloads_remained] = __ldg(&a_data[srcA + left_to_finish_1]);
        if (left_to_finish_2 < kn) lbb_r[load_unroll_factor_2 + knloads_remained] = __ldg(&b_data[srcB + left_to_finish_2]);
      }

      is_loaded = true;
    }

    /* Do actual multiplication. */
    if (c < cmax  && r < rmax){
      for (int l = 0; l < k; l++){
        for (int i = 0; i < M; i++){
          for (int j = 0; j < N; j++){
//            if (M * r + i < m && N * c + j < n)          // ^ref1
            myc[N * i + j] += buff_l[l * m + M * r + i] * buff_r[l * n + N * c + j];
          }
        }
      }
    }


    if (run == nrun - 1 || param_stack_s[psp + 2] != param_stack_s[psp + 2 + npar]){

      int srcC = param_stack_s[psp + 2];

      if (M > 1 || N > 1){
        if (need_sync) syncthreads();
        /* decompress results to buffer */
        if (c < cmax  && r < rmax){
          for (int i = 0; i < M; i++){
            for (int j = 0; j < N; j++){
              if (M * r + i < m && N * c + j < n){
                buff[(N * c + j) * m + M * r + i] = myc[N * i + j];
                myc[N * i + j] = 0.0;
              }
            }
          }
        }
        if (need_sync) syncthreads();

        /* Add results to global C block. */
#pragma unroll
        for (int i = tidx; i < mn; i += threads){
          atomicAdd (&c_data[srcC + i], buff[i]);
        }
      } else {
        /* Add results to global C block. */
#pragma unroll
        for (int i = tidx; i < mn; i += threads){
          atomicAdd (&c_data[srcC + i], myc[0]);
        }
        myc[0] = 0.0;
      }
    }
    if (need_sync) syncthreads ();
  }
}
