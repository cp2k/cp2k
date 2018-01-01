/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2018  CP2K developers group                         *
 *****************************************************************************/

/*****************************************************************************
 *  Authors: Peter Messmer <pmessmer@nvidia.com>,                            *
 *           Nikolay Markovskiy <nmarkovskiy@nvidia.com>                     *
 *           Ole Schuett <ole.schuett@mat.ethz.ch>                           *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include "cusmm_common.h"



namespace ns_cusmm_dnt_largeDB1 {
//**************************************************************************//
__device__ static inline void load_gmem_into_smem(const double* __restrict__ from, double* dest,
                                           const int length, const int threads){
  if (length < threads){ // are there enough threads to load in one step?
    if (threadIdx.x < length)
      dest[threadIdx.x] = __ldg(&from[threadIdx.x]);
  } else {
    for (int i = threadIdx.x; i < length; i += threads)
      dest[i] = __ldg(&from[i]);
  }
}


//**************************************************************************//
__device__ static inline void load_gmem_into_regs(const double* __restrict__ from, double* dest,
                                           const int length, const int threads){
  const int NR = (length + threads - 1) / threads;

  if (length < threads){ // are there enough threads to load in one step?
    if (threadIdx.x < length)
      dest[0] = __ldg(&from[threadIdx.x]);
  } else {
    int i = threadIdx.x;
    for (int ri = 0; ri < NR; ri++){  //loop with fixed bounds
      if (i < length)
        dest[ri] = __ldg(&from[i]);
      i += threads;
    }
  }
}


//**************************************************************************//
__device__ static inline void load_regs_into_smem(double* from, double* dest,
                                           const int length, const int threads){
   const int NR = (length + threads - 1) / threads;

  if (length < threads){ // are there enough threads to load in one step?
    if (threadIdx.x < length)
      dest[threadIdx.x] = from[0];
  } else {
    int i = threadIdx.x;
    for (int ri = 0; ri < NR; ri++){  //loop with fixed bounds
      if (i < length)
        dest[i] = from[ri];
      i += threads;
    }
  }
}


//**************************************************************************//
__device__ static inline void multiply(double* buff_a, double* buff_b, double* buff_c,
                                const int w, const int m, const int n,
                                const int M, const int N){
    // There might be more threads than needed for the calculation.
    // Only the first cmax*rmax threads participate in the calculation.

  const int cmax = (n + N - 1) / N; // max tile-column
  const int rmax = (m + M - 1) / M; //  max tile-row
  const int c = threadIdx.x / rmax; // this thread's tile-column
  const int r = threadIdx.x - c * rmax; // this thread's tile-row

  if (c < cmax && r < rmax) // is this thread participating?
    for (int l = 0; l < w; l++)
      for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++)
          buff_c[M * i + j] += buff_a[l * m + M * r + j] * buff_b[l * n + N * c + i];
}


//**************************************************************************//
__device__ static inline void store_results_into_smem(double* from, double* dest,
                                               const int t, const int v,
                                               const int m, const int n,
                                               const int M, const int N){
  const int rmax = (m + M - 1) / M; //  max tile-row
  const int c = threadIdx.x / rmax; // this thread's tile-column
  const int r = threadIdx.x - c * rmax; // this thread's tile-row

  int ctmp = c * N - t;
  if (ctmp >= -(N - 1) && ctmp < v)
    for (int i = 0; i < N; i++)
      if (ctmp + i >= 0 && ctmp + i < v)
        for (int j = 0; j < M; j++)
          if (M * r + j < m){
            dest[(ctmp + i) * m + M * r + j] = from[M * i + j];
            from[M * i + j] = 0.0; // reset result tile
          }
}

} //end of namespace

//**************************************************************************//
template < int m, int n, int k, int M, int N, int w, int v, int threads, int grouping, int minblocks >
__global__
__launch_bounds__(threads, minblocks)
void cusmm_dnt_largeDB1(const int* param_stack, const int stack_size,
  const double* a_data, const double* b_data, double* c_data){

  using namespace ns_cusmm_dnt_largeDB1;

  const int npar = 3;
  const int nrun = (((blockIdx.x + 1) * grouping) > stack_size)? stack_size % grouping: grouping;

  // registers to store input slabs during double buffering
  // If there are too few thread, each thread has to store
  // multiple elements of the input slabs in it's registers.
  const int mya_size = (w * m + threads - 1) / threads;
  const int myb_size = (w * n + threads - 1) / threads;
  const int buff_tmp  = MAX(        (w - 1) * m + ((m + M - 1) / M) * M,
                            m * w + (w - 1) * n + ((n + N - 1) / N) * N);
  const int buff_size = MAX(buff_tmp, v * m);

  // registers to buffer and store thread's result tile
  double mya[mya_size];
  double myb[myb_size];
  double myc[M * N];

  //__shared__ int    param_stack_s[npar * grouping]; // by putting this here: regs= 92-->128 (why?)
  __shared__ double buff[buff_size];
  __shared__ int    param_stack_s[npar * grouping];

  double* buff_l = buff;
  double* buff_r = &(buff[m * w]);

  /* set the partial sum to zero */
  for (int i = 0; i < M * N; i++)
      myc[i] = 0.0;

  /* load and pack stack data for current block into smem */
  int psp = blockIdx.x * npar * grouping;
#pragma unroll 3
  for (int i = threadIdx.x; i < nrun; i += threads){
    param_stack_s[i * npar    ] = __ldg(&param_stack[psp + i * npar    ]) - 1;
    param_stack_s[i * npar + 1] = __ldg(&param_stack[psp + i * npar + 1]) - 1;
    param_stack_s[i * npar + 2] = __ldg(&param_stack[psp + i * npar + 2]) - 1;
  }

  // in each run we process one stack entry
  for (int run = 0; run < nrun; run++){
    psp = run * npar;

    syncthreads();

    // get the offsets for the a-block and the b-block from the stack
    int srcA = param_stack_s[psp];
    int srcB = param_stack_s[psp + 1];

    // start off double buffering by loading the first
    // input slab directly from global into shared memory
    load_gmem_into_smem(&a_data[srcA], buff_l, m * w, threads);
    load_gmem_into_smem(&b_data[srcB], buff_r, n * w, threads);

    // this is the acutall double buffering loop
    for (int t = 0; t < (k / w - 1) * w ; t += w){
      syncthreads();
      // load next input slab from global memory into registers
      srcA += m * w;
      srcB += n * w;
      load_gmem_into_regs(&a_data[srcA], mya, m * w, threads);
      load_gmem_into_regs(&b_data[srcB], myb, n * w, threads);
      // multiply previous slab, which is stored in shared memory,
      // and accumulate the results in the registers myc
      multiply(buff_l, buff_r, myc, w, m, n, M, N);
      syncthreads();
      // copy next slab from registers to shared memory
      load_regs_into_smem(mya, buff_l, m * w, threads);
      load_regs_into_smem(myb, buff_r, n * w, threads);
    }

    syncthreads();

    // If the input slab witdh w is not a divisor of k,
    // a smaller tail-slab of width wa has to be process
    const int wa = k - (k / w) * w;
    if (wa != 0){ // is there a tail-slab?
      // load tail-slab into registers
      srcA += m * w;
      srcB += n * w;
      load_gmem_into_regs(&a_data[srcA], mya, m * wa, threads);
      load_gmem_into_regs(&b_data[srcB], myb, n * wa, threads);
    }

    // multiply last regular slab, which the loop left in shared memory
    multiply(buff_l, buff_r, myc, w, m, n, M, N);
    syncthreads();

    if (wa != 0){ // is there a tail-slab?
      // copy tail-slab from register into shared mem
      load_regs_into_smem(mya, buff_l, m * wa, threads);
      load_regs_into_smem(myb, buff_r, n * wa, threads);
      syncthreads();
      // multiply the tail-slab
      multiply(buff_l, buff_r, myc, wa, m, n, M, N);
      syncthreads();
    }


    // multiplication for this run done
    // do we have to flush the result tile?
    if (run == nrun - 1
      || param_stack_s[psp + 2] != param_stack_s[psp + 2 + npar]){
      int srcC = param_stack_s[psp + 2];

      syncthreads();

      // results are written in output-slabs of width v
      for (int t = 0; t < (n / v) * v; t += v) {
        // copy output slab from registers to shared memory
        store_results_into_smem(myc, buff, t, v, m, n, M, N);
        syncthreads();
        // Add our results to the accumulator in global memory
        for (int i = threadIdx.x; i < m * v; i += threads)
          atomicAdd(&c_data[srcC + i], buff[i]);
        srcC += m * v;
        syncthreads();
      }

      // If the output slab witdh v is not a divisor of n,
      // a smaller tail-slab of width va has to be process
      const int va = n - (n / v) * v;
      if (va != 0){  // is there a tail-slab?
        int t = (n / v) * v;
        store_results_into_smem(myc, buff, t, va, m, n, M, N);
        syncthreads();
        for (int i = threadIdx.x; i < m * va; i += threads)
          atomicAdd(&c_data[srcC + i], buff[i]);
        syncthreads();
      }

    }
  }
}

//EOF
