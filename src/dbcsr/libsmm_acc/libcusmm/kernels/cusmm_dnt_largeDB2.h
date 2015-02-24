/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015  CP2K developers group                         *
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


namespace ns_cusmm_dnt_largeDB2 {

//**************************************************************************//
__device__ static inline void load_gmem_into_regs(double *from, double *dest,
                                           const int length, const int blockdim)
{
    const int NR = (length + blockdim - 1) / blockdim;

    if (length <= blockdim) { // are there enough threads to load in one step?
        if (threadIdx.x < length)
            dest[0] = __ldg(from + threadIdx.x);
    }  else if (length % blockdim == 0) {
        int i = threadIdx.x;
        for (int ri = 0; ri < NR; ri++) {  //loop with fixed bounds
            dest[ri] = __ldg(from + i);
            i += blockdim;
        }
   } else {
        int i = threadIdx.x;
        for (int ri = 0; ri < NR-1; ri++) {  //loop with fixed bounds
            dest[ri] = __ldg(from + i);
            i += blockdim;
        }
        if (i < length)
            dest[NR-1] = __ldg(from + i);
    }
}


//**************************************************************************//
__device__ static inline void load_regs_into_smem(double *from, double *dest,
                                           const int length, const int blockdim)
{
   const int NR = (length + blockdim - 1) / blockdim;

   if (length <= blockdim) { // are there enough threads to load in one step?
       if (threadIdx.x < length)
           dest[threadIdx.x] = from[0];
   } else if (length % blockdim == 0) {
        int i = threadIdx.x;
        for (int ri = 0; ri < NR; ri++) {  //loop with fixed bounds
            dest[i] = from[ri];
            i += blockdim;
        }
   } else {
        int i = threadIdx.x;
        for (int ri = 0; ri < NR-1; ri++) {  //loop with fixed bounds
            dest[i] = from[ri];
            i += blockdim;
        }
        if (i < length)
            dest[i] = from[NR-1];
    }
}


//**************************************************************************//
__device__ static inline void multiply(double *buff_a, double *buff_b, double *buff_c,
                                const int w, const int m, const int n,
                                const int M, const int N, const int blockdim)
{
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
                    buff_c[M * i + j] +=
                        buff_a[l * m + M * r + j] * buff_b[l * n + N * c + i];
}


//**************************************************************************//
__device__ static inline void store_results_into_smem(double *from, double *dest,
                                               const int t, const int v,
                                               const int m, const int n,
                                               const int M, const int N,
                                               const int blockdim)
{
    const int rmax = (m + M - 1) / M; //  max tile-row
    const int c = threadIdx.x / rmax; // this thread's tile-column
    const int r = threadIdx.x - c * rmax; // this thread's tile-row

    const int ctmp = c * N - t;
    if (ctmp >= -(N - 1) && ctmp < v)
        for (int i = 0; i < N; i++)
            if (ctmp + i >= 0 && ctmp + i < v)
                for (int j = 0; j < M; j++)
                    if (M * r + j < m) {
                        dest[(ctmp + i) * m + M * r + j] = from[M * i + j];
                        from[M * i + j] = 0.0; // reset result tile
                    }

}


//**************************************************************************//
__device__ static inline void writeback_results(double *from, double *dest,
                                                double *buff,
                                           const int m, const int n,
                                           const int M, const int N,
                                           const int v, const int blockdim)
{
   // results are written in output-slabs of width v
      for (int t = 0; t < (n / v) * v; t += v) {
          // copy output slab from registers to shared memory
          store_results_into_smem(from, buff, t, v, m, n, M, N, blockdim);
          syncthreads();
          // Add our results to the accumulator in global memory
          for (int i = threadIdx.x; i < m * v; i += blockdim)
              atomicAdd(dest + i, buff[i]);
          dest += m * v;
          syncthreads();
      }

      // If the output slab witdh v is not a divisor of n,
      // a smaller tail-slab of width va has to be process
      const int va = n - (n / v) * v;
      if (va != 0) {  // is there a tail-slab?
          int t = (n / v) * v;
          store_results_into_smem(from, buff, t, va, m, n, M, N, blockdim);
          syncthreads();
          for (int i = threadIdx.x; i < m * va; i += blockdim)
              atomicAdd(dest + i, buff[i]);
          syncthreads();
      }

}

} //end of namespace


//**************************************************************************//
template < int m, int n, int k, int M, int N, int w, int v, int blockdim, int grouping, int minblocks >
__global__
__launch_bounds__(blockdim,minblocks)
void cusmm_dnt_largeDB2(const int *__restrict__ param_stack, int careful, int nruns, double *a_data, double *b_data, double *c_data){

	using namespace ns_cusmm_dnt_largeDB2;

    // registers to store thread's result tile
    double myc[N * M];

    // registers to store input slabs during double buffering
    // If there are too few thread, each thread has to store
    // multiple elements of the input slabs in it's registers.
    const int mya_size = (w*m + blockdim-1)/blockdim;
    const int myb_size = (w*n + blockdim-1)/blockdim;
    double mya[mya_size];
    double myb[myb_size];

     // initialize the thread's result tile to zero
    for (int i = 0; i < N * M; i++)
        myc[i] = 0.0;

    // buffer needs to hold input and output slabs (not both simultaneously).
    const int buff_size = MAX(m*w + w*n, v * m);
    __shared__ double buff[buff_size];

    // conveniece pointers
    double *buff_l = buff;
    double *buff_r = &(buff[m * w]);

    // first stack entry to be processed by this thread-block
    int psp = 7 * (blockIdx.x * grouping);

    // grouping is the number of stack entries process by each thread-block
    // careful is the number of launched thread-blocks.
    // nruns is the number of stack entries process by the last thread-block
    int nrun = (blockIdx.x == careful) ? nruns : grouping;

    // all stack entries relavant for this thread-block are loaded at once
    // allows to look ahead and and flush result tile only when really needed
    __shared__ int param_stack_s[4 * grouping];

    // load parameter stack, might read beyond
    for (int i = threadIdx.x; i < 7 * nrun; i += blockdim) {
        int p_tmp = __ldg(param_stack + psp + i);
        if (i % 7 > 2)
            param_stack_s[(i / 7) * 4 + i % 7 - 3] = p_tmp - 1;
    }

    psp = 0;

    syncthreads();

    // get the offsets for the a-block and the b-block from the stack
    int srcA = param_stack_s[psp];
    int srcB = param_stack_s[psp + 1];
    // start off double buffering by loading the first data
    load_gmem_into_regs(a_data + srcA, mya, m * w, blockdim);
    load_gmem_into_regs(b_data + srcB, myb, n * w, blockdim);

    syncthreads();
    const int wa = k - (k / w) * w;
    // in each run we process one stack entry
    for (int run = 0; run < nrun; run++) {

        // load the first slab for multiplication into the smem
        load_regs_into_smem(mya, buff_l, m * w, blockdim);
        load_regs_into_smem(myb, buff_r, n * w, blockdim);
        syncthreads();

        // this is the actual double buffering loop
        for (int t = 0; t < (k/w -1) * w ; t += w) {
            // load next input slab from global memory into registers
            srcA += m * w;
            srcB += n * w;
            load_gmem_into_regs(a_data + srcA, mya, m * w, blockdim);
            load_gmem_into_regs(b_data + srcB, myb, n * w, blockdim);
            // multiply previous slab, which is stored in shared memory,
            // and accumulate the results in the registers myc
            multiply(buff_l, buff_r, myc, w, m, n, M, N, blockdim);
            syncthreads();
            // copy next slab from registers to shared memory
            load_regs_into_smem(mya, buff_l, m * w, blockdim);
            load_regs_into_smem(myb, buff_r, n * w, blockdim);
            syncthreads();
        }

        if (wa != 0) { // is there a tail-slab?
            // If the input slab witdh w is not a divisor of k,
            // a smaller tail-slab of width wa has to be process
            // load tail-slab into registers
            srcA += m * w;
            srcB += n * w;
            load_gmem_into_regs(a_data + srcA, mya, m * wa, blockdim);
            load_gmem_into_regs(b_data + srcB, myb, n * wa, blockdim);
            // multiply last regular slab, which the loop left in shared memory
            multiply(buff_l, buff_r, myc, w, m, n, M, N, blockdim);
            syncthreads();
            // copy tail-slab from register into shared mem
            load_regs_into_smem(mya, buff_l, m * wa, blockdim);
            load_regs_into_smem(myb, buff_r, n * wa, blockdim);
            syncthreads();
        }

        psp = 4 * run + 4;

        if(run < nrun-1){
            // get the offsets for the a-block and the b-block from the stack
            srcA = param_stack_s[psp];
            srcB = param_stack_s[psp + 1];
            // load the data for the next iteration of the loop
            load_gmem_into_regs(a_data + srcA, mya, m * w, blockdim);
            load_gmem_into_regs(b_data + srcB, myb, n * w, blockdim);
        }

        if (wa != 0) { // is there a tail-slab?
            // multiply the tail-slab
            multiply(buff_l, buff_r, myc, wa, m, n, M, N, blockdim);
        }else{
            // multiply last regular slab, which the loop left in shared memory
            multiply(buff_l, buff_r, myc, w, m, n, M, N, blockdim);
        }
        syncthreads();

        // multiplication for this run done
        // do we have to flush the result tile?
        if(run == nrun-1 || (param_stack_s[psp - 1] != param_stack_s[psp + 3])) {
            int c_loc = param_stack_s[psp - 2];
            writeback_results(myc, c_data+c_loc, buff, m,n,M,N,v,blockdim);
        }
    }
}

//EOF
