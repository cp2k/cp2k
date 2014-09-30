/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2014 the CP2K developers group                      *
 *  Authors: Andreas Gloess <andreas.gloess@chem.uzh.ch>                     *
 *****************************************************************************/

#if defined (__ACC)

#include "clsmm_common.h"

//**************************************************************************//
__kernel void clsmm_dnt_largeDB_16_23_23_12_23_96_2_3_12_10 (
                __global int    *param_stack,
                         int    careful,
                         int    nruns,
                __global double *a_data,
                __global double *b_data,
                __global double *c_data)
{
    const int m = 23;
    const int n = 23;
    const int k = 23;
    const int M = 2;
    const int N = 3;
    const int w = 10;
    const int v = 12;
    const int blockdim = 96;
    const int grouping = 16;
    const int minblocks = 12;
    // registers to store thread's result tile
    //double myc[N * M];
    __private double myc[3 * 2];

    // registers to store input slabs during double buffering
    // If there are too few thread, each thread has to store
    // multiple elements of the input slabs in it's registers.
    const int mya_size = (w * m + blockdim - 1) / blockdim;
    const int myb_size = (w * n + blockdim - 1) / blockdim;
    //double mya[mya_size];
    //double myb[myb_size];
    __private double mya[(10 * 23 + 96 - 1) / 96];
    __private double myb[(10 * 23 + 96 - 1) / 96];

     // initialize the thread's result tile to zero
    for (int i = 0; i < N * M; i++)
        myc[i] = 0.0;

    // buffer needs to hold input and output slabs (not both simultaneously).
    //const int buff_size = MAX(m * w + w * n, v * m);
    //__local double buff[buff_size];
    const int buff_size = 460;
    __local double buff[460];

    // conveniece pointers
    // double *buff_l = buff;
    //double *buff_r = &(buff[m * w]);
    __local double *buff_l = buff;
    __local double *buff_r = &(buff[23 * 10]);

    // first stack entry to be processed by this thread-block
    int psp = 7 * (get_group_id(0) * grouping);

    // grouping is the number of stack entries process by each thread-block
    // careful is the number of launched thread-blocks.
    // nruns is the number of stack entries process by the last thread-block
    int nrun = (get_group_id(0) == careful) ? nruns : grouping;

    // all stack entries relavant for this thread-block are loaded at once
    // allows to look ahead and and flush result tile only when really needed
    //__local int param_stack_s[4 * grouping];
    __local int param_stack_s[4 * 16];

    // load parameter stack, might read beyond
    for (int i = get_local_id(0); i < 7 * nrun; i += blockdim) {
        //int p_tmp = __ldg(param_stack + psp + i);
        int p_tmp = *(param_stack + psp + i);
        if (i % 7 > 2)
            param_stack_s[(i / 7) * 4 + i % 7 - 3] = p_tmp - 1;
    }

    // in each run we process one stack entry
    for (int run = 0; run < nrun; run++) {
        psp = run * 4;

        barrier(CLK_LOCAL_MEM_FENCE);

        // get the offsets for the a-block and the b-block from the stack
        int srcA = param_stack_s[psp];
        int srcB = param_stack_s[psp + 1];

        // start off double buffering by loading the first
        // input slab directly from global into shared memory
        load_gmem_into_smem(a_data + srcA, buff_l, m * w, blockdim);
        load_gmem_into_smem(b_data + srcB, buff_r, n * w, blockdim);

        // this is the actual double buffering loop
        for (int t = 0; t < (k / w -1) * w ; t += w) {
            barrier(CLK_LOCAL_MEM_FENCE);
            // load next input slab from global memory into registers
            srcA += m * w;
            srcB += n * w;
            load_gmem_into_regs(a_data + srcA, mya, m * w, blockdim);
            load_gmem_into_regs(b_data + srcB, myb, n * w, blockdim);
            // multiply previous slab, which is stored in shared memory,
            // and accumulate the results in the registers myc
            multiply(buff_l, buff_r, myc, w, m, n, M, N, blockdim);
            barrier(CLK_LOCAL_MEM_FENCE);
            // copy next slab from registers to shared memory
            load_regs_into_smem(mya, buff_l, m * w, blockdim);
            load_regs_into_smem(myb, buff_r, n * w, blockdim);
        }

        barrier(CLK_LOCAL_MEM_FENCE);

        // If the input slab witdh w is not a divisor of k,
        // a smaller tail-slab of width wa has to be process
        const int wa = k - (k / w) * w;
        if (wa != 0) { // is there a tail-slab?
            // load tail-slab into registers
            srcA += m * w;
            srcB += n * w;
            load_gmem_into_regs(a_data + srcA, mya, m * wa, blockdim);
            load_gmem_into_regs(b_data + srcB, myb, n * wa, blockdim);
        }

        // multiply last regular slab, which the loop left in shared memory
        multiply(buff_l, buff_r, myc, w, m, n, M, N, blockdim);
        barrier(CLK_LOCAL_MEM_FENCE);

        if (wa != 0) { // is there a tail-slab?
            // copy tail-slab from register into shared mem
            load_regs_into_smem(mya, buff_l, m * wa, blockdim);
            load_regs_into_smem(myb, buff_r, n * wa, blockdim);
            barrier(CLK_LOCAL_MEM_FENCE);
            // multiply the tail-slab
            multiply(buff_l, buff_r, myc, wa, m, n, M, N, blockdim);
            barrier(CLK_LOCAL_MEM_FENCE);
        }


        // multiplication for this run done
        // do we have to flush the result tile?
        if (run == nrun - 1
            || param_stack_s[psp + 3] != param_stack_s[psp + 3 + 4]) {
            int c_loc = param_stack_s[psp + 2];

            barrier(CLK_LOCAL_MEM_FENCE);

            // results are written in output-slabs of width v
            for (int t = 0; t < (n / v) * v; t += v) {
                // copy output slab from registers to shared memory
                store_results_into_smem(myc, buff, t, v, m, n, M, N, blockdim);
                barrier(CLK_LOCAL_MEM_FENCE);
                // Add our results to the accumulator in global memory
                for (int i = get_local_id(0); i < m * v; i += blockdim)
                    AddAtomic(&c_data[c_loc + i], buff[i]);
                c_loc += m * v;
                barrier(CLK_LOCAL_MEM_FENCE);
            }

            // If the output slab witdh v is not a divisor of n,
            // a smaller tail-slab of width va has to be process
            const int va = n - (n / v) * v;
            if (va != 0) {  // is there a tail-slab?
                int t = (n / v) * v;
                store_results_into_smem(myc, buff, t, va, m, n, M, N, blockdim);
                barrier(CLK_LOCAL_MEM_FENCE);
                for (int i = get_local_id(0); i < m * va; i += blockdim)
                    AddAtomic(&c_data[c_loc + i], buff[i]);
                barrier(CLK_LOCAL_MEM_FENCE);

            }
        }
    }
}


#endif

//EOF
