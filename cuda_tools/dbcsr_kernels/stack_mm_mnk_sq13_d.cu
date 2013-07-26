/******************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations
 *  Copyright (C) 2000 - 2013 the CP2K developers group
 *****************************************************************************/


#include "dbcsr_kernel.h"
#include "stack_mm_mnk_sq13_d.h"


__global__ void
//__launch_bounds__(32, 16)
stack_mm_mnk_sq13_d (const int *__restrict__ param_stack,
                    const int careful, const int nruns,
                    const int m, const int n, const int k,
                    const int liter,
                    const double *__restrict__ a_data,
                    const double *__restrict__ b_data,
                    double *__restrict__ c_data, int *__restrict__ c_locks,
                    int lock_offset)
{

  int psp;

  __shared__ double buff_a[13*13];
  __shared__ double buff_b[13*13];

  int col    = threadIdx.x / 13;
  int row    = threadIdx.x - col * 13;

  int nrun = GROUPING;

  if (blockIdx.x == careful)
    nrun = nruns;

  /* Set the partial sum to zero */
  double myc = 0.0;

  for (int run = 0; run < nrun; run++) {
      psp = 7 * (blockIdx.x * GROUPING + run);
      int a_base = param_stack[psp+3]-1;
      int b_base = param_stack[psp+4]-1;

      // load a matrix into smem
      if(threadIdx.x < 13*13){
         //buff_a[threadIdx.x] = a_data[a_base + threadIdx.x];
         buff_a[threadIdx.x] = __ldg(a_data+a_base + threadIdx.x);
         //buff_b[threadIdx.x] = b_data[b_base + threadIdx.x];
         buff_b[threadIdx.x] = __ldg(b_data+b_base + threadIdx.x);
      }

      syncthreads();

      if(threadIdx.x < 13*13){
        for(int k=0; k<13; k++){
          myc += buff_a[row + k*13] * buff_b[k + col*13];
        }

        if (run == nrun - 1
          || param_stack[psp + 6] - 1 != param_stack[psp + 6 + 7] - 1) {
          /* Add results to global C block. */
         int  c_loc = param_stack[psp + 5] - 1;
         atomicAdd (&c_data[c_loc + threadIdx.x], myc);

          myc = 0.0l;
        }
      }
      syncthreads();
  }
}


