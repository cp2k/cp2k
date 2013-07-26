/******************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations
 *  Copyright (C) 2000 - 2013 the CP2K developers group
 *****************************************************************************/

// DBCSR_KERNEL datatype=dbcsr_type_real_8, homogeneous_only=True, m=8, n=8, k=8

#include "dbcsr_kernel.h"
#include "stack_mm_mnk_sq8_d.h"

__global__ void
//__launch_bounds__(32, 16)
stack_mm_mnk_sq8_d (const int *__restrict__ param_stack,
                    const int careful, const int nruns,
                    const int m, const int n, const int k,
                    const double *__restrict__ a_data,
                    const double *__restrict__ b_data,
                    double *__restrict__ c_data)
{
  int psp;

  __shared__ double buff_a1[8*8];
  __shared__ double buff_b1[8*8];

  __shared__ double buff_a2[8*8];
  __shared__ double buff_b2[8*8];

  int col    = threadIdx.x / 8;
  int row    = threadIdx.x - col * 8;

  int nrun = GROUPING;

  //if (blockIdx.x == gridDim.x-1)
  if (blockIdx.x == careful)
    nrun = nruns;

  /* Set the partial sum to zero */
  double myc = 0.0;

  for (int run = 0; run < nrun; run++)
    {
      psp = 7 * (blockIdx.x * GROUPING + run);

      // load a matrix into smem
      if(threadIdx.x >= 64){
         //int a_base = param_stack[psp+3]-1;
         //int b_base = param_stack[psp+4]-1;
         //buff_a1[threadIdx.x-64] = a_data[a_base + threadIdx.x-64];
         //buff_b1[threadIdx.x-64] = b_data[b_base + threadIdx.x-64];
//         buff_a1[threadIdx.x] = a_data[a_base + threadIdx.x];
//         buff_b1[threadIdx.x] = b_data[b_base + threadIdx.x];

         int a_base = param_stack[psp+3]-1;
         int b_base = param_stack[psp+4]-1;
         buff_a1[threadIdx.x-64] = __ldg(a_data+a_base + threadIdx.x-64);
         buff_b1[threadIdx.x-64] = __ldg(b_data+b_base + threadIdx.x-64);

      }

      syncthreads();

      if(threadIdx.x < 64){
        for(int k=0; k<8; k++){
          myc += buff_a1[row + k*8] * buff_b1[k + col*8];
        }

        if (run == nrun - 1
          || param_stack[psp + 6] - 1 != param_stack[psp + 6 + 7] - 1) {
         int  c_loc = param_stack[psp + 5] - 1;
         atomicAdd (&c_data[c_loc + threadIdx.x], myc);

          myc = 0.0l;
        }
      }
//      syncthreads();

      run++;
      psp = 7 * (blockIdx.x * GROUPING + run);
      if(threadIdx.x >= 64) {
         if(run < nrun){
           int a_base = param_stack[psp+3]-1;
           int b_base = param_stack[psp+4]-1;
           buff_a2[threadIdx.x-64] = a_data[a_base + threadIdx.x-64];
           buff_b2[threadIdx.x-64] = b_data[b_base + threadIdx.x-64];
         //  int a_base = param_stack[psp+3]-1;
         //  int b_base = param_stack[psp+4]-1;
         //  buff_a2[threadIdx.x-64] = __ldg(a_data+a_base + threadIdx.x-64);
         //  buff_b2[threadIdx.x-64] = __ldg(b_data+b_base + threadIdx.x-64);
         }
      }

      syncthreads();

      if(threadIdx.x < 64 && run < nrun){
        for(int k=0; k<8; k++){
          myc += buff_a2[row + k*8] * buff_b2[k + col*8];
        }

        if (run == nrun - 1
          || param_stack[psp + 6] - 1 != param_stack[psp + 6 + 7] - 1) {
         int  c_loc = param_stack[psp + 5] - 1;
         atomicAdd (&c_data[c_loc + threadIdx.x], myc);

          myc = 0.0l;
        }
      }
   }

}


//==============================================================================
int launch_stack_mm_mnk_sq8_d(int *param_stack, int stack_size, cudaStream_t stream,
     int m_max, int n_max, int k_max, double *a_data, double *b_data, double *c_data){

      int careful = (stack_size / GROUPING);
      int nruns = stack_size - careful * GROUPING;
      int shared_size = 0;

      stack_mm_mnk_sq8_d <<< ((stack_size + GROUPING - 1) / GROUPING),
      128, shared_size, stream >>> (param_stack, careful, nruns,
                                              m_max, n_max, k_max,
                                              a_data,
                                              b_data,
                                              c_data);
      return(0);
}
