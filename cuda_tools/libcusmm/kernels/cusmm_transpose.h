/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2013 the CP2K developers group                      *
 *  Authors: Peter Messmer <pmessmer@nvidia.com>,                            *
 *           Nikolay Markovskiy <nmarkovskiy@nvidia.com>                     *
 *****************************************************************************/
template < int m, int n>
__global__ void transpose_d(int *trs_stack, int nblks, double* mat){
 __shared__ double buf[m*n];
 int offset = trs_stack[blockIdx.x];
 for(int i=threadIdx.x; i < m*n; i+=blockDim.x){
     buf[i] = mat[offset + i]; 
 }
 syncthreads();

 for(int i=threadIdx.x; i < m*n; i+=blockDim.x){
     int r_out = i % n;
     int c_out = i / n;
     int idx = r_out * m + c_out;
     mat[offset + i] = buf[idx];
 }

}

//EOF
