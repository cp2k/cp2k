/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2013 the CP2K developers group                      *
 *  Authors: Peter Messmer <pmessmer@nvidia.com>,                            *
 *           Nikolay Markovskiy <nmarkovskiy@nvidia.com>                     *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include "cusmm_common.h"


template < int m, int n, int k, int M, int N, int w, int v, int blockdim, int grouping>
__global__ void
cusmm_dnt_largeDB(const int* __restrict__ param_stack, int careful, int nruns,
     double* a_data, double* b_data, double* c_data){

  const int buf_sz = (m*w + w*n < v*m) ? v*m : m*w + w*n;


  const int NR =( max(w*n, w*m)+ blockdim - 1 )/blockdim;
  const int cmax = (n % N == 0)? n/N: (n+N-1)/N;
  const int rmax = (m % M == 0)? m/M: (m+M-1)/M;

  const int c = threadIdx.x / rmax;
  const int r = threadIdx.x - c * rmax < rmax ? threadIdx.x - c * rmax : rmax ;


  double myc[N*M];
  double mya[3];
  double myb[3];

  double *buff_l, *buff_r;

  int psp;

  __shared__ double buff[buf_sz];
  __shared__ int param_stack_s[4*grouping];
// __shared__ int c_loc;

  buff_l = buff;
  buff_r = &(buff[m*w]);

  int nrun = grouping;

  psp = 7 * (blockIdx.x * grouping);
  //if (blockIdx.x == gridDim.x-1)
  if (blockIdx.x == careful)
    nrun = nruns;


  for (int i = 0; i < N*M; i++)
    myc[i] = 0.0;


  for(int i = threadIdx.x; i < 7*nrun; i += blockdim){
    int p_tmp = __ldg(param_stack + psp + i);

    if (i % 7 > 2)
      param_stack_s[(i / 7)*4 + i % 7 - 3] = p_tmp-1;
  }


  for (int run = 0; run < nrun; run++)
  {
    psp = run*4;

    syncthreads ();
    int srcA = param_stack_s[psp];
    int srcB = param_stack_s[psp+1];

/*
    if(threadIdx.x == 0) {
       c_loc = param_stack_s[psp + 2];
       if( run < nrun-1 && param_stack_s[psp+3] == param_stack_s[psp+3+4]) c_loc = -1;
    }
*/

    if((m*n+31)/(N*M) < m*w) {
        for(int i=threadIdx.x; i < m*w; i += blockdim){
            buff[i] = __ldg(a_data+srcA+i);
        }
    } else {
        if(threadIdx.x < m*w) buff[threadIdx.x] = __ldg(a_data+srcA+threadIdx.x);
        //if(threadIdx.x < m*w) buff[threadIdx.x] = a_data[srcA+threadIdx.x];
    }

    if(w == k) {
        for(int i=threadIdx.x; i < w*n; i+= blockdim){
          buff_r[i] = __ldg(b_data+srcB+i);
        }
    } else if((m*n+31)/(N*M) < w*n) {
        for(int i=threadIdx.x; i < w*n; i += blockdim){
           buff_r[i] = __ldg(b_data+srcB+i);
        }
    } else {

        if(threadIdx.x < w*n){
          buff_r[threadIdx.x] = __ldg(b_data+srcB + threadIdx.x);
          //buff_r[threadIdx.x] = b_data[srcB + threadIdx.x];
        }
    }


    for(int t=0; t < (k/w)*w-w; t += w) {

      syncthreads ();

      srcA += m*w;
      srcB += w*n;

      if((m*n+31)/(N*M) < m*w) {
        int i = threadIdx.x;
        for(int ri=0; ri<NR; ri++){ 
          if(i < m*w)  mya[ri] = __ldg(a_data+srcA+i);
          i+=blockdim;
        }
      } else {
        if(threadIdx.x < m*w) mya[0] = __ldg(a_data+srcA+threadIdx.x);
        //if(threadIdx.x < m*w) mya[0] = a_data[srcA+threadIdx.x];
      }

      if(w == k) {
        for(int i=threadIdx.x, ri=0; i < w*n; i+= blockdim, ri++){
          myb[ri] = __ldg(b_data+srcB+i);
        }
      } else if((m*n+31)/(N*M) < w*n) {
        //for(int i=threadIdx.x, ri=0; i < w*n; i += blockdim, ri++){
        int i = threadIdx.x;
        for(int ri=0; ri<NR;  ri++){
           if(i < w*n) myb[ri] = __ldg(b_data+srcB + i);
           i+=blockdim;
        }
      } else {
        if(threadIdx.x < w*n) myb[0] = __ldg(b_data+srcB+threadIdx.x);
        //if(threadIdx.x < w*n) myb[0] = b_data[srcB+threadIdx.x];
      }

      if (c < cmax && r < rmax)
       {
         for (int l=0 ; l<w; l++) {
           for(int i=0; i<N;i++){
            for(int j=0; j<M; j++){
             myc[M*i+j] += buff_l[l*m + M*r + j] * buff_r[l*n +N*c +i];
            }
           }
         }
      }
      syncthreads();

      if((m*n+31)/(N*M) < m*w) {
      //  for(int i=threadIdx.x, ri=0; i < m*w; i += blockdim, ri++){
      //      buff[i] = mya[ri];
      //  }
         int i= threadIdx.x;
         for(int ri=0; ri<NR; ri++){
           if(i < m*w)  buff[i] = mya[ri];
           i+= blockdim;
         }
      } else {
        if(threadIdx.x < m*w) buff[threadIdx.x] = mya[0];
      }

      if(w == k) {
        for(int i=threadIdx.x, ri=0; i < w*n; i+= blockdim, ri++){
          buff_r[i] = myb[ri];
        }
      } else if((m*n+31)/(N*M) < w*n) {
        int i = threadIdx.x;
        for(int ri=0; ri < NR; ri++){
           if(i < w*n) buff_r[i] = myb[ri];
           i+=blockdim;
        }
      } else {

        if(threadIdx.x < w*n){
          buff_r[threadIdx.x] = myb[0];
        }
      }


    }


    syncthreads();

    const int wa = k - (k/w)*w;

    if(wa != 0){
      srcA += m*w;
      //srcB += w;
      srcB += w*n;

      if((m*n+31)/(N*M) < m*wa) {
//        for(int i=threadIdx.x, ri=0; i < m*wa; i += blockdim, ri++){
//            mya[ri] = a_data[srcA+i];
//        }
       int i = threadIdx.x;
        for(int ri=0; ri<NR; ri++){
          if(i < m*wa)  mya[ri] = __ldg(a_data+srcA+i);
          i+=blockdim;
        }
      } else {

        if(threadIdx.x < m*wa) mya[0] = __ldg(a_data+srcA+threadIdx.x);
        //if(threadIdx.x < m*wa) mya[0] = a_data[srcA+threadIdx.x];
      }

      if(w == k) {
        for(int i=threadIdx.x, ri=0; i < wa*n; i+= blockdim, ri++){

          myb[ri] = __ldg(b_data+srcB+i);

        }
      } else if((m*n+31)/(N*M) < wa*n) {
        int i = threadIdx.x;
        for(int ri=0; ri<NR; ri++){
           if(i < wa*n) myb[ri] = __ldg(b_data+srcB + i);
           i+=blockdim;
        }
      } else {
        if(threadIdx.x < wa*n) myb[0] = __ldg(b_data+srcB+threadIdx.x);
        //if(threadIdx.x < wa*n) myb[0] = b_data[srcB+threadIdx.x];
      }

    }


    if (c < cmax && r < rmax)
       {
         for (int l=0 ; l<w; l++) {
           for(int i=0; i<N;i++){
            for(int j=0; j<M; j++){
             //myc[M*i+j] += buff_l[l*m + M*r + j] * buff_r[(N*c+i)*w +l];
             myc[M*i+j] += buff_l[l*m + M*r + j] * buff_r[l*n + N*c + i];
            }
           }
         }
       }
    syncthreads();


    if(wa != 0) {

      if((m*n+31)/(N*M) < m*wa) {
        int i= threadIdx.x;
         for(int ri=0; ri<NR; ri++){
           if(i < m*wa)  buff[i] = mya[ri];
           i+= blockdim;
         }

      } else {
        if(threadIdx.x < m*wa) buff[threadIdx.x] = mya[0];
      }

     if((m*n+31)/(N*M) < m*wa) {
       int i = threadIdx.x;
        for(int ri=0; ri < NR; ri++){
           if(i < wa*n) buff_r[i] = myb[ri];
           i+=blockdim;
        }
      } else {
        if(threadIdx.x < wa*n){
          buff_r[threadIdx.x] = myb[0];
        }
      }

      syncthreads ();

      if (c < cmax && r < rmax)
        {
          for (int l = 0; l < wa; l++) {
           for(int i=0; i<N;i++)
            for(int j=0; j<M; j++){
             //myc[M*i+j] += buff_l[l*m + M*r + j] * buff_r[(N*c+i)*wa +l];
             myc[M*i+j] += buff_l[l*m + M*r + j] * buff_r[l*n + N*c + i];
            }
          }
        }
    }

 
    syncthreads();

//    if (run == nrun - 1 || c_loc != -1 )
//    {


    if (run == nrun - 1
        || param_stack_s[psp + 3] != param_stack_s[psp + 3 + 4] )
    {
      int c_loc = param_stack_s[psp + 2];

      syncthreads();

/*
      for(int i=0; i<N; i++)
        for(int j=0; j<M; j++){
          if(M*r+j < m && N*c+i < n){
             atomicAdd (&c_data[c_loc + i], myc[M*i+j]);
          }
        }

*/
      for(int t=0; t < (n/v)*v; t+= v){
        int ctmp = c*N - t;
        if ( ctmp >= -(N-1) && ctmp < v)
          {
            for (int i = 0; i < N; i++)

             if(ctmp+i >=0 && ctmp+i< v){
               for (int j = 0; j < M; j++)
                 if (M*r+j < m) {
                    buff[(ctmp+i) * m + M*r+j] = myc[M*i+j];
                    myc[M*i+j]=0.0;
                 }
             }

          }

         syncthreads();


         for(int i = threadIdx.x; i < m*v; i += blockdim){
             atomicAdd (&c_data[c_loc + i], buff[i]);
         }
         c_loc += m*v;
         syncthreads();
      }
      const int va = n - (n/v)*v;

      if(va != 0) {
        int t= (n/v)*v;
        int ctmp = c*N - t;

        if ( ctmp >= -(N-1) && ctmp < va)
          {
            for (int i = 0; i < N; i++)
             if (ctmp+i >= 0 && N*c+i < n){
               for (int j = 0; j < M; j++)
                 if (M*r+j < m) {
                    buff[(ctmp+i) * m + M*r+j] = myc[M*i+j];
                    myc[M*i+j]=0.0;
                 }
             }

          }

         syncthreads();


         for(int i = threadIdx.x; i < m*va; i += blockdim){
             atomicAdd (&c_data[c_loc + i], buff[i]);
          }
         syncthreads();

       }

    }
  }

}
