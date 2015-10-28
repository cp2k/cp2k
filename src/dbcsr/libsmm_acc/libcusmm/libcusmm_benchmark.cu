/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015  CP2K developers group                         *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include "libcusmm_benchmark.h"

//===========================================================================
static const int n_a = 10000;
static const int n_b = 10000;
static const int n_c = 1000;

//===========================================================================
// Allocate memory and cuda events
void libcusmm_benchmark_init(libcusmm_benchmark_t** handle, int max_m, int max_n, int max_k){
    libcusmm_benchmark_t* h = (libcusmm_benchmark_t*) malloc(sizeof(libcusmm_benchmark_t));
    *handle = h;

    h->max_m = max_m;
    h->max_n = max_n;
    h->max_k = max_k;
    h->max_stack = 16005;

    h->mat_a = (double*) malloc(n_a * max_m * max_k * sizeof(double));
    h->mat_b = (double*) malloc(n_b * max_k * max_n * sizeof(double));
    h->mat_c = (double*) malloc(n_c * max_m * max_n * sizeof(double));
    h->stack = (int*) malloc(h->max_stack * 7 * sizeof(int));

    cudaMalloc(&h->d_mat_a, n_a * max_m * max_k * sizeof(double));
    cudaMalloc(&h->d_mat_b, n_b * max_k * max_n * sizeof(double));
    cudaMalloc(&h->d_mat_c, n_c * max_m * max_n * sizeof(double));
    cudaMalloc(&h->d_stack, h->max_stack * 7 * sizeof(int));

    cudaEventCreate(&h->t_start);
    cudaEventCreate(&h->t_stop);

    cudaError_t cudaError = cudaGetLastError();
    if (cudaError != cudaSuccess){
      printf("libcusmm_benchmark_init: %s\n", cudaGetErrorString(cudaError));
      exit(1);
    }
}

//===========================================================================
// Free memory and cuda events
void libcusmm_benchmark_finalize(libcusmm_benchmark_t* handle){
    cudaEventDestroy(handle->t_stop);
    cudaEventDestroy(handle->t_start);
    cudaFree(handle->d_stack);
    cudaFree(handle->d_mat_c);
    cudaFree(handle->d_mat_b);
    cudaFree(handle->d_mat_a);
    free(handle->stack);
    free(handle->mat_c);
    free(handle->mat_b);
    free(handle->mat_a);
    free(handle);
    cudaError_t cudaError = cudaGetLastError();
    if (cudaError != cudaSuccess){
      printf("libcusmm_benchmark_init: %s\n", cudaGetErrorString(cudaError));
      exit(1);
    }
}

//===========================================================================
// initialize matrix
static void matInit(double* mat, int mat_n, int x, int y, int seed){

 double *m = mat;

 for(int n=0; n<mat_n; n++){
   for(int j=0; j<y; j++) {
     for(int i=0; i<x; i++, m++) {
     *m = (double) j*x + i + n + seed;
     //printf("matrix [%d, %d]=%g\n", i, j, *m);
     }
   }
 }

}


//===========================================================================
// initialize the task list ("stack" in CP2K lingo)
// for each of the result matrices we have a random number
static void stackInit(int *stack, int n_stack, int n_c, double* mat_c,
               int n_a, double * mat_a, int n_b, double* mat_b,
               int mat_m, int mat_n, int mat_k){

  // on average, we have n_avg matrix prodcuts contributing to a result mat_c
  int n_avg = max(1, n_stack / n_c);

  int n_imbalance = max(1, n_avg-4);

  int c = 0;
  int n_top = 0;
  int p = 0;

 int* s = stack;
  while( p < n_stack ){
    if(c >= n_c) c = n_c-1;

    n_top += n_avg + (rand() % (2*n_imbalance) - n_imbalance); 
    if(n_top > n_stack) n_top = n_stack;

    for(;p < n_top; p++){
     int a = rand() % n_a;
     int b = rand() % n_b;

     *s++ =  mat_m;                        // matrix size m
     *s++ =  mat_n;                        // matrix size n
     *s++ =  mat_k;                        // matrix size k
     *s++ =  a * mat_m * mat_k + 1;        // factor 1
     *s++ =  b * mat_k * mat_n + 1;        // factor 2
     *s++ =  c * mat_m * mat_n + 1;        // result
     *s++ =  c * mat_m * mat_n + 1;        // just an identifier..
    }
    c++;
 }
}


//===========================================================================
static void stackCalc(int* stack, int n_stack, double* mat_c, double *mat_a, double* mat_b,
               int mat_m, int mat_n, int mat_k){

  for(int s=0; s<n_stack; s++){
     int a_base = stack[7*s + 3]-1;
     int b_base = stack[7*s + 4]-1;
     int c_base = stack[7*s + 5]-1;

     for(int n=0; n<mat_n; n++){
       for(int m=0; m<mat_m; m++){
         double res = 0.;
         for(int k=0; k<mat_k; k++){
          int a_ind = k * mat_m + m;


//         // initialize with non-transpose matrix
//         int b_ind = n * mat_k + k;
//         res += mat_a[a_base + a_ind] * mat_b[b_base + b_ind];

          // initialize with transpose matrix
          int b_ind = k * mat_n + n;
          res += mat_a[a_base + a_ind] * mat_b[b_base + b_ind];
         }
         int c_ind = n * mat_m +  m;
         mat_c[c_base + c_ind] += res;
       }
     }
  }

}


//===========================================================================
static double checkSum(double* mat_c, int n_c, int mat_m, int mat_n){
   double res = 0;
   for(int i=0; i<n_c * mat_m * mat_n; i++){
     res += mat_c[i];
   }
   return res;
}



//===========================================================================
//Removes special symbols so that the output is usefull for awk and gnuplot.
static void clean_string(char* str_in, char* str_out){
    for(int i=0; i<1000 ; i++){
        if(str_in[i] == '=' || str_in[i] == ',' || str_in[i] == '(' || str_in[i] == ')'){
            str_out[i] = ' ';
         }else{
             str_out[i] = str_in[i];
         }
         if(str_in[i] == 0)
             break;
    }
}


//===========================================================================
int libcusmm_benchmark(libcusmm_benchmark_t* handle,
                       int mat_m, int mat_n, int mat_k,
                       int nkernels, KernelLauncher* launchers, char ** kernel_descr, bool tune_mode){
  int stack_n, n_iter;

  if(tune_mode){
      // for larger matrices we get enough statistics from fewer iterations
      n_iter = max(3, 1250/(mat_m * mat_n * mat_k));
      stack_n = 16005;
  }else{
      n_iter = 1;
      stack_n = 10  ;
  }

  if(mat_m>handle->max_m ||mat_n>handle->max_n || mat_k>handle->max_k|| stack_n>handle->max_stack){
      printf("libcusmm_benchmark: got handle with too few resources");
      exit(1);
  }

 const int stream = 0;

 int error_counter = 0;
 int best_kernel = -1;
 double best_gflops = 0.0;
 double sumCPU, sumGPU;
 float t_duration;
 char descr[1000], msg_prefix[100]="";
 cudaError_t cudaError;

 memset(handle->mat_c, 0, n_c * mat_m * mat_n * sizeof(double));
 matInit(handle->mat_a, n_a, mat_m, mat_k, 42);
 matInit(handle->mat_b, n_b, mat_k, mat_n, 24);

 if(tune_mode)
     printf("Initializing ...\n");
 stackInit(handle->stack, stack_n, n_c, handle->mat_c, n_a, handle->mat_a, n_b, handle->mat_b, mat_m, mat_n, mat_k);

 // Actually, we would have to calculate the stack n_iter times.
 // We cheat by simply scaling the results of a single stack calulcation.
 stackCalc(handle->stack, stack_n, handle->mat_c, handle->mat_a, handle->mat_b, mat_m, mat_n, mat_k);
 for(int i=0 ; i < n_c*mat_m*mat_n ; i++)
     handle->mat_c[i] *= n_iter;

 sumCPU =  checkSum(handle->mat_c, n_c, mat_m, mat_n);

 cudaMemcpy(handle->d_mat_a, handle->mat_a, n_a * mat_m * mat_k * sizeof(double), cudaMemcpyHostToDevice);
 cudaMemcpy(handle->d_mat_b, handle->mat_b, n_b * mat_k * mat_n * sizeof(double), cudaMemcpyHostToDevice);
 cudaMemcpy(handle->d_stack, handle->stack, stack_n * 7 * sizeof(int), cudaMemcpyHostToDevice);
 //d_mat_c get's zeroed after warmup run

 for(int ikern=0; ikern < nkernels; ikern++){
    //warmup run
    launchers[ikern](handle->d_stack, stack_n, stream, mat_m, mat_n, mat_k, handle->d_mat_a, handle->d_mat_b, handle->d_mat_c);
    cudaMemset(handle->d_mat_c, stream, n_c * mat_m * mat_n * sizeof(double));

    cudaEventRecord(handle->t_start, stream);

    for(int i=0; i<n_iter; i++)
        launchers[ikern](handle->d_stack, stack_n, stream, mat_m, mat_n, mat_k, handle->d_mat_a, handle->d_mat_b, handle->d_mat_c);

    cudaEventRecord(handle->t_stop, stream);
    cudaEventSynchronize(handle->t_stop);
    cudaEventElapsedTime(&t_duration, handle->t_start, handle->t_stop);

    cudaMemcpy(handle->mat_c, handle->d_mat_c, n_c * mat_m * mat_n * sizeof(double), cudaMemcpyDeviceToHost);

    clean_string(kernel_descr[ikern], descr);

    if(tune_mode)
        sprintf(msg_prefix, "params %d / %d",ikern+1, nkernels);

    cudaError = cudaGetLastError();
    if (cudaError != cudaSuccess){
      printf("%sERROR %s cuda_error: %s\n", msg_prefix, descr, cudaGetErrorString(cudaError));
      error_counter++;
      continue;
    }

    sumGPU =  checkSum(handle->mat_c, n_c, mat_m, mat_n);
    if(sumGPU != sumCPU){
        printf("%sERROR %s checksum_diff: %g\n",msg_prefix, descr, sumGPU-sumCPU);
        error_counter++;
        continue;
    }

    if(tune_mode){
       double gflops = ((double) n_iter *stack_n * mat_m * mat_n * mat_k * 2 / (1e9))/(t_duration * 1e-3);
       printf("%sOK %s GFlop/s %g\n", msg_prefix, descr, gflops);
       if(best_gflops < gflops){
           best_gflops = gflops;
           best_kernel = ikern;
       }
    }else{
       printf("%sOK %s\n", msg_prefix, descr);
    }
 }

 if(tune_mode){
    printf("\n\n");
    if(best_kernel > -1){
        printf("WINNER: %d %s , # %g GFlop/s \n", best_kernel+1, kernel_descr[best_kernel], best_gflops);
    }else{
       printf("WINNER: None\n");
    }
    printf("Number of errors: %d\n", error_counter);
    cudaDeviceReset();
 }

 // cleanup
 return(error_counter);
}

//EOF
