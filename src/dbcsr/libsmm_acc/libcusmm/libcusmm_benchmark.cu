/*****************************************************************************
*  CP2K: A general program to perform molecular dynamics simulations         *
*  Copyright (C) 2000 - 2014 the CP2K developers group                       *
*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include "libcusmm_benchmark.h"


//===========================================================================
// allocate matrix
static void matAlloc(double** mat, int mat_n, int x, int y){

 double *m = (double*) malloc(mat_n * x * y * sizeof(double));
 *mat = m;
 memset(m, 0, mat_n * x * y * sizeof(double));
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
static void stackInit(int **stack, int n_stack, int n_c, double* mat_c,
               int n_a, double * mat_a, int n_b, double* mat_b,
               int mat_m, int mat_n, int mat_k){


  int* s = (int*) malloc(n_stack * 7 * sizeof(int));
  *stack = s;

  if(n_stack < n_c){
    printf("Error: n_stack < n_c\n");
    exit(1);
  }

  int n_avg = n_stack / n_c;  // on average, we have n_avg matrix prodcuts contributing to 
                              // a result mat_c

  int n_imbalance = max(1, n_avg-4);

  int c = 0;
  int n_top = 0;
  int p = 0;

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
int libcusmm_benchmark(int mat_m, int mat_n, int mat_k, int nkernels, KernelLauncher* launchers, char ** kernel_descr, bool tune_mode){
const int stack_n = 16005;
//const int stack_n = 14336;

  // for larger matrices we get enough statistics from fewer iterations:
 int n_iter = max(3, 1250/(mat_m * mat_n * mat_k));

 const int stream = 0;
 const int n_a = 10000;
 const int n_b = 10000;
//const int n_a = 1;
//const int n_b = 1;
 const int n_c = 1000;
// const int n_c = 1;

// const int n_a = 100;
// const int n_b = 100;
// const int n_a = 1;
// const int n_b = 1;
// const int n_c = 100;

 double* mat_a;
 double* mat_b;
 double* mat_c;
 int* stack;

 double* d_mat_a;
 double* d_mat_b;
 double* d_mat_c;
 int* d_stack;

 int error_counter = 0;
 int best_kernel = -1;
 double best_gflops = 0.0;
 double sumCPU, sumGPU;
 float t_duration;
 char descr[1000], msg_prefix[100]="";
 cudaError_t cudaError;

 matAlloc(&mat_a, n_a, mat_m, mat_k);
 matAlloc(&mat_b, n_b, mat_k, mat_n);
 matAlloc(&mat_c, n_c, mat_m, mat_n);

 matInit(mat_a, n_a, mat_m, mat_k, 42);
 matInit(mat_b, n_b, mat_k, mat_n, 24);

 if(tune_mode)
     printf("Initializing ...\n");
 stackInit(&stack, stack_n, n_c, mat_c, n_a, mat_a, n_b, mat_b, mat_m, mat_n, mat_k);

 for(int i=0; i<n_iter; i++)
   stackCalc(stack, stack_n, mat_c, mat_a, mat_b, mat_m, mat_n, mat_k);

 sumCPU =  checkSum(mat_c, n_c, mat_m, mat_n);

 cudaMalloc(&d_mat_a, n_a * mat_m * mat_k * sizeof(double));
 cudaMemcpy(d_mat_a, mat_a, n_a * mat_m * mat_k * sizeof(double), cudaMemcpyHostToDevice);

 cudaMalloc(&d_mat_b, n_b * mat_k * mat_n * sizeof(double));
 cudaMemcpy(d_mat_b, mat_b, n_b * mat_k * mat_n * sizeof(double), cudaMemcpyHostToDevice);

 cudaMalloc(&d_mat_c, n_c * mat_m * mat_n * sizeof(double));
 cudaMemset(d_mat_c, 0, n_c * mat_m * mat_n * sizeof(double));

 cudaMalloc(&d_stack, stack_n * 7 * sizeof(int));
 cudaMemcpy(d_stack, stack, stack_n * 7 * sizeof(int), cudaMemcpyHostToDevice);

 cudaEvent_t t_start, t_stop;
 cudaEventCreate(&t_start);
 cudaEventCreate(&t_stop);

 for(int ikern=0; ikern < nkernels; ikern++){
    //warmup run
    launchers[ikern](d_stack, stack_n, stream, mat_m, mat_n, mat_k, d_mat_a, d_mat_b, d_mat_c);

    cudaMemset(d_mat_c, stream, n_c * mat_m * mat_n * sizeof(double));
    cudaEventRecord(t_start, stream);

    for(int i=0; i<n_iter; i++)
        launchers[ikern](d_stack, stack_n, stream, mat_m, mat_n, mat_k, d_mat_a, d_mat_b, d_mat_c);

    cudaEventRecord(t_stop, stream);
    cudaEventSynchronize(t_stop);
    cudaEventElapsedTime(&t_duration, t_start, t_stop);

    memset(mat_c, stream, n_c * mat_m * mat_n *sizeof(double));
    cudaMemcpy(mat_c, d_mat_c, n_c * mat_m * mat_n * sizeof(double), cudaMemcpyDeviceToHost);

    clean_string(kernel_descr[ikern], descr);

    if(tune_mode)
        sprintf(msg_prefix, "params %d / %d",ikern+1, nkernels);

    cudaError = cudaGetLastError();
    if (cudaError != cudaSuccess){
      printf("%sERROR %s cuda_error: %s\n", msg_prefix, descr, cudaGetErrorString(cudaError));
      error_counter++;
      continue;
    }

    sumGPU =  checkSum(mat_c, n_c, mat_m, mat_n);
    if(sumGPU != sumCPU){
        printf("%sERROR %s checksum_diff: %g\n",msg_prefix, descr, sumGPU-sumCPU);
        error_counter++;
        continue;
    }

    double gflops = ((double) n_iter *stack_n * mat_m * mat_n * mat_k * 2 / (1e9))/(t_duration * 1e-3);
    printf("%sOK %s GFlop/s %g\n", msg_prefix, descr, gflops);
    if(best_gflops < gflops){
        best_gflops = gflops;
        best_kernel = ikern;
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
 }

 // cleanup
 cudaEventDestroy(t_stop);
 cudaEventDestroy(t_start);
 cudaFree(d_stack);
 cudaFree(d_mat_c);
 cudaFree(d_mat_b);
 cudaFree(d_mat_a);
 cudaDeviceReset();
 free(stack);
 free(mat_c);
 free(mat_b);
 free(mat_a);
 return(error_counter);
}

//EOF
