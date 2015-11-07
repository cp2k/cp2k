/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015  CP2K developers group                         *
 *****************************************************************************/

/******************************************************************************
 *  Author: Andreas Gloess
 *
 *****************************************************************************/
#if defined ( __PW_CUDA )

// global dependencies
#include <cuda_runtime.h>
#include <cufft.h>
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>

// debug flag
#define CHECK 1
#define VERBOSE 0

// --- CODE -------------------------------------------------------------------


extern void cufft_error_check (cufftResult_t cufftError, int line) {
  int         pid;
  size_t      free, total;
  cudaError_t cErr2;

  cErr2 = cudaGetLastError();
  if (cufftError != CUFFT_SUCCESS || cErr2 != cudaSuccess) {
    pid = getpid();
    printf("%d CUDA FFT Error line: %d \n", pid, line);
    switch (cufftError) {
      case CUFFT_INVALID_PLAN:   printf("%d CUDA FFT1 Error (CUFFT_INVALID_PLAN)\n", pid); break;
      case CUFFT_ALLOC_FAILED:   printf("%d CUDA FFT1 Error (CUFFT_ALLOC_FAILED)\n", pid); break;
      case CUFFT_INVALID_VALUE:  printf("%d CUDA FFT1 Error (CUFFT_INVALID_VALUE)\n", pid); break;
      case CUFFT_INTERNAL_ERROR: printf("%d CUDA FFT1 Error (CUFFT_INTERNAL_ERROR)\n", pid); break;
      case CUFFT_EXEC_FAILED:    printf("%d CUDA FFT1 Error (CUFFT_EXEC_FAILED)\n", pid); break;
      case CUFFT_INVALID_SIZE:   printf("%d CUDA FFT1 Error (CUFFT_INVALID_SIZE)\n", pid); break;
      default: printf("%d CUDA FFT1 Error (--unimplemented--) %d %d\n", pid, cufftError, cErr2); break;
    }
    printf("%d CUDA FFT2 Error %s \n", pid, cudaGetErrorString(cErr2));
    cudaMemGetInfo(&free,&total);
    printf("%d Free: %zu , Total: %zu\n", pid, free, total);
    fflush(stdout);
    exit(-1);
  }
}

#endif
