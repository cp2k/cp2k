/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2014 the CP2K developers group                      *
 *****************************************************************************/

#include <cuda_runtime.h>
#include <stdio.h>
#include <math.h>
#include "acc_cuda_error.h"
#include "../include/acc.h"

#ifdef __CUDA_PROFILING
#include <nvToolsExtCudaRt.h>
#endif

static const int verbose_print = 0;


/****************************************************************************/
extern "C" int acc_stream_priority_range(int* least, int* greatest){
  *least = -1;
  *greatest = -1;

#ifndef __HAS_NO_CUDA_STREAM_PRIORITIES
  cudaError_t cErr = cudaDeviceGetStreamPriorityRange(least, greatest);
  if (cuda_error_check(cErr)) return -1;
  if (cuda_error_check(cudaGetLastError())) return -1;
#endif

  return 0;
}


/****************************************************************************/
extern "C" int acc_stream_create(void** stream_p, char* name, int priority){
  cudaError_t cErr;
  *stream_p = malloc(sizeof(cudaStream_t));

  cudaStream_t* custream = (cudaStream_t*) *stream_p;

#ifndef __HAS_NO_CUDA_STREAM_PRIORITIES
  if(priority > 0){
      unsigned int flags = cudaStreamNonBlocking;
      cErr =  cudaStreamCreateWithPriority(custream, flags, priority);
  }else
#endif
      cErr = cudaStreamCreate(custream);


  if (verbose_print) printf("cuda_stream_create: %p -> %d \n", *stream_p, *custream);
  if (cuda_error_check(cErr)) return -1;
  if (cuda_error_check(cudaGetLastError())) return -1;

#ifdef __CUDA_PROFILING
  nvtxNameCudaStreamA(*custream, name);
#endif

    return 0;
}


/****************************************************************************/
extern "C" int acc_stream_destroy(void* stream){
    cudaStream_t* custream = (cudaStream_t*) stream;

    if(verbose_print) printf("cuda_stream_destroy called\n");
    cudaError_t cErr = cudaStreamDestroy(*custream);
    free(custream);
    if (cuda_error_check (cErr)) return -1;
    if (cuda_error_check(cudaGetLastError ()))return -1;
    return 0;
}

/****************************************************************************/
extern "C" int acc_stream_sync(void* stream)
{
  cudaStream_t* custream = (cudaStream_t*) stream;
  cudaError_t cErr = cudaStreamSynchronize(*custream);
  if (cuda_error_check (cErr)) return -1;
  return 0;
}

//EOF
