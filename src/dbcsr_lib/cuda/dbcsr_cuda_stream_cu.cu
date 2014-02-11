/******************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations
 *  Copyright (C) 2000 - 2013 CP2K developers group
 *****************************************************************************/

#include <cuda_runtime.h>
#include <stdio.h>

#include "dbcsr_cuda.h"
#include <math.h>

#ifdef __CUDA_PROFILING
#include <nvToolsExtCudaRt.h>
#endif

  static const int verbose_print = 0;

extern "C" int cuda_stream_create(cudaStream_t** stream_p){
  *stream_p = (cudaStream_t*) malloc(sizeof(cudaStream_t));
  cudaError_t cErr = cudaStreamCreate(*stream_p);
  if(verbose_print) printf("cuda_stream_create: %p -> %d \n", *stream_p, **stream_p);
  if (cuda_error_check(cErr)) return 1;
  if (cuda_error_check(cudaGetLastError())) return 1;
  return 0;
}

#ifndef __HAS_NO_CUDA_STREAM_PRIORITIES
extern "C" int cuda_stream_create_with_priority(cudaStream_t** stream_p, int priority){
  if(verbose_print) printf("cuda_stream_create_with_priority called\n");
  *stream_p = (cudaStream_t*) malloc(sizeof(cudaStream_t));
  unsigned int flags = cudaStreamNonBlocking;
  cudaError_t cErr =  cudaStreamCreateWithPriority(*stream_p, flags, priority);
  if (cuda_error_check(cErr)) return 1;
  if (cuda_error_check(cudaGetLastError())) return 1;
  return 0;
}

extern "C" int cuda_stream_priority_range(int* least, int* greatest){
  cudaError_t cErr = cudaDeviceGetStreamPriorityRange(least, greatest);
  if (cuda_error_check(cErr)) return 1;
  if (cuda_error_check(cudaGetLastError())) return 1;
  return 0;
}
#endif

extern "C" int cuda_stream_destroy(cudaStream_t* stream){
    if(verbose_print) printf("cuda_stream_destroy called\n");
    cudaError_t cErr = cudaStreamDestroy(*stream);
    free(stream);
    if (cuda_error_check (cErr)) return 1;
    if (cuda_error_check(cudaGetLastError ()))return 1;
    return 0;
}

extern "C" int cuda_stream_sync(cudaStream_t* stream)
{
  cudaError_t cErr;
  cErr = cudaStreamSynchronize(*stream);
  if (cuda_error_check (cErr))
    return 1;
  return 0;
}

#ifdef __CUDA_PROFILING
extern "C" void cuda_stream_set_name(cudaStream_t* stream_p, char* name){
  nvtxNameCudaStreamA(*stream_p, name);
}
#endif
