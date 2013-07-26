/******************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations
 *  Copyright (C) 2000 - 2013 CP2K developers group
 *****************************************************************************/

#include <cuda_runtime.h>
#include <stdio.h>

#include "dbcsr_cuda.h"
#include <math.h>

  static const int verbose_print = 0;

extern "C" int cuda_stream_create(cudaStream_t** stream_p){
  if(verbose_print) printf("cuda_stream_create called\n");
  *stream_p = (cudaStream_t*) malloc(sizeof(cudaStream_t));
  cudaError_t cErr = cudaStreamCreate(*stream_p);
  if (cuda_error_check(cErr)) return 1;
  if (cuda_error_check(cudaGetLastError())) return 1;
  return 0;
}

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
