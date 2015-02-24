/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015  CP2K developers group                         *
 *****************************************************************************/

#include <cuda_runtime.h>
#include <stdio.h>
#include <math.h>
#include "acc_cuda_error.h"

/****************************************************************************/
int cuda_error_check (cudaError_t cudaError){
  if (cudaError != cudaSuccess){
      printf ("CUDA Error: %s\n", cudaGetErrorString (cudaError));
      return -1;
    }
  return 0;
};

//EOF
