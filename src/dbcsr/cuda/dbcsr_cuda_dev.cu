/******************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations
 *  Copyright (C) 2000 - 2013  Urban Borstnik and the CP2K developers group
 *****************************************************************************/

#include <cuda_runtime.h>
#include <stdio.h>
#include <sm_11_atomic_functions.h>

#include "dbcsr_cuda.h"


struct cudaDeviceProp devProperties;
#pragma omp threadprivate(devProperties)


int cuda_error_check (cudaError_t cudaError){
  if (cudaError != cudaSuccess){
      printf ("CUDA Error: %s\n", cudaGetErrorString (cudaError));
      return 1;
    }
  return 0;
};


extern "C" int
dc_device_sync_cu ()
{
  cudaError_t cErr;

  cErr = cudaDeviceSynchronize ();
  if (cuda_error_check (cErr))
    return 1;
  return 0;
}


extern "C" int
dc_set_device_cu (int device_id)
{
  cudaError_t cErr;
  int myDevice;

  cErr = cudaSetDevice (device_id);
  if (cuda_error_check (cErr))
    return 1;

  cErr = cudaGetDevice (&myDevice);
  if (cuda_error_check (cErr))
    return 1;

  if (myDevice != device_id)
    return 1;

  cErr = cudaGetDeviceProperties (&devProperties, myDevice);
  if (cuda_error_check (cErr))
    return 1;

  return 0;
}

extern "C" int
dc_get_ndevices_cu (int *n_devices)
{
  cudaError_t cErr;

  cErr = cudaGetDeviceCount (n_devices);
  if (cuda_error_check (cErr))
    return 1;
  return 0;
}

extern "C" int dc_set_shared_mem_config(int config){
  if(config==0)
      return cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeDefault);
  if(config==1)
      return cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeFourByte);
  if(config==2)
      return cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);
  return(-1);
}

