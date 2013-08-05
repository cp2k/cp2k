/******************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations
 *  Copyright (C) 2000 - 2013 the CP2K developers group
 *****************************************************************************/

#include <cuda_runtime.h>
#include "error_cuda.h"

//==============================================================================
extern "C" int cp_set_device_cu (int device_id)
{
  cudaError_t cErr;
  int myDevice;

  cErr = cudaSetDevice (device_id);
  cuda_error_check2(cErr, __LINE__);

  cErr = cudaGetDevice (&myDevice);
  cuda_error_check2(cErr, __LINE__);

  if (myDevice != device_id)
    return 1;

  //cErr = cudaGetDeviceProperties (&devProperties, myDevice);
  //cuda_error_check2(cErr, __LINE__);

  return 0;
}


//==============================================================================
extern "C" int cp_get_ndevices_cu (int *n_devices)
{
  cudaError_t cErr;

  cErr = cudaGetDeviceCount (n_devices);
  cuda_error_check2(cErr, __LINE__);
  return 0;
}


//==============================================================================
extern "C" int cp_device_sync_cu ()
{
  cudaError_t cErr;

  cErr = cudaDeviceSynchronize ();
  cuda_error_check2(cErr, __LINE__);
  return 0;
}


//==============================================================================
extern "C" int cp_device_reset_cu ()
{
  cudaError_t cErr;

  cErr = cudaDeviceReset();
  cuda_error_check2(cErr, __LINE__);
  return 0; 
}

//EOF
