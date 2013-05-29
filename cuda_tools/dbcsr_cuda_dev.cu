/******************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations
 *  Copyright (C) 2000 - 2013  Urban Borstnik and the CP2K developers group
 *****************************************************************************/

#include <cuda_runtime.h>
#include <stdio.h>
#include <sm_11_atomic_functions.h>

#include "dbcsr_cuda.h"

static cudaStream_t *streams = 0;
static int nStreams = 0;

struct cudaDeviceProp devProperties;
#pragma omp threadprivate(devProperties)

//static const int verbose_print = 0;

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
dc_device_reset_cu ()
{
  cudaError_t cErr;

  cErr = cudaDeviceReset();
  if (cuda_error_check (cErr))
    return 1;
  return 0;
}

extern "C" int
dc_stream_sync_cu (int stream_id)
{
  cudaError_t cErr;
  cudaStream_t stream;

  stream = (cudaStream_t) dc_get_stream (stream_id);
  cErr = cudaStreamSynchronize (stream);
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

cudaStream_t
dc_get_stream (int stream_id)
{
  if (nStreams == 0)
    return (cudaStream_t) 0;
  else
    return streams[stream_id];
}


extern "C" int
dc_create_streams (int n_streams)
{
  cudaError_t cErr;
  int i;

  nStreams = n_streams;
  streams = (cudaStream_t *) malloc (sizeof (cudaStream_t) * (n_streams + 1));
  if (streams == NULL)
    return 2;
  for (i = 1; i <= n_streams; i++)
    {
      cErr = cudaStreamCreate (&(streams[i]));
      if (cuda_error_check (cErr))
	{
	  free ((void *) streams);
	  return 1;
	}
    }
  streams[0] = (cudaStream_t) 0;
  return 0;
}

extern "C" int
dc_destroy_streams ()
{
  cudaError_t cErr;
  int i;

  for (i = 1; i <= nStreams; i++)
    {
      cErr = cudaStreamDestroy (streams[i]);
      if (cuda_error_check (cErr))
	{
	  free ((void *) streams);
	  return 1;
	}
    }

  free ((void *) streams);
  return 0;
}
