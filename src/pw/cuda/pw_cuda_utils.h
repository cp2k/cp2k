/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015  CP2K developers group                         *
 *****************************************************************************/

#ifndef PW_CUDA_UTILS_H
#define PW_CUDA_UTILS_H

extern void pw_cuda_error_check (cudaError_t cudaError, int line);

// STREAMS INIT/GET/RELEASE
extern void pw_cuda_device_streams_alloc (cudaStream_t **streams);
extern void pw_cuda_get_streams (cudaStream_t **streams);
extern void pw_cuda_device_streams_release (cudaStream_t **streams);

// EVENTS INIT/GET/RELEASE
extern void pw_cuda_device_events_alloc (cudaEvent_t **events);
extern void pw_cuda_get_events (cudaEvent_t **events);
extern void pw_cuda_device_events_release (cudaEvent_t **events);

// MEMORY ALLOC/RELEASE
extern void pw_cuda_device_mem_alloc (int **ptr, int n);
extern void pw_cuda_device_mem_alloc (float **ptr, int n);
extern void pw_cuda_device_mem_alloc (double **ptr, int n);
extern void pw_cuda_device_mem_free (int **ptr);
extern void pw_cuda_device_mem_free (float **ptr);
extern void pw_cuda_device_mem_free (double **ptr);

// DEVICE INIT/RELEASE
extern "C" void pw_cuda_device_mem_init (int memory);
extern "C" void pw_cuda_device_mem_release ();

#endif
