/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015  CP2K developers group                         *
 *****************************************************************************/

#if defined ( __PW_CUDA )

/* This file contains memory management routines for device memory. */ 
// Revised: Sept. 2012, Andreas Gloess
// Author:  Benjamin G Levine

// global dependencies
#include <cuda_runtime.h>
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>

// local dependencies
#include "fft_cuda.h"

// debug flag
#define CHECK 1
#define VERBOSE 0

// --- CODE --------------------------------------------------------------------
static const int    nstreams      = 3;
static const int    nevents       = 2;

cudaError_t          cErr;
static cudaStream_t *cuda_streams;
static cudaEvent_t  *cuda_events;
static int           is_configured = 0;

extern void pw_cuda_error_check (cudaError_t cudaError, int line) {
  int         pid;
  size_t      free, total;
  cudaError_t cErr2;

  cErr2 = cudaGetLastError();
  if (cudaError != cudaSuccess || cErr2 != cudaSuccess) {
    pid = getpid();
    printf("%d CUDA RT Error line %d\n", pid, line);
    printf("%d CUDA RT1 Error: %s\n", pid, cudaGetErrorString(cudaError));
    printf("%d CUDA RT2 Error: %s\n", pid, cudaGetErrorString(cErr2));
    cudaMemGetInfo(&free,&total);
    printf("%d Free: %zu , Total: %zu\n", pid, free, total);
    fflush(stdout);
    exit(-1);
  }
}

// STREAMS INIT/GET/RELEASE
void pw_cuda_device_streams_alloc (cudaStream_t **streams) {
  cudaStream_t *cuda_streams_ptr;
  cuda_streams_ptr = (cudaStream_t *) malloc(nstreams * sizeof(cudaStream_t));
  for (int i = 0; i < nstreams; i++) {
    cErr = cudaStreamCreateWithFlags(&cuda_streams_ptr[i], cudaStreamNonBlocking);
    if (CHECK) pw_cuda_error_check (cErr, __LINE__);
  }
  *streams = cuda_streams_ptr;
}

extern void pw_cuda_get_streams (cudaStream_t **streams) {
  *streams = cuda_streams;
}

void pw_cuda_device_streams_release (cudaStream_t **streams) {
  cudaStream_t *cuda_streams_ptr;
  cuda_streams_ptr = *streams;
  for (int i = 0; i < nstreams; i++) {
     cErr = cudaStreamDestroy(cuda_streams_ptr[i]);
     if (CHECK) pw_cuda_error_check (cErr, __LINE__);
  }
  free(cuda_streams_ptr);
  cuda_streams_ptr = NULL;
}

// EVENTS INIT/GET/RELEASE
void pw_cuda_device_events_alloc (cudaEvent_t **events) {
  cudaEvent_t *cuda_events_ptr;
  cuda_events_ptr = (cudaEvent_t *) malloc(nevents * sizeof(cudaEvent_t));
  for (int i =0; i < nevents; i++) {
    cErr = cudaEventCreateWithFlags(&cuda_events_ptr[i], cudaEventDisableTiming);
    //cErr = cudaEventCreateWithFlags(&(cuda_events_ptr[i]), cudaEventDefault);
    //cErr = cudaEventCreateWithFlags(&(cuda_events_ptr[i]), cudaEventBlockingSync);
    if (CHECK) pw_cuda_error_check (cErr, __LINE__);
  }
  *events = cuda_events_ptr;
}

extern void pw_cuda_get_events (cudaEvent_t **events) {
  *events = cuda_events;
}

void pw_cuda_device_events_release (cudaEvent_t **events) {
  cudaEvent_t *cuda_events_ptr;
  cuda_events_ptr = *events;
  for (int i = 0; i < nevents; i++) {
    cErr = cudaEventDestroy(cuda_events_ptr[i]);
    if (CHECK) pw_cuda_error_check (cErr, __LINE__);
  }
  free(cuda_events_ptr);
  cuda_events_ptr = NULL;
}

// MEMORY ALLOC/RELEASE
extern void pw_cuda_device_mem_alloc (int **ptr, int n) {
  cErr = cudaMalloc((void **) ptr, (size_t) sizeof(int)*n);
  if (CHECK) pw_cuda_error_check (cErr, __LINE__);
}

extern void pw_cuda_device_mem_alloc (float **ptr, int n) {
  cErr = cudaMalloc((void **) ptr, (size_t) sizeof(float)*n);
  if (CHECK) pw_cuda_error_check (cErr, __LINE__);
}

extern void pw_cuda_device_mem_alloc (double **ptr, int n) {
  cErr = cudaMalloc((void **) ptr, (size_t) sizeof(double)*n);
  if (CHECK) pw_cuda_error_check (cErr, __LINE__);
}

extern void pw_cuda_device_mem_free (int **ptr) {
  cErr = cudaFree((void *) *ptr);
  if (CHECK) pw_cuda_error_check (cErr, __LINE__);
  *ptr = NULL;
}

extern void pw_cuda_device_mem_free (float **ptr) {
  cErr = cudaFree((void *) *ptr);
  if (CHECK) pw_cuda_error_check (cErr, __LINE__);
  *ptr = NULL;
}

extern void pw_cuda_device_mem_free (double **ptr) {
  cErr = cudaFree((void *) *ptr);
  if (CHECK) pw_cuda_error_check (cErr, __LINE__);
  *ptr = NULL;
}

// INIT/RELEASE
extern "C" void pw_cuda_init (int memory) {
  if ( is_configured == 0 ) {
    pw_cuda_device_streams_alloc (&cuda_streams);
    pw_cuda_device_events_alloc (&cuda_events);
    is_configured = 1;
  }
}

extern "C" void pw_cuda_finalize () {
  if ( is_configured == 1 ) {
    fftcu_release_();
    pw_cuda_device_streams_release (&cuda_streams);
    pw_cuda_device_events_release (&cuda_events);
    is_configured = 0;
  }
}

#endif
