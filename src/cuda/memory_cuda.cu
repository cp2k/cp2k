#if defined ( __PW_CUDA ) || defined ( __CUBLASDP )

/* This file contains memory management routines for device memory. */ 
// Revised: Sept. 2012, Andreas Gloess
// Author:  Benjamin G Levine

// global dependencies
#include <cuda_runtime.h>
#include <stdio.h>

// local dependencies
#include "error_cuda.h"

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


// STREAMS INIT/GET/RELEASE
extern void cuda_device_streams_alloc_cu_ (cudaStream_t **streams) {
  cudaStream_t *cuda_streams_ptr;
  cuda_streams_ptr = (cudaStream_t *) malloc(nstreams * sizeof(cudaStream_t));
  for (int i = 0; i < nstreams; i++) {
    cErr = cudaStreamCreateWithFlags(&cuda_streams_ptr[i], cudaStreamNonBlocking);
    if (CHECK) cuda_error_check2 (cErr, __LINE__);
  }
  *streams = cuda_streams_ptr;
}

extern void cuda_get_streams_cu_ (cudaStream_t **streams) {
  *streams = cuda_streams;
}

extern void cuda_device_streams_release_cu_ (cudaStream_t **streams) {
  cudaStream_t *cuda_streams_ptr;
  cuda_streams_ptr = *streams;
  for (int i = 0; i < nstreams; i++) {
     cErr = cudaStreamDestroy(cuda_streams_ptr[i]);
     if (CHECK) cuda_error_check2 (cErr, __LINE__);
  }
  cuda_streams_ptr = NULL;
  free(cuda_streams_ptr);
}

// EVENTS INIT/GET/RELEASE
extern void cuda_device_events_alloc_cu_ (cudaEvent_t **events) {
  cudaEvent_t *cuda_events_ptr;
  cuda_events_ptr = (cudaEvent_t *) malloc(nevents * sizeof(cudaEvent_t));
  for (int i =0; i < nevents; i++) {
    cErr = cudaEventCreateWithFlags(&cuda_events_ptr[i], cudaEventDisableTiming);
    //cErr = cudaEventCreateWithFlags(&(cuda_events_ptr[i]), cudaEventDefault);
    //cErr = cudaEventCreateWithFlags(&(cuda_events_ptr[i]), cudaEventBlockingSync);
    if (CHECK) cuda_error_check2 (cErr, __LINE__);
  }
  *events = cuda_events_ptr;
}

extern void cuda_get_events_cu_ (cudaEvent_t **events) {
  *events = cuda_events;
}

extern void cuda_device_events_release_cu_ (cudaEvent_t **events) {
  cudaEvent_t *cuda_events_ptr;
  cuda_events_ptr = *events;
  for (int i = 0; i < nevents; i++) {
    cErr = cudaEventDestroy(cuda_events_ptr[i]);
    if (CHECK) cuda_error_check2 (cErr, __LINE__);
  }
  cuda_events_ptr = NULL;
  free(cuda_events_ptr);
}

// MEMORY ALLOC/RELEASE
extern void cuda_device_mem_alloc_cu_ (int **ptr, int n) {
  cErr = cudaMalloc((void **) ptr, (size_t) sizeof(int)*n);
  if (CHECK) cuda_error_check2 (cErr, __LINE__);
}

extern void cuda_device_mem_alloc_cu_ (float **ptr, int n) {
  cErr = cudaMalloc((void **) ptr, (size_t) sizeof(float)*n);
  if (CHECK) cuda_error_check2 (cErr, __LINE__);
}

extern void cuda_device_mem_alloc_cu_ (double **ptr, int n) {
  cErr = cudaMalloc((void **) ptr, (size_t) sizeof(double)*n);
  if (CHECK) cuda_error_check2 (cErr, __LINE__);
}

extern void cuda_device_mem_free_cu_ (int **ptr) {
  cErr = cudaFree((void *) *ptr);
  if (CHECK) cuda_error_check2 (cErr, __LINE__);
  *ptr = NULL;
}

extern void cuda_device_mem_free_cu_ (float **ptr) {
  cErr = cudaFree((void *) *ptr);
  if (CHECK) cuda_error_check2 (cErr, __LINE__);
  *ptr = NULL;
}

extern void cuda_device_mem_free_cu_ (double **ptr) {
  cErr = cudaFree((void *) *ptr);
  if (CHECK) cuda_error_check2 (cErr, __LINE__);
  *ptr = NULL;
}

// DEVICE INIT/RELEASE
extern "C" void cuda_device_mem_init_cu_ (int memory) {
  if ( is_configured == 0 ) {
    cuda_device_streams_alloc_cu_ (&cuda_streams);
    cuda_device_events_alloc_cu_ (&cuda_events);
    is_configured = 1;
  }
}

extern "C" void cuda_device_mem_release_cu_ () {
  if ( is_configured == 1 ) {
    cuda_device_streams_release_cu_ (&cuda_streams);
    cuda_device_events_release_cu_ (&cuda_events);
    is_configured = 0;
  }
}

#endif
