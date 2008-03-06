#if defined ( __FFTCU ) || defined ( __CUDA )

/* This file contains memory management routines for device memory.  The goal of these routines is eliminate the need to allocate device memory more than once during a run.  Basically, CP2K allocates a large chunk of memory at the beginning of the job using cuda_device_mem_init_cu_.  The size of this chunk is specified in the input file.  The below routines cuda_device_mem_alloc_cu_ and cuda_device_mem_free_cu_ are used to allocate and free memory from this chunk.  The bookkeeping is done by the static variables below.  Finally, the chunk can be freed at the end of the job using cuda_device_mem_release_cu_ */

// Author:  Benjamin G Levine

#include <cuda_runtime.h>
#include <cufft.h>
#include <stdio.h>

static const int max_allocs=20;
static int ioffset[max_allocs];
static int isize[max_allocs];
static int n_allocs;
static int length;
static float* device_memory;

extern "C" void cuda_device_mem_init_cu_ (int *memory) {

  
  n_allocs=0;
  length=(*memory)*256;
  cudaMalloc((void**)&device_memory, sizeof(float)*(length));
}

extern void cuda_device_mem_alloc_cu_ (float **ptr, int n) {

  if ( n_allocs==0 && (n) <= length ) {
    isize[n_allocs] = n;
    ioffset[n_allocs] = 0;
    *ptr = device_memory;
    n_allocs++;
  }
  else {
    if ((n) <= (length - (ioffset[n_allocs-1] + isize[n_allocs-1] ))) { 
      isize[n_allocs] = (n);
      ioffset[n_allocs] = ioffset[n_allocs-1] + isize[n_allocs-1];
      *ptr = device_memory + ioffset[n_allocs];
      n_allocs++;
    }
    else {
      *ptr=NULL;
      printf("NOT ENOUGH GPU DEVICE MEMORY!\n");
      fflush(stdout);
    }
  }
}

extern void cuda_device_mem_alloc_cu_ (int **ptr, int n) {

  if ( n_allocs==0 && (n) <= length ) {
    isize[n_allocs] = n;
    ioffset[n_allocs] = 0;
    *ptr = (int*)device_memory;
    n_allocs++;
  }
  else {
    if ((n) <= (length - (ioffset[n_allocs-1] + isize[n_allocs-1] ))) { 
      isize[n_allocs] = (n);
      ioffset[n_allocs] = ioffset[n_allocs-1] + isize[n_allocs-1];
      *ptr = (int*)device_memory + ioffset[n_allocs];
      n_allocs++;
    }
    else {
      *ptr=NULL;
      printf("NOT ENOUGH GPU DEVICE MEMORY!\n");
      fflush(stdout);
    }
  }
}

extern void cuda_device_mem_free_cu_ (float **ptr) {
  int ioffsettmp, zfound, i;
  
  zfound=0;
  ioffsettmp = *ptr - device_memory;
  for (i=0; i < n_allocs; i++) {
    if (zfound) {
      ioffset[i-1] = ioffset[i];
      isize[i-1] = isize[i];
    }
    if (ioffsettmp==ioffset[i]) zfound=1;
  }
  n_allocs--;
  *ptr=NULL;
}

extern void cuda_device_mem_free_cu_ (int **ptr) {
  int ioffsettmp, zfound, i;
  
  zfound=0;
  ioffsettmp = (float*)(*ptr) - device_memory;
  for (i=0; i < n_allocs; i++) {
    if (zfound) {
      ioffset[i-1] = ioffset[i];
      isize[i-1] = isize[i];
    }
    if (ioffsettmp==ioffset[i]) zfound=1;
  }
  n_allocs--;
  *ptr=NULL;
}

extern "C" void cuda_device_mem_release_cu_ () {

  cudaFree(device_memory);

}

#endif
