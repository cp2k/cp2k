/******************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations
 *  Copyright (C) 2000 - 2011  Urban Borstnik and the CP2K developers group
 *****************************************************************************/

#include <cuda_runtime.h>
#include <stdio.h>

#include "dbcsr_cuda.h"
#include <math.h>


static const int verbose_print = 0;


__global__ void zeroMem (char *mem, int len) {
  int offset;
  offset = blockIdx.x*blockDim.x + threadIdx.x;
  if (offset < len)
    *(mem+offset) = 0;
}


extern "C" int dc_dev_mem_alloc (void **dev_mem, size_t n) {
  cudaError_t cErr;

  cErr = cudaMalloc((void**) dev_mem, (size_t) n);
  if (cuda_error_check (cErr)) return 1;
  if (cuda_error_check (cudaGetLastError())) return 1;
  if (dev_mem == NULL) return 2;
  if (verbose_print) printf("Device allocation address %p\n", *dev_mem);

  return 0;
}


extern "C" int dc_dev_mem_realloc (void **dev_mem, size_t n, size_t old_n, int *memory_crunch) {
  cudaError_t cErr;
  void *new_dev_mem;
  size_t count;

  *memory_crunch = 0;
  cErr = cudaMalloc((void**) &new_dev_mem, (size_t) n);
  if (cuda_error_check (cErr)) return 1;
  if (cuda_error_check (cudaGetLastError())) return 1;
  if (dev_mem == NULL) return 2;

  count = MIN(old_n, n);
  if (count > 0) {
    cErr = cudaMemcpy (new_dev_mem, *dev_mem, count, cudaMemcpyDeviceToDevice);
    if (cuda_error_check (cErr)) return 1;
    if (cuda_error_check (cudaGetLastError())) return 1;
  }

  cErr = cudaFree((void *) *dev_mem);
  if (cuda_error_check (cErr)) return 1;
  if (cuda_error_check (cudaGetLastError())) return 1;

  dev_mem = &new_dev_mem;
  return 0;
}

extern "C" int dc_dev_mem_dealloc (void *dev_mem) {
  cudaError_t cErr;

  if (verbose_print)  printf("Device deallocation address %p\n", dev_mem);
  cErr = cudaFree((void *) dev_mem);
  if (cuda_error_check (cErr)) return 1;
  if (cuda_error_check (cudaGetLastError())) return 1;

  return 0;
}

extern "C" int dc_host_mem_alloc (void **host_mem, size_t n, int wc, int port) {
  cudaError_t cErr;
  unsigned int flag;

  flag = cudaHostAllocDefault;
  if (wc) flag |= cudaHostAllocWriteCombined;
  if (port) flag |= cudaHostAllocPortable;
  cErr = cudaHostAlloc((void**) host_mem, (size_t) n, flag);
  if (cuda_error_check (cErr)) return 1;
  if (cuda_error_check (cudaGetLastError())) return 1;
  if (host_mem == NULL) return 2;
  if (verbose_print) printf("Host pinned allocation address %p\n", *host_mem);

  return 0;
}

extern "C" int dc_host_mem_dealloc (void *host_mem) {
  cudaError_t cErr;

  if (verbose_print) printf("Host pinned deallocation address %p\n", host_mem);
  cErr = cudaFreeHost((void *) host_mem);
  if (cuda_error_check (cErr)) return 1;
  if (cuda_error_check (cudaGetLastError())) return 1;

  return 0;
}


extern "C" int dc_memcpy_h2d_cu (const void *host_mem, void *dev_mem, size_t count) {
  cudaError_t cErr;

  if (verbose_print) {
    printf("Copy from host address %p\n", host_mem);
    printf("Copy to device address %p\n", dev_mem);
    printf("h2d %f\n", *((double *) host_mem));
  }
  cErr = cudaMemcpy(dev_mem, host_mem, count, cudaMemcpyHostToDevice);
  if (cuda_error_check (cErr)) return 1;
  if (cuda_error_check (cudaGetLastError())) return 1;

  return 0;
}


extern "C" int dc_memcpy_d2h_cu (const void *dev_mem, void *host_mem, size_t count) {
  cudaError_t cErr;

  if (verbose_print) {
    printf("Copy from device address %p\n", dev_mem);
    printf("Copy to host address %p\n", host_mem);
  }
  cErr = cudaMemcpy(host_mem, dev_mem, count, cudaMemcpyDeviceToHost);
  if (cuda_error_check (cErr)) return 1;
  if (cuda_error_check (cudaGetLastError())) return 1;
  if (verbose_print) printf("d2h %f\n", *((double *) host_mem));

  return 0;
}


extern "C" int dc_memzero_cu (void *dev_mem, size_t offset, size_t length) {
  cudaError_t cErr;

  cErr = cudaMemset ((void *) (((char *) dev_mem)+offset), (int) 0, length);
  if (cuda_error_check (cErr)) return 1;
  if (cuda_error_check (cudaGetLastError())) return 1;

  /*  struct cudaDeviceProp devProperties;
  int myDevice, nt, nb, ws, maxt;

  cErr = cudaGetDevice(&myDevice);
  if (cuda_error_check (cErr)) return 1;

  cErr = cudaGetDeviceProperties(&devProperties, myDevice);
  if (cuda_error_check (cErr)) return 1;

  ws = devProperties.warpSize;
  maxt = devProperties.maxThreadsPerBlock;
  printf("count %d, ws %d, maxt %d", (int) count, ws, maxt);

  nt = (int) sqrt(count);
  nt = ((int) (nt + ws-1)/ws) * ws;
  nt = MAX(MIN(nt, maxt), ws);

  printf("nt", nt);

  nb = (count+nt-1) / nt;
  printf("nb", nb);

   
  zeroMem <<< nb, nt >>> ((char *) dev_mem, (int) count);
  if (cuda_error_check (cudaGetLastError())) return 1;*/
  return 0;
}
