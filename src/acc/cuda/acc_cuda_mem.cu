/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015  CP2K developers group                         *
 *****************************************************************************/

#include <cuda_runtime.h>
#include <stdio.h>
#include <math.h>
#include "acc_cuda_error.h"
#include "../include/acc.h"

static const int verbose_print = 0;


/****************************************************************************/
extern "C" int acc_dev_mem_allocate(void **dev_mem, size_t n){
  cudaError_t cErr;

  cErr = cudaMalloc ((void **) dev_mem, (size_t) n);
  if (cuda_error_check (cErr))
    return -1;
  if (cuda_error_check (cudaGetLastError ()))
    return -1;
  if (dev_mem == NULL)
    return -2;
  if (verbose_print)
    printf ("Device allocation address %p, size %ld\n", *dev_mem, (long) n);

  return 0;
}


/****************************************************************************/
extern "C" int acc_dev_mem_deallocate(void *dev_mem){
  cudaError_t cErr;

  if (verbose_print)
    printf ("Device deallocation address %p\n", dev_mem);
  cErr = cudaFree ((void *) dev_mem);
  if (cuda_error_check (cErr))
    return -1;
  if (cuda_error_check (cudaGetLastError ()))
    return -1;

  return 0;
}


/****************************************************************************/
extern "C" int acc_host_mem_allocate(void **host_mem, size_t n, void *stream){
  cudaError_t cErr;
  unsigned int flag;

  flag = cudaHostAllocDefault;
  cErr = cudaHostAlloc ((void **) host_mem, (size_t) n, flag);
  if (cuda_error_check (cErr))
    return -1;
  if (cuda_error_check (cudaGetLastError ()))
    return -1;
  if (host_mem == NULL)
    return -2;
  if (verbose_print)
    printf ("Allocating %d bytes of host pinned memory at %p\n",n,  *host_mem);

  return 0;
}


/****************************************************************************/
extern "C" int acc_host_mem_deallocate(void *host_mem, void *stream){
  cudaError_t cErr;

  if (verbose_print)
    printf ("Host pinned deallocation address %p\n", host_mem);
  cErr = cudaFreeHost ((void *) host_mem);
  if (cuda_error_check (cErr))
    return -1;
  if (cuda_error_check (cudaGetLastError ()))
    return -1;

  return 0;
}


/****************************************************************************/
extern "C" int acc_memcpy_h2d(const void *host_mem, void *dev_mem, size_t count, void* stream){
  cudaError_t cErr;
  cudaStream_t* custream = (cudaStream_t*) stream;
  if (verbose_print)
      printf ("Copyint %d bytes from host address %p to device address %p \n",count, host_mem, dev_mem);

  cErr = cudaMemcpyAsync (dev_mem, host_mem, count, cudaMemcpyHostToDevice, *custream);

  if (cuda_error_check (cErr))
    return -1;
  if (cuda_error_check (cudaGetLastError ()))
    return -1;

  return 0;
}


/****************************************************************************/
extern "C" int acc_memcpy_d2h(const void *dev_mem, void *host_mem, size_t count, void* stream){
  cudaError_t cErr;
  cudaStream_t* custream = (cudaStream_t*) stream;
  if (verbose_print)
      printf ("Copying %d bytes from device address %p to host address %p\n", count, dev_mem, host_mem);

  cErr = cudaMemcpyAsync (host_mem, dev_mem, count, cudaMemcpyDeviceToHost, *custream);

  if (cuda_error_check (cErr))
    return -1;
  if (cuda_error_check (cudaGetLastError ()))
    return -1;
  if (verbose_print)
    printf ("d2h %f\n", *((double *) host_mem));

  return 0;
}


/****************************************************************************/
extern "C" int acc_memcpy_d2d(const void *devmem_src, void *devmem_dst, size_t count, void* stream){
  cudaError_t cErr;
  cudaStream_t* custream = (cudaStream_t*) stream;
  if (verbose_print)
      printf ("Coping %d bytes from device address %p to device address %p \n", count, devmem_src, devmem_dst);


  if(stream == NULL){
      cErr = cudaMemcpy (devmem_dst, devmem_src, count, cudaMemcpyDeviceToDevice);
  }else{
      cErr = cudaMemcpyAsync (devmem_dst, devmem_src, count, cudaMemcpyDeviceToDevice, *custream);
  }

  if (cuda_error_check (cErr))
    return -1;
  if (cuda_error_check (cudaGetLastError ()))
    return -1;

  return 0;
}


/****************************************************************************/
extern "C" int acc_memset_zero(void *dev_mem, size_t offset, size_t length, void* stream){
  cudaError_t cErr;
  cudaStream_t* custream = (cudaStream_t*) stream;
  if(stream == NULL){
      cErr = cudaMemset ((void *) (((char *) dev_mem) + offset), (int) 0, length);
  }else{
      cErr = cudaMemsetAsync ((void *) (((char *) dev_mem) + offset), (int) 0, length, *custream);
  }

  if (verbose_print)
    printf ("Zero at device address %p, offset %d, len %d\n",
     dev_mem, (int) offset, (int) length);
  if (cuda_error_check (cErr))
    return -1;
  if (cuda_error_check (cudaGetLastError ()))
    return -1;

  return 0;
}


/****************************************************************************/
extern "C" int acc_dev_mem_info(size_t* free, size_t* avail){
  cudaError_t cErr;
  cErr = cudaMemGetInfo (free, avail);
  if (cuda_error_check (cErr))
    return 1;
  if (cuda_error_check (cudaGetLastError ()))
    return 1;
  return 0;
}

//EOF
