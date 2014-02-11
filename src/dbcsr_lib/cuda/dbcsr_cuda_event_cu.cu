/******************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations
 *  Copyright (C) 2000 - 2013  Urban Borstnik and the CP2K developers group
 *****************************************************************************/

#include <cuda_runtime.h>
#include <stdio.h>

#include "dbcsr_cuda.h"
#include <math.h>


  static const int verbose_print = 0;


extern "C" int cuda_event_create(cudaEvent_t** event_p){
  *event_p = (cudaEvent_t*) malloc(sizeof(cudaEvent_t));
  cudaError_t cErr = cudaEventCreate(*event_p);
  if(verbose_print) printf("cuda_event_created:  %p -> %d\n", *event_p, **event_p);
  if (cuda_error_check(cErr)) return 1;
  if (cuda_error_check(cudaGetLastError())) return 1;
  return 0;
}

extern "C" int cuda_event_destroy(cudaEvent_t* event){
    if(verbose_print) printf("cuda_event_destroy called\n");
    cudaError_t cErr = cudaEventDestroy(*event);
    free(event);
    if (cuda_error_check (cErr)) return 1;
    if (cuda_error_check(cudaGetLastError ()))return 1;
    return 0;
}

extern "C" int cuda_event_record(cudaEvent_t* event, cudaStream_t* stream){
    if(verbose_print) printf("cuda_event_record: %p -> %d,  %p -> %d\n", event, * event,  stream, *stream);
    cudaError_t cErr = cudaEventRecord (*event, *stream);
    if (cuda_error_check (cErr)) return 1;
    //if (cuda_error_check(cudaGetLastError ()))return 1;
    return 0;
}

extern "C" int cuda_event_query(cudaEvent_t* event){
    if(verbose_print) printf("cuda_event_query called\n");
    cudaError_t cErr = cudaEventQuery(*event);
    //if(cuda_error_check(cudaGetLastError ())) return -1;
    if(cErr==cudaSuccess) return 0;
    if(cErr==cudaErrorNotReady) return 1;
    return -2;
}

extern "C" int cuda_stream_wait_event(cudaStream_t* stream, cudaEvent_t* event){
    if(verbose_print) printf("cuda_stream_wait_event called\n");
    // flags: Parameters for the operation (must be 0)
    cudaError_t cErr = cudaStreamWaitEvent(*stream, *event, 0);
    if (cuda_error_check (cErr)) return 1;
    //if (cuda_error_check(cudaGetLastError ()))return 1;
    return 0;
}

extern "C" int cuda_event_synchronize(cudaEvent_t* event){
    if(verbose_print) printf("cuda_event_synchronize called\n");
    cudaError_t cErr = cudaEventSynchronize(*event);
    if (cuda_error_check (cErr)) return 1;
    if (cuda_error_check(cudaGetLastError ()))return 1;
    return 0;
}

//
// extern "C" int
// dc_dev_mem_realloc (void **dev_mem, size_t n, size_t old_n,
// 		    int *memory_crunch)
// {
//   cudaError_t cErr;
//   void *new_dev_mem;
//   size_t count;
// 
//   *memory_crunch = 0;
//   cErr = cudaMalloc ((void **) &new_dev_mem, (size_t) n);
//   if (cuda_error_check (cErr))
//     return 1;
//   if (cuda_error_check (cudaGetLastError ()))
//     return 1;
//   if (new_dev_mem == NULL)
//     return 2;
//   if (verbose_print)
//     printf ("Device allocation address %p, size %ld\n", new_dev_mem,
// 	    (long) n);
//   count = MIN (old_n, n);
//   if (count > 0)
//     {
//       if (verbose_print)
// 	printf ("Copy %d bytes.\n", (int) count);
//       cErr =
// 	cudaMemcpy (new_dev_mem, *dev_mem, count, cudaMemcpyDeviceToDevice);
//       if (cuda_error_check (cErr))
// 	return 1;
//       if (cuda_error_check (cudaGetLastError ()))
// 	return 1;
//     }
// 
//   cErr = cudaFree ((void *) *dev_mem);
//   if (cuda_error_check (cErr))
//     return 1;
//   if (cuda_error_check (cudaGetLastError ()))
//     return 1;
// 
//   *dev_mem = new_dev_mem;
//   return 0;
// }
// 
// extern "C" int
// dc_dev_mem_dealloc (void *dev_mem)
// {
//   cudaError_t cErr;
// 
//   if (verbose_print)
//     printf ("Device deallocation address %p\n", dev_mem);
//   cErr = cudaFree ((void *) dev_mem);
//   if (cuda_error_check (cErr))
//     return 1;
//   if (cuda_error_check (cudaGetLastError ()))
//     return 1;
// 
//   return 0;
// }
// 
// extern "C" int
// dc_host_mem_alloc (void **host_mem, size_t n, int wc, int port)
// {
//   cudaError_t cErr;
//   unsigned int flag;
// 
//   flag = cudaHostAllocDefault;
//   if (wc)
//     flag |= cudaHostAllocWriteCombined;
//   if (port)
//     flag |= cudaHostAllocPortable;
//   cErr = cudaHostAlloc ((void **) host_mem, (size_t) n, flag);
//   if (cuda_error_check (cErr))
//     return 1;
//   if (cuda_error_check (cudaGetLastError ()))
//     return 1;
//   if (host_mem == NULL)
//     return 2;
//   if (verbose_print)
//     printf ("Host pinned allocation address %p\n", *host_mem);
// 
//   return 0;
// }
// 
// extern "C" int
// dc_host_mem_dealloc (void *host_mem)
// {
//   cudaError_t cErr;
// 
//   if (verbose_print)
//     printf ("Host pinned deallocation address %p\n", host_mem);
//   cErr = cudaFreeHost ((void *) host_mem);
//   if (cuda_error_check (cErr))
//     return 1;
//   if (cuda_error_check (cudaGetLastError ()))
//     return 1;
// 
//   return 0;
// }
// 
// 
// extern "C" int
// dc_memcpy_h2d_cu (const void *host_mem, void *dev_mem, size_t count,
// 		  int async_type, int stream_id)
// {
//   cudaError_t cErr;
// 
//   if (verbose_print)
//     {
//       printf ("Copy from host address %p\n", host_mem);
//       printf ("Copy to device address %p\n", dev_mem);
//       printf ("h2d %f\n", *((double *) host_mem));
//       printf ("Async? %d\n", async_type);
//     }
// 
//   switch (async_type)
//     {
//     case 0:
//       /* Synchronous */
//       cErr = cudaMemcpy (dev_mem, host_mem, count, cudaMemcpyHostToDevice);
//       break;
//     case 1:
//       /* Asynchronous */
//       cErr =
// 	cudaMemcpyAsync (dev_mem, host_mem, count, cudaMemcpyHostToDevice,
// 			 (cudaStream_t) dc_get_stream (stream_id));
//       break;
//     case 2:
//       cErr =
// 	cudaMemcpyAsync (dev_mem, host_mem, count, cudaMemcpyHostToDevice,
// 			 (cudaStream_t) dc_get_stream (stream_id));
//       /* Try async if sync is unsuccessful. */
//       if (cuda_error_check (cErr))
// 	{
// 	  if (verbose_print)
// 	    printf ("Async unsuccessful, trying sync.\n");
// 	  cErr =
// 	    cudaMemcpy (dev_mem, host_mem, count, cudaMemcpyHostToDevice);
// 	}
//       break;
//     }
//   if (cuda_error_check (cErr))
//     return 1;
//   if (cuda_error_check (cudaGetLastError ()))
//     return 1;
// 
//   return 0;
// }
// 
// 
// extern "C" int
// dc_memcpy_d2h_cu (const void *dev_mem, void *host_mem, size_t count,
// 		  int async_type, int stream_id)
// {
//   cudaError_t cErr;
// 
//   if (verbose_print)
//     {
//       printf ("Copy from device address %p\n", dev_mem);
//       printf ("Copy to host address %p\n", host_mem);
//       printf ("Async? %d\n", async_type);
//     }
//   switch (async_type)
//     {
//     case 0:
//       /* Synchronous */
//       cErr = cudaMemcpy (host_mem, dev_mem, count, cudaMemcpyDeviceToHost);
//       break;
//     case 1:
//       /* Asynchronous */
//       cErr =
// 	cudaMemcpyAsync (host_mem, dev_mem, count, cudaMemcpyDeviceToHost,
// 			 (cudaStream_t) dc_get_stream (stream_id));
//       break;
//     case 2:
//       cErr =
// 	cudaMemcpyAsync (host_mem, dev_mem, count, cudaMemcpyDeviceToHost,
// 			 (cudaStream_t) dc_get_stream (stream_id));
//       /* Try async if sync is unsuccessful. */
//       if (cuda_error_check (cErr))
// 	{
// 	  if (verbose_print)
// 	    printf ("Async unsuccessful, trying sync.\n");
// 	  cErr =
// 	    cudaMemcpy (host_mem, dev_mem, count, cudaMemcpyDeviceToHost);
// 	}
//       break;
//     }
//   if (cuda_error_check (cErr))
//     return 1;
//   if (cuda_error_check (cudaGetLastError ()))
//     return 1;
//   if (verbose_print)
//     printf ("d2h %f\n", *((double *) host_mem));
// 
//   return 0;
// }
// 
// 
// extern "C" int
// dc_memzero_cu (void *dev_mem, size_t offset, size_t length)
// {
//   cudaError_t cErr;
// 
//   cErr = cudaMemset ((void *) (((char *) dev_mem) + offset), (int) 0, length);
//   if (verbose_print)
//     printf ("Zero at device address %p, offset %d, len %d\n",
// 	    dev_mem, (int) offset, (int) length);
//   if (cuda_error_check (cErr))
//     return 1;
//   if (cuda_error_check (cudaGetLastError ()))
//     return 1;
// 
//   /*  struct cudaDeviceProp devProperties;
//      int myDevice, nt, nb, ws, maxt;
// 
//      cErr = cudaGetDevice(&myDevice);
//      if (cuda_error_check (cErr)) return 1;
// 
//      cErr = cudaGetDeviceProperties(&devProperties, myDevice);
//      if (cuda_error_check (cErr)) return 1;
// 
//      ws = devProperties.warpSize;
//      maxt = devProperties.maxThreadsPerBlock;
//      printf("count %d, ws %d, maxt %d", (int) count, ws, maxt);
// 
//      nt = (int) sqrt(count);
//      nt = ((int) (nt + ws-1)/ws) * ws;
//      nt = MAX(MIN(nt, maxt), ws);
// 
//      printf("nt", nt);
// 
//      nb = (count+nt-1) / nt;
//      printf("nb", nb);
// 
// 
//      zeroMem <<< nb, nt >>> ((char *) dev_mem, (int) count);
//      if (cuda_error_check (cudaGetLastError())) return 1; */
//   return 0;
// }
// 
// extern "C" int
// dc_dev_mem_info_cu (size_t * free, size_t * avail)
// {
//   cudaError_t cErr;
//   cErr = cudaMemGetInfo (free, avail);
//   if (cuda_error_check (cErr))
//     return 1;
//   if (cuda_error_check (cudaGetLastError ()))
//     return 1;
//   return 0;
// }
