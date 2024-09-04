/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2024 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#ifndef OFFLOAD_RUNTIME_H
#define OFFLOAD_RUNTIME_H

#if defined(__OFFLOAD_OPENCL) && !defined(__DBCSR_ACC)
#undef __OFFLOAD_OPENCL
#endif
/* TODO: implement support or missing features */
#if defined(__OFFLOAD_OPENCL)
#if !defined(__NO_OFFLOAD_GRID)
#define __NO_OFFLOAD_GRID
#endif
#if !defined(__NO_OFFLOAD_PW)
#define __NO_OFFLOAD_PW
#endif
#endif

#if defined(__OFFLOAD_CUDA) || defined(__OFFLOAD_HIP) ||                       \
    defined(__OFFLOAD_OPENCL)
#define __OFFLOAD

#include <stdio.h>
#include <stdlib.h>

#if defined(__OFFLOAD_CUDA)
#include <cuda_runtime.h>
#elif defined(__OFFLOAD_HIP)
#include <hip/hip_runtime.h>
#include <hip/hip_version.h>
#elif defined(__OFFLOAD_OPENCL)
/* No relative path aka double-quote to avoid PACKAGE deps (-Iexts/dbcsr). */
#include <src/acc/opencl/acc_opencl.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__OFFLOAD_CUDA)
typedef cudaStream_t offloadStream_t;
typedef cudaEvent_t offloadEvent_t;
typedef cudaError_t offloadError_t;
#elif defined(__OFFLOAD_HIP)
typedef hipStream_t offloadStream_t;
typedef hipEvent_t offloadEvent_t;
typedef hipError_t offloadError_t;
#elif defined(__OFFLOAD_OPENCL)
typedef void *offloadStream_t;
typedef void *offloadEvent_t;
typedef int offloadError_t;
#endif

#if defined(__OFFLOAD_CUDA)
#define offloadSuccess cudaSuccess
#elif defined(__OFFLOAD_HIP)
#define offloadSuccess hipSuccess
#elif defined(__OFFLOAD_OPENCL)
#define offloadSuccess EXIT_SUCCESS
#endif

/*******************************************************************************
 * \brief Check given Cuda status and upon failure abort with a nice message.
 * \author Ole Schuett
 ******************************************************************************/
#if !defined(OFFLOAD_CHECK)
#define OFFLOAD_CHECK(CMD)                                                     \
  do {                                                                         \
    const offloadError_t error = (CMD);                                        \
    if (error != offloadSuccess) {                                             \
      const char *const name = offloadGetErrorName(error);                     \
      if (NULL != name && '\0' != *name) {                                     \
        fprintf(stderr, "ERROR: \"%s\" at %s:%d\n", name, __FILE__, __LINE__); \
      } else {                                                                 \
        fprintf(stderr, "ERROR %i: %s:%d\n", (int)error, __FILE__, __LINE__);  \
      }                                                                        \
      abort();                                                                 \
    }                                                                          \
  } while (0)
#endif

/*******************************************************************************
 * \brief Wrapper around cudaGetErrorName.
 ******************************************************************************/
static inline const char *offloadGetErrorName(offloadError_t error) {
#if defined(__OFFLOAD_CUDA)
  return cudaGetErrorName(error);
#elif defined(__OFFLOAD_HIP)
  return hipGetErrorName(error);
#elif defined(__OFFLOAD_OPENCL)
  (void)error; /* mark used */
  return "";   /* TODO */
#endif
}

/*******************************************************************************
 * \brief Wrapper around cudaGetLastError.
 ******************************************************************************/
static inline offloadError_t offloadGetLastError(void) {
#if defined(__OFFLOAD_CUDA)
  return cudaGetLastError();
#elif defined(__OFFLOAD_HIP)
  return hipGetLastError();
#elif defined(__OFFLOAD_OPENCL)
  return offloadSuccess; /* TODO */
#endif
}

/*******************************************************************************
 * \brief Wrapper around cudaMemsetAsync.
 ******************************************************************************/
static inline void offloadMemsetAsync(void *const ptr, const int val,
                                      const size_t size,
                                      offloadStream_t stream) {
#if defined(__OFFLOAD_CUDA)
  OFFLOAD_CHECK(cudaMemsetAsync(ptr, val, size, stream));
#elif defined(__OFFLOAD_HIP)
  OFFLOAD_CHECK(hipMemsetAsync(ptr, val, size, stream));
#elif defined(__OFFLOAD_OPENCL)
  OFFLOAD_CHECK(
      c_dbcsr_acc_opencl_memset(ptr, val, 0 /*offset*/, size, stream));
#endif
}

/*******************************************************************************
 * \brief Wrapper around cudaMemset.
 ******************************************************************************/
static inline void offloadMemset(void *ptr, const int val, size_t size) {
#if defined(__OFFLOAD_CUDA)
  OFFLOAD_CHECK(cudaMemset(ptr, val, size));
#elif defined(__OFFLOAD_HIP)
  OFFLOAD_CHECK(hipMemset(ptr, val, size));
#elif defined(__OFFLOAD_OPENCL)
  offloadMemsetAsync(ptr, val, size, NULL /*stream*/);
#endif
}

/*******************************************************************************
 * \brief Wrapper around cudaMemcpyAsync(...,cudaMemcpyHostToDevice,...).
 ******************************************************************************/
static inline void offloadMemcpyAsyncHtoD(void *const ptr1, const void *ptr2,
                                          const size_t size,
                                          offloadStream_t stream) {
#if defined(__OFFLOAD_CUDA)
  OFFLOAD_CHECK(
      cudaMemcpyAsync(ptr1, ptr2, size, cudaMemcpyHostToDevice, stream));
#elif defined(__OFFLOAD_HIP)
#if defined(__OFFLOAD_UNIFIED_MEMORY)
  if (ptr1 == ptr2) {
    return;
  }
#endif
  OFFLOAD_CHECK(
      hipMemcpyAsync(ptr1, ptr2, size, hipMemcpyHostToDevice, stream));
#elif defined(__OFFLOAD_OPENCL)
  OFFLOAD_CHECK(c_dbcsr_acc_memcpy_h2d(ptr2, ptr1, size, stream));
#endif
}

/*******************************************************************************
 * \brief Wrapper around cudaMemcpyAsync(...,cudaMemcpyDeviceToHost,...).
 ******************************************************************************/
static inline void offloadMemcpyAsyncDtoH(void *const ptr1, const void *ptr2,
                                          const size_t size,
                                          const offloadStream_t stream) {
#if defined(__OFFLOAD_CUDA)
  OFFLOAD_CHECK(
      cudaMemcpyAsync(ptr1, ptr2, size, cudaMemcpyDeviceToHost, stream));
#elif defined(__OFFLOAD_HIP)
#if defined(__OFFLOAD_UNIFIED_MEMORY)
  if (ptr1 == ptr2) {
    return;
  }
#endif
  OFFLOAD_CHECK(
      hipMemcpyAsync(ptr1, ptr2, size, hipMemcpyDeviceToHost, stream));
#elif defined(__OFFLOAD_OPENCL)
  OFFLOAD_CHECK(c_dbcsr_acc_memcpy_d2h(ptr2, ptr1, size, stream));
#endif
}

/*******************************************************************************
 * \brief Wrapper around cudaMemcpyAsync(...,cudaMemcpyDeviceToDevice).
 ******************************************************************************/
static inline void offloadMemcpyAsyncDtoD(void *ptr1, const void *ptr2,
                                          const size_t size,
                                          const offloadStream_t stream) {
#if defined(__OFFLOAD_CUDA)
  OFFLOAD_CHECK(
      cudaMemcpyAsync(ptr1, ptr2, size, cudaMemcpyDeviceToDevice, stream));
#elif defined(__OFFLOAD_HIP)
  OFFLOAD_CHECK(
      hipMemcpyAsync(ptr1, ptr2, size, hipMemcpyDeviceToDevice, stream));
#elif defined(__OFFLOAD_OPENCL)
  OFFLOAD_CHECK(c_dbcsr_acc_memcpy_d2d(ptr2, ptr1, size, stream));
#endif
}

/*******************************************************************************
 * \brief Wrapper around cudaMemcpy(...,cudaMemcpyHostToDevice).
 ******************************************************************************/
static inline void offloadMemcpyHtoD(void *ptr_device, const void *ptr_host,
                                     const size_t size) {
#if defined(__OFFLOAD_CUDA)
  OFFLOAD_CHECK(cudaMemcpy(ptr_device, ptr_host, size, cudaMemcpyHostToDevice));
#elif defined(__OFFLOAD_HIP)
#if defined(__OFFLOAD_UNIFIED_MEMORY)
  if (ptr_device == ptr_host) {
    return;
  }
#endif
  OFFLOAD_CHECK(hipMemcpy(ptr_device, ptr_host, size, hipMemcpyHostToDevice));
#elif defined(__OFFLOAD_OPENCL)
  offloadMemcpyAsyncHtoD(ptr_device, ptr_host, size, NULL /*stream*/);
#endif
}

/*******************************************************************************
 * \brief Wrapper around cudaMemcpy(...,cudaMemcpyDeviceToHost).
 ******************************************************************************/
static inline void offloadMemcpyDtoH(void *ptr_device, const void *ptr_host,
                                     const size_t size) {
#if defined(__OFFLOAD_CUDA)
  OFFLOAD_CHECK(cudaMemcpy(ptr_device, ptr_host, size, cudaMemcpyDeviceToHost));
#elif defined(__OFFLOAD_HIP)
#if defined(__OFFLOAD_UNIFIED_MEMORY)
  if (ptr_device == ptr_host) {
    return;
  }
#endif
  OFFLOAD_CHECK(hipMemcpy(ptr_device, ptr_host, size, hipMemcpyDeviceToHost));
#elif defined(__OFFLOAD_OPENCL)
  offloadMemcpyAsyncDtoH(ptr_device, ptr_host, size, NULL /*stream*/);
#endif
}

/*******************************************************************************
 * \brief Wrapper around cudaMemcpyToSymbol.
 ******************************************************************************/
static inline void offloadMemcpyToSymbol(const void *symbol, const void *src,
                                         const size_t count) {
#if defined(__OFFLOAD_CUDA)
  OFFLOAD_CHECK(
      cudaMemcpyToSymbol(symbol, src, count, 0, cudaMemcpyHostToDevice));
#elif defined(__OFFLOAD_HIP)
  OFFLOAD_CHECK(
      hipMemcpyToSymbol(symbol, src, count, 0, hipMemcpyHostToDevice));
#elif defined(__OFFLOAD_OPENCL)
  assert(NULL == symbol || NULL == src || 0 == count); /* TODO */
#endif
}

/*******************************************************************************
 * \brief Wrapper around cudaEventCreate.
 ******************************************************************************/
static inline void offloadEventCreate(offloadEvent_t *event) {
#if defined(__OFFLOAD_CUDA)
  OFFLOAD_CHECK(cudaEventCreate(event));
#elif defined(__OFFLOAD_HIP)
  OFFLOAD_CHECK(hipEventCreate(event));
#elif defined(__OFFLOAD_OPENCL)
  OFFLOAD_CHECK(c_dbcsr_acc_event_create(event));
#endif
}

/*******************************************************************************
 * \brief Wrapper around cudaEventDestroy.
 ******************************************************************************/
static inline void offloadEventDestroy(offloadEvent_t event) {
#if defined(__OFFLOAD_CUDA)
  OFFLOAD_CHECK(cudaEventDestroy(event));
#elif defined(__OFFLOAD_HIP)
  OFFLOAD_CHECK(hipEventDestroy(event));
#elif defined(__OFFLOAD_OPENCL)
  OFFLOAD_CHECK(c_dbcsr_acc_event_destroy(event));
#endif
}

/*******************************************************************************
 * \brief Wrapper around cudaStreamCreate.
 ******************************************************************************/
static inline void offloadStreamCreate(offloadStream_t *stream) {
#if defined(__OFFLOAD_CUDA)
  OFFLOAD_CHECK(cudaStreamCreate(stream));
#elif defined(__OFFLOAD_HIP)
  OFFLOAD_CHECK(hipStreamCreate(stream));
#elif defined(__OFFLOAD_OPENCL)
  int least = 0;
  OFFLOAD_CHECK(c_dbcsr_acc_stream_priority_range(&least, NULL /*greatest*/));
  OFFLOAD_CHECK(c_dbcsr_acc_stream_create(stream, "Offload Stream", least));
#endif
}

/*******************************************************************************
 * \brief Wrapper around cudaStreamDestroy.
 ******************************************************************************/
static inline void offloadStreamDestroy(offloadStream_t stream) {
#if defined(__OFFLOAD_CUDA)
  OFFLOAD_CHECK(cudaStreamDestroy(stream));
#elif defined(__OFFLOAD_HIP)
  OFFLOAD_CHECK(hipStreamDestroy(stream));
#elif defined(__OFFLOAD_OPENCL)
  OFFLOAD_CHECK(c_dbcsr_acc_stream_destroy(stream));
#endif
}

/*******************************************************************************
 * \brief Wrapper around cudaEventSynchronize.
 ******************************************************************************/
static inline void offloadEventSynchronize(offloadEvent_t event) {
#if defined(__OFFLOAD_CUDA)
  OFFLOAD_CHECK(cudaEventSynchronize(event));
#elif defined(__OFFLOAD_HIP)
  OFFLOAD_CHECK(hipEventSynchronize(event));
#elif defined(__OFFLOAD_OPENCL)
  OFFLOAD_CHECK(c_dbcsr_acc_event_synchronize(event));
#endif
}

/*******************************************************************************
 * \brief Wrapper around cudaStreamSynchronize.
 ******************************************************************************/
static inline void offloadStreamSynchronize(offloadStream_t stream) {
#if defined(__OFFLOAD_CUDA)
  OFFLOAD_CHECK(cudaStreamSynchronize(stream));
#elif defined(__OFFLOAD_HIP)
  OFFLOAD_CHECK(hipStreamSynchronize(stream));
#elif defined(__OFFLOAD_OPENCL)
  OFFLOAD_CHECK(c_dbcsr_acc_stream_sync(stream));
#endif
}

/*******************************************************************************
 * \brief Wrapper around cudaEventRecord.
 ******************************************************************************/
static inline void offloadEventRecord(offloadEvent_t event,
                                      offloadStream_t stream) {
#if defined(__OFFLOAD_CUDA)
  OFFLOAD_CHECK(cudaEventRecord(event, stream));
#elif defined(__OFFLOAD_HIP)
  OFFLOAD_CHECK(hipEventRecord(event, stream));
#elif defined(__OFFLOAD_OPENCL)
  OFFLOAD_CHECK(c_dbcsr_acc_event_record(event, stream));
#endif
}

/*******************************************************************************
 * \brief Wrapper around cudaMallocHost.
 ******************************************************************************/
static inline void offloadMallocHost(void **ptr, size_t size) {
#if defined(__OFFLOAD_CUDA)
  OFFLOAD_CHECK(cudaMallocHost(ptr, size));
#elif defined(__OFFLOAD_HIP)
#if !defined(__OFFLOAD_UNIFIED_MEMORY)
  OFFLOAD_CHECK(hipHostMalloc(ptr, size, hipHostMallocDefault)); // inconsistent
#else
  assert(NULL != ptr);
  *ptr = NULL;
#endif
#elif defined(__OFFLOAD_OPENCL)
  OFFLOAD_CHECK(c_dbcsr_acc_host_mem_allocate(ptr, size, NULL /*stream*/));
#else
  assert(NULL != ptr);
  *ptr = malloc(size);
#endif
}

/*******************************************************************************
 * \brief Wrapper around cudaMalloc.
 ******************************************************************************/
static inline void offloadMalloc(void **ptr, size_t size) {
#if defined(__OFFLOAD_CUDA)
  OFFLOAD_CHECK(cudaMalloc(ptr, size));
#elif defined(__OFFLOAD_HIP)
  OFFLOAD_CHECK(hipMalloc(ptr, size));
#elif defined(__OFFLOAD_OPENCL)
  OFFLOAD_CHECK(c_dbcsr_acc_dev_mem_allocate(ptr, size));
#else
  assert(NULL != ptr);
  *ptr = NULL;
#endif
}

/*******************************************************************************
 * \brief Wrapper around cudaFree.
 ******************************************************************************/
static inline void offloadFree(void *ptr) {
#if defined(__OFFLOAD_CUDA)
  OFFLOAD_CHECK(cudaFree(ptr));
#elif defined(__OFFLOAD_HIP)
  OFFLOAD_CHECK(hipFree(ptr));
#elif defined(__OFFLOAD_OPENCL)
  OFFLOAD_CHECK(c_dbcsr_acc_dev_mem_deallocate(ptr));
#else
  assert(NULL == ptr);
#endif
}

/*******************************************************************************
 * \brief Wrapper around cudaFreeHost.
 ******************************************************************************/
static inline void offloadFreeHost(void *ptr) {
#if defined(__OFFLOAD_CUDA)
  OFFLOAD_CHECK(cudaFreeHost(ptr));
#elif defined(__OFFLOAD_HIP)
#if !defined(__OFFLOAD_UNIFIED_MEMORY)
  OFFLOAD_CHECK(hipHostFree(ptr)); // inconsistent
#endif
#elif defined(__OFFLOAD_OPENCL)
  OFFLOAD_CHECK(c_dbcsr_acc_host_mem_deallocate(ptr, NULL /*stream*/));
#else
  free(ptr);
#endif
}

/*******************************************************************************
 * \brief Wrapper around cudaStreamWaitEvent.
 ******************************************************************************/
static inline void offloadStreamWaitEvent(offloadStream_t stream,
                                          offloadEvent_t event, const int val) {
#if defined(__OFFLOAD_CUDA)
  OFFLOAD_CHECK(cudaStreamWaitEvent(stream, event, val));
#elif defined(__OFFLOAD_HIP)
  OFFLOAD_CHECK(hipStreamWaitEvent(stream, event, val));
#elif defined(__OFFLOAD_OPENCL)
  assert(0 == val); /* TODO */
  OFFLOAD_CHECK(c_dbcsr_acc_stream_wait_event(stream, event));
#endif
}

/*******************************************************************************
 * \brief Wrapper around cudaDeviceSynchronize.
 ******************************************************************************/
static inline void offloadDeviceSynchronize(void) {
#if defined(__OFFLOAD_CUDA)
  OFFLOAD_CHECK(cudaDeviceSynchronize());
#elif defined(__OFFLOAD_HIP)
  OFFLOAD_CHECK(hipDeviceSynchronize());
#elif defined(__OFFLOAD_OPENCL)
  OFFLOAD_CHECK(c_dbcsr_acc_device_synchronize());
#endif
}

/*******************************************************************************
 * \brief Wrapper around cudaDeviceSetLimit(cudaLimitMallocHeapSize,...).
 ******************************************************************************/
static inline void offloadEnsureMallocHeapSize(const size_t required_size) {
#if defined(__OFFLOAD_CUDA)
  size_t current_size;
  OFFLOAD_CHECK(cudaDeviceGetLimit(&current_size, cudaLimitMallocHeapSize));
  if (current_size < required_size) {
    OFFLOAD_CHECK(cudaDeviceSetLimit(cudaLimitMallocHeapSize, required_size));
  }
#elif defined(__OFFLOAD_HIP) && (HIP_VERSION >= 50300000)
  size_t current_size;
  OFFLOAD_CHECK(hipDeviceGetLimit(&current_size, hipLimitMallocHeapSize));
  if (current_size < required_size) {
    OFFLOAD_CHECK(hipDeviceSetLimit(hipLimitMallocHeapSize, required_size));
  }
#elif defined(__OFFLOAD_OPENCL)
  assert(0 == required_size); /* TODO */
#else
  (void)required_size; /* mark used */
#endif
}

#ifdef __cplusplus
}
#endif

#endif // defined(__OFFLOAD_CUDA) || defined(__OFFLOAD_HIP) ||
       // defined(__OFFLOAD_OPENCL)
#endif
