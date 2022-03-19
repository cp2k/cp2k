/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#ifndef OFFLOAD_RUNTIME_H
#define OFFLOAD_RUNTIME_H

#if defined(__GRID_CUDA) || defined(__DBM_CUDA) || defined(__PW_CUDA)
#define __OFFLOAD_CUDA
#elif defined(__GRID_HIP) || defined(__DBM_HIP) || defined(__PW_HIP)
#define __OFFLOAD_HIP
#endif

#if (defined(__OFFLOAD_CUDA) || defined(__OFFLOAD_HIP))

#include <stdio.h>
#include <stdlib.h>

#if defined(__OFFLOAD_CUDA)
#include <cuda_runtime.h>
#elif defined(__OFFLOAD_HIP)
#include <hip/hip_runtime.h>
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
#endif

#if defined(__OFFLOAD_CUDA)
#define offloadSuccess cudaSuccess
#elif defined(__OFFLOAD_HIP)
#define offloadSuccess hipSuccess
#endif

/*******************************************************************************
 * \brief Check given Cuda status and upon failure abort with a nice message.
 * \author Ole Schuett
 ******************************************************************************/
#define OFFLOAD_CHECK(status)                                                  \
  if (status != offloadSuccess) {                                              \
    fprintf(stderr, "ERROR: %s %s %d\n", offloadGetErrorName(status),          \
            __FILE__, __LINE__);                                               \
    abort();                                                                   \
  }

/*******************************************************************************
 * \brief Wrapper around cudaGetErrorName.
 ******************************************************************************/
static inline const char *offloadGetErrorName(offloadError_t error) {
#if defined(__OFFLOAD_CUDA)
  return cudaGetErrorName(error);
#elif defined(__OFFLOAD_HIP)
  return hipGetErrorName(error);
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
  OFFLOAD_CHECK(
      hipMemcpyAsync(ptr1, ptr2, size, hipMemcpyHostToDevice, stream));
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
  OFFLOAD_CHECK(
      hipMemcpyAsync(ptr1, ptr2, size, hipMemcpyDeviceToHost, stream));
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
  OFFLOAD_CHECK(hipMemcpy(ptr_device, ptr_host, size, hipMemcpyHostToDevice));
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
  OFFLOAD_CHECK(hipMemcpy(ptr_device, ptr_host, size, hipMemcpyDeviceToHost));
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
#endif
}

/*******************************************************************************
 * \brief Wrapper around cudaMallocHost.
 ******************************************************************************/
static inline void offloadMallocHost(void **ptr, size_t size) {
#if defined(__OFFLOAD_CUDA)
  OFFLOAD_CHECK(cudaMallocHost(ptr, size));
#elif defined(__OFFLOAD_HIP)
  OFFLOAD_CHECK(hipHostMalloc(ptr, size, hipHostMallocDefault)); // inconsistent
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
#endif
}

/*******************************************************************************
 * \brief Wrapper around cudaFreeHost.
 ******************************************************************************/
static inline void offloadFreeHost(void *ptr) {
#if defined(__OFFLOAD_CUDA)
  OFFLOAD_CHECK(cudaFreeHost(ptr));
#elif defined(__OFFLOAD_HIP)
  OFFLOAD_CHECK(hipHostFree(ptr)); // inconsistent
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
#endif
}

#ifdef __cplusplus
}
#endif

#endif // #if (defined(__OFFLOAD_CUDA) || defined(__OFFLOAD_HIP))

#endif
