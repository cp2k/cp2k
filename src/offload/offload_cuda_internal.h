/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#ifndef OFFLOAD_CUDA_INTERNAL_H
#define OFFLOAD_CUDA_INTERNAL_H
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif
typedef cudaStream_t offloadStream_t;
typedef cudaEvent_t offloadEvent_t;

/*******************************************************************************
 * \brief Check given Cuda status and upon failure abort with a nice message.
 * \author Ole Schuett
 ******************************************************************************/
#define OFFLOAD_CHECK(status)                                                  \
  if (status != cudaSuccess) {                                                 \
    fprintf(stderr, "ERROR: %s %s %d\n", cudaGetErrorString(status), __FILE__, \
            __LINE__);                                                         \
    abort();                                                                   \
  }

static inline void offloadMemsetAsync(void *const ptr__, const int val__,
                                      const size_t size__,
                                      offloadStream_t stream__) {
  OFFLOAD_CHECK(cudaMemsetAsync(ptr__, val__, size__, stream__));
}

static inline void offloadMemcpyAsyncHtoD(void *const ptr1__,
                                          const void *ptr2__,
                                          const size_t size__,
                                          offloadStream_t stream__) {
  OFFLOAD_CHECK(cudaMemcpyAsync(ptr1__, ptr2__, size__, cudaMemcpyHostToDevice,
                                stream__));
}

static inline void offloadMemcpyAsyncDtoH(void *const ptr1__,
                                          const void *ptr2__,
                                          const size_t size__,
                                          const offloadStream_t stream__) {
  OFFLOAD_CHECK(cudaMemcpyAsync(ptr1__, ptr2__, size__, cudaMemcpyDeviceToHost,
                                stream__));
}

static inline void offloadEventCreate(offloadEvent_t *event__) {
  OFFLOAD_CHECK(cudaEventCreate(event__));
}

static inline void offloadEventDestroy(offloadEvent_t event__) {
  OFFLOAD_CHECK(cudaEventDestroy(event__));
}

static inline void offloadStreamCreate(offloadStream_t *stream__) {
  OFFLOAD_CHECK(cudaStreamCreate(stream__));
}

static inline void offloadStreamDestroy(offloadStream_t stream__) {
  OFFLOAD_CHECK(cudaStreamDestroy(stream__));
}

static inline void offloadEventSynchronize(offloadEvent_t event__) {
  OFFLOAD_CHECK(cudaEventSynchronize(event__));
}

static inline void offloadStreamSynchronize(offloadStream_t stream__) {
  OFFLOAD_CHECK(cudaStreamSynchronize(stream__));
}

static inline void offloadEventRecord(offloadEvent_t event__,
                                      offloadStream_t stream__) {
  OFFLOAD_CHECK(cudaEventRecord(event__, stream__));
}

static inline void offloadMallocHost(void **ptr__, size_t size__) {
  OFFLOAD_CHECK(cudaMallocHost(ptr__, size__));
}

static inline void offloadMalloc(void **ptr__, size_t size__) {
  OFFLOAD_CHECK(cudaMalloc(ptr__, size__));
}

static inline void offloadFree(void *ptr__) { OFFLOAD_CHECK(cudaFree(ptr__)); }

static inline void offloadFreeHost(void *ptr__) {
  OFFLOAD_CHECK(cudaFreeHost(ptr__));
}

static inline void offloadStreamWaitEvent(offloadStream_t stream__,
                                          offloadEvent_t event__,
                                          const int val__) {
  cudaStreamWaitEvent(stream__, event__, val__);
}

#ifdef __cplusplus
}
#endif

#endif
