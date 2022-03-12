/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "offload_buffer.h"
#include "offload_library.h"
#include "offload_runtime.h"
/*******************************************************************************
 * \brief Allocates a buffer of given length, ie. number of elements.
 * \author Ole Schuett
 ******************************************************************************/
void offload_create_buffer(const int length, offload_buffer **buffer) {

  const size_t requested_size = length * sizeof(double);

  if (*buffer != NULL) {
    if ((*buffer)->size >= requested_size) {
      return; // reuse existing buffer
    } else {
      offload_free_buffer(*buffer);
    }
  }

  (*buffer) = malloc(sizeof(offload_buffer));
  (*buffer)->size = requested_size;

#if defined(__OFFLOAD_CUDA)
  // With size 0 cudaMallocHost doesn't null the pointer and cudaFreeHost fails.
  (*buffer)->host_buffer = NULL;
  offload_set_device();
  OFFLOAD_CHECK(
      cudaMallocHost((void **)&(*buffer)->host_buffer, requested_size));
  OFFLOAD_CHECK(cudaMalloc((void **)&(*buffer)->device_buffer, requested_size));
#elif defined(__OFFLOAD_HIP)
  // With size 0 cudaMallocHost doesn't null the pointer and cudaFreeHost fails.
  (*buffer)->host_buffer = NULL;
  offload_set_device();
  OFFLOAD_CHECK(hipHostMalloc((void **)&(*buffer)->host_buffer, requested_size,
                              hipHostMallocDefault));
  OFFLOAD_CHECK(hipMalloc((void **)&(*buffer)->device_buffer, requested_size));
#else
  (*buffer)->host_buffer = malloc(requested_size);
  (*buffer)->device_buffer = NULL;
#endif
  return;
}

/*******************************************************************************
 * \brief Deallocate given buffer.
 * \author Ole Schuett
 ******************************************************************************/
void offload_free_buffer(offload_buffer *buffer) {

  if (buffer == NULL)
    return;

#ifdef __OFFLOAD_CUDA
  OFFLOAD_CHECK(cudaFreeHost(buffer->host_buffer));
  OFFLOAD_CHECK(cudaFree(buffer->device_buffer));
#elif defined(__OFFLOAD_HIP)
  OFFLOAD_CHECK(hipHostFree(buffer->host_buffer));
  OFFLOAD_CHECK(hipFree(buffer->device_buffer));
#else
  free(buffer->host_buffer);
#endif
  free(buffer);
  return;
}

/*******************************************************************************
 * \brief Returns a pointer to the host buffer.
 * \author Ole Schuett
 ******************************************************************************/
double *offload_get_buffer_host_pointer(offload_buffer *buffer) {
  assert(buffer != NULL);
  return buffer->host_buffer;
}

// EOF
