/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2021 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "offload_buffer.h"
#include "offload_library.h"

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

#ifdef __OFFLOAD_CUDA
  // With size 0 cudaMallocHost doesn't null the pointer and cudaFreeHost fails.
  (*buffer)->host_buffer = NULL;
  OFFLOAD_CHECK(cudaSetDevice(offload_get_device_id()));
  OFFLOAD_CHECK(
      cudaMallocHost((void **)&(*buffer)->host_buffer, requested_size));
  OFFLOAD_CHECK(cudaMalloc((void **)&(*buffer)->device_buffer, requested_size));
#else
  (*buffer)->host_buffer = malloc(requested_size);
  (*buffer)->device_buffer = NULL;
#endif
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
#else
  free(buffer->host_buffer);
#endif

  free(buffer);
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
