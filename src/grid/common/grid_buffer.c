/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "grid_buffer.h"

#ifdef __GRID_CUDA
#include <cuda_runtime.h>

/*******************************************************************************
 * \brief Check given Cuda status and upon failure abort with a nice message.
 * \author Ole Schuett
 ******************************************************************************/
#define CHECK(status)                                                          \
  if (status != cudaSuccess) {                                                 \
    fprintf(stderr, "ERROR: %s %s %d\n", cudaGetErrorString(status), __FILE__, \
            __LINE__);                                                         \
    abort();                                                                   \
  }

#endif // __GRID_CUDA

/*******************************************************************************
 * \brief Allocates a buffer of given length, ie. number of elements.
 * \author Ole Schuett
 ******************************************************************************/
void grid_create_buffer(const int length, grid_buffer **buffer) {

  const size_t requested_size = length * sizeof(double);

  if (*buffer != NULL) {
    if ((*buffer)->size >= requested_size) {
      return; // reuse existing buffer
    } else {
      grid_free_buffer(*buffer);
    }
  }

  (*buffer) = malloc(sizeof(grid_buffer));
  (*buffer)->size = requested_size;

#ifdef __GRID_CUDA
  // With size 0 cudaMallocHost doesn't null the pointer and cudaFreeHost fails.
  (*buffer)->host_buffer = NULL;
  CHECK(cudaMallocHost((void **)&(*buffer)->host_buffer, requested_size));
  CHECK(cudaMalloc((void **)&(*buffer)->device_buffer, requested_size));
#else
  (*buffer)->host_buffer = malloc(requested_size);
  (*buffer)->device_buffer = NULL;
#endif
}

/*******************************************************************************
 * \brief Deallocate given buffer.
 * \author Ole Schuett
 ******************************************************************************/
void grid_free_buffer(grid_buffer *buffer) {

  if (buffer == NULL)
    return;

#ifdef __GRID_CUDA
  CHECK(cudaFreeHost(buffer->host_buffer));
  CHECK(cudaFree(buffer->device_buffer));
#else
  free(buffer->host_buffer);
#endif

  free(buffer);
}

/*******************************************************************************
 * \brief Returns a pointer to the host buffer.
 * \author Ole Schuett
 ******************************************************************************/
double *grid_buffer_get_host_pointer(grid_buffer *buffer) {
  assert(buffer != NULL);
  return buffer->host_buffer;
}

// EOF
