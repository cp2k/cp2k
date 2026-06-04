/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#include "offload_buffer.h"
#include "offload_library.h"
#include "offload_mempool.h"
#include "offload_runtime.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#if 1
#define OFFLOAD_BUFFER_MEMPOOL
#endif

/*******************************************************************************
 * \brief Internal routine to deallocate given buffer.
 * \author Hans Pabst
 ******************************************************************************/
static void offload_free_buffer_internal(offload_buffer *buffer) {
  assert(NULL != buffer);
#if defined(OFFLOAD_BUFFER_MEMPOOL)
  offload_mempool_host_free(buffer->host_buffer);
  offload_mempool_device_free(buffer->device_buffer);
#elif defined(__OFFLOAD)
  offloadFreeHost(buffer->host_buffer);
  offloadFree(buffer->device_buffer);
#else
  free(buffer->host_buffer);
  assert(NULL == buffer->device_buffer);
#endif
}

/*******************************************************************************
 * \brief Allocates a buffer (NULL or valid) with the given number of elements.
 * \author Ole Schuett and Hans Pabst
 ******************************************************************************/
void offload_create_buffer(const int length, offload_buffer **buffer) {
  const size_t requested_size = sizeof(double) * length;

  if (*buffer != NULL) {
    if (requested_size <= (*buffer)->size) {
      return; // reuse existing buffer
    } else {
      offload_free_buffer_internal(*buffer);
    }
  } else {
    (*buffer) = malloc(sizeof(offload_buffer));
    assert(NULL != *buffer);
  }

  (*buffer)->size = requested_size;
  (*buffer)->host_buffer = NULL;
  (*buffer)->device_buffer = NULL;
#if defined(OFFLOAD_BUFFER_MEMPOOL)
#if !defined(__OFFLOAD_UNIFIED_MEMORY)
  (*buffer)->host_buffer = offload_mempool_host_malloc(requested_size);
#endif
  (*buffer)->device_buffer = offload_mempool_device_malloc(requested_size);
#elif defined(__OFFLOAD)
  offload_activate_chosen_device();
#if !defined(__OFFLOAD_UNIFIED_MEMORY)
  offloadMallocHost((void **)&(*buffer)->host_buffer, requested_size);
#endif
  offloadMalloc((void **)&(*buffer)->device_buffer, requested_size);
#else
  (*buffer)->host_buffer = malloc(requested_size);
  (*buffer)->device_buffer = NULL;
#endif
  if (NULL == (*buffer)->host_buffer) { /* unified memory */
    (*buffer)->host_buffer = (*buffer)->device_buffer;
  }
}

/*******************************************************************************
 * \brief Deallocate given buffer.
 * \author Ole Schuett
 ******************************************************************************/
void offload_free_buffer(offload_buffer *buffer) {
  if (NULL != buffer) {
    offload_free_buffer_internal(buffer);
    free(buffer);
  }
}

/*******************************************************************************
 * \brief Returns a pointer to the host buffer (Fortran API).
 * \author Ole Schuett
 ******************************************************************************/
double *offload_get_buffer_host_pointer(offload_buffer *buffer) {
  assert(NULL != buffer);
  return buffer->host_buffer;
}

// EOF
