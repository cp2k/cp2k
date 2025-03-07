/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2025 CP2K developers group <https://cp2k.org>              */
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
 * \brief Allocates a buffer of given length, ie., number of elements.
 * \author Ole Schuett
 ******************************************************************************/
void offload_create_buffer(const int length, offload_buffer **buffer) {
  const size_t requested_size = sizeof(double) * length;

  if (*buffer != NULL) {
    if ((*buffer)->size >= requested_size) {
      return; // reuse existing buffer
    } else {
      offload_free_buffer(*buffer);
    }
  }

  (*buffer) = malloc(sizeof(offload_buffer));
  (*buffer)->size = requested_size;
  (*buffer)->host_buffer = NULL;
  (*buffer)->device_buffer = NULL;
#if defined(__OFFLOAD)
  offload_activate_chosen_device();
  offloadMallocHost((void **)&(*buffer)->host_buffer, requested_size);
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
  if (NULL == buffer)
    return;
#if defined(__OFFLOAD)
  offloadFreeHost(buffer->host_buffer);
  offloadFree(buffer->device_buffer);
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
  assert(NULL != buffer);
  return buffer->host_buffer;
}

// EOF
