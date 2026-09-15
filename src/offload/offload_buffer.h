/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#ifndef OFFLOAD_BUFFER_H
#define OFFLOAD_BUFFER_H

#include <stddef.h>

/*******************************************************************************
 * \brief Internal representation of a buffer.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  size_t size;
  double *host_buffer;
  double *device_buffer;
} offload_buffer;

/*******************************************************************************
 * \brief Allocate a buffer (NULL or valid) with the given number of elements.
 * \author Ole Schuett
 ******************************************************************************/
void offload_create_buffer(const int length, offload_buffer **buffer);

/*******************************************************************************
 * \brief Deallocate given buffer.
 * \author Ole Schuett
 ******************************************************************************/
void offload_free_buffer(offload_buffer *buffer);

/*******************************************************************************
 * \brief Return a pointer to the host buffer (Fortran API).
 * \author Ole Schuett
 ******************************************************************************/
double *offload_get_buffer_host_pointer(offload_buffer *buffer);

/*******************************************************************************
 * \brief Return a pointer to the device buffer (Fortran API).
 * \author Ole Schuett
 ******************************************************************************/
double *offload_get_buffer_device_pointer(offload_buffer *buffer);

/*******************************************************************************
 * \brief Copy data from host to device (Fortran API).
 * \author Ole Schuett
 ******************************************************************************/
void offload_buffer_h2d(offload_buffer *buffer, const int length);

/*******************************************************************************
 * \brief Copy data from device to host (Fortran API).
 * \author Ole Schuett
 ******************************************************************************/
void offload_buffer_d2h(offload_buffer *buffer, const int length);

#endif

// EOF
