/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/
#ifndef GRID_BUFFER_H
#define GRID_BUFFER_H

#include <stddef.h>

/*******************************************************************************
 * \brief Internal representation of a buffer.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  size_t size;
  double *host_buffer;
  double *device_buffer;
} grid_buffer;

/*******************************************************************************
 * \brief Allocates a buffer of given length, ie. number of elements.
 * \author Ole Schuett
 ******************************************************************************/
void grid_create_buffer(const int length, grid_buffer **buffer);

/*******************************************************************************
 * \brief Deallocate given buffer.
 * \author Ole Schuett
 ******************************************************************************/
void grid_free_buffer(grid_buffer *buffer);

/*******************************************************************************
 * \brief Returns a pointer to the host buffer.
 * \author Ole Schuett
 ******************************************************************************/
double *grid_buffer_get_host_pointer(grid_buffer *buffer);

#endif

// EOF
