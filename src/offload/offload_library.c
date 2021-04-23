/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2021 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "offload_library.h"

static int current_device_id = -1;

/*******************************************************************************
 * \brief Returns the number of available devices.
 * \author Ole Schuett
 ******************************************************************************/
int offload_get_device_count() {
  int count = 0;
#ifdef __OFFLOAD_CUDA
  OFFLOAD_CHECK(cudaGetDeviceCount(&count));
#endif
  return count;
}

/*******************************************************************************
 * \brief Selects the device to be used.
 * \author Ole Schuett
 ******************************************************************************/
void offload_set_device_id(int device_id) { current_device_id = device_id; }

/*******************************************************************************
 * \brief Returns the device to be used.
 * \author Ole Schuett
 ******************************************************************************/
int offload_get_device_id() { return current_device_id; }

// EOF
