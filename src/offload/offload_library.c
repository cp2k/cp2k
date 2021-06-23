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
int offload_get_device_count(void) {
  int count = 0;
#ifdef __OFFLOAD_CUDA
  OFFLOAD_CHECK(cudaGetDeviceCount(&count));
#elif defined(__OFFLOAD_HIP)
  OFFLOAD_CHECK(hipGetDeviceCount(&count));
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
int offload_get_device_id(void) { return current_device_id; }

/*******************************************************************************
 * \brief Activates the device selected via offload_set_device_id()
 * \author Ole Schuett
 ******************************************************************************/
void offload_set_device(void) {
#ifdef __OFFLOAD_CUDA
  OFFLOAD_CHECK(cudaSetDevice(current_device_id));
#elif defined(__OFFLOAD_HIP)
  OFFLOAD_CHECK(hipSetDevice(current_device_id));
#endif
}

// EOF
