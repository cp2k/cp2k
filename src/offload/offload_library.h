/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2021 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/
#ifndef OFFLOAD_LIBRARY_H
#define OFFLOAD_LIBRARY_H

#if defined(__GRID_CUDA) || defined(__PW_CUDA)
#define __OFFLOAD_CUDA
#endif

#ifdef __OFFLOAD_CUDA
#include <cuda_runtime.h>

/*******************************************************************************
 * \brief Checks given Cuda status and upon failure abort with a nice message.
 * \author Ole Schuett
 ******************************************************************************/
#define OFFLOAD_CHECK(status)                                                  \
  if (status != cudaSuccess) {                                                 \
    fprintf(stderr, "ERROR: %s %s %d\n", cudaGetErrorString(status), __FILE__, \
            __LINE__);                                                         \
    abort();                                                                   \
  }

#endif

#ifdef __cplusplus
extern "C" {
#endif

/*******************************************************************************
 * \brief Returns the number of available devices.
 * \author Ole Schuett
 ******************************************************************************/
int offload_get_device_count();

/*******************************************************************************
 * \brief Selects the device to be used.
 * \author Ole Schuett
 ******************************************************************************/
void offload_set_device_id(int device_id);

/*******************************************************************************
 * \brief Returns the device to be used.
 * \author Ole Schuett
 ******************************************************************************/
int offload_get_device_id();

#ifdef __cplusplus
}
#endif

#endif

// EOF
