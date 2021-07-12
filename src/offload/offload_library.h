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
#elif defined(__GRID_HIP)
#define __OFFLOAD_HIP
#endif

#if defined(__OFFLOAD_CUDA)
#include <cuda.h>
#include <cuda_runtime.h>
#elif defined(__OFFLOAD_HIP)
#include <hip/hip_runtime_api.h>
#endif

#if defined(__OFFLOAD_CUDA)
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

#if defined(__OFFLOAD_HIP)
/*******************************************************************************
 * \brief Checks given rocm status and upon failure abort with a nice message.
 ******************************************************************************/
#define OFFLOAD_CHECK(status)                                                  \
  if (status != hipSuccess) {                                                  \
    fprintf(stderr, "ERROR: %s %s %d\n", hipGetErrorString(status), __FILE__,  \
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
int offload_get_device_count(void);

/*******************************************************************************
 * \brief Selects the device to be used.
 * \author Ole Schuett
 ******************************************************************************/
void offload_set_device_id(int device_id);

/*******************************************************************************
 * \brief Returns the device to be used.
 * \author Ole Schuett
 ******************************************************************************/
int offload_get_device_id(void);

/*******************************************************************************
 * \brief Activates the device selected via offload_set_device_id()
 * \author Ole Schuett
 ******************************************************************************/
void offload_set_device(void);

/*******************************************************************************
 * \brief Starts a timing range.
 * \author Ole Schuett
 ******************************************************************************/
void offload_timeset(const char *message);

/*******************************************************************************
 * \brief Ends a timing range.
 * \author Ole Schuett
 ******************************************************************************/
void offload_timestop(void);

/*******************************************************************************
 * \brief Gets free and total device memory.
 * \author Ole Schuett
 ******************************************************************************/
void offload_mem_info(size_t *free, size_t *total);

#ifdef __cplusplus
}
#endif

#endif

// EOF
