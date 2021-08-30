/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2021 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "offload_library.h"

#if defined(__OFFLOAD_PROFILING)
#if defined(__OFFLOAD_CUDA)
#include <nvToolsExt.h>
#elif defined(__OFFLOAD_HIP) && defined(__HIP_PLATFORM_AMD__)
#include <roctracer/roctx.h>
#endif
#endif

static int current_device_id = -1;

const uint32_t colormap[] = {0xFFFFFF00,  // Yellow
                             0xFFFF00FF,  // Fuchsia
                             0xFFFF0000,  // Red
                             0xFFC0C0C0,  // Silver
                             0xFF808080,  // Gray
                             0xFF808000,  // Olive
                             0xFF800080,  // Purple
                             0xFF800000,  // Maroon
                             0xFF00FFFF,  // Aqua
                             0xFF00FF00,  // Lime
                             0xFF008080,  // Teal
                             0xFF008000,  // Green
                             0xFF0000FF,  // Blue
                             0xFF000080}; // Navy

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

/*******************************************************************************
 * \brief Starts a timing range.
 * \author Ole Schuett
 ******************************************************************************/
void offload_timeset(const char *message) {
#if defined(__OFFLOAD_PROFILING)
#if defined(__OFFLOAD_CUDA)
  // colors are picked based on a (very simple) hash value of the message
  int hash = 0;
  for (size_t i = 0; i < strlen(message); i++) {
    hash += i * message[i] * message[i];
  }
  nvtxEventAttributes_t eventAttrib = {0};
  eventAttrib.version = NVTX_VERSION;
  eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
  eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII;
  eventAttrib.message.ascii = message;
  eventAttrib.colorType = NVTX_COLOR_ARGB;
  eventAttrib.color = colormap[hash % 14];
  eventAttrib.payloadType = NVTX_PAYLOAD_TYPE_INT64;
  eventAttrib.payload.llValue = 123;
  eventAttrib.category = 42;
  nvtxRangePushEx(&eventAttrib);
#elif defined(__OFFLOAD_HIP) && defined(__HIP_PLATFORM_AMD__)
  roctxRangePushA(message);
#endif
#endif
  (void)message; // mark argument as used
}

/*******************************************************************************
 * \brief Ends a timing range.
 * \author Ole Schuett
 ******************************************************************************/
void offload_timestop(void) {
#if defined(__OFFLOAD_PROFILING)
#if defined(__OFFLOAD_CUDA)
  nvtxRangePop();
#elif defined(__OFFLOAD_HIP) && defined(__HIP_PLATFORM_AMD__)
  roctxRangePop();
#endif
#endif
}

/*******************************************************************************
 * \brief Gets free and total device memory.
 * \author Ole Schuett
 ******************************************************************************/
void offload_mem_info(size_t *free, size_t *total) {
#if defined(__OFFLOAD_CUDA)
  OFFLOAD_CHECK(cudaMemGetInfo(free, total));
#elif defined(__OFFLOAD_HIP)
  OFFLOAD_CHECK(hipMemGetInfo(free, total));
#else
  *free = 0;
  *total = 0;
#endif
}

// EOF
