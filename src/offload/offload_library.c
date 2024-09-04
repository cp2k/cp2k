/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2024 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "offload_library.h"
#include "offload_runtime.h"

#if defined(__OFFLOAD_CUDA)
#include <cuda.h>
#elif defined(__OFFLOAD_HIP)
#include <hip/hip_runtime_api.h>
#endif

#if defined(__OFFLOAD_PROFILING)
#if defined(__OFFLOAD_CUDA)
#include <nvToolsExt.h>
#elif defined(__OFFLOAD_HIP) && defined(__HIP_PLATFORM_AMD__)
#include <roctracer/roctx.h>
#endif
#endif

static int chosen_device_id = -1;

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
 * \brief Initialize runtime.
 * \author Rocco Meli
 ******************************************************************************/
void offload_init(void) {
#if defined(__OFFLOAD_CUDA)
  CUresult error = cuInit(0);
  if (error != CUDA_SUCCESS) {
    fprintf(stderr, "ERROR: %s %d %s %d\n", "cuInit failed with error: ", error,
            __FILE__, __LINE__);
    abort();
  }
#elif defined(__OFFLOAD_HIP)
  OFFLOAD_CHECK(hipInit(0));
#elif defined(__OFFLOAD_OPENCL)
  OFFLOAD_CHECK(c_dbcsr_acc_init());
#endif
}

/*******************************************************************************
 * \brief Returns the number of available devices.
 * \author Ole Schuett
 ******************************************************************************/
int offload_get_device_count(void) {
  int count = 0;
#if defined(__OFFLOAD_CUDA)
  OFFLOAD_CHECK(cudaGetDeviceCount(&count));
#elif defined(__OFFLOAD_HIP)
  OFFLOAD_CHECK(hipGetDeviceCount(&count));
#elif defined(__OFFLOAD_OPENCL)
  OFFLOAD_CHECK(c_dbcsr_acc_get_ndevices(&count));
#endif
  return count;
}

/*******************************************************************************
 * \brief Selects the chosen device to be used.
 * \author Ole Schuett
 ******************************************************************************/
void offload_set_chosen_device(int device_id) { chosen_device_id = device_id; }

/*******************************************************************************
 * \brief Returns the chosen device.
 * \author Ole Schuett
 ******************************************************************************/
int offload_get_chosen_device(void) { return chosen_device_id; }

/*******************************************************************************
 * \brief Activates the device selected via offload_set_chosen_device()
 * \author Ole Schuett
 ******************************************************************************/
void offload_activate_chosen_device(void) {
#if defined(__OFFLOAD_CUDA)
  OFFLOAD_CHECK(cudaSetDevice(chosen_device_id));
#elif defined(__OFFLOAD_HIP)
  OFFLOAD_CHECK(hipSetDevice(chosen_device_id));
#elif defined(__OFFLOAD_OPENCL)
  OFFLOAD_CHECK(c_dbcsr_acc_set_active_device(chosen_device_id));
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
#elif defined(__OFFLOAD_OPENCL)
  OFFLOAD_CHECK(c_dbcsr_acc_dev_mem_info(free, total));
#else
  *free = *total = 0;
#endif
}

int offload_host_malloc(void **ptr__, const size_t size__) {
#if defined(__OFFLOAD)
  offloadMallocHost(ptr__, size__); /* checked */
  return offloadSuccess;
#else
  *ptr__ = malloc(size__);
  return EXIT_SUCCESS;
#endif
}

int offload_host_free(void *ptr__) {
#if defined(__OFFLOAD)
  offloadFreeHost(ptr__); /* checked */
  return offloadSuccess;
#else
  free(ptr__);
  return EXIT_SUCCESS;
#endif
}

// EOF
