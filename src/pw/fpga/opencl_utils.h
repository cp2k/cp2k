/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2021 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

/*******************************************************************************
 * \author Arjun Ramaswami
 ******************************************************************************/

#ifndef OPENCL_UTILS_H
#define OPENCL_UTILS_H

#if defined(__PW_FPGA)

extern void queue_cleanup();
extern void cleanup();
// Search for a platform that contains the search string
// Returns platform id if found
// Return NULL if none found
cl_platform_id findPlatform(char *platform_name);

// Search for a device based on the platform
// Return array of device ids
cl_device_id *getDevices(cl_platform_id pid, cl_device_type device_type,
                         cl_uint *num_devices);

// OpenCL program created for all devices of the context with the same binary
cl_program getProgramWithBinary(cl_context context, const cl_device_id *devices,
                                unsigned num_devices, int N[3],
                                char *data_path);

void openCLContextCallBackFxn(const char *errinfo, const void *private_info,
                              size_t cb, void *user_data);

void *alignedMalloc(size_t size);

void printError(cl_int error);

void _checkError(const char *file, int line, const char *func, cl_int err,
                 const char *msg, ...);

#define checkError(status, ...)                                                \
  _checkError(__FILE__, __LINE__, __FUNCTION__, status, __VA_ARGS__)

#endif

#endif // OPENCL_UTILS_H
