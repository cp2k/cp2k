/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2014 the CP2K developers group                      *
 *****************************************************************************/

#if defined (__ACC) && defined (__OPENCL)

#include <CL/cl.h>
#include <string.h>
#include <stdio.h>

// defines error check functions and 'cl_error'
#include "acc_opencl_error.h"

// defines 'acc_opencl_my_device' and some default lenghts
#include "acc_opencl_dev.h"
uint acc_opencl_ndevices;
acc_opencl_dev_type *acc_opencl_devices;
acc_opencl_dev_type *acc_opencl_my_device;

// defines the ACC interface
#include "../include/acc.h"

static const int verbose_print = 0;

/****************************************************************************/
#ifdef __cplusplus
extern "C" {
#endif
int acc_get_ndevices (int *n_devices){
  // debug info
  if (verbose_print) fprintf(stdout, "Entering: acc_get_ndevices.\n");

  if (acc_opencl_ndevices_configured) {
    // just get it from the global variable
    *n_devices = acc_opencl_ndevices;
  }
  else {
    // declarations
    int nplats_avail, ndevs_avail, ndevs;
    uint i, j, k;
    cl_uint plat_count, dev_count;
    cl_platform_id *platforms;
    cl_device_id *devices;
    cl_uint device_max_compute_units;
    cl_bitfield device_type;
    size_t vendor_name_size, device_name_size;
    char *vendor_name = NULL, *device_name = NULL;
    char dev_type[MAX_DEV_TYPE_LEN];
  
    // initialization
    nplats_avail = 0;
    ndevs_avail = 0;
    ndevs = 0;
  
    // get platforms
    cl_error = clGetPlatformIDs(0, NULL, &nplats_avail);
    if (acc_opencl_error_check(cl_error, __LINE__))
      return -1;
    platforms = (cl_platform_id *) malloc(nplats_avail * sizeof(cl_platform_id));
    cl_error = clGetPlatformIDs(nplats_avail, platforms, NULL);
    if (acc_opencl_error_check(cl_error, __LINE__))
      return -1;
  
    // print some informations
    for (i = 0; i < 80; i++) fprintf(stdout, "-"); fprintf(stdout, "\n");
    if (nplats_avail > 1) {
      fprintf(stdout, " OPENCL| Number of platforms:                 %d --- IDs: [0-%d]\n", nplats_avail, nplats_avail-1);
    } else {
      fprintf(stdout, " OPENCL| Number of platforms:                 %d --- IDs: [0]\n", nplats_avail);
    }
  
    // loop over platforms and count devices
    for (i = 0; i < nplats_avail; i++) {
      // get vendor name
      cl_error = clGetPlatformInfo(platforms[i], CL_PLATFORM_VENDOR, 0, NULL, &vendor_name_size);
      if (acc_opencl_error_check(cl_error, __LINE__))
        return -1;
      vendor_name = (char *) malloc(vendor_name_size * sizeof(char));
      cl_error = clGetPlatformInfo(platforms[i], CL_PLATFORM_VENDOR, vendor_name_size, vendor_name, NULL);
      if (acc_opencl_error_check(cl_error, __LINE__))
        return -1;
  
      // get number of devices for this platform and increase overall counter
      cl_error = clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 0, NULL, &dev_count);
      if (acc_opencl_error_check(cl_error, __LINE__))
        return -1;
      devices = (cl_device_id *) malloc(dev_count * sizeof(cl_device_id));
      cl_error = clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, dev_count, devices, NULL);
      if (acc_opencl_error_check(cl_error, __LINE__))
        return -1;
      ndevs_avail += (int) dev_count;
  
      // print some informations
      fprintf(stdout, " OPENCL| "); for (j = 0; j < 71; j++) fprintf(stdout, "."); fprintf(stdout, "\n");
      fprintf(stdout, " OPENCL| Platform [%d] vendor:                 %s\n", i, vendor_name);
      if (dev_count > 1) {
        fprintf(stdout, " OPENCL| Number of devices for platform[%d]:   %d --- IDs: [0-%d]\n", i, dev_count, dev_count-1);
      } else {
        fprintf(stdout, " OPENCL| Number of devices for platform[%d]:   %d --- IDs: [0]\n", i, dev_count);
      }

      // free memory
      free(vendor_name);
  
      // loop over devices
      for (j = 0; j < dev_count; j++) {
          // get device name 
          cl_error = clGetDeviceInfo(devices[j], CL_DEVICE_NAME, 0, NULL, &device_name_size);
          if (acc_opencl_error_check(cl_error, __LINE__))
            return -1;
          device_name = (char *) malloc(device_name_size * sizeof(char));
          cl_error = clGetDeviceInfo(devices[j], CL_DEVICE_NAME, device_name_size, device_name, NULL);
          if (acc_opencl_error_check(cl_error, __LINE__))
            return -1;

          // get device type and choose only GPU, CPU and/or ACC
          cl_error = clGetDeviceInfo(devices[j], CL_DEVICE_TYPE, sizeof(device_type), &device_type, NULL);
          if (acc_opencl_error_check(cl_error, __LINE__))
            return -1;
          switch(device_type) {
            case(CL_DEVICE_TYPE_DEFAULT):     strcpy(dev_type, "DEF"); break;
            case(CL_DEVICE_TYPE_GPU):         strcpy(dev_type, "GPU"); ndevs += 1; break;
            case(CL_DEVICE_TYPE_CPU):         strcpy(dev_type, "CPU"); ndevs += 1; break;
            case(CL_DEVICE_TYPE_ACCELERATOR): strcpy(dev_type, "ACC"); ndevs += 1; break;
            case(CL_DEVICE_TYPE_ALL):         strcpy(dev_type, "ALL"); break;
            default: fprintf(stderr, "OPENCL Error: Uninplemented device type.\n"); exit(-1);
          }

          // get device max compute units 
          cl_error = clGetDeviceInfo(devices[j], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(device_max_compute_units), &device_max_compute_units, NULL);
          if (acc_opencl_error_check(cl_error, __LINE__))
            return -1;

          // print some informations
          fprintf(stdout, " OPENCL| Device [%d] name:                     %s\n", j, device_name);
          fprintf(stdout, " OPENCL| Device [%d] type:                     %s\n", j, dev_type);
          fprintf(stdout, " OPENCL| Device [%d] max compute units:        %d\n", j, device_max_compute_units);
          fprintf(stdout, " OPENCL| Unified device number:                %d\n", (ndevs - 1));

          // free memory
          free(device_name);
      }
      // free memory
      free(devices);
    }
  
    // print some informations
    fprintf(stdout, " OPENCL| "); for (j = 0; j < 71; j++) fprintf(stdout, "."); fprintf(stdout, "\n");
    fprintf(stdout, " OPENCL| Number of available devices:          %d\n", ndevs_avail);
    fprintf(stdout, " OPENCL| Number of selected devices:           %d\n", ndevs);
    for (i = 0; i < 80; i++) fprintf(stdout, "-"); fprintf(stdout, "\n");
  
    // allocate space for selected devices
    acc_opencl_devices = (acc_opencl_dev_type *) malloc(ndevs * sizeof(acc_opencl_dev_type));
  
    // second loop over platforms
    k = 0;
    for (i = 0; i < nplats_avail; i++) {
      // get number of devices for this platform
      cl_error = clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 0, NULL, &dev_count);
      if (acc_opencl_error_check(cl_error, __LINE__))
        return -1;
      devices = (cl_device_id *) malloc(dev_count * sizeof(cl_device_id));
      cl_error = clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, dev_count, devices, NULL);
      if (acc_opencl_error_check(cl_error, __LINE__))
        return -1;
  
      // loop over devices
      for (j = 0; j < dev_count; j++) {
        // get device type and choose only GPU, CPU and/or ACC
        cl_error = clGetDeviceInfo(devices[j], CL_DEVICE_TYPE, sizeof(device_type), &device_type, NULL);
        if (acc_opencl_error_check(cl_error, __LINE__))
          return -1;
        switch(device_type) {
          case(CL_DEVICE_TYPE_GPU):
          case(CL_DEVICE_TYPE_CPU):
          case(CL_DEVICE_TYPE_ACCELERATOR):
            acc_opencl_devices[k].platform_id = platforms[i];
            acc_opencl_devices[k].device_id   = devices[j];
            k += 1;
            break;
          default:
            fprintf(stdout, " OPENCL| Skipping platform [%d] device [%d].\n", i, j);
            break;
        }
      }
      // free memory
      free(devices);
    }
  
    // check that both loops count equivalent
    if (k != ndevs) return -1;
  
    // free memory
    free(platforms);
  
    // assign number of useable devices
    acc_opencl_ndevices = ndevs;
    *n_devices = acc_opencl_ndevices;
  
    //set configuration flag
    acc_opencl_ndevices_configured = 1;
  }

  // debug info
  if (verbose_print) fprintf(stdout, "Leaving: acc_get_ndevices.\n");

  // assign return value
  return 0;
}
#ifdef __cplusplus
}
#endif


/****************************************************************************/
#ifdef __cplusplus
extern "C" {
#endif
int acc_set_active_device (int device_id){
  // debug info
  if (verbose_print) fprintf(stdout, "Entering: acc_set_active_device.\n");

  // declarations
  uint i;

  if (acc_opencl_set_device_configured) {
    // reset to (new) 'device_id'
    (*acc_opencl_my_device).platform_id = acc_opencl_devices[device_id].platform_id;
    (*acc_opencl_my_device).device_id = acc_opencl_devices[device_id].device_id;
    (*acc_opencl_my_device).ctx = acc_opencl_devices[device_id].ctx;
  }
  else {
    // allocate space for selected devices
    acc_opencl_my_device = (acc_opencl_dev_type *) malloc(sizeof(acc_opencl_dev_type));
  
  
    // initialize device number 'device_id'
    for (i = 0; i < 80; i++) fprintf(stdout, "="); fprintf(stdout, "\n");
    fprintf(stdout, " OpenCL| My device number: %d \n", device_id);
    for (i = 0; i < 80; i++) fprintf(stdout, "="); fprintf(stdout, "\n");
    (*acc_opencl_my_device).platform_id = acc_opencl_devices[device_id].platform_id;
    (*acc_opencl_my_device).device_id = acc_opencl_devices[device_id].device_id;
  
    // create context with platform_id and device_id
    cl_context_properties ctx_properties[3] = {CL_CONTEXT_PLATFORM, (cl_context_properties) (*acc_opencl_my_device).platform_id, 0};
    acc_opencl_devices[device_id].ctx = clCreateContext(ctx_properties, (cl_uint) 1, &(*acc_opencl_my_device).device_id, NULL, NULL, &cl_error);
    if (acc_opencl_error_check(cl_error, __LINE__))
      return -1;
    (*acc_opencl_my_device).ctx = acc_opencl_devices[device_id].ctx;
  
    // just configure once
    if (!acc_opencl_set_device_configured) acc_opencl_set_device_configured = 1;
  }

  // debug info
  if (verbose_print) fprintf(stdout, "Leaving: acc_set_active_device.\n");

  // assign return value
  return 0;
}
#ifdef __cplusplus
}
#endif

#endif
//EOF
