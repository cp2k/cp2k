/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015  CP2K developers group                         *
 *****************************************************************************/

#ifndef ACC_OPENCL_DEV_H
#define ACC_OPENCL_DEV_H

#if defined (__ACC) && defined (__OPENCL)

// maximum information line lenght
#define MAX_DEV_TYPE_LEN 3

// struct definitions
typedef struct {
   cl_platform_id   platform_id;
   cl_device_id     device_id;
   cl_context       ctx;
} acc_opencl_dev_type;

// global (per MPI) device information
extern uint acc_opencl_ndevices;
extern acc_opencl_dev_type *acc_opencl_devices;
extern acc_opencl_dev_type *acc_opencl_my_device;

// global configuration information
static uint acc_opencl_ndevices_configured = 0;
static uint acc_opencl_set_device_configured = 0;

#endif
#endif
//EOF
