/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2014 the CP2K developers group                      *
 *****************************************************************************/

#ifndef ACC_OPENCL_STREAM_H
#define ACC_OPENCL_STREAM_H

#if defined (__ACC) && defined (__OPENCL)

// struct definitions
typedef struct {
   acc_opencl_dev_type  device;
   cl_command_queue     queue;
} acc_opencl_stream_type;

#endif

#endif
//EOF
