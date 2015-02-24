/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015  CP2K developers group                         *
 *****************************************************************************/

#ifndef ACC_OPENCL_ERROR_H
#define ACC_OPENCL_ERROR_H

#if defined (__ACC) && defined (__OPENCL)
// define global opencl error variable and type
typedef cl_int cl_error_type;
cl_error_type cl_error;

// define custom OpenCL error check function
int acc_opencl_error_check (cl_int cl_error, int line);

#endif

#endif
//EOF
