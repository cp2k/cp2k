/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015 the CP2K developers group                      *
 *****************************************************************************/

#if defined (__ACC) && defined (__OPENCL)

#include <CL/cl.h>
#include <stdio.h>

// defines error check functions and 'cl_error'
#include "acc_opencl_error.h"

/****************************************************************************/
int acc_opencl_error_check (cl_int cl_error, int line){
  int pid;

  if (cl_error != CL_SUCCESS) {
    pid = getpid();
    fprintf(stderr, "%d OPENCL RT Error line: %d, ERROR_CODE: %d\n", pid, line, cl_error);
    fflush(stdout);
    fflush(stderr);
    return -1;
  }
  return 0;
}

#endif
//EOF
