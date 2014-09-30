/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2014 the CP2K developers group                      *
 *****************************************************************************/


/* 
 *
 * NOTE: In OpenCL streams are called queues and the related device.ctx and
 *       device.platform is used in combination with it. Therefore we need 
 *       a struct 'acc_opencl_queue' which combines this information.
 *
 *       For convenience the routine names are called 'xxx_stream_xxx' to
 *       match the ACC interface.
 */

#if defined (__ACC) && defined (__OPENCL)

#include <CL/cl.h>
#include <string.h>
#include <stdio.h>

// defines error check functions and 'cl_error'
#include "acc_opencl_error.h"

// defines 'acc_opencl_my_device' and some default lenghts
#include "acc_opencl_dev.h"

// defines 'acc_opencl_stream_type' struct
#include "acc_opencl_stream.h"

// defines the ACC interface
#include "../include/acc.h"

static const int verbose_print = 0;


/****************************************************************************/
#ifdef __cplusplus
extern "C" {
#endif
int acc_stream_priority_range (int* least, int* greatest){
  // debug info
  if (verbose_print){
    fprintf(stdout, "\n +++ STREAM PRIORITY RANGE SETUP +++ \n");
    fprintf(stdout, " ---> Entering: acc_stream_priority_range.\n");
  }

  // NOTE: This functionality is not available in OpenCL.
  *least = -1;
  *greatest = -1;

  // debug info
  if (verbose_print){
    fprintf(stdout, " ---> Leaving: acc_stream_priority_range.\n");
  }

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
// NOTE: 'priority' and 'name' are ignored.
int acc_stream_create (void** stream_p, char* name, int priority){
  // debug info
  if (verbose_print){
    fprintf(stdout, "\n +++ STREAM CREATION +++ \n");
    fprintf(stdout, " ---> Entering: acc_stream_create.\n");
  }

  // get memory on pointer
  *stream_p = (void *) malloc(sizeof(acc_opencl_stream_type));

  // local queue pointer 
  acc_opencl_stream_type *clstream = (acc_opencl_stream_type *) *stream_p;
  (*clstream).device = *acc_opencl_my_device;

  // create a command queue
  cl_command_queue_properties queue_properties = 0;
  (*clstream).queue = (cl_command_queue) clCreateCommandQueue(
                                           (*acc_opencl_my_device).ctx,       // cl_context                  context
                                           (*acc_opencl_my_device).device_id, // cl_device_id                device
                                           queue_properties,                  // cl_command_queue_properties properties
                                           &cl_error);                        // cl_int                      *errcode_ret
  if (acc_opencl_error_check(cl_error, __LINE__))
    return -1;

  // debug info
  if (verbose_print){
    fprintf(stdout, " +    STREAM address:  HEX=%p INT=%ld\n", &((*clstream).queue), (uintptr_t) &((*clstream).queue));
    fprintf(stdout, " +    STREAM value:  %u\n", (*clstream).queue);
    fprintf(stdout, " ---> Leaving: acc_stream_create.\n");
  }

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
int acc_stream_destroy (void* stream){
  // debug info
  if (verbose_print){
    fprintf(stdout, "\n +++ STREAM DESTRUCTION +++ \n");
    fprintf(stdout, " ---> Entering: acc_stream_destroy.\n");
  }

  // local queue pointer 
  acc_opencl_stream_type *clstream = (acc_opencl_stream_type *) stream;

  // release the command queue
  cl_error = clReleaseCommandQueue((*clstream).queue);
  if (acc_opencl_error_check(cl_error, __LINE__))
    return -1;
  // free the struct acc_opencl_queue 'stream'
  free(clstream);

  // debug info
  if (verbose_print){
    fprintf(stdout, " -    STREAM value:  %u\n", (*clstream).queue);
    fprintf(stdout, " ---> Leaving: acc_stream_destroy.\n");
  }

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
int acc_stream_sync (void* stream){
  // debug info
  if (verbose_print){
    fprintf(stdout, "\n +++ STREAM SYNCHRONIZATION +++ \n");
    fprintf(stdout, " ---> Entering: acc_stream_sync.\n");
  }

  // local queue pointer 
  acc_opencl_stream_type *clstream = (acc_opencl_stream_type *) stream;

  // synchronize the command queue
  // A ' clEnqueueBarrier is probably enough
  cl_error = clFlush((*clstream).queue);
//ToDo: Flush sends all commands in a queue to the device but does not
//      guarantee that they will be processed while return to host.
//  cl_error = clFinish((*clstream).queue);

  if (acc_opencl_error_check(cl_error, __LINE__))
    return -1;

  // debug info
  if (verbose_print){
    fprintf(stdout, "      STREAM value:  %u\n", (*clstream).queue);
    fprintf(stdout, " ---> Leaving: acc_stream_sync.\n");
  }

  // assign return value
  return 0;
}
#ifdef __cplusplus
}
#endif

#endif
//EOF
