/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015  CP2K developers group                         *
 *****************************************************************************/

#if defined (__ACC) && defined (__OPENCL)

#include <CL/cl.h>
#include <string.h>
#include <stdio.h>

// defines error check functions and 'cl_error'
#include "acc_opencl_error.h"

// defines 'acc_opencl_my_device' and some default lenghts
#include "acc_opencl_dev.h"

// defines 'acc_opencl_stream_type'
#include "acc_opencl_stream.h"

// defines the ACC interface
#include "../include/acc.h"

// debug flag
static const int verbose_print = 0;

#ifdef __cplusplus
extern "C" {
#endif

/****************************************************************************/
int acc_event_create (void** event_p){
  // debug info
  if (verbose_print){
    fprintf(stdout, "\n ... EVENT CREATION ... \n");
    fprintf(stdout, " ---> Entering: acc_event_create.\n");
  }

  // local event object pointer
  *event_p = malloc(sizeof(cl_event));
  cl_event *clevent = (cl_event *) *event_p;

  // get a device event object
  *clevent = clCreateUserEvent((*acc_opencl_my_device).ctx, &cl_error);
  if (acc_opencl_error_check(cl_error, __LINE__))
    return -1;

  // debug info
  if (verbose_print){
    fprintf(stdout, " ---> Leaving: acc_event_create.\n");
  }

  cl_error = clSetUserEventStatus(*clevent, CL_COMPLETE);

  // assign return value
  return 0;
}


/****************************************************************************/
int acc_event_destroy (void* event){
  // debug info
  if (verbose_print){
    fprintf(stdout, "\n ... EVENT DESTRUCTION ... \n");
    fprintf(stdout, " ---> Entering: acc_event_destroy.\n");
  }

  // local event object pointer
  cl_event *clevent = (cl_event *) event;

  // release event object
  cl_error = clReleaseEvent(*clevent);
  if (acc_opencl_error_check(cl_error, __LINE__))
    return -1;
  free(clevent);

  // debug info
  if (verbose_print){
    fprintf(stdout, " ---> Leaving: acc_event_destroy.\n");
  }

  // assign return value
  return 0;
}


/****************************************************************************/
int acc_event_record (void* event, void* stream){
  // debug info
  if (verbose_print){
    fprintf(stdout, "\n ... EVENT RECORDING ... \n");
    fprintf(stdout, " ---> Entering: acc_event_record.\n");
  }

  // local event and queue pointers
  cl_event *clevent = (cl_event *) event;
  acc_opencl_stream_type *clstream = (acc_opencl_stream_type *) stream;

  // set a marker 'event' to listen on to queue 'stream' 
  cl_error = clEnqueueMarker((*clstream).queue, clevent);
  if (acc_opencl_error_check(cl_error, __LINE__))
    return -1;

  // debug info
  if (verbose_print){
    fprintf(stdout, " ---> Leaving: acc_event_record.\n");
  }

  // assign return value
  return 0;
}


/****************************************************************************/
int acc_event_query (void *event, int *has_occured){
  //declarations
  cl_int param_value;

  // debug info
  if (verbose_print){
    fprintf(stdout, "\n ... EVENT QUERYING ... \n");
    fprintf(stdout, " ---> Entering: acc_event_query.\n");
  }

  // local event pointer
  cl_event *clevent = (cl_event *) event;

  // get event status
  cl_error = clGetEventInfo(
              *clevent,                          // cl_event      event
              CL_EVENT_COMMAND_EXECUTION_STATUS, // cl_event_info param_name
              (size_t) sizeof(cl_int),           // size_t        param_value_size
              &param_value,                      // void          *param_value
              NULL);                             // size_t        *param_value_size_ret
  if (acc_opencl_error_check(cl_error, __LINE__))
    return -1;

  // check event status
  if (param_value == CL_COMPLETE){
    *has_occured = 1;
  } else {
    *has_occured = 0;
  }

  // debug info
  if (verbose_print) fprintf(stdout, "Leaving: acc_event_query.\n");
  if (verbose_print){
    fprintf(stdout, "      Result: %d\n", *has_occured);
    fprintf(stdout, " ---> Entering: acc_event_query.\n");
  }

  // assign return value
  return 0;
}


/****************************************************************************/
int acc_stream_wait_event (void* stream, void* event){
  // debug info
  if (verbose_print){
    fprintf(stdout, "\n ... STREAM waiting for EVENT ... \n");
    fprintf(stdout, " ---> Entering: acc_stream_wait_event.\n");
  }

  // local event and queue pointers
  cl_event *clevent = (cl_event *) event;
  acc_opencl_stream_type *clstream = (acc_opencl_stream_type *) stream;

  // wait for an event on a stream
  cl_error = clEnqueueWaitForEvents(  // cl_int
               (*clstream).queue,     // cl_command_queue command_queue
               (cl_uint) 1,           // cl_uint          num_events
               clevent);              // cnst cl_event    *event_list
  if (acc_opencl_error_check(cl_error, __LINE__))
    return -1;

  // debug info
  if (verbose_print){
    fprintf(stdout, " ---> Leaving: acc_stream_wait_event.\n");
  }

  // assign return value
  return 0;
}


/****************************************************************************/
int acc_event_synchronize (void* event){
  // debug info
  if (verbose_print){
    fprintf(stdout, "\n ... EVENT SYNCHRONIZATION ... \n");
    fprintf(stdout, " ---> Entering: acc_event_synchronize.\n");
  }

  // local event and queue pointers
  cl_event *clevent = (cl_event *) event;

  // wait for an event ( !!! need to share the same ctx !!! )
  cl_error = clWaitForEvents((cl_uint) 1, clevent);
  if (acc_opencl_error_check(cl_error, __LINE__))
    return -1;

  // debug info
  if (verbose_print){
    fprintf(stdout, " ---> Leaving: acc_event_synchronize.\n");
  }

  // assign return value
  return 0;
}


#ifdef __cplusplus
}
#endif

#endif
//EOF
