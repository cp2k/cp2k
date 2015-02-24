/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015  CP2K developers group                         *
 *****************************************************************************/

#include <cuda_runtime.h>
#include <stdio.h>
#include <math.h>
#include "acc_cuda_error.h"
#include "../include/acc.h"

static const int verbose_print = 0;

/****************************************************************************/
extern "C" int acc_event_create(void** event_p){
  *event_p = malloc(sizeof(cudaEvent_t));
  cudaEvent_t* cuevent = (cudaEvent_t*) *event_p;

  cudaError_t cErr = cudaEventCreate(cuevent);
  if(verbose_print) printf("cuda_event_created:  %p -> %d\n", *event_p, *cuevent);
  if (cuda_error_check(cErr)) return -1;
  if (cuda_error_check(cudaGetLastError())) return -1;
  return 0;
}


/****************************************************************************/
extern "C" int acc_event_destroy(void* event){
    cudaEvent_t* cuevent = (cudaEvent_t*) event;

    if(verbose_print) printf("cuda_event_destroy called\n");
    cudaError_t cErr = cudaEventDestroy(*cuevent);
    free(cuevent);
    if (cuda_error_check (cErr)) return -1;
    if (cuda_error_check(cudaGetLastError ())) return -1;
    return 0;
}


/****************************************************************************/
extern "C" int acc_event_record(void* event, void* stream){
    cudaEvent_t* cuevent = (cudaEvent_t*) event;
    cudaStream_t* custream = (cudaStream_t*) stream;

    if(verbose_print) printf("cuda_event_record: %p -> %d,  %p -> %d\n", cuevent, *cuevent,  custream, *custream);
    cudaError_t cErr = cudaEventRecord (*cuevent, *custream);
    if (cuda_error_check (cErr)) return -1;
    //if (cuda_error_check(cudaGetLastError ())) return -1;
    return 0;
}


/****************************************************************************/
extern "C" int acc_event_query(void* event, int* has_occured){
    if(verbose_print) printf("cuda_event_query called\n");

    cudaEvent_t* cuevent = (cudaEvent_t*) event;
    cudaError_t cErr = cudaEventQuery(*cuevent);
    //if(cuda_error_check(cudaGetLastError ())) return -1;
    if(cErr==cudaSuccess){
         *has_occured = 1;
         return 0;
    }

    if(cErr==cudaErrorNotReady){
        *has_occured = 0;
        return 0;
    }

    return -1; // something went wrong
}


/****************************************************************************/
extern "C" int acc_stream_wait_event(void* stream, void* event){
    if(verbose_print) printf("cuda_stream_wait_event called\n");

    cudaEvent_t* cuevent = (cudaEvent_t*) event;
    cudaStream_t* custream = (cudaStream_t*) stream;

    // flags: Parameters for the operation (must be 0)
    cudaError_t cErr = cudaStreamWaitEvent(*custream, *cuevent, 0);
    if (cuda_error_check (cErr)) return -1;
    //if (cuda_error_check(cudaGetLastError ())) return -1;
    return 0;
}


/****************************************************************************/
extern "C" int acc_event_synchronize(void* event){
    if(verbose_print) printf("cuda_event_synchronize called\n");
    cudaEvent_t* cuevent = (cudaEvent_t*) event;
    cudaError_t cErr = cudaEventSynchronize(*cuevent);
    if (cuda_error_check (cErr)) return -1;
    if (cuda_error_check(cudaGetLastError ())) return -1;
    return 0;
}

//EOF
