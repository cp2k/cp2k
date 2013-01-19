/******************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations
 *  Copyright (C) 2000 - 2013  Urban Borstnik and the CP2K developers group
 *****************************************************************************/

#ifndef MAX
#define MAX(a,b) ((a > b) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b) ((a < b) ? (a) : (b))
#endif


int cuda_error_check (cudaError_t cudaError);

#define GROUPING 16

/* to get this really threadprivate we need nvcc with --compiler-options -fopenmp 
 * however, if all devices have the same properties, the only sane thing, this will be OK more generally */
extern struct cudaDeviceProp devProperties;
#pragma omp threadprivate(devProperties)

cudaStream_t dc_get_stream (int stream_id);

/* The following are only required or used by the 23 squared case */
#define SQUARESIZE 23
#define BLOCK 24
#define TDIM 3
#define NUMTHREADS23SQ ((BLOCK*BLOCK)/(TDIM*TDIM))
#define BLOCKSPLIT23SQ 2
#define TOTALTHREADS23SQ ((NUMTHREADS23SQ)*(BLOCKSPLIT23SQ))
