/******************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations
 *  Copyright (C) 2000 - 2013  Urban Borstnik and the CP2K developers group
 *****************************************************************************/
#ifndef DBCSR_CUDA_H
#define DBCSR_CUDA_H

#ifndef MAX
#define MAX(a,b) ((a > b) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b) ((a < b) ? (a) : (b))
#endif


// function is defined in dbcsr_cuda_dev.cu
int cuda_error_check (cudaError_t cudaError);




#define GROUPING 16

/* to get this really threadprivate we need nvcc with --compiler-options -fopenmp 
 * however, if all devices have the same properties, the only sane thing, this will be OK more generally */
extern struct cudaDeviceProp devProperties;
#pragma omp threadprivate(devProperties)


#endif
