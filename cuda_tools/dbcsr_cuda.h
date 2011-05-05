/******************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations
 *  Copyright (C) 2000 - 2011  Urban Borstnik and the CP2K developers group
 *****************************************************************************/

#ifndef MAX
#define MAX(a,b) ((a > b) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b) ((a < b) ? (a) : (b))
#endif


int cuda_error_check (cudaError_t cudaError);
