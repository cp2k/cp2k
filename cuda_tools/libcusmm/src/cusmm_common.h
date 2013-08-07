/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2013 the CP2K developers group                      *
 *****************************************************************************/

#ifndef CUSMM_COMMON_H
#define CUSMM_COMMON_H


/******************************************************************************
 * There is no nativ support for atomicAdd on doubles in Cuda 5.0. However the*
 * following implementation is provided in the CUDA C Programing guide.       *
 ******************************************************************************/
static __device__ double atomicAdd(double *address, double val) {
    unsigned long long int *address_as_ull =
        (unsigned long long int *) address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                                             __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}


/******************************************************************************
 ******************************************************************************/
static __device__ __forceinline__ double
dbl_shfl (double r, int lane)
{
  int hi = __shfl (__double2hiint (r), lane);
  int lo = __shfl (__double2loint (r), lane);
  return __hiloint2double (hi, lo);
}

#endif
