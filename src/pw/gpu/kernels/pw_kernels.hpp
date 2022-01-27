/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#ifndef PW_KERNELS_HPP
#define PW_KERNELS_HPP

#ifndef __PW_CUDA_HIP_KERNELS
#error "this header file can only be included in pw_hip_z.cc or pw_cuda_z.cu"
#endif

/*******************************************************************************
 * \brief   Performs a (double precision complex) gather and scale on the GPU.
 * \version 0.01
 ******************************************************************************/
template <typename T>
__global__ void pw_gather_z(const T scale, const int ngpts,
                            const int *__restrict__ const ghatmap,
                            const T *__restrict__ const c,
                            T *__restrict__ pwcc) {

  const int igpt = blockIdx.x * blockDim.x + threadIdx.x;

  if (igpt >= ngpts)
    return;

  pwcc[2 * igpt] = scale * c[2 * ghatmap[igpt]];
  pwcc[2 * igpt + 1] = scale * c[2 * ghatmap[igpt] + 1];
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) scatter and scale on the GPU.
 * \version 0.01
 ******************************************************************************/
template <typename T>
__global__ void pw_scatter_z(const T scale, const int nmaps, const int ngpts,
                             const int *ghatmap, const T *pwcc, T *c) {
  const int igpt = blockIdx.x * blockDim.x + threadIdx.x;

  if (igpt >= ngpts)
    return;

  c[2 * ghatmap[igpt]] = scale * pwcc[2 * igpt];
  c[2 * ghatmap[igpt] + 1] = scale * pwcc[2 * igpt + 1];
  if (nmaps == 2) {
    c[2 * ghatmap[igpt + ngpts]] = scale * pwcc[2 * igpt];
    c[2 * ghatmap[igpt + ngpts] + 1] = -scale * pwcc[2 * igpt + 1];
  }
}

/* I really do not think this kernel will be faster than memset and a "blas"
 * dcopy */
template <typename T>
__global__ void real_to_complex_gpu(const int length__,
                                    const T *__restrict__ src__,
                                    T *const __restrict__ dst__) {
  const int ind = blockIdx.x * blockDim.x + threadIdx.x;

  if (ind >= length__)
    return;

  dst__[2 * ind] = src__[ind];
  dst__[2 * ind + 1] = 0.0;
}

/* I really do not think this kernel will be faster than memset and a "blas"
 * dcopy */
template <typename T>
__global__ void complex_to_real_gpu(const int length__,
                                    const T *__restrict__ src__,
                                    T *const __restrict__ dst__) {
  const int ind = blockIdx.x * blockDim.x + threadIdx.x;

  if (ind >= length__)
    return;

  dst__[ind] = src__[2 * ind];
}

#endif
