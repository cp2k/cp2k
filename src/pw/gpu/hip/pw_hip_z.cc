/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

// global dependencies
#include <cstdio>
#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>
#include <hipblas.h>
#include <hipfft.h>

#define __PW_CUDA_HIP_KERNELS
#include "../kernels/pw_kernels.hpp"

void gpu_scatter(hipStream_t &stream__, const double scale__, int num_points__,
                 const int nmaps__, const int *map_index__,
                 const double *dataIn__, double *dataOut__) {
  dim3 blocksPerGrid, threadsPerBlock;

  blocksPerGrid.x = num_points__ / 512 + ((num_points__ % 512) != 0);
  threadsPerBlock.x = 512;

  pw_scatter_z<double><<<blocksPerGrid, threadsPerBlock, 0, stream__>>>(
      scale__, nmaps__, num_points__, map_index__, dataIn__, dataOut__);
}

void gpu_gather(hipStream_t &stream__, const double scale__, int num_points__,
                const int *map_index__, const double *dataIn__,
                double *dataOut__) {
  dim3 blocksPerGrid, threadsPerBlock;

  blocksPerGrid.x = num_points__ / 512 + ((num_points__ % 512) != 0);
  threadsPerBlock.x = 512;

  pw_gather_z<double><<<blocksPerGrid, threadsPerBlock, 0, stream__>>>(
      scale__, num_points__, map_index__, dataIn__, dataOut__);
}
