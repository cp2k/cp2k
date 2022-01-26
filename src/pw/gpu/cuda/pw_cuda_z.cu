/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include "../../../offload/offload_library.h"
#include <cstdio>
#include <cstdlib>
#include <cuda.h>
#include <cuda_runtime.h>

#define __PW_CUDA_HIP_KERNELS

#include "../kernels/pw_kernels.hpp"

void gpu_scatter(cudaStream_t &stream__, const double scale__, int num_points__,
                 const int nmaps__, const int *map_index__,
                 const double *dataIn__, double *dataOut__) {
  dim3 blocksPerGrid, threadsPerBlock;

  blocksPerGrid.x = num_points__ / 512 + ((num_points__ % 512) != 0);
  threadsPerBlock.x = 512;

  pw_scatter_z<double><<<blocksPerGrid, threadsPerBlock, 0, stream__>>>(
      scale__, nmaps__, num_points__, map_index__, dataIn__, dataOut__);
}

void gpu_gather(cudaStream_t &stream__, const double scale__, int num_points__,
                const int *map_index__, const double *dataIn__,
                double *dataOut__) {
  dim3 blocksPerGrid, threadsPerBlock;

  blocksPerGrid.x = num_points__ / 512 + ((num_points__ % 512) != 0);
  threadsPerBlock.x = 512;

  pw_gather_z<double><<<blocksPerGrid, threadsPerBlock, 0, stream__>>>(
      scale__, num_points__, map_index__, dataIn__, dataOut__);
}

void real_to_complex(cudaStream_t &stream__, const int length__,
                     const double *src__, double *const dst__) {
  dim3 blocksPerGrid, threadsPerBlock;
  blocksPerGrid.x = length__ / 512 + ((length__ % 512) != 0);
  threadsPerBlock.x = 512;
  real_to_complex_gpu<double>
      <<<blocksPerGrid, threadsPerBlock, 0, stream__>>>(length__, src__, dst__);
}

void complex_to_real(cudaStream_t &stream__, const int length__,
                     const double *src__, double *const dst__) {
  dim3 blocksPerGrid, threadsPerBlock;
  blocksPerGrid.x = length__ / 512 + ((length__ % 512) != 0);
  threadsPerBlock.x = 512;
  complex_to_real_gpu<double>
      <<<blocksPerGrid, threadsPerBlock, 0, stream__>>>(length__, src__, dst__);
}
