/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include "pw_gpu_kernels.h"

/*******************************************************************************
 * \brief   Performs a out-of-place copy of a double precision vector (first
 *          half filled) into a double precision complex vector on the GPU.
 *          It requires a global double precision vector 'zout' of length '2n'.
 *          [memory (shared):  none Byte
 *           memory (private): 4 Byte
 *           memory (global):  16*n Byte]
 *          n - size of double precision input vector
 * \author  Andreas Gloess
 ******************************************************************************/
__global__ void pw_real_to_complex(const double *din, double *zout,
                                   const int ngpts) {
  const int igpt = blockIdx.x * blockDim.x + threadIdx.x;
  if (igpt < ngpts) {
    zout[2 * igpt] = din[igpt];
    zout[2 * igpt + 1] = 0.0;
  }
}

/*******************************************************************************
 * \brief Launcher for pw_real_to_complex kernel.
 * \author Ole Schütt
 ******************************************************************************/
void pw_gpu_launch_real_to_complex(const double *din, double *zout,
                                   const int ngpts, cudaStream_t stream) {
  const int threadsPerBlock = 1024;
  const int numBlocks = (ngpts + threadsPerBlock - 1) / threadsPerBlock;
  pw_real_to_complex<<<numBlocks, threadsPerBlock, 0, stream>>>(din, zout,
                                                                ngpts);
}

/*******************************************************************************
 * \brief   Performs a out-of-place copy of a double precision complex vector
 *          (real part) into a double precision vector on the GPU.
 *          It requires a global double precision vector 'dout' of length 'n'.
 *          [memory (shared):  none Byte
 *           memory (private): 4 Byte
 *           memory (global):  16*n Byte]
 *          n - size of double precision output vector
 * \author  Andreas Gloess
 ******************************************************************************/
__global__ void pw_complex_to_real(const double *zin, double *dout,
                                   const int ngpts) {
  const int igpt = blockIdx.x * blockDim.x + threadIdx.x;
  if (igpt < ngpts) {
    dout[igpt] = zin[2 * igpt];
  }
}

/*******************************************************************************
 * \brief Launcher for pw_complex_to_real kernel.
 * \author Ole Schütt
 ******************************************************************************/
void pw_gpu_launch_complex_to_real(const double *zin, double *dout,
                                   const int ngpts, offloadStream_t stream) {
  const int threadsPerBlock = 1024;
  const int numBlocks = (ngpts + threadsPerBlock - 1) / threadsPerBlock;
  pw_complex_to_real<<<numBlocks, threadsPerBlock, 0, stream>>>(zin, dout,
                                                                ngpts);
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) gather and scale on the GPU.
 * \author  Andreas Gloess
 ******************************************************************************/
__global__ void pw_gather_z(double *pwcc, const double *c, const double scale,
                            const int ngpts, const int *ghatmap) {
  const int igpt = blockIdx.x * blockDim.x + threadIdx.x;
  if (igpt < ngpts) {
    pwcc[2 * igpt] = scale * c[2 * ghatmap[igpt]];
    pwcc[2 * igpt + 1] = scale * c[2 * ghatmap[igpt] + 1];
  }
}

/*******************************************************************************
 * \brief Launcher for pw_gather_z kernel.
 * \author Ole Schütt
 ******************************************************************************/
void pw_gpu_launch_gather_z(double *pwcc, const double *c, const double scale,
                            const int ngpts, const int *ghatmap,
                            offloadStream_t stream) {
  const int threadsPerBlock = 32;
  const int numBlocks = (ngpts + threadsPerBlock - 1) / threadsPerBlock;
  pw_gather_z<<<numBlocks, threadsPerBlock, 0, stream>>>(pwcc, c, scale, ngpts,
                                                         ghatmap);
}

/*******************************************************************************
 * \brief   Performs a (double precision complex) scatter and scale on the GPU.
 * \author  Andreas Gloess
 ******************************************************************************/
__global__ void pw_scatter_z(double *c, const double *pwcc, const double scale,
                             const int ngpts, const int nmaps,
                             const int *ghatmap) {
  const int igpt = blockIdx.x * blockDim.x + threadIdx.x;
  if (igpt < ngpts) {
    c[2 * ghatmap[igpt]] = scale * pwcc[2 * igpt];
    c[2 * ghatmap[igpt] + 1] = scale * pwcc[2 * igpt + 1];
    if (nmaps == 2) {
      c[2 * ghatmap[igpt + ngpts]] = scale * pwcc[2 * igpt];
      c[2 * ghatmap[igpt + ngpts] + 1] = -scale * pwcc[2 * igpt + 1];
    }
  }
}

/*******************************************************************************
 * \brief Launcher for pw_scatter_z kernel.
 * \author Ole Schütt
 ******************************************************************************/
void pw_gpu_launch_scatter_z(double *c, const double *pwcc, const double scale,
                             const int ngpts, const int nmaps,
                             const int *ghatmap, offloadStream_t stream) {
  const int threadsPerBlock = 32;
  const int numBlocks = (ngpts + threadsPerBlock - 1) / threadsPerBlock;
  pw_scatter_z<<<numBlocks, threadsPerBlock, 0, stream>>>(c, pwcc, scale, ngpts,
                                                          nmaps, ghatmap);
}

// EOF
