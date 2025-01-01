/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2025 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#ifndef PW_GPU_KERNELS_H
#define PW_GPU_KERNELS_H

#include "../../offload/offload_runtime.h"
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)

#ifdef __cplusplus
extern "C" {
#endif

/*******************************************************************************
 * \brief Launcher for pw_real_to_complex kernel.
 * \author Ole Schuett
 ******************************************************************************/
void pw_gpu_launch_real_to_complex(const double *din, double *zout,
                                   const int ngpts, offloadStream_t stream);

/*******************************************************************************
 * \brief Launcher for pw_complex_to_real kernel.
 * \author Ole Schuett
 ******************************************************************************/
void pw_gpu_launch_complex_to_real(const double *zin, double *dout,
                                   const int ngpts, offloadStream_t stream);

/*******************************************************************************
 * \brief Launcher for pw_gather kernel.
 * \author Ole Schuett
 ******************************************************************************/
void pw_gpu_launch_gather(double *pwcc, const double *c, const double scale,
                          const int ngpts, const int *ghatmap,
                          offloadStream_t stream);

/*******************************************************************************
 * \brief Launcher for pw_scatter kernel.
 * \author Ole Schuett
 ******************************************************************************/
void pw_gpu_launch_scatter(double *c, const double *pwcc, const double scale,
                           const int ngpts, const int nmaps, const int *ghatmap,
                           offloadStream_t stream);

#ifdef __cplusplus
}
#endif

#endif // defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
#endif

// EOF
