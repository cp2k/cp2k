/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#ifndef PW_GPU_KERNELS_H
#define PW_GPU_KERNELS_H

#include "../../offload/offload_operations.h"

#ifdef __cplusplus
extern "C" {
#endif

/*******************************************************************************
 * \brief Launcher for pw_real_to_complex kernel.
 * \author Ole Schütt
 ******************************************************************************/
void pw_gpu_launch_real_to_complex(const double *din, double *zout,
                                   const int ngpts, offloadStream_t stream);

/*******************************************************************************
 * \brief Launcher for pw_complex_to_real kernel.
 * \author Ole Schütt
 ******************************************************************************/
void pw_gpu_launch_complex_to_real(const double *zin, double *dout,
                                   const int ngpts, offloadStream_t stream);

/*******************************************************************************
 * \brief Launcher for pw_gather_z kernel.
 * \author Ole Schütt
 ******************************************************************************/
void pw_gpu_launch_gather_z(double *pwcc, const double *c, const double scale,
                            const int ngpts, const int *ghatmap,
                            offloadStream_t stream);

/*******************************************************************************
 * \brief Launcher for pw_scatter_z kernel.
 * \author Ole Schütt
 ******************************************************************************/
void pw_gpu_launch_scatter_z(double *c, const double *pwcc, const double scale,
                             const int ngpts, const int nmaps,
                             const int *ghatmap, offloadStream_t stream);

#ifdef __cplusplus
}
#endif

#endif

// EOF
