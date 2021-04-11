/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2021 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/
#ifndef GRID_GPU_INTEGRATE_H
#define GRID_GPU_INTEGRATE_H

#ifdef __GRID_CUDA

#include "grid_gpu_task_list.h"
#include <cuda_runtime.h>

#ifdef __cplusplus
extern "C" {
#endif

/*******************************************************************************
 * \brief Launches the Cuda kernel that integrates all tasks of one grid level.
 * \author Ole Schuett
 ******************************************************************************/
void grid_gpu_integrate_one_grid_level(
    const grid_gpu_task_list *task_list, const int first_task,
    const int last_task, const bool compute_tau, const grid_gpu_layout *layout,
    const cudaStream_t stream, const double *pab_blocks_dev,
    const double *grid_dev, double *hab_blocks_dev, double *forces_dev,
    double *virial_dev, int *lp_diff);

#ifdef __cplusplus
}
#endif

#endif // __GRID_CUDA
#endif
// EOF
