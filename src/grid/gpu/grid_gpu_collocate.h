/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#ifndef GRID_GPU_COLLOCATE_H
#define GRID_GPU_COLLOCATE_H

#include "grid_gpu_task_list.h"
#include <cuda_runtime.h>

#ifdef __cplusplus
extern "C" {
#endif

/*******************************************************************************
 * \brief Launches the Cuda kernel that collocates all tasks of one grid level.
 * \author Ole Schuett
 ******************************************************************************/
void grid_gpu_collocate_one_grid_level(
    const grid_gpu_task_list *task_list, const int first_task,
    const int last_task, const enum grid_func func,
    const grid_gpu_layout *layout, const cudaStream_t stream,
    const double *pab_blocks_dev, double *grid_dev, int *lp_diff);

#ifdef __cplusplus
}
#endif
#endif
// EOF
