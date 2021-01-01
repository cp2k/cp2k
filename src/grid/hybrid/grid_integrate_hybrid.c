/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2021 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#ifdef __GRID_CUDA

#include <assert.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "../common/grid_common.h"
#include "../common/grid_library.h"
#include "../cpu/coefficients.h"
#include "../cpu/collocation_integration.h"
#include "../cpu/cpu_private_header.h"
#include "../cpu/grid_context_cpu.h"
#include "../cpu/grid_prepare_pab_dgemm.h"
#include "../cpu/tensor_local.h"
#include "../cpu/utils.h"
#include "grid_hybrid_task_list.h"

void grid_cpu_integrate_task_list(
    void *const ptr, const bool orthorhombic, const bool compute_tau,
    const int natoms, const int nlevels, const int npts_global[nlevels][3],
    const int npts_local[nlevels][3], const int shift_local[nlevels][3],
    const int border_width[nlevels][3], const double dh[nlevels][3][3],
    const double dh_inv[nlevels][3][3], const grid_buffer *const pab_blocks,
    const double *const grid[nlevels], grid_buffer *hab_blocks,
    double forces[natoms][3], double virial[3][3]);

void grid_hybrid_integrate_task_list(
    void *const ptr, const bool orthorhombic, const bool compute_tau,
    const int natoms, const int nlevels, const int npts_global[nlevels][3],
    const int npts_local[nlevels][3], const int shift_local[nlevels][3],
    const int border_width[nlevels][3], const double dh[nlevels][3][3],
    const double dh_inv[nlevels][3][3], const grid_buffer *const pab_blocks,
    const double *const grid[nlevels], grid_buffer *hab_blocks,
    double forces[natoms][3], double virial[3][3]) {
  grid_cpu_integrate_task_list(ptr, orthorhombic, compute_tau, natoms, nlevels,
                               npts_global, npts_local, shift_local,
                               border_width, dh, dh_inv, pab_blocks, grid,
                               hab_blocks, forces, virial);
}

#endif
