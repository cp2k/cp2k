/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2021 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#ifdef __GRID_CUDA

#include <algorithm>
#include <assert.h>
#include <cuda.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define GRID_DO_COLLOCATE 1
#include "../common/grid_common.h"
#include "../common/grid_prepare_pab.h"
#include "grid_gpu_collint.h"
#include "grid_gpu_collocate.h"

/*******************************************************************************
 * \brief Collocate a single grid point with distance d{xyz} from center.
 * \author Ole Schuett
 ******************************************************************************/
__device__ static void cxyz_to_gridpoint(const double dx, const double dy,
                                         const double dz, const double zetp,
                                         const int lp, const cxyz_store *cxyz,
                                         double *gridpoint) {

  // Squared distance of point from center.
  const double r2 = dx * dx + dy * dy + dz * dz;
  const double gaussian = exp(-zetp * r2);

  // accumulate into register
  double gridpoint_reg = 0.0;

  // Manually unrolled loops based on terms in coset_inv.
  // For lp > 6 the register usage increases and the kernel becomes slower.
  gridpoint_reg += cxyz[0];

  if (lp >= 1) {
    gridpoint_reg += cxyz[1] * dx;
    gridpoint_reg += cxyz[2] * dy;
    gridpoint_reg += cxyz[3] * dz;
    if (lp >= 2) {
      const double dx2 = dx * dx;
      const double dy2 = dy * dy;
      const double dz2 = dz * dz;
      gridpoint_reg += cxyz[4] * dx2;
      gridpoint_reg += cxyz[5] * dx * dy;
      gridpoint_reg += cxyz[6] * dx * dz;
      gridpoint_reg += cxyz[7] * dy2;
      gridpoint_reg += cxyz[8] * dy * dz;
      gridpoint_reg += cxyz[9] * dz2;
      if (lp >= 3) {
        const double dx3 = dx2 * dx;
        const double dy3 = dy2 * dy;
        const double dz3 = dz2 * dz;
        gridpoint_reg += cxyz[10] * dx3;
        gridpoint_reg += cxyz[11] * dx2 * dy;
        gridpoint_reg += cxyz[12] * dx2 * dz;
        gridpoint_reg += cxyz[13] * dx * dy2;
        gridpoint_reg += cxyz[14] * dx * dy * dz;
        gridpoint_reg += cxyz[15] * dx * dz2;
        gridpoint_reg += cxyz[16] * dy3;
        gridpoint_reg += cxyz[17] * dy2 * dz;
        gridpoint_reg += cxyz[18] * dy * dz2;
        gridpoint_reg += cxyz[19] * dz3;
        if (lp >= 4) {
          const double dx4 = dx3 * dx;
          const double dy4 = dy3 * dy;
          const double dz4 = dz3 * dz;
          gridpoint_reg += cxyz[20] * dx4;
          gridpoint_reg += cxyz[21] * dx3 * dy;
          gridpoint_reg += cxyz[22] * dx3 * dz;
          gridpoint_reg += cxyz[23] * dx2 * dy2;
          gridpoint_reg += cxyz[24] * dx2 * dy * dz;
          gridpoint_reg += cxyz[25] * dx2 * dz2;
          gridpoint_reg += cxyz[26] * dx * dy3;
          gridpoint_reg += cxyz[27] * dx * dy2 * dz;
          gridpoint_reg += cxyz[28] * dx * dy * dz2;
          gridpoint_reg += cxyz[29] * dx * dz3;
          gridpoint_reg += cxyz[30] * dy4;
          gridpoint_reg += cxyz[31] * dy3 * dz;
          gridpoint_reg += cxyz[32] * dy2 * dz2;
          gridpoint_reg += cxyz[33] * dy * dz3;
          gridpoint_reg += cxyz[34] * dz4;
          if (lp >= 5) {
            const double dx5 = dx4 * dx;
            const double dy5 = dy4 * dy;
            const double dz5 = dz4 * dz;
            gridpoint_reg += cxyz[35] * dx5;
            gridpoint_reg += cxyz[36] * dx4 * dy;
            gridpoint_reg += cxyz[37] * dx4 * dz;
            gridpoint_reg += cxyz[38] * dx3 * dy2;
            gridpoint_reg += cxyz[39] * dx3 * dy * dz;
            gridpoint_reg += cxyz[40] * dx3 * dz2;
            gridpoint_reg += cxyz[41] * dx2 * dy3;
            gridpoint_reg += cxyz[42] * dx2 * dy2 * dz;
            gridpoint_reg += cxyz[43] * dx2 * dy * dz2;
            gridpoint_reg += cxyz[44] * dx2 * dz3;
            gridpoint_reg += cxyz[45] * dx * dy4;
            gridpoint_reg += cxyz[46] * dx * dy3 * dz;
            gridpoint_reg += cxyz[47] * dx * dy2 * dz2;
            gridpoint_reg += cxyz[48] * dx * dy * dz3;
            gridpoint_reg += cxyz[49] * dx * dz4;
            gridpoint_reg += cxyz[50] * dy5;
            gridpoint_reg += cxyz[51] * dy4 * dz;
            gridpoint_reg += cxyz[52] * dy3 * dz2;
            gridpoint_reg += cxyz[53] * dy2 * dz3;
            gridpoint_reg += cxyz[54] * dy * dz4;
            gridpoint_reg += cxyz[55] * dz5;
            if (lp >= 6) {
              const double dx6 = dx5 * dx;
              const double dy6 = dy5 * dy;
              const double dz6 = dz5 * dz;
              gridpoint_reg += cxyz[56] * dx6;
              gridpoint_reg += cxyz[57] * dx5 * dy;
              gridpoint_reg += cxyz[58] * dx5 * dz;
              gridpoint_reg += cxyz[59] * dx4 * dy2;
              gridpoint_reg += cxyz[60] * dx4 * dy * dz;
              gridpoint_reg += cxyz[61] * dx4 * dz2;
              gridpoint_reg += cxyz[62] * dx3 * dy3;
              gridpoint_reg += cxyz[63] * dx3 * dy2 * dz;
              gridpoint_reg += cxyz[64] * dx3 * dy * dz2;
              gridpoint_reg += cxyz[65] * dx3 * dz3;
              gridpoint_reg += cxyz[66] * dx2 * dy4;
              gridpoint_reg += cxyz[67] * dx2 * dy3 * dz;
              gridpoint_reg += cxyz[68] * dx2 * dy2 * dz2;
              gridpoint_reg += cxyz[69] * dx2 * dy * dz3;
              gridpoint_reg += cxyz[70] * dx2 * dz4;
              gridpoint_reg += cxyz[71] * dx * dy5;
              gridpoint_reg += cxyz[72] * dx * dy4 * dz;
              gridpoint_reg += cxyz[73] * dx * dy3 * dz2;
              gridpoint_reg += cxyz[74] * dx * dy2 * dz3;
              gridpoint_reg += cxyz[75] * dx * dy * dz4;
              gridpoint_reg += cxyz[76] * dx * dz5;
              gridpoint_reg += cxyz[77] * dy6;
              gridpoint_reg += cxyz[78] * dy5 * dz;
              gridpoint_reg += cxyz[79] * dy4 * dz2;
              gridpoint_reg += cxyz[80] * dy3 * dz3;
              gridpoint_reg += cxyz[81] * dy2 * dz4;
              gridpoint_reg += cxyz[82] * dy * dz5;
              gridpoint_reg += cxyz[83] * dz6;
            }
          }
        }
      }
    }
  }

  // Handle higher values of lp.
  if (lp >= 7) {
    for (int i = 84; i < ncoset(lp); i++) {
      double val = cxyz[i];
      const orbital a = coset_inv[i];
      for (int j = 0; j < a.l[0]; j++) {
        val *= dx;
      }
      for (int j = 0; j < a.l[1]; j++) {
        val *= dy;
      }
      for (int j = 0; j < a.l[2]; j++) {
        val *= dz;
      }
      gridpoint_reg += val;
    }
  }

  atomicAddDouble(gridpoint, gridpoint_reg * gaussian);
}

/*******************************************************************************
 * \brief Collocates coefficients C_xyz onto the grid.
 * \author Ole Schuett
 ******************************************************************************/
__device__ static void cxyz_to_grid(const kernel_params *params,
                                    const smem_task *task, const double *cxyz,
                                    double *grid) {

  if (task->use_orthorhombic_kernel) {
    ortho_cxyz_to_grid(params, task, cxyz, grid);
  } else {
    general_cxyz_to_grid(params, task, cxyz, grid);
  }
}

/*******************************************************************************
 * \brief Adds given value to matrix element cab[idx(b)][idx(a)].
 * \author Ole Schuett
 ******************************************************************************/
__device__ static inline void prep_term(const orbital a, const orbital b,
                                        const double value, const int n,
                                        double *cab) {
  atomicAddDouble(&cab[idx(b) * n + idx(a)], value);
}

/*******************************************************************************
 * \brief Decontracts the subblock, going from spherical to cartesian harmonics.
 * \author Ole Schuett
 ******************************************************************************/
template <bool IS_FUNC_AB>
__device__ static void block_to_cab(const kernel_params *params,
                                    const smem_task *task, double *cab) {

  // The spherical index runs over angular momentum and then over contractions.
  // The carthesian index runs over exponents and then over angular momentum.

  // Decontract block, apply prepare_pab, and store in cab.
  // This is a double matrix product. Since the pab block can be quite large the
  // two products are fused to conserve shared memory.
  for (int i = threadIdx.z; i < task->nsgf_setb; i += blockDim.z) {
    for (int j = threadIdx.y; j < task->nsgf_seta; j += blockDim.y) {
      double block_val;
      if (task->block_transposed) {
        block_val = task->pab_block[j * task->nsgfb + i] * task->off_diag_twice;
      } else {
        block_val = task->pab_block[i * task->nsgfa + j] * task->off_diag_twice;
      }

      if (IS_FUNC_AB) {
        // fast path for common case
        const int jco_start = task->first_cosetb + threadIdx.x;
        for (int jco = jco_start; jco < task->ncosetb; jco += blockDim.x) {
          const double sphib = task->sphib[i * task->maxcob + jco];
          for (int ico = task->first_coseta; ico < task->ncoseta; ico++) {
            const double sphia = task->sphia[j * task->maxcoa + ico];
            const double pab_val = block_val * sphia * sphib;
            atomicAddDouble(&cab[jco * task->ncoseta + ico], pab_val);
          }
        }
      } else {
        // Since prepare_pab is a register hog we use it only when really needed
        const int jco_start = task->first_cosetb + threadIdx.x;
        for (int jco = jco_start; jco < task->ncosetb; jco += blockDim.x) {
          const orbital b = coset_inv[jco];
          for (int ico = task->first_coseta; ico < task->ncoseta; ico++) {
            const orbital a = coset_inv[ico];
            const double sphia = task->sphia[j * task->maxcoa + idx(a)];
            const double sphib = task->sphib[i * task->maxcob + idx(b)];
            const double pab_val = block_val * sphia * sphib;
            prepare_pab(params->func, a, b, task->zeta, task->zetb, pab_val,
                        task->n1, cab);
          }
        }
      }
    }
  }
  __syncthreads(); // because of concurrent writes to cab
}

/*******************************************************************************
 * \brief Cuda kernel for collocating all tasks of one grid level.
 * \author Ole Schuett
 ******************************************************************************/
template <bool IS_FUNC_AB>
__device__ static void collocate_kernel(const kernel_params *params) {

  // Copy task from global to shared memory and precompute some stuff.
  __shared__ smem_task task;
  load_task(params, &task);

  // Check if radius is below the resolution of the grid.
  if (2.0 * task.radius < task.dh_max) {
    return; // nothing to do
  }

  // Allot dynamic shared memory.
  extern __shared__ double shared_memory[];
  double *smem_cab = &shared_memory[params->smem_cab_offset];
  double *smem_alpha = &shared_memory[params->smem_alpha_offset];
  double *smem_cxyz = &shared_memory[params->smem_cxyz_offset];

  zero_cab(&task, smem_cab);
  block_to_cab<IS_FUNC_AB>(params, &task, smem_cab);

  compute_alpha(params, &task, smem_alpha);
  cab_to_cxyz(params, &task, smem_alpha, smem_cab, smem_cxyz);

  cxyz_to_grid(params, &task, smem_cxyz, params->grid);
}

/*******************************************************************************
 * \brief Specialized Cuda kernel that can only collocate GRID_FUNC_AB.
 * \author Ole Schuett
 ******************************************************************************/
__global__ static void collocate_kernel_density(const kernel_params params) {
  collocate_kernel<true>(&params);
}

/*******************************************************************************
 * \brief Cuda kernel that can collocate any function, ie. GRID_FUNC_*.
 * \author Ole Schuett
 ******************************************************************************/
__global__ static void collocate_kernel_anyfunc(const kernel_params params) {
  collocate_kernel<false>(&params);
}

/*******************************************************************************
 * \brief Launches the Cuda kernel that collocates all tasks of one grid level.
 * \author Ole Schuett
 ******************************************************************************/
void grid_gpu_collocate_one_grid_level(
    const grid_gpu_task_list *task_list, const int first_task,
    const int last_task, const enum grid_func func,
    const grid_gpu_layout *layout, const cudaStream_t stream,
    const double *pab_blocks_dev, double *grid_dev, int *lp_diff) {

  // Compute max angular momentum.
  const prepare_ldiffs ldiffs = prepare_get_ldiffs(func);
  *lp_diff = ldiffs.la_max_diff + ldiffs.lb_max_diff; // for reporting stats
  const int la_max = task_list->lmax + ldiffs.la_max_diff;
  const int lb_max = task_list->lmax + ldiffs.lb_max_diff;
  const int lp_max = la_max + lb_max;

  const int ntasks = last_task - first_task + 1;
  if (ntasks == 0) {
    return; // Nothing to do and lp_diff already set.
  }

  init_constant_memory();

  // Compute required shared memory.
  // TODO: Currently, cab's indicies run over 0...ncoset[lmax],
  //       however only ncoset(lmin)...ncoset(lmax) are actually needed.
  const int cab_len = ncoset(lb_max) * ncoset(la_max);
  const int alpha_len = 3 * (lb_max + 1) * (la_max + 1) * (lp_max + 1);
  const int cxyz_len = ncoset(lp_max);
  const size_t smem_per_block =
      (cab_len + alpha_len + cxyz_len) * sizeof(double);

  if (smem_per_block > 48 * 1024) {
    fprintf(stderr, "ERROR: Not enough shared memory in grid_gpu_collocate.\n");
    fprintf(stderr, "cab_len: %i, ", cab_len);
    fprintf(stderr, "alpha_len: %i, ", alpha_len);
    fprintf(stderr, "cxyz_len: %i, ", cxyz_len);
    fprintf(stderr, "total smem_per_block: %f kb\n\n", smem_per_block / 1024.0);
    abort();
  }

  // kernel parameters
  kernel_params params;
  params.smem_cab_offset = 0;
  params.smem_alpha_offset = cab_len;
  params.smem_cxyz_offset = params.smem_alpha_offset + alpha_len;
  params.first_task = first_task;
  params.func = func;
  params.grid = grid_dev;
  params.la_min_diff = ldiffs.la_min_diff;
  params.lb_min_diff = ldiffs.lb_min_diff;
  params.la_max_diff = ldiffs.la_max_diff;
  params.lb_max_diff = ldiffs.lb_max_diff;
  params.tasks = task_list->tasks_dev;
  params.pab_blocks = pab_blocks_dev;
  memcpy(params.dh, layout->dh, 9 * sizeof(double));
  memcpy(params.dh_inv, layout->dh_inv, 9 * sizeof(double));
  memcpy(params.npts_global, layout->npts_global, 3 * sizeof(int));
  memcpy(params.npts_local, layout->npts_local, 3 * sizeof(int));
  memcpy(params.shift_local, layout->shift_local, 3 * sizeof(int));

  // Launch !
  const int nblocks = ntasks;
  const dim3 threads_per_block(4, 4, 4);

  if (func == GRID_FUNC_AB) {
    collocate_kernel_density<<<nblocks, threads_per_block, smem_per_block,
                               stream>>>(params);
  } else {
    collocate_kernel_anyfunc<<<nblocks, threads_per_block, smem_per_block,
                               stream>>>(params);
  }
}

#endif // __GRID_CUDA
// EOF
