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
  fill_smem_task(params, &task);

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
    const int last_task, const enum grid_func func, const int npts_global[3],
    const int npts_local[3], const int shift_local[3],
    const int border_width[3], const double dh[3][3], const double dh_inv[3][3],
    const cudaStream_t stream, const double *pab_blocks_dev, double *grid_dev,
    int *lp_diff) {

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
  params.orthorhombic = task_list->orthorhombic;
  params.func = func;
  params.grid = grid_dev;
  params.la_min_diff = ldiffs.la_min_diff;
  params.lb_min_diff = ldiffs.lb_min_diff;
  params.la_max_diff = ldiffs.la_max_diff;
  params.lb_max_diff = ldiffs.lb_max_diff;
  params.tasks = task_list->tasks_dev;
  params.atom_kinds = task_list->atom_kinds_dev;
  params.basis_sets = task_list->basis_sets_dev;
  params.block_offsets = task_list->block_offsets_dev;
  params.atom_positions = task_list->atom_positions_dev;
  params.pab_blocks = pab_blocks_dev;
  memcpy(params.dh, dh, 9 * sizeof(double));
  memcpy(params.dh_inv, dh_inv, 9 * sizeof(double));
  memcpy(params.npts_global, npts_global, 3 * sizeof(int));
  memcpy(params.npts_local, npts_local, 3 * sizeof(int));
  memcpy(params.shift_local, shift_local, 3 * sizeof(int));
  memcpy(params.border_width, border_width, 3 * sizeof(int));

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
