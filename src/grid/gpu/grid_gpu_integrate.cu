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

#define GRID_DO_COLLOCATE 0
#include "../common/grid_common.h"
#include "../common/grid_process_vab.h"
#include "grid_gpu_collint.h"
#include "grid_gpu_integrate.h"

/*******************************************************************************
 * \brief Contracts the subblock, going from cartesian harmonics to spherical.
 * \author Ole Schuett
 ******************************************************************************/
template <bool COMPUTE_TAU>
__device__ static void store_hab(const kernel_params *params,
                                 const smem_task *task, const double *cab) {

  // The spherical index runs over angular momentum and then over contractions.
  // The carthesian index runs over exponents and then over angular momentum.

  // This is a double matrix product. Since the block can be quite large the
  // two products are fused to conserve shared memory.
  for (int i = threadIdx.x; i < task->nsgf_setb; i += blockDim.x) {
    for (int j = threadIdx.y; j < task->nsgf_seta; j += blockDim.y) {
      const int jco_start = task->first_cosetb + threadIdx.z;
      for (int jco = jco_start; jco < task->ncosetb; jco += blockDim.z) {
        const orbital b = coset_inv[jco];
        double block_val = 0.0;
        const double sphib = task->sphib[i * task->maxcob + jco];
        for (int ico = task->first_coseta; ico < task->ncoseta; ico++) {
          const orbital a = coset_inv[ico];
          const double hab =
              get_hab(a, b, task->zeta, task->zetb, task->n1, cab, COMPUTE_TAU);
          const double sphia = task->sphia[j * task->maxcoa + ico];
          block_val += hab * sphia * sphib;
        }
        if (task->block_transposed) {
          atomicAddDouble(&task->hab_block[j * task->nsgfb + i], block_val);
        } else {
          atomicAddDouble(&task->hab_block[i * task->nsgfa + j], block_val);
        }
      }
    }
  }
  __syncthreads(); // Not needed, but coalesced threads are nice.
}

/*******************************************************************************
 * \brief Adds contributions from cab to forces and virial.
 * \author Ole Schuett
 ******************************************************************************/
template <bool COMPUTE_TAU>
__device__ static void store_forces_and_virial(const kernel_params *params,
                                               const smem_task *task,
                                               const double *cab) {
  for (int i = threadIdx.x; i < task->nsgf_setb; i += blockDim.x) {
    for (int j = threadIdx.y; j < task->nsgf_seta; j += blockDim.y) {
      double block_val;
      if (task->block_transposed) {
        block_val = task->pab_block[j * task->nsgfb + i] * task->off_diag_twice;
      } else {
        block_val = task->pab_block[i * task->nsgfa + j] * task->off_diag_twice;
      }
      const int jco_start = task->first_cosetb + threadIdx.z;
      for (int jco = jco_start; jco < task->ncosetb; jco += blockDim.z) {
        const double sphib = task->sphib[i * task->maxcob + jco];
        for (int ico = task->first_coseta; ico < task->ncoseta; ico++) {
          const double sphia = task->sphia[j * task->maxcoa + ico];
          const double pabval = block_val * sphia * sphib;
          const orbital b = coset_inv[jco];
          const orbital a = coset_inv[ico];
          for (int k = 0; k < 3; k++) {
            const double force_a = get_force_a(a, b, k, task->zeta, task->zetb,
                                               task->n1, cab, COMPUTE_TAU);
            coalescedAtomicAdd(&task->forces_a[k], force_a * pabval);
            const double force_b =
                get_force_b(a, b, k, task->zeta, task->zetb, task->rab,
                            task->n1, cab, COMPUTE_TAU);
            coalescedAtomicAdd(&task->forces_b[k], force_b * pabval);
          }
          if (params->virial != NULL) {
            for (int k = 0; k < 3; k++) {
              for (int l = 0; l < 3; l++) {
                const double virial_a =
                    get_virial_a(a, b, k, l, task->zeta, task->zetb, task->n1,
                                 cab, COMPUTE_TAU);
                const double virial_b =
                    get_virial_b(a, b, k, l, task->zeta, task->zetb, task->rab,
                                 task->n1, cab, COMPUTE_TAU);
                const double virial = pabval * (virial_a + virial_b);
                coalescedAtomicAdd(&params->virial[k * 3 + l], virial);
              }
            }
          }
        }
      }
    }
  }
  __syncthreads(); // Not needed, but coalesced threads are nice.
}

/*******************************************************************************
 * \brief Initializes the cxyz matrix with zeros.
 * \author Ole Schuett
 ******************************************************************************/
__device__ static void zero_cxyz(const smem_task *task, double *cxyz) {
  if (threadIdx.z == 0 && threadIdx.y == 0) {
    for (int i = threadIdx.x; i < ncoset(task->lp); i += blockDim.x) {
      cxyz[i] = 0.0;
    }
  }
  __syncthreads(); // because of concurrent writes to cxyz
}

/*******************************************************************************
 * \brief Cuda kernel for integrating all tasks of one grid level.
 * \author Ole Schuett
 ******************************************************************************/
template <bool COMPUTE_TAU, bool CALCULATE_FORCES>
__device__ static void integrate_kernel(const kernel_params *params) {

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

  zero_cxyz(&task, smem_cxyz);
  cxyz_to_grid(params, &task, smem_cxyz, params->grid);

  zero_cab(&task, smem_cab);
  compute_alpha(params, &task, smem_alpha);
  cab_to_cxyz(params, &task, smem_alpha, smem_cab, smem_cxyz);

  store_hab<COMPUTE_TAU>(params, &task, smem_cab);

  if (CALCULATE_FORCES) {
    store_forces_and_virial<COMPUTE_TAU>(params, &task, smem_cab);
  }
}

/*******************************************************************************
 * \brief Specialized Cuda kernel for compute_tau=false & calculate_forces=false
 * \author Ole Schuett
 ******************************************************************************/
__global__ static void grid_integrate_density(const kernel_params params) {
  integrate_kernel<false, false>(&params);
}

/*******************************************************************************
 * \brief Specialized Cuda kernel for compute_tau=true & calculate_forces=false.
 * \author Ole Schuett
 ******************************************************************************/
__global__ static void grid_integrate_tau(const kernel_params params) {
  integrate_kernel<true, false>(&params);
}

/*******************************************************************************
 * \brief Specialized Cuda kernel for compute_tau=false & calculate_forces=true.
 * \author Ole Schuett
 ******************************************************************************/
__global__ static void
grid_integrate_density_forces(const kernel_params params) {
  integrate_kernel<false, true>(&params);
}

/*******************************************************************************
 * \brief Specialized Cuda kernel for compute_tau=true & calculate_forces=true.
 * \author Ole Schuett
 ******************************************************************************/
__global__ static void grid_integrate_tau_forces(const kernel_params params) {
  integrate_kernel<true, true>(&params);
}

/*******************************************************************************
 * \brief Launches the Cuda kernel that integrates all tasks of one grid level.
 * \author Ole Schuett
 ******************************************************************************/
void grid_gpu_integrate_one_grid_level(
    const grid_gpu_task_list *task_list, const int first_task,
    const int last_task, const bool orthorhombic, const bool compute_tau,
    const int npts_global[3], const int npts_local[3], const int shift_local[3],
    const int border_width[3], const double dh[3][3], const double dh_inv[3][3],
    const cudaStream_t stream, const double *pab_blocks_dev,
    const double *grid_dev, double *hab_blocks_dev, double *forces_dev,
    double *virial_dev, int *lp_diff) {

  // Compute max angular momentum.
  const bool calculate_forces = (forces_dev != NULL);
  const bool calculate_virial = (virial_dev != NULL);
  assert(!calculate_virial || calculate_forces);
  const process_ldiffs ldiffs =
      process_get_ldiffs(calculate_forces, calculate_virial, compute_tau);
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
    fprintf(stderr, "ERROR: Not enough shared memory in grid_gpu_integrate.\n");
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
  params.orthorhombic = orthorhombic;
  params.grid = grid_dev;
  params.tasks = task_list->tasks_dev;
  params.atom_kinds = task_list->atom_kinds_dev;
  params.basis_sets = task_list->basis_sets_dev;
  params.block_offsets = task_list->block_offsets_dev;
  params.atom_positions = task_list->atom_positions_dev;
  params.pab_blocks = pab_blocks_dev;
  params.hab_blocks = hab_blocks_dev;
  params.forces = forces_dev;
  params.virial = virial_dev;
  params.la_min_diff = ldiffs.la_min_diff;
  params.lb_min_diff = ldiffs.lb_min_diff;
  params.la_max_diff = ldiffs.la_max_diff;
  params.lb_max_diff = ldiffs.lb_max_diff;
  memcpy(params.dh, dh, 9 * sizeof(double));
  memcpy(params.dh_inv, dh_inv, 9 * sizeof(double));
  memcpy(params.npts_global, npts_global, 3 * sizeof(int));
  memcpy(params.npts_local, npts_local, 3 * sizeof(int));
  memcpy(params.shift_local, shift_local, 3 * sizeof(int));
  memcpy(params.border_width, border_width, 3 * sizeof(int));

  // Launch !
  const int nblocks = ntasks;
  const dim3 threads_per_block(4, 8, 8);

  if (!compute_tau && !calculate_forces) {
    grid_integrate_density<<<nblocks, threads_per_block, smem_per_block,
                             stream>>>(params);
  } else if (compute_tau && !calculate_forces) {
    grid_integrate_tau<<<nblocks, threads_per_block, smem_per_block, stream>>>(
        params);
  } else if (!compute_tau && calculate_forces) {
    grid_integrate_density_forces<<<nblocks, threads_per_block, smem_per_block,
                                    stream>>>(params);
  } else if (compute_tau && calculate_forces) {
    grid_integrate_tau_forces<<<nblocks, threads_per_block, smem_per_block,
                                stream>>>(params);
  }
}

#endif // __GRID_CUDA
// EOF
