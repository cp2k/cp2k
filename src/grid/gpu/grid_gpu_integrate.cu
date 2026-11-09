/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

/*
 * inspirations from the gpu backend
 * Authors :
 - Mathieu Taillefumier (ETH Zurich / CSCS)
 - Advanced Micro Devices, Inc.
 - Ole Schuett
*/

#include <algorithm>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "grid_gpu_context.h"
#include "grid_gpu_internal_header.h"
#include "grid_gpu_process_vab.h"
#include "grid_gpu_separable.h"

#if defined(_OMP_H)
#error "OpenMP should not be used in .cu files to accommodate HIP."
#endif

namespace rocm_backend {

// tasks with lp >= this value go to integrate_separable_kernel, the others to
// integrate_kernel. The standalone driver shows the separable kernel 20-37%
// faster than integrate_kernel for lp = 0..3 as well (V100), so all tasks go
// to it and integrate_kernel is only launched as a fallback.
constexpr int integrate_separable_lp_threshold = 0;

__inline__ __device__ double sep_row_reduce(double val) {
  for (int offset = CP2K_SEP_INTEGRATE_ROW_LANES / 2; offset > 0;
       offset >>= 1) {
#if defined(__CUDACC__)
    val +=
        __shfl_down_sync(0xffffffff, val, offset, CP2K_SEP_INTEGRATE_ROW_LANES);
#else
    val += __shfl_down(val, offset, CP2K_SEP_INTEGRATE_ROW_LANES);
#endif
  }
  return val;
}

// do a block reduction for a block of size 64 and return the final sum to
// thread_id = 0
template <typename T>
__device__ __inline__ T block_reduce_64(T *table, const T val_, const int tid) {
  // AMD GPU have warp size of 64 while nvidia GPUs have warpSize of 32 so the
  // first step is common to both platforms.
  table[tid] = val_;
  __syncthreads();

  T val = 0.0;

  if (tid < 32) {
    val = table[tid] + table[tid + 32];

    for (int offset = 16; offset > 0; offset >>= 1) {
#if defined(__CUDACC__)
      val += __shfl_down_sync(0xffffffff, val, offset);
#else
      val += __shfl_down(val, offset);
#endif
    }
  }

  // prevents threads >= 32 (which never touch table[]) from racing ahead and
  // overwriting it for the next call before the reduction above has read it.
  __syncthreads();
  return val;
}

template <typename T, bool COMPUTE_TAU>
__global__ __launch_bounds__(64) void compute_hab(const kernel_params dev_) {
  // Copy task from global to shared memory and precompute some stuff.
  extern __shared__ T shared_memory[];
  // T *smem_cab = &shared_memory[dev_.smem_cab_offset];
  const int number_of_tasks = dev_.num_tasks_per_block_dev[block_index()];

  if (number_of_tasks == 0)
    return;

  T *smem_alpha = &shared_memory[0];
  const int tid = thread_global_index();
  const int offset = dev_.sorted_blocks_offset_dev[block_index()];

  T *__restrict__ smem_cab = nullptr;

  smem_cab = allocate_workspace<T>(dev_);

  for (int tk = 0; tk < number_of_tasks; tk++) {
    __shared__ smem_task<T> task;
    const int task_id = dev_.task_sorted_by_blocks_dev[offset + tk];
    if (dev_.tasks[task_id].skip_task)
      continue;

    // all warps need to be synchronized here before modifying the task
    // information.
    __syncthreads();
    fill_smem_task_coef(dev_, task_id, task);

    T *__restrict__ coef_ = reinterpret_cast<T *>(__builtin_assume_aligned(
        &dev_.buffers_dev.coef[dev_.tasks[task_id].coef_offset], 32));

    __syncthreads();
    compute_alpha(task, smem_alpha);
    __syncthreads();
    cxyz_to_cab(task, smem_alpha, coef_, smem_cab);
    __syncthreads();

    for (int i = tid / 8; i < task.nsgf_setb; i += 8) {
      for (int j = tid % 8; j < task.nsgf_seta; j += 8) {
        T tmp = 0.0;
        for (int jco = task.first_cosetb; jco < task.ncosetb; jco++) {
          const T sphib = task.sphib[i * task.maxcob + jco];
          const auto &b = coset_inv[jco];
          for (int ico = task.first_coseta; ico < task.ncoseta; ico++) {
            const auto &a = coset_inv[ico];
            const T hab = get_hab<COMPUTE_TAU, T>(a, b, task.zeta, task.zetb,
                                                  task.n1, smem_cab);
            const T sphia_times_sphib =
                task.sphia[j * task.maxcoa + ico] * sphib;
            tmp += hab * sphia_times_sphib;
          }
        }
        if (task.block_transposed) {
          task.hab_block[j * task.nsgfb + i] += tmp;
        } else {
          task.hab_block[i * task.nsgfa + j] += tmp;
        }
      }
    }
  }
}

/*******************************************************************************
 * \brief Same as compute_hab, with the cartesian to spherical transformation
 *        done in two steps,
 *
 *   H[jco][ico] = get_hab(a, b, cab)                     (ncob ncoa)
 *   U[i][ico]   = sum_jco sphib[i][jco] H[jco][ico]      (nsgfb ncob ncoa)
 *   hab[i][j]  += sum_ico U[i][ico] sphia[j][ico]        (nsgfb nsgfa ncoa)
 *
 *        instead of evaluating get_hab and both sphi factors in the innermost
 *        loop of a nsgfb nsgfa ncob ncoa loop. The work arrays H [ld][ld] and
 *        U [max_nsgf_set][ld] live in the dynamic shared memory, after alpha,
 *        at work_offset (in doubles), see sphi_work_size.
 ******************************************************************************/
template <typename T, bool COMPUTE_TAU>
__global__ __launch_bounds__(64) void compute_hab_two_steps(
    const kernel_params dev_, const int work_offset, const int ld,
    const int max_nsgf_set) {
  extern __shared__ T shared_memory[];
  const int number_of_tasks = dev_.num_tasks_per_block_dev[block_index()];

  if (number_of_tasks == 0)
    return;

  T *smem_alpha = &shared_memory[0];
  T *smem_h = &shared_memory[work_offset]; // [ld][ld], H[jco][ico]
  T *smem_u = smem_h + ld * ld;            // [max_nsgf_set][ld]
  const int tid = thread_global_index();
  const int nthreads = blockDim.x * blockDim.y * blockDim.z;
  const int offset = dev_.sorted_blocks_offset_dev[block_index()];

  T *__restrict__ smem_cab = allocate_workspace<T>(dev_);

  for (int tk = 0; tk < number_of_tasks; tk++) {
    __shared__ smem_task<T> task;
    const int task_id = dev_.task_sorted_by_blocks_dev[offset + tk];
    if (dev_.tasks[task_id].skip_task)
      continue;

    __syncthreads();
    fill_smem_task_coef(dev_, task_id, task);

    T *__restrict__ coef_ = reinterpret_cast<T *>(__builtin_assume_aligned(
        &dev_.buffers_dev.coef[dev_.tasks[task_id].coef_offset], 32));

    __syncthreads();
    cxyz_to_cab_factorized(task, coef_, smem_cab);
    __syncthreads();

    const int nco_a = task.ncoseta - task.first_coseta;
    const int nco_b = task.ncosetb - task.first_cosetb;

    // H[jco][ico] = get_hab(a, b)
    for (int k = tid; k < nco_b * nco_a; k += nthreads) {
      const int jj = k / nco_a;
      const int jco = task.first_cosetb + jj;
      const int ico = task.first_coseta + (k - jj * nco_a);
      smem_h[jco * ld + ico] =
          get_hab<COMPUTE_TAU, T>(coset_inv[ico], coset_inv[jco], task.zeta,
                                  task.zetb, task.n1, smem_cab);
    }
    __syncthreads();

    // U[i][ico] = sum_jco sphib[i][jco] H[jco][ico]
    for (int k = tid; k < task.nsgf_setb * nco_a; k += nthreads) {
      const int i = k / nco_a;
      const int ico = task.first_coseta + (k - i * nco_a);
      T val = 0.0;
      for (int jco = task.first_cosetb; jco < task.ncosetb; jco++)
        val += task.sphib[i * task.maxcob + jco] * smem_h[jco * ld + ico];
      smem_u[i * ld + ico] = val;
    }
    __syncthreads();

    // hab[i][j] += sum_ico U[i][ico] sphia[j][ico]
    for (int k = tid; k < task.nsgf_setb * task.nsgf_seta; k += nthreads) {
      const int i = k / task.nsgf_seta;
      const int j = k - i * task.nsgf_seta;
      T val = 0.0;
      for (int ico = task.first_coseta; ico < task.ncoseta; ico++)
        val += smem_u[i * ld + ico] * task.sphia[j * task.maxcoa + ico];
      if (task.block_transposed) {
        task.hab_block[j * task.nsgfb + i] += val;
      } else {
        task.hab_block[i * task.nsgfa + j] += val;
      }
    }
  }
}

template <typename T, typename T3, bool COMPUTE_TAU>
__global__
__launch_bounds__(64) void compute_hab_forces(const kernel_params dev_) {
  // Copy task from global to shared memory and precompute some stuff.
  extern __shared__ T shared_memory[];
  const int number_of_tasks = dev_.num_tasks_per_block_dev[block_index()];

  if (number_of_tasks == 0)
    return;

  T *smem_alpha = &shared_memory[0];
  const int tid = thread_global_index();
  const int offset = dev_.sorted_blocks_offset_dev[block_index()];

  T *__restrict__ smem_cab = allocate_workspace<T>(dev_);

  T fa[3], fb[3];
  T virial[9];

  fa[0] = 0.0;
  fa[1] = 0.0;
  fa[2] = 0.0;
  fb[0] = 0.0;
  fb[1] = 0.0;
  fb[2] = 0.0;

  virial[0] = 0.0;
  virial[1] = 0.0;
  virial[2] = 0.0;
  virial[3] = 0.0;
  virial[4] = 0.0;
  virial[5] = 0.0;
  virial[6] = 0.0;
  virial[7] = 0.0;
  virial[8] = 0.0;

  for (int tk = 0; tk < number_of_tasks; tk++) {
    __shared__ smem_task<T> task;
    const int task_id = dev_.task_sorted_by_blocks_dev[offset + tk];
    if (dev_.tasks[task_id].skip_task)
      continue;
    fill_smem_task_coef(dev_, task_id, task);

    T *__restrict__ coef_ =
        &dev_.buffers_dev.coef[dev_.tasks[task_id].coef_offset];
    __syncthreads();
    compute_alpha(task, smem_alpha);
    __syncthreads();
    cxyz_to_cab(task, smem_alpha, coef_, smem_cab);
    __syncthreads();

    for (int i = tid / 8; i < task.nsgf_setb; i += 8) {
      for (int j = tid % 8; j < task.nsgf_seta; j += 8) {
        T tmp = 0.0;
        T block_value = 0.0;
        if (task.block_transposed) {
          block_value =
              task.pab_block[j * task.nsgfb + i] * task.off_diag_twice;
        } else {
          block_value =
              task.pab_block[i * task.nsgfa + j] * task.off_diag_twice;
        }
        for (int jco = task.first_cosetb; jco < task.ncosetb; jco++) {
          const T sphib = task.sphib[i * task.maxcob + jco];
          const auto &b = coset_inv[jco];
          for (int ico = task.first_coseta; ico < task.ncoseta; ico++) {
            const auto &a = coset_inv[ico];
            const T hab = get_hab<COMPUTE_TAU, T>(a, b, task.zeta, task.zetb,
                                                  task.n1, smem_cab);
            T sphia_times_sphib = task.sphia[j * task.maxcoa + ico] * sphib;
            tmp += hab * sphia_times_sphib;

            sphia_times_sphib *= block_value;
            fa[0] += sphia_times_sphib *
                     get_force_a<COMPUTE_TAU, T>(a, b, 0, task.zeta, task.zetb,
                                                 task.n1, smem_cab);
            fa[1] += sphia_times_sphib *
                     get_force_a<COMPUTE_TAU, T>(a, b, 1, task.zeta, task.zetb,
                                                 task.n1, smem_cab);
            fa[2] += sphia_times_sphib *
                     get_force_a<COMPUTE_TAU, T>(a, b, 2, task.zeta, task.zetb,
                                                 task.n1, smem_cab);

            fb[0] += sphia_times_sphib *
                     get_force_b<COMPUTE_TAU, T>(a, b, 0, task.zeta, task.zetb,
                                                 task.rab, task.n1, smem_cab);
            fb[1] += sphia_times_sphib *
                     get_force_b<COMPUTE_TAU, T>(a, b, 1, task.zeta, task.zetb,
                                                 task.rab, task.n1, smem_cab);
            fb[2] += sphia_times_sphib *
                     get_force_b<COMPUTE_TAU, T>(a, b, 2, task.zeta, task.zetb,
                                                 task.rab, task.n1, smem_cab);

            if (dev_.buffers_dev.virial != nullptr) {
              virial[0] +=
                  sphia_times_sphib *
                  (get_virial_a<COMPUTE_TAU, T>(a, b, 0, 0, task.zeta,
                                                task.zetb, task.n1, smem_cab) +
                   get_virial_b<COMPUTE_TAU, T>(a, b, 0, 0, task.zeta,
                                                task.zetb, task.rab, task.n1,
                                                smem_cab));
              virial[1] +=
                  sphia_times_sphib *
                  (get_virial_a<COMPUTE_TAU, T>(a, b, 0, 1, task.zeta,
                                                task.zetb, task.n1, smem_cab) +
                   get_virial_b<COMPUTE_TAU, T>(a, b, 0, 1, task.zeta,
                                                task.zetb, task.rab, task.n1,
                                                smem_cab));
              virial[2] +=
                  sphia_times_sphib *
                  (get_virial_a<COMPUTE_TAU, T>(a, b, 0, 2, task.zeta,
                                                task.zetb, task.n1, smem_cab) +
                   get_virial_b<COMPUTE_TAU, T>(a, b, 0, 2, task.zeta,
                                                task.zetb, task.rab, task.n1,
                                                smem_cab));
              virial[3] +=
                  sphia_times_sphib *
                  (get_virial_a<COMPUTE_TAU, T>(a, b, 1, 0, task.zeta,
                                                task.zetb, task.n1, smem_cab) +
                   get_virial_b<COMPUTE_TAU, T>(a, b, 1, 0, task.zeta,
                                                task.zetb, task.rab, task.n1,
                                                smem_cab));
              virial[4] +=
                  sphia_times_sphib *
                  (get_virial_a<COMPUTE_TAU, T>(a, b, 1, 1, task.zeta,
                                                task.zetb, task.n1, smem_cab) +
                   get_virial_b<COMPUTE_TAU, T>(a, b, 1, 1, task.zeta,
                                                task.zetb, task.rab, task.n1,
                                                smem_cab));
              virial[5] +=
                  sphia_times_sphib *
                  (get_virial_a<COMPUTE_TAU, T>(a, b, 1, 2, task.zeta,
                                                task.zetb, task.n1, smem_cab) +
                   get_virial_b<COMPUTE_TAU, T>(a, b, 1, 2, task.zeta,
                                                task.zetb, task.rab, task.n1,
                                                smem_cab));
              virial[6] +=
                  sphia_times_sphib *
                  (get_virial_a<COMPUTE_TAU, T>(a, b, 2, 0, task.zeta,
                                                task.zetb, task.n1, smem_cab) +
                   get_virial_b<COMPUTE_TAU, T>(a, b, 2, 0, task.zeta,
                                                task.zetb, task.rab, task.n1,
                                                smem_cab));
              virial[7] +=
                  sphia_times_sphib *
                  (get_virial_a<COMPUTE_TAU, T>(a, b, 2, 1, task.zeta,
                                                task.zetb, task.n1, smem_cab) +
                   get_virial_b<COMPUTE_TAU, T>(a, b, 2, 1, task.zeta,
                                                task.zetb, task.rab, task.n1,
                                                smem_cab));
              virial[8] +=
                  sphia_times_sphib *
                  (get_virial_a<COMPUTE_TAU, T>(a, b, 2, 2, task.zeta,
                                                task.zetb, task.n1, smem_cab) +
                   get_virial_b<COMPUTE_TAU, T>(a, b, 2, 2, task.zeta,
                                                task.zetb, task.rab, task.n1,
                                                smem_cab));
            }
          }
        }

        if (task.block_transposed) {
          task.hab_block[j * task.nsgfb + i] += tmp;
        } else {
          task.hab_block[i * task.nsgfa + j] += tmp;
        }
      }
    }
    __syncthreads();
  }

  // theoretically not needed
  __syncthreads();

  const int task_id = dev_.task_sorted_by_blocks_dev[offset];
  const auto &glb_task = dev_.tasks[task_id];
  const int iatom = glb_task.iatom;
  const int jatom = glb_task.jatom;
  T *forces_a = &dev_.buffers_dev.forces[3 * iatom];
  T *forces_b = &dev_.buffers_dev.forces[3 * jatom];

  T *sum = (T *)shared_memory;
  if (dev_.buffers_dev.virial != nullptr) {

    for (int i = 0; i < 9; i++) {
      virial[i] = block_reduce_64<T>(sum, virial[i], tid);

      if (tid == 0)
        atomicAdd(dev_.buffers_dev.virial + i, virial[i]);
    }
  }

  for (int i = 0; i < 3; i++) {
    fa[i] = block_reduce_64<T>(sum, fa[i], tid);

    if (tid == 0)
      atomicAdd(forces_a + i, fa[i]);

    fb[i] = block_reduce_64<T>(sum, fb[i], tid);

    if (tid == 0)
      atomicAdd(forces_b + i, fb[i]);
  }
}

/*******************************************************************************
 * \brief Same as compute_hab_forces with the spherical to cartesian
 *        transformations done as matrix products. The forces and the virial
 *        only depend on the cartesian pair (a, b), so they are contracted with
 *        the cartesian density
 *
 *   P[jco][ico] = sum_ij sphib[i][jco] pab[i][j] sphia[j][ico]
 *
 *        and evaluated once per pair instead of once per (i, j, a, b):
 *
 *   W[i][ico]   = sum_j pab[i][j] sphia[j][ico]              (nsgfb nsgfa ncoa)
 *   for each (jco, ico):
 *     p          = sum_i sphib[i][jco] W[i][ico]             (nsgfb ncob ncoa)
 *     fa, fb, virial += p * get_force / get_virial (a, b)
 *     H[jco][ico] = get_hab(a, b)
 *   U[i][ico]   = sum_jco sphib[i][jco] H[jco][ico]          (nsgfb ncob ncoa)
 *   hab[i][j]  += sum_ico U[i][ico] sphia[j][ico]            (nsgfb nsgfa ncoa)
 *
 *        P is never stored. The work arrays H [ld][ld] and W / U
 *        [max_nsgf_set][ld] (U overwrites W) live in the dynamic shared memory
 *        after alpha, at work_offset (in doubles), see sphi_work_size.
 ******************************************************************************/
template <typename T, bool COMPUTE_TAU>
__global__ __launch_bounds__(64) void compute_hab_forces_two_steps(
    const kernel_params dev_, const int work_offset, const int ld,
    const int max_nsgf_set) {
  extern __shared__ T shared_memory[];
  const int number_of_tasks = dev_.num_tasks_per_block_dev[block_index()];

  if (number_of_tasks == 0)
    return;

  T *smem_alpha = &shared_memory[0];
  T *smem_h = &shared_memory[work_offset]; // [ld][ld], H[jco][ico]
  T *smem_w = smem_h + ld * ld;            // [max_nsgf_set][ld], W then U
  const int tid = thread_global_index();
  const int nthreads = blockDim.x * blockDim.y * blockDim.z;
  const int offset = dev_.sorted_blocks_offset_dev[block_index()];
  const bool compute_virial = (dev_.buffers_dev.virial != nullptr);

  T *__restrict__ smem_cab = allocate_workspace<T>(dev_);

  T fa[3] = {0.0, 0.0, 0.0}, fb[3] = {0.0, 0.0, 0.0};
  T virial[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  for (int tk = 0; tk < number_of_tasks; tk++) {
    __shared__ smem_task<T> task;
    const int task_id = dev_.task_sorted_by_blocks_dev[offset + tk];
    if (dev_.tasks[task_id].skip_task)
      continue;

    // the previous task may still be reading task, W / U and H
    __syncthreads();
    fill_smem_task_coef(dev_, task_id, task);

    T *__restrict__ coef_ = reinterpret_cast<T *>(__builtin_assume_aligned(
        &dev_.buffers_dev.coef[dev_.tasks[task_id].coef_offset], 32));

    __syncthreads();
    cxyz_to_cab_factorized(task, coef_, smem_cab);

    const int nco_a = task.ncoseta - task.first_coseta;
    const int nco_b = task.ncosetb - task.first_cosetb;

    // W[i][ico] = sum_j pab[i][j] sphia[j][ico]
    for (int k = tid; k < task.nsgf_setb * nco_a; k += nthreads) {
      const int i = k / nco_a;
      const int ico = task.first_coseta + (k - i * nco_a);
      T val = 0.0;
      for (int j = 0; j < task.nsgf_seta; j++) {
        const T block_value = task.block_transposed
                                  ? task.pab_block[j * task.nsgfb + i]
                                  : task.pab_block[i * task.nsgfa + j];
        val += block_value * task.sphia[j * task.maxcoa + ico];
      }
      smem_w[i * ld + ico] = val * task.off_diag_twice;
    }
    __syncthreads(); // W and cab are complete

    // one pass over the cartesian pairs: density, forces, virial and H
    for (int k = tid; k < nco_b * nco_a; k += nthreads) {
      const int jj = k / nco_a;
      const int jco = task.first_cosetb + jj;
      const int ico = task.first_coseta + (k - jj * nco_a);
      const auto &a = coset_inv[ico];
      const auto &b = coset_inv[jco];

      T p = 0.0;
      for (int i = 0; i < task.nsgf_setb; i++)
        p += task.sphib[i * task.maxcob + jco] * smem_w[i * ld + ico];

      for (int d = 0; d < 3; d++) {
        fa[d] += p * get_force_a<COMPUTE_TAU, T>(a, b, d, task.zeta, task.zetb,
                                                 task.n1, smem_cab);
        fb[d] += p * get_force_b<COMPUTE_TAU, T>(a, b, d, task.zeta, task.zetb,
                                                 task.rab, task.n1, smem_cab);
      }
      if (compute_virial) {
        for (int d1 = 0; d1 < 3; d1++)
          for (int d2 = 0; d2 < 3; d2++)
            virial[3 * d1 + d2] +=
                p *
                (get_virial_a<COMPUTE_TAU, T>(a, b, d1, d2, task.zeta,
                                              task.zetb, task.n1, smem_cab) +
                 get_virial_b<COMPUTE_TAU, T>(a, b, d1, d2, task.zeta,
                                              task.zetb, task.rab, task.n1,
                                              smem_cab));
      }

      smem_h[jco * ld + ico] = get_hab<COMPUTE_TAU, T>(
          a, b, task.zeta, task.zetb, task.n1, smem_cab);
    }
    __syncthreads(); // H is complete, W is no longer read

    // U[i][ico] = sum_jco sphib[i][jco] H[jco][ico], overwrites W
    for (int k = tid; k < task.nsgf_setb * nco_a; k += nthreads) {
      const int i = k / nco_a;
      const int ico = task.first_coseta + (k - i * nco_a);
      T val = 0.0;
      for (int jco = task.first_cosetb; jco < task.ncosetb; jco++)
        val += task.sphib[i * task.maxcob + jco] * smem_h[jco * ld + ico];
      smem_w[i * ld + ico] = val;
    }
    __syncthreads();

    // hab[i][j] += sum_ico U[i][ico] sphia[j][ico]
    for (int k = tid; k < task.nsgf_setb * task.nsgf_seta; k += nthreads) {
      const int i = k / task.nsgf_seta;
      const int j = k - i * task.nsgf_seta;
      T val = 0.0;
      for (int ico = task.first_coseta; ico < task.ncoseta; ico++)
        val += smem_w[i * ld + ico] * task.sphia[j * task.maxcoa + ico];
      if (task.block_transposed) {
        task.hab_block[j * task.nsgfb + i] += val;
      } else {
        task.hab_block[i * task.nsgfa + j] += val;
      }
    }
  }

  // all tasks of a block belong to the same atom pair
  __syncthreads();

  const int task_id = dev_.task_sorted_by_blocks_dev[offset];
  const auto &glb_task = dev_.tasks[task_id];
  T *forces_a = &dev_.buffers_dev.forces[3 * glb_task.iatom];
  T *forces_b = &dev_.buffers_dev.forces[3 * glb_task.jatom];

  T *sum = (T *)shared_memory;
  if (compute_virial) {
    for (int i = 0; i < 9; i++) {
      virial[i] = block_reduce_64<T>(sum, virial[i], tid);
      if (tid == 0)
        atomicAdd(dev_.buffers_dev.virial + i, virial[i]);
    }
  }

  for (int i = 0; i < 3; i++) {
    fa[i] = block_reduce_64<T>(sum, fa[i], tid);
    if (tid == 0)
      atomicAdd(forces_a + i, fa[i]);
    fb[i] = block_reduce_64<T>(sum, fb[i], tid);
    if (tid == 0)
      atomicAdd(forces_b + i, fb[i]);
  }
}

/*******************************************************************************
 * Cuda kernel for calculating the coefficients of a potential (density,
 etc...) for a given pair of gaussian basis sets

We compute the discretized version of the following integral

$$\int _\infty ^\infty V_{ijk} P^\alpha_iP^beta_j P ^ \gamma _ k Exp(- \eta
|r_{ijk} - r_c|^2)$$

where in practice the summation is over a finite domain. The discrete form has
this shape

$$
\sum_{ijk < domain} V_{ijk} (x_i - x_c)^\alpha (y_j - y_c)^\beta (z_k - z_c) ^
\gamma Exp(- \eta |r_{ijk} - r_c|^2)
$$

where $0 \le \alpha + \beta + \gamma \le lmax$

It is formely the same operation than collocate (from a discrete point of view)
but the implementation differ because of technical reasons.

So most of the code is the same except the core of the routine that is
specialized to the integration.

******************************************************************************/
template <typename T, typename T3, bool distributed__, bool orthogonal_>
__global__
__launch_bounds__(64) void integrate_kernel(const kernel_params dev_) {
  const int lbatch = 20;
  const task_info &glb_task = dev_.tasks[dev_.first_task + block_index()];
  if (glb_task.skip_task)
    return;

  // handled by integrate_separable_kernel
  if (glb_task.la_max + dev_.la_max_diff + glb_task.lb_max + dev_.lb_max_diff >=
      dev_.separable_lp_min)
    return;

  const int tid = thread_global_index();

  __shared__ T dh_[9];

  if (tid < 9)
    dh_[tid] = dev_.dh_[tid];

  __shared__ smem_task_reduced<T, T3> task;
  fill_smem_task_reduced(dev_, dev_.first_task + block_index(), task);

  if (tid == 0) {
    setup_task_cube_center<T, T3, distributed__>(dev_, task);
  }
  __syncthreads();

  __shared__ T accumulator[20][64];

  // we use a multi pass algorithm here because shared memory usage (or
  // register) would become too high for high angular momentum

  const short int size_loop =
      (ncoset(task.lp) / lbatch + ((ncoset(task.lp) % lbatch) != 0)) * lbatch;
  const short int length = ncoset(task.lp);
  for (int ico = 0; ico < size_loop; ico += lbatch) {
#pragma unroll lbatch
    for (int i = 0; i < lbatch; i++)
      accumulator[i][tid] = 0.0;

    __syncthreads();

    for (int z = threadIdx.z; z < task.cube_size.z; z += blockDim.z) {
      int z2 = wrap_grid_index(z + task.cube_center.z, dev_.grid_full_size_.z);

      if (distributed__) {
        // known at compile time. Will be stripped away
        if (task.apply_border_mask) {
          /* check if the point is within the window */
          if ((z2 < task.window_shift.z) || (z2 > task.window_size.z)) {
            continue;
          }
        }
      }

      /* compute the coordinates of the point in atomic coordinates */
      int ymin = 0;
      int ymax = task.cube_size.y - 1;
      T kremain = 0.0;

      if (orthogonal_ && !task.apply_border_mask) {
        kremain = calculate_ymix_ymax_boundaries(task, z, ymin, ymax);
      }

      for (int y = ymin + threadIdx.y; y <= ymax; y += blockDim.y) {
        int y2 =
            wrap_grid_index(y + task.cube_center.y, dev_.grid_full_size_.y);

        if (distributed__) {
          /* check if the point is within the window */
          if (task.apply_border_mask) {
            if ((y2 < task.window_shift.y) || (y2 > task.window_size.y)) {
              continue;
            }
          }
        }

        int xmin = 0;
        int xmax = task.cube_size.x - 1;

        if (orthogonal_ && !task.apply_border_mask) {
          calculate_xmin_xmax_boundaries<T, T3>(task, y, kremain, xmin, xmax);
        }

        for (int x = xmin + threadIdx.x; x <= xmax; x += blockDim.x) {
          int x2 =
              wrap_grid_index(x + task.cube_center.x, dev_.grid_full_size_.x);

          if (distributed__) {
            /* check if the point is within the window */
            if (task.apply_border_mask) {
              if ((x2 < task.window_shift.x) || (x2 > task.window_size.x)) {
                continue;
              }
            }
          }

          // I make no distinction between orthogonal and non orthogonal
          // cases
          T3 r3;

          r3 = compute_coordinates(dh_, (x + task.lb_cube.x + task.roffset.x),
                                   (y + task.lb_cube.y + task.roffset.y),
                                   (z + task.lb_cube.z + task.roffset.z));
          // check if the point is inside the sphere or not. Note that it does
          // not apply for the orthogonal case when the full sphere is inside
          // the region of interest.
          const T r3x2 = r3.x * r3.x;
          const T r3y2 = r3.y * r3.y;
          const T r3z2 = r3.z * r3.z;

          if (distributed__) {
            // check if the point is inside the sphere or not. Note that it does
            // not apply for the orthorhombic case when the full sphere is
            // inside the region of interest.

            if (((task.radius * task.radius) <= (r3x2 + r3y2 + r3z2)) &&
                (!orthogonal_ || task.apply_border_mask))
              continue;
          } else {
            // we do not need to do this test for the orthorhombic case
            if ((!orthogonal_) &&
                ((task.radius * task.radius) <= (r3x2 + r3y2 + r3z2)))
              continue;
          }

          // read the next point assuming that reading is non blocking until
          // the register is actually needed for computation. This is true on
          // NVIDIA hardware

          const int grid_index =
              (z2 * dev_.grid_local_size_.y + y2) * dev_.grid_local_size_.x +
              x2;
          T grid_value = __ldg(&dev_.buffers_dev.grid[grid_index]);

          const T r3xy = r3.x * r3.y;
          const T r3xz = r3.x * r3.z;
          const T r3yz = r3.y * r3.z;

          grid_value *= exp(-(r3x2 + r3y2 + r3z2) * task.zetp);

          switch (ico / lbatch) {
          case 0: {
            accumulator[0][tid] += grid_value;

            if (task.lp >= 1) {
              accumulator[1][tid] += grid_value * r3.x;
              accumulator[2][tid] += grid_value * r3.y;
              accumulator[3][tid] += grid_value * r3.z;
            }

            if (task.lp >= 2) {
              accumulator[4][tid] += grid_value * r3x2;
              accumulator[5][tid] += grid_value * r3xy;
              accumulator[6][tid] += grid_value * r3xz;
              accumulator[7][tid] += grid_value * r3y2;
              accumulator[8][tid] += grid_value * r3yz;
              accumulator[9][tid] += grid_value * r3z2;
            }
            if (task.lp >= 3) {
              T tmp = grid_value * r3x2;
              accumulator[10][tid] += tmp * r3.x;
              accumulator[11][tid] += tmp * r3.y;
              accumulator[12][tid] += tmp * r3.z;
              tmp = grid_value * r3.x;
              accumulator[13][tid] += tmp * r3y2;
              accumulator[14][tid] += tmp * r3yz;
              accumulator[15][tid] += tmp * r3z2;
              tmp = grid_value * r3y2;
              accumulator[16][tid] += tmp * r3.y;
              accumulator[17][tid] += tmp * r3.z;
              tmp = grid_value * r3z2;
              accumulator[18][tid] += tmp * r3.y;
              accumulator[19][tid] += tmp * r3.z;
            }
          } break;
          case 1: {
            if (task.lp >= 4) {
              T tmp = grid_value * r3x2;
              accumulator[0][tid] += tmp * r3x2;
              accumulator[1][tid] += tmp * r3xy;
              accumulator[2][tid] += tmp * r3xz;
              accumulator[3][tid] += tmp * r3y2;
              accumulator[4][tid] += tmp * r3yz;
              accumulator[5][tid] += tmp * r3z2;
              tmp = grid_value * r3y2;
              accumulator[6][tid] += tmp * r3xy;
              accumulator[7][tid] += tmp * r3xz;
              tmp = grid_value * r3z2;
              accumulator[8][tid] += tmp * r3xy;
              accumulator[9][tid] += tmp * r3xz;
            }
            if (task.lp >= 4) {
              T tmp = grid_value * r3y2;
              accumulator[10][tid] += tmp * r3y2;
              accumulator[11][tid] += tmp * r3yz;
              accumulator[12][tid] += tmp * r3z2;
              accumulator[13][tid] += grid_value * r3yz * r3z2;
              accumulator[14][tid] += grid_value * r3z2 * r3z2;
            }
            if (task.lp >= 5) {
              T tmp = grid_value * r3x2 * r3x2;
              accumulator[15][tid] += tmp * r3.x;
              accumulator[16][tid] += tmp * r3.y;
              accumulator[17][tid] += tmp * r3.z;
              tmp = grid_value * r3x2;
              accumulator[18][tid] += tmp * r3.x * r3y2;
              accumulator[19][tid] += tmp * r3xy * r3.z;
            }
          } break;
          case 2: {
            T tmp = grid_value * r3x2;
            accumulator[0][tid] += tmp * r3.x * r3z2;
            accumulator[1][tid] += tmp * r3y2 * r3.y;
            accumulator[2][tid] += tmp * r3y2 * r3.z;
            accumulator[3][tid] += tmp * r3.y * r3z2;
            accumulator[4][tid] += tmp * r3z2 * r3.z;
            tmp = grid_value * r3.x * r3y2;
            accumulator[5][tid] += tmp * r3y2;
            accumulator[6][tid] += tmp * r3yz;
            accumulator[7][tid] += tmp * r3z2;
            tmp = grid_value * r3.x * r3z2;
            accumulator[8][tid] += tmp * r3yz;
            accumulator[9][tid] += tmp * r3z2;
            tmp = grid_value * r3y2 * r3.y;
            accumulator[10][tid] += tmp * r3y2;
            accumulator[11][tid] += tmp * r3yz;
            accumulator[12][tid] += tmp * r3z2;
            accumulator[13][tid] += grid_value * r3y2 * r3z2 * r3.z;
            accumulator[14][tid] += grid_value * r3.y * r3z2 * r3z2;
            accumulator[15][tid] += grid_value * r3z2 * r3z2 * r3.z;
            if (task.lp >= 6) {
              tmp = grid_value * r3x2 * r3x2;
              accumulator[16][tid] += tmp * r3x2; // x^6
              accumulator[17][tid] += tmp * r3xy; // x^5 y
              accumulator[18][tid] += tmp * r3xz; // x^5 z
              accumulator[19][tid] += tmp * r3y2; // x^4 y^2
            }
          } break;
          case 3: {
            T tmp = grid_value * r3x2;
            accumulator[0][tid] += tmp * r3x2 * r3yz;  // x^4 y z
            accumulator[1][tid] += tmp * r3x2 * r3z2;  // x^4 z^2
            accumulator[2][tid] += tmp * r3y2 * r3xy;  // x^3 y^3
            accumulator[3][tid] += tmp * r3y2 * r3xz;  // x^3 y^2 z
            accumulator[4][tid] += tmp * r3xy * r3z2;  // x^3 y z^2
            accumulator[5][tid] += tmp * r3z2 * r3xz;  // x^3 z^3
            accumulator[6][tid] += tmp * r3y2 * r3y2;  // x^2 y^4
            accumulator[7][tid] += tmp * r3y2 * r3yz;  // x^3 y^2 z
            accumulator[8][tid] += tmp * r3y2 * r3z2;  // x^2 y^2 z^2
            accumulator[9][tid] += tmp * r3z2 * r3yz;  // x^2 y z^3
            accumulator[10][tid] += tmp * r3z2 * r3z2; // x^2 z^4
            tmp = grid_value * r3y2 * r3y2;
            accumulator[11][tid] += tmp * r3xy; // x y^5
            accumulator[12][tid] += tmp * r3xz; // x y^4 z
            accumulator[13][tid] +=
                grid_value * r3y2 * r3xy * r3z2; // x y^3 z^2
            accumulator[14][tid] +=
                grid_value * r3y2 * r3z2 * r3xz; // x y^2 z^3
            accumulator[15][tid] += grid_value * r3xy * r3z2 * r3z2; // x y z^4
            accumulator[16][tid] += grid_value * r3z2 * r3z2 * r3xz; // x z^5
            accumulator[17][tid] += tmp * r3y2;                      // y^6
            accumulator[18][tid] += tmp * r3yz;                      // y^5 z
            accumulator[19][tid] += tmp * r3z2;                      // y^4 z^2
          } break;
          default:
            for (int ic = 0; (ic < lbatch) && ((ic + ico) < length); ic++) {
              auto &co = coset_inv[ic + ico];
              T tmp = 1.0;
              for (int po = 0; po < (co.l[2] >> 1); po++)
                tmp *= r3z2;
              if (co.l[2] & 0x1)
                tmp *= r3.z;
              for (int po = 0; po < (co.l[1] >> 1); po++)
                tmp *= r3y2;
              if (co.l[1] & 0x1)
                tmp *= r3.y;
              for (int po = 0; po < (co.l[0] >> 1); po++)
                tmp *= r3x2;
              if (co.l[0] & 0x1)
                tmp *= r3.x;
              accumulator[ic][tid] += tmp * grid_value;
            }
            break;
          }
        }
      }
    }

    const int max_i = min(length - ico, lbatch);

    __syncthreads();

    // we know there is only 1 wavefront in each block lbatch threads could
    // reduce the values saved by all threads of the warp and save results do a
    // shuffle_down by hand

    if (tid < 32) {
      for (int i = 0; i < max_i; i++) {
        // Load local value
        T val = accumulator[i][tid] + accumulator[i][tid + 32];

        // Warp-level reduction (32 lanes) After this loop, lane 0 of warp 0
        // holds the result All threads should execute this operation so
        // __activemask() is incorrect here.

        // The ifdef statement is needed because hip does not necessarily
        // support the instruction either. Reverting back to __shfl_down is
        // required.
        //
        for (int offset = 16; offset > 0; offset >>= 1) {
#if defined(__CUDACC__)
          val += __shfl_down_sync(0xffffffff, val, offset);
#else
          val += __shfl_down(val, offset);
#endif
        }

        if (tid == 0)
          accumulator[i][0] = val;
      }
#if defined(__CUDACC__)
      __syncwarp();
#endif
    }

#if !defined(__CUDACC__)
    __syncthreads();
#endif

    if (tid < min(length - ico, lbatch)) {
      const size_t coef_offset =
          dev_.tasks[dev_.first_task + block_index()].coef_offset;
      dev_.buffers_dev.coef[coef_offset + tid + ico] = accumulator[tid][0];
    }
    __syncthreads();
  }
}

/*******************************************************************************
 * Separable (plane streaming) integrate kernel.
 *
 * Computes the same coefficients as integrate_kernel, but with a cost per grid
 * point linear in lp instead of proportional to ncoset(lp), and reading the
 * cube only once whatever lp is.
 *
 * It uses the same underlying method than the DGEMM backend which is a tensor
 * contraction. Speedup comes from different areas:
 *
 * Axis aligned cell (axis_aligned_ = true): r3.x, r3.y, r3.z only depend on
 * x, y and z respectively, so that
 *
 *   x^a y^b z^c exp(-zetp r^2) = [x^a e^{-zetp x^2}] [y^b e^{-zetp y^2}]
 *                                [z^c e^{-zetp z^2}]
 *
 * Non orthogonal cell (axis_aligned_ = false): the moments are accumulated in
 * lattice monomials u^i v^j w^k (u, v, w are the lattice coordinates of the
 * point relative to the gaussian center) with the full gaussian evaluated per
 * point, and transformed to cartesian monomials at the end with the matrix
 * built by build_lattice_to_cartesian_matrix. It is the same transformation
 * as grid_transform_coef_jik_to_yxz in the dgemm backend, precomputed once per
 * grid level.
 *
 * For each z plane of the cube:
 *
 *   A[y][a]  = sum_x V[z][y][x] px[a][x]                       (step A)
 *   B[a][b]  = sum_y A[y][a] py[b][y]              a + b <= lp (step B)
 *   acc[coset(a, b, c)] += B[a][b] pz[c]       c <= lp - a - b (step C)
 *
 * B only lives in registers and the reduction over z is a running sum in acc,
 * so the shared memory does not depend on the z extent of the cube. The set of
 * grid points visited is the same as for integrate_kernel (window masks, y / x
 * trimming in the orthogonal case, sphere test).
 *
 * Shared memory layout: mc = max cube size, ls = lp_max + 1:
 *   ex[mc], py[ls][mc], pz[ls], A[mc][ls], acc[ncoset(lp_max)]
 * pz only holds the values of the current plane, ex is only used for an axis
 * aligned cell.
 ******************************************************************************/

template <bool distributed__, bool axis_aligned_, int lmax_>
__global__ __launch_bounds__(64) void integrate_separable_kernel(
    const kernel_params dev_, const int mc, const int lp_max,
    const double *__restrict__ lattice_to_cartesian) {
  // an axis aligned cell is orthogonal, otherwise the cell is treated as non
  // orthogonal (see integrate_one_grid_level).
  constexpr bool orthogonal_ = axis_aligned_;

  // the per thread accumulators have lmax_ + 1 entries, lmax_ >= lp_max (see
  // separable_lmax_bucket)
  static_assert(lmax_ <= CP2K_SEP_MAX_L, "lmax_ larger than CP2K_SEP_MAX_L");

  const task_info &glb_task = dev_.tasks[dev_.first_task + block_index()];
  if (glb_task.skip_task)
    return;

  // handled by integrate_kernel
  if (glb_task.la_max + dev_.la_max_diff + glb_task.lb_max + dev_.lb_max_diff <
      dev_.separable_lp_min)
    return;

  const int tid = thread_global_index();
  const int ls = lp_max + 1;

  extern __shared__ double sep_smem[];
  // exponetials along the x direction (note that it is the lattice direction
  // not necessarily the cartesian coordinate)
  double *ex = sep_smem;
  // y coordinate and its powers (note that it is the lattice direction not
  // necessarily the cartesian coordinate)
  double *py = ex + mc;
  // z coordinate and its powers common to all threads and single point (note
  // that it is the lattice direction not necessarily the cartesian coordinate)
  double *pz = py + ls * mc;
  double *A = pz + ls;
  double *acc = A + mc * ls;

  __shared__ smem_task_reduced<double, double3> task;
  fill_smem_task_reduced(dev_, dev_.first_task + block_index(), task);
  if (tid == 0) {
    setup_task_cube_center<double, double3, distributed__>(dev_, task);
  }
  __syncthreads();

  const int lp = task.lp;
  const int ncoef = ncoset(lp);
  const double zetp = task.zetp;
  const double r2max = task.radius * task.radius;
  const bool trim = orthogonal_ && !task.apply_border_mask;

  // 1D tables. Gaussian along x for an axis aligned cell (the x coordinates
  // are computed on the fly) and powers of y in py.
  for (int x = tid; x < task.cube_size.x; x += 64) {
    const double u = x + task.lb_cube.x + task.roffset.x;
    if (axis_aligned_) {
      const double r = u * dev_.dh_[0];
      ex[x] = exp(-zetp * r * r);
    }
  }

  // p[beta][y]
  for (int y = tid; y < task.cube_size.y; y += 64) {
    const double v = y + task.lb_cube.y + task.roffset.y;
    const double r = axis_aligned_ ? v * dev_.dh_[4] : v;
    double p = axis_aligned_ ? exp(-zetp * r * r) : 1.0;
    for (int b = 0; b <= lp; b++, p *= r)
      py[b * mc + y] = p;
  }
  for (int i = tid; i < ncoef; i += 64)
    acc[i] = 0.0;
  __syncthreads();

  const int npairs = (lp + 1) * (lp + 2) / 2;
  const int row_group = tid / CP2K_SEP_INTEGRATE_ROW_LANES;
  const int lane = tid % CP2K_SEP_INTEGRATE_ROW_LANES;
  constexpr int nrow_groups = 64 / CP2K_SEP_INTEGRATE_ROW_LANES;

  for (int z = 0; z < task.cube_size.z; z++) {
    // the plane decisions below are uniform over the block
    const int z2 =
        wrap_grid_index(z + task.cube_center.z, dev_.grid_full_size_.z);
    if (distributed__ && task.apply_border_mask &&
        ((z2 < task.window_shift.z) || (z2 > task.window_size.z)))
      continue;

    int ymin = 0, ymax = task.cube_size.y - 1;
    double kremain = 0.0;
    if (trim)
      kremain = calculate_ymix_ymax_boundaries(task, z, ymin, ymax);
    if (ymin > ymax)
      continue;

    const double w = z + task.lb_cube.z + task.roffset.z;
    const double rz = axis_aligned_ ? w * dev_.dh_[8] : w;

    // pz for this plane. Made visible to step C by the barrier after step A,
    // and only overwritten after the barrier that ends step C.

    // Calculate 1, z, z^2, .., z^lp recursively

    if (tid == 0) {
      double p = axis_aligned_ ? exp(-zetp * rz * rz) : 1.0;
      for (int c = 0; c <= lp; c++, p *= rz)
        pz[c] = p;
    }

    // step A: one group of CP2K_SEP_INTEGRATE_ROW_LANES lanes per row y. The
    // trip count of the y0 loop is uniform over the block so that every lane
    // takes part in the shuffle reduction, even when its row is past ymax.
    for (int y0 = ymin; y0 <= ymax; y0 += nrow_groups) {
      const int y = y0 + row_group;
      const bool row_valid = (y <= ymax);
      double part[lmax_ + 1];
#pragma unroll
      for (int a = 0; a <= lmax_; a++)
        part[a] = 0.0;

      const int y2 =
          wrap_grid_index(y + task.cube_center.y, dev_.grid_full_size_.y);
      const bool row_masked =
          distributed__ && task.apply_border_mask &&
          ((y2 < task.window_shift.y) || (y2 > task.window_size.y));

      if (row_valid && !row_masked) {
        int xmin = 0, xmax = task.cube_size.x - 1;
        if (trim)
          calculate_xmin_xmax_boundaries<double, double3>(task, y, kremain,
                                                          xmin, xmax);
        const double v = y + task.lb_cube.y + task.roffset.y;
        const double ry = v * dev_.dh_[4];
        const double r2yz = ry * ry + rz * rz; // only used when axis aligned
        const double *__restrict__ row =
            &dev_.buffers_dev.grid[(z2 * dev_.grid_local_size_.y + y2) *
                                   dev_.grid_local_size_.x];

        // periodic index of x, updated incrementally
        const int nx_full = dev_.grid_full_size_.x;
        int x2 = wrap_grid_index(xmin + lane + task.cube_center.x, nx_full);
        for (int x = xmin + lane; x <= xmax; x += CP2K_SEP_INTEGRATE_ROW_LANES,
                 x2 = sep_advance_wrapped_index<CP2K_SEP_INTEGRATE_ROW_LANES>(
                     x2, nx_full)) {
          if (distributed__ && task.apply_border_mask &&
              ((x2 < task.window_shift.x) || (x2 > task.window_size.x)))
            continue;

          const double u = x + task.lb_cube.x + task.roffset.x;
          const double r = axis_aligned_ ? u * dev_.dh_[0] : u;
          double r2;
          if (axis_aligned_) {
            r2 = r * r + r2yz;
          } else {
            const double3 r3 = compute_coordinates(dev_.dh_, r, v, w);
            r2 = r3.x * r3.x + r3.y * r3.y + r3.z * r3.z;
          }

          // same sphere test as integrate_kernel
          const bool outside = (r2max <= r2);
          if (distributed__) {
            if (outside && (!orthogonal_ || task.apply_border_mask))
              continue;
          } else if (!orthogonal_ && outside) {
            continue;
          }

          // Calculate x_i^\alpha \exp (-x_i^2) V_{ijk}

          double p = __ldg(&row[x2]);
          p *= axis_aligned_ ? ex[x] : exp(-zetp * r2);
#pragma unroll
          for (int a = 0; a <= lmax_; a++) {
            if (a <= lp) {
              part[a] += p;
              p *= r;
            }
          }
        }
      }

      // now calculate B_{j k \alpha} = \sum_{i} V_{ijk} x_i ^ \alpha \exp(
      // -x_i^2)
#pragma unroll
      for (int a = 0; a <= lmax_; a++) {
        if (a <= lp) {
          const double val = sep_row_reduce(part[a]);
          if (row_valid && lane == 0)
            A[(y - ymin) * ls + a] = val;
        }
      }
    }
    __syncthreads();

    // steps B and C: each thread owns the (a, b, *) coefficients of its pairs,
    // so no atomics are needed.

    // finally calculate
    //
    // C_{\alpha\beta\gamma} = \sum_{jk} B_{jk\alpha} p^y_{j\beta}
    // p^z_{k,\gamma}
    //
    // with \alpha+\beta+\gamma<=l. The pair list contains the \alpha,\beta
    // pairs that are required. $\gamma$ can be derived from the constraint over
    // angular momentum.

    for (int ip = tid; ip < npairs; ip += 64) {
      const int a = sep_pair_ab[ip].x;
      const int b = sep_pair_ab[ip].y;
      double bab = 0.0;
      for (int y = ymin; y <= ymax; y++)
        bab += A[(y - ymin) * ls + a] * py[b * mc + y];
      for (int c = 0; c <= lp - a - b; c++)
        acc[coset(a, b, c)] += bab * pz[c];
    }
    __syncthreads();
  }

  // finally store the results directly when the lattice is orthorhombic or
  // rotate the coefficients when the lattice is generic.
  //
  // the dgemm backend uses the same trick.
  //

  double *coef = &dev_.buffers_dev.coef[glb_task.coef_offset];
  if (axis_aligned_) {
    for (int i = tid; i < ncoef; i += 64)
      coef[i] = acc[i];
  } else {
    // lattice -> cartesian monomials, one block per total angular momentum
    for (int i = tid; i < ncoef; i += 64) {
      const orbital &co = coset_inv[i];
      const int l = co.l[0] + co.l[1] + co.l[2];
      const int base = (l > 0) ? ncoset(l - 1) : 0;
      const int nl = ncoset(l) - base;
      const double *__restrict__ m =
          lattice_to_cartesian + sep_m_offset[l] + (i - base) * nl;
      double val = 0.0;
      for (int j = 0; j < nl; j++)
        val += m[j] * acc[base + j];
      coef[i] = val;
    }
  }
}

static size_t integrate_separable_smem_size(const int mc, const int lp_max) {
  const int ls = lp_max + 1;
  return sizeof(double) * (mc + ls * mc + ls + mc * ls + ncoset(lp_max));
}

/*******************************************************************************
 * \brief Launches integrate_separable_kernel for one value of lmax_.
 ******************************************************************************/
template <int lmax_>
static void launch_integrate_separable(const separable_launch_config &cfg,
                                       const bool distributed,
                                       const bool axis_aligned,
                                       const kernel_params &params,
                                       const double *lattice_to_cartesian) {
  if (distributed) {
    if (axis_aligned)
      integrate_separable_kernel<true, true, lmax_>
          <<<cfg.nblocks, 64, cfg.smem, cfg.stream>>>(
              params, cfg.max_cube_size, cfg.lp_max, lattice_to_cartesian);
    else
      integrate_separable_kernel<true, false, lmax_>
          <<<cfg.nblocks, 64, cfg.smem, cfg.stream>>>(
              params, cfg.max_cube_size, cfg.lp_max, lattice_to_cartesian);
  } else {
    if (axis_aligned)
      integrate_separable_kernel<false, true, lmax_>
          <<<cfg.nblocks, 64, cfg.smem, cfg.stream>>>(
              params, cfg.max_cube_size, cfg.lp_max, lattice_to_cartesian);
    else
      integrate_separable_kernel<false, false, lmax_>
          <<<cfg.nblocks, 64, cfg.smem, cfg.stream>>>(
              params, cfg.max_cube_size, cfg.lp_max, lattice_to_cartesian);
  }
}

/*******************************************************************************
 * \brief Launches the Cuda kernel that integrates all tasks of one grid level.
 ******************************************************************************/
void context_info::integrate_one_grid_level(const int level, int *lp_diff) {
  if (number_of_tasks_per_level_[level] == 0)
    return;
  assert(!calculate_virial || calculate_forces);

  // Compute max angular momentum.
  const ldiffs_value ldiffs =
      process_get_ldiffs(calculate_forces, calculate_virial, compute_tau);

  smem_parameters smem_params(ldiffs, lmax());

  *lp_diff = smem_params.lp_diff();
  init_constant_memory();

  kernel_params params = set_kernel_parameters(level, smem_params);

  // Tasks with lp >= integrate_separable_lp_threshold go to the separable
  // kernel when it applies:
  // - axis aligned cells (cartesian factorization) and non orthogonal cells
  //   (lattice coordinates + transformation). Rotated orthogonal cells keep
  //   integrate_kernel since their cube is set up by the orthogonal code path.
  // - the angular momentum fits the tables and the shared memory the budget.
  auto &grid = grid_[level];
  const int lp_max = smem_params.lp_max();
  const int max_cube_size = grid.max_cube_size();
  const bool axis_aligned = grid.is_axis_aligned();

  // fall back to the old integrate if shared memory is not large enough or
  // anything else happens
  const bool use_separable =
      (lp_max >= integrate_separable_lp_threshold) &&
      (lp_max <= CP2K_SEP_MAX_L) && (max_cube_size > 0) &&
      (axis_aligned || !grid.is_orthogonal()) &&
      (integrate_separable_smem_size(max_cube_size, lp_max) <=
       separable_max_smem);

  if (use_separable)
    params.separable_lp_min = integrate_separable_lp_threshold;

  /* WARNING : if you change the block size please be aware of that the number
   * of warps is hardcoded when we do the block reduction in the integrate
   * kernel. The number of warps should be explicitly indicated in the
   * templating parameters or simply adjust the switch statements inside the
   * integrate kernels */

  // integrate_kernel handles the tasks the separable kernel does not take
  const bool use_integrate_kernel =
      !use_separable || (integrate_separable_lp_threshold > 0);
  const dim3 threads_per_block(4, 4, 4);
  if (!use_integrate_kernel) {
    // nothing to do
  } else if (grid_[level].is_distributed()) {
    if (grid_[level].is_orthogonal()) {
      integrate_kernel<double, double3, true, true>
          <<<number_of_tasks_per_level_[level], threads_per_block, 0,
             level_streams[level]>>>(params);
    } else {
      integrate_kernel<double, double3, true, false>
          <<<number_of_tasks_per_level_[level], threads_per_block, 0,
             level_streams[level]>>>(params);
    }
  } else {
    if (grid_[level].is_orthogonal()) {
      integrate_kernel<double, double3, false, true>
          <<<number_of_tasks_per_level_[level], threads_per_block, 0,
             level_streams[level]>>>(params);
    } else {
      integrate_kernel<double, double3, false, false>
          <<<number_of_tasks_per_level_[level], threads_per_block, 0,
             level_streams[level]>>>(params);
    }
  }

  if (!use_separable)
    return;

  // both kernels write the coefficients of disjoint sets of tasks.
  init_separable_constant_memory();
  const double *lattice_to_cartesian =
      axis_aligned ? nullptr : get_lattice_to_cartesian(grid, lp_max);
  const size_t smem = integrate_separable_smem_size(max_cube_size, lp_max);
  const int nblocks_level = number_of_tasks_per_level_[level];
  const separable_launch_config cfg{nblocks_level, smem, level_streams[level],
                                    max_cube_size, lp_max};
  switch (separable_lmax_bucket(lp_max)) {
  case 4:
    launch_integrate_separable<4>(cfg, grid.is_distributed(), axis_aligned,
                                  params, lattice_to_cartesian);
    break;
  case 8:
    launch_integrate_separable<8>(cfg, grid.is_distributed(), axis_aligned,
                                  params, lattice_to_cartesian);
    break;
  case 12:
    launch_integrate_separable<12>(cfg, grid.is_distributed(), axis_aligned,
                                   params, lattice_to_cartesian);
    break;
  default:
    launch_integrate_separable<CP2K_SEP_MAX_L>(
        cfg, grid.is_distributed(), axis_aligned, params, lattice_to_cartesian);
  }
}

void context_info::compute_hab_coefficients() {

  assert(!calculate_virial || calculate_forces);

  // Compute max angular momentum.
  const ldiffs_value ldiffs =
      process_get_ldiffs(calculate_forces, calculate_virial, compute_tau);

  smem_parameters smem_params(ldiffs, lmax());
  init_constant_memory();

  kernel_params params = set_kernel_parameters(-1, smem_params);

  /* WARNING if you change the block size. The number
   * of warps is hardcoded when we do the block reduction in the integrate
   * kernel. */

  const dim3 threads_per_block(4, 4, 4);

  // two step sphi transformation when its work arrays fit in the shared memory
  const int ld = ncoset(lmax());
  const size_t smem_two_steps =
      smem_params.smem_per_block() +
      sizeof(double) * sphi_work_size(lmax(), max_nsgf_set());
  const int work_offset = smem_params.smem_per_block() / sizeof(double);
  const bool two_step = (smem_two_steps <= sphi_two_step_max_smem);

  if (!compute_tau && !calculate_forces) {
    if (two_step)
      compute_hab_two_steps<double, false>
          <<<this->nblocks, threads_per_block, smem_two_steps,
             this->main_stream>>>(params, work_offset, ld, max_nsgf_set());
    else
      compute_hab<double, false>
          <<<this->nblocks, threads_per_block, smem_params.smem_per_block(),
             this->main_stream>>>(params);
    return;
  }

  if (compute_tau && !calculate_forces) {
    if (two_step)
      compute_hab_two_steps<double, true>
          <<<this->nblocks, threads_per_block, smem_two_steps,
             this->main_stream>>>(params, work_offset, ld, max_nsgf_set());
    else
      compute_hab<double, true>
          <<<this->nblocks, threads_per_block, smem_params.smem_per_block(),
             this->main_stream>>>(params);
  }

  if (!compute_tau && calculate_forces) {
    if (two_step)
      compute_hab_forces_two_steps<double, false>
          <<<this->nblocks, threads_per_block, smem_two_steps,
             this->main_stream>>>(params, work_offset, ld, max_nsgf_set());
    else
      compute_hab_forces<double, double3, false>
          <<<this->nblocks, threads_per_block, smem_params.smem_per_block(),
             this->main_stream>>>(params);
    return;
  }

  if (compute_tau && calculate_forces) {
    if (two_step)
      compute_hab_forces_two_steps<double, true>
          <<<this->nblocks, threads_per_block, smem_two_steps,
             this->main_stream>>>(params, work_offset, ld, max_nsgf_set());
    else
      compute_hab_forces<double, double3, true>
          <<<this->nblocks, threads_per_block, smem_params.smem_per_block(),
             this->main_stream>>>(params);
  }
}
}; // namespace rocm_backend
