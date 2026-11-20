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

#include "../../base/cp_assert.h"
#include <algorithm>
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

namespace grid_gpu_backend {

__inline__ __device__ double sep_row_reduce(double val) {
  for (int offset = GRID_SEP_INTEGRATE_ROW_LANES / 2; offset > 0;
       offset >>= 1) {
#if defined(__CUDACC__)
    val +=
        __shfl_down_sync(0xffffffff, val, offset, GRID_SEP_INTEGRATE_ROW_LANES);
#else
    val += __shfl_down(val, offset, GRID_SEP_INTEGRATE_ROW_LANES);
#endif
  }
  return val;
}

// do a block reduction for a block of size 64 and return the final sum to
// thread_id = 0. table is 64 elements of the dynamic shared memory, see
// block_reduce_64_smem.
constexpr size_t block_reduce_64_smem = 64 * sizeof(double);

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

// Reduces the virial of the block and stores it in virial_buffer[9][nblocks].
// The last block to finish then sums the per block values in a fixed order,
// so that the virial is deterministic and needs no atomics. Must be called by
// every block of the launch, including the ones without tasks (virial = 0).
template <typename T>
__device__ void reduce_virial(const kernel_params &dev_, T *sum,
                              const T virial[9], const int tid) {
  T *buf = dev_.buffers_dev.virial_buffer;
  const unsigned int nb = gridDim.x;
  const unsigned int b = block_index();
  unsigned int *retirement_count =
      (unsigned int *)(dev_.buffers_dev.virial_buffer + 9 * nb);

  for (int i = 0; i < 9; i++) {
    const T v = block_reduce_64<T>(sum, virial[i], tid);
    if (tid == 0)
      buf[i * nb + b] = v;
  }

  // make the partial sums visible before taking a ticket
  __shared__ bool last;
  __threadfence();
  if (tid == 0)
    last = (atomicInc(retirement_count, gridDim.x) == (gridDim.x - 1));
  __syncthreads();
  if (!last)
    return;

  if (dev_.buffers_dev.virial != nullptr) {
    for (int i = 0; i < 9; i++) {
      T s = 0.0;
      for (unsigned int k = tid; k < nb; k += 64)
        s += buf[i * nb + k]; // L2, never a stale L1 line
      s = block_reduce_64<T>(sum, s, tid);
      if (tid == 0)
        dev_.buffers_dev.virial[i] = s;
    }
  }
  if (tid == 0)
    *retirement_count = 0; // ready for the next launch
}
template <typename T>
__device__ void reduce_forces(const kernel_params &dev_, T *sum, T *fa, T *fb,
                              const int tid) {
  T *forces_a_tmp = &dev_.buffers_dev.forces_a[3 * block_index()];
  T *forces_b_tmp = &dev_.buffers_dev.forces_b[3 * block_index()];

  for (int i = 0; i < 3; i++) {
    fa[i] = block_reduce_64<T>(sum, fa[i], tid);
    if (tid == 0) {
      forces_a_tmp[i] = fa[i];
    }
    fb[i] = block_reduce_64<T>(sum, fb[i], tid);
    if (tid == 0) {
      forces_b_tmp[i] = fb[i];
    }
  }
}

template <typename T, bool COMPUTE_TAU>
__global__ __launch_bounds__(64) void compute_hab(const kernel_params dev_) {
  const int number_of_tasks = dev_.num_tasks_per_block_dev[block_index()];

  if (number_of_tasks == 0)
    return;

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
    cxyz_to_cab_factorized(task, coef_, smem_cab);
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
 *        U [max_nsgf_set][ld] live in the dynamic shared memory, see
 *        sphi_work_size.
 ******************************************************************************/
template <typename T, bool COMPUTE_TAU>
__global__ __launch_bounds__(64) void compute_hab_two_steps(
    const kernel_params dev_, const int ld, const int max_nsgf_set) {
  extern __shared__ T shared_memory[];
  const int number_of_tasks = dev_.num_tasks_per_block_dev[block_index()];

  if (number_of_tasks == 0)
    return;

  T *smem_h = &shared_memory[0]; // [ld][ld], H[jco][ico]
  T *smem_u = smem_h + ld * ld;  // [max_nsgf_set][ld]
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

/*******************************************************************************
 * \brief Adds the contribution of the cartesian pair (a, b), weighted by the
 *        cartesian density p, to the forces and the virial.
 ******************************************************************************/
template <typename T, bool COMPUTE_TAU>
__device__ __inline__ void
accumulate_pair_forces(const smem_task<T> &task, const orbital &a,
                       const orbital &b, const T p, const bool compute_virial,
                       const T *__restrict__ cab, T *fa, T *fb, T *virial) {
  for (int d = 0; d < 3; d++) {
    fa[d] += p * get_force_a<COMPUTE_TAU, T>(a, b, d, task.zeta, task.zetb,
                                             task.n1, cab);
    fb[d] += p * get_force_b<COMPUTE_TAU, T>(a, b, d, task.zeta, task.zetb,
                                             task.rab, task.n1, cab);
  }
  if (compute_virial) {
    for (int d1 = 0; d1 < 3; d1++)
      for (int d2 = 0; d2 < 3; d2++)
        virial[3 * d1 + d2] +=
            p *
            (get_virial_a<COMPUTE_TAU, T>(a, b, d1, d2, task.zeta, task.zetb,
                                          task.n1, cab) +
             get_virial_b<COMPUTE_TAU, T>(a, b, d1, d2, task.zeta, task.zetb,
                                          task.rab, task.n1, cab));
  }
}

/*******************************************************************************
 * \brief Fallback of compute_hab_forces_two_steps when its work arrays do not
 *        fit in the shared memory. hab is computed as in compute_hab. The
 *        forces and the virial only depend on the cartesian pair (a, b), so
 *        they are evaluated once per pair with the cartesian density
 *
 *          p[jco][ico] = sum_i sphib[i][jco] sum_j pab[i][j] sphia[j][ico]
 *
 *        computed on the fly (nsgfb (nsgfa + 1) FMA per pair), instead of
 *        once per (i, j, a, b) in the hab loop.
 ******************************************************************************/
template <typename T, typename T3, bool COMPUTE_TAU>
__global__
__launch_bounds__(64) void compute_hab_forces(const kernel_params dev_) {
  // only used by the final block reductions, block_reduce_64_smem bytes
  extern __shared__ T shared_memory[];
  const int number_of_tasks = dev_.num_tasks_per_block_dev[block_index()];

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
    fill_smem_task_coef(dev_, task_id, task);

    T *__restrict__ coef_ =
        &dev_.buffers_dev.coef[dev_.tasks[task_id].coef_offset];
    __syncthreads();
    cxyz_to_cab_factorized(task, coef_, smem_cab);
    __syncthreads();

    // hab[i][j] += sum_{jco, ico} sphib[i][jco] get_hab(a, b) sphia[j][ico]
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
            tmp += hab * task.sphia[j * task.maxcoa + ico] * sphib;
          }
        }
        if (task.block_transposed) {
          task.hab_block[j * task.nsgfb + i] += tmp;
        } else {
          task.hab_block[i * task.nsgfa + j] += tmp;
        }
      }
    }

    // forces and virial, once per cartesian pair
    const int nco_a = task.ncoseta - task.first_coseta;
    const int nco_b = task.ncosetb - task.first_cosetb;
    for (int k = tid; k < nco_b * nco_a; k += nthreads) {
      const int jj = k / nco_a;
      const int jco = task.first_cosetb + jj;
      const int ico = task.first_coseta + (k - jj * nco_a);

      T p = 0.0;
      for (int i = 0; i < task.nsgf_setb; i++) {
        T w = 0.0;
        for (int j = 0; j < task.nsgf_seta; j++) {
          const T block_value = task.block_transposed
                                    ? task.pab_block[j * task.nsgfb + i]
                                    : task.pab_block[i * task.nsgfa + j];
          w += block_value * task.sphia[j * task.maxcoa + ico];
        }
        p += task.sphib[i * task.maxcob + jco] * w;
      }
      accumulate_pair_forces<T, COMPUTE_TAU>(
          task, coset_inv[ico], coset_inv[jco], p * task.off_diag_twice,
          compute_virial, smem_cab, fa, fb, virial);
    }
    __syncthreads();
  }

  // theoretically not needed
  __syncthreads();

  T *sum = (T *)shared_memory;

  if (dev_.buffers_dev.virial != nullptr) {
    reduce_virial(dev_, sum, virial, tid);
  }

  if (number_of_tasks) {
    reduce_forces(dev_, sum, fa, fb, tid);
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
 *        [max_nsgf_set][ld] (U overwrites W) live in the dynamic shared
 *        memory, see sphi_work_size.
 ******************************************************************************/
template <typename T, bool COMPUTE_TAU>
__global__ __launch_bounds__(64) void compute_hab_forces_two_steps(
    const kernel_params dev_, const int ld, const int max_nsgf_set) {
  extern __shared__ T shared_memory[];
  const int number_of_tasks = dev_.num_tasks_per_block_dev[block_index()];

  T *smem_h = &shared_memory[0]; // [ld][ld], H[jco][ico]
  T *smem_w = smem_h + ld * ld;  // [max_nsgf_set][ld], W then U
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

      accumulate_pair_forces<T, COMPUTE_TAU>(task, a, b, p, compute_virial,
                                             smem_cab, fa, fb, virial);

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

  T *sum = (T *)shared_memory;
  if (compute_virial) {
    reduce_virial(dev_, sum, virial, tid);
  }

  if (number_of_tasks) {
    reduce_forces(dev_, sum, fa, fb, tid);
  }
}

/*******************************************************************************
 * Separable (plane streaming) integrate kernel.
 *
 * Computes the coefficients of the potential V for each task
 *
 *   coef[coset(a, b, c)] = sum_{ijk} V_{ijk} (x_i - x_c)^a (y_j - y_c)^b
 *                          (z_k - z_c)^c exp(-zetp |r_{ijk} - r_c|^2)
 *
 * for a + b + c <= lp, with a cost per grid point linear in lp instead of
 * proportional to ncoset(lp), and reading the cube only once whatever lp is.
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
 * grid points visited is the same as for collocate_separable_kernel (window
 * masks, y / x trimming in the orthogonal case, sphere test).
 *
 * The tasks are split between two instances (see sep_skip_task):
 *
 * - chunked_ = false: tasks with a cube edge of at most mc points, the whole
 *   plane in one chunk (yrows = mc), which the compiler reduces to a single
 *   pass, and the gaussian along x tabulated in ex.
 * - chunked_ = true: larger cubes, the rows of a plane by chunks of yrows rows
 *   (see integrate_separable_yrows), B accumulated over the chunks in the
 *   shared memory and step C once per plane. The gaussian along x is evaluated
 *   on the fly, so that the shared memory does not depend on the cube size.
 *
 * Shared memory layout (doubles), ls = lp_max + 1:
 *   ex[mc], pz[ls], acc[ncoset(lp_max)], py[ls][yrows], A[yrows][ls],
 *   B[ls (ls + 1) / 2] (chunked_ only)
 * ex is only allocated for an axis aligned cell and chunked_ = false. pz only
 * holds the values of the current plane. py is filled once per task, or for
 * each chunk when chunked_.
 ******************************************************************************/

template <bool distributed__, bool axis_aligned_, int lmax_, bool chunked_>
__global__ __launch_bounds__(64) void integrate_separable_kernel(
    const kernel_params dev_, const int mc, const int lp_max, const int yrows,
    const double *__restrict__ lattice_to_cartesian) {
  // an axis aligned cell is orthogonal, otherwise the cell is treated as non
  // orthogonal (see integrate_one_grid_level).
  constexpr bool orthogonal_ = axis_aligned_;

  // the per thread accumulators have lmax_ + 1 entries, lmax_ >= lp_max (see
  // separable_lmax_bucket)
  static_assert(lmax_ <= GRID_SEP_MAX_L, "lmax_ larger than GRID_SEP_MAX_L");
  // the gaussian along x is tabulated
  constexpr bool tables_ = axis_aligned_ && !chunked_;

  const task_info &glb_task = dev_.tasks[dev_.first_task + block_index()];
  if (sep_skip_task<chunked_>(glb_task, mc))
    return;

  const int tid = thread_global_index();
  const int ls = lp_max + 1;

  extern __shared__ double sep_smem[];
  // exponetials along the x direction (note that it is the lattice direction
  // not necessarily the cartesian coordinate)
  double *ex = sep_smem; // tables_ only
  // z coordinate and its powers common to all threads and single point (note
  // that it is the lattice direction not necessarily the cartesian coordinate)
  double *pz = ex + (tables_ ? mc : 0);
  double *acc = pz + ls;
  // y coordinate and its powers (note that it is the lattice direction not
  // necessarily the cartesian coordinate), py[b * yrows + row]
  double *py = acc + ncoset(lp_max);
  double *A = py + ls * yrows; // [yrows][ls]
  // B[a][b] of the current plane, accumulated over the chunks of rows
  double *B = A + yrows * ls; // chunked_ only

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
    if (tables_) {
      const double r = u * dev_.dh_[0];
      ex[x] = exp(-zetp * r * r);
    }
  }

  // py[beta][y - ybase] for the rows y0 <= y <= y1
  auto fill_py = [&](const int y0, const int y1, const int ybase) {
    for (int y = y0 + tid; y <= y1; y += 64) {
      const double v = y + task.lb_cube.y + task.roffset.y;
      const double r = axis_aligned_ ? v * dev_.dh_[4] : v;
      double p = axis_aligned_ ? exp(-zetp * r * r) : 1.0;
      for (int b = 0; b <= lp; b++, p *= r)
        py[b * yrows + y - ybase] = p;
    }
  };

  if (!chunked_)
    fill_py(0, task.cube_size.y - 1, 0);
  for (int i = tid; i < ncoef; i += 64)
    acc[i] = 0.0;
  __syncthreads();

  const int npairs = (lp + 1) * (lp + 2) / 2;
  const int row_group = tid / GRID_SEP_INTEGRATE_ROW_LANES;
  const int lane = tid % GRID_SEP_INTEGRATE_ROW_LANES;
  constexpr int nrow_groups = 64 / GRID_SEP_INTEGRATE_ROW_LANES;

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

    // B is owned by the threads of step B, pair ip by thread ip % 64, so that
    // it needs no barrier of its own.
    if (chunked_)
      for (int ip = tid; ip < npairs; ip += 64)
        B[ip] = 0.0;

    // the trip count of the chunk loop is uniform over the block, a single
    // pass over the whole plane unless chunked_
    for (int yc = ymin; yc <= ymax; yc = chunked_ ? yc + yrows : ymax + 1) {
      const int yc1 = chunked_ ? min(ymax, yc + yrows - 1) : ymax;
      // py is read in step B, after the barrier that ends step A
      const int ybase = chunked_ ? yc : 0;
      if (chunked_)
        fill_py(yc, yc1, yc);

      // step A: one group of GRID_SEP_INTEGRATE_ROW_LANES lanes per row y. The
      // trip count of the y0 loop is uniform over the block so that every lane
      // takes part in the shuffle reduction, even when its row is past yc1.
      for (int y0 = yc; y0 <= yc1; y0 += nrow_groups) {
        const int y = y0 + row_group;
        const bool row_valid = (y <= yc1);
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
          for (int x = xmin + lane; x <= xmax;
               x += GRID_SEP_INTEGRATE_ROW_LANES,
                   x2 = sep_advance_wrapped_index<GRID_SEP_INTEGRATE_ROW_LANES>(
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

            // same sphere test as collocate_separable_kernel
            const bool outside = (r2max <= r2);
            if (distributed__) {
              if (outside && (!orthogonal_ || task.apply_border_mask))
                continue;
            } else if (!orthogonal_ && outside) {
              continue;
            }

            // Calculate x_i^\alpha \exp (-x_i^2) V_{ijk}

            double p = __ldg(&row[x2]);
            p *= tables_ ? ex[x] : exp(-zetp * (axis_aligned_ ? r * r : r2));
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
              A[(y - yc) * ls + a] = val;
          }
        }
      }
      __syncthreads();

      // steps B and C: each thread owns the (a, b, *) coefficients of its
      // pairs, so no atomics are needed.

      // finally calculate
      //
      // C_{\alpha\beta\gamma} = \sum_{jk} B_{jk\alpha} p^y_{j\beta}
      // p^z_{k,\gamma}
      //
      // with \alpha+\beta+\gamma<=l. The pair list contains the \alpha,\beta
      // pairs that are required. $\gamma$ can be derived from the constraint
      // over angular momentum.

      for (int ip = tid; ip < npairs; ip += 64) {
        const int a = sep_pair_ab[ip].x;
        const int b = sep_pair_ab[ip].y;
        double bab = 0.0;
        for (int y = yc; y <= yc1; y++)
          bab += A[(y - yc) * ls + a] * py[b * yrows + y - ybase];
        if (chunked_) {
          B[ip] += bab;
        } else {
          for (int c = 0; c <= lp - a - b; c++)
            acc[coset(a, b, c)] += bab * pz[c];
        }
      }
      // A and py are overwritten by the next chunk
      if (chunked_ && (yc1 < ymax))
        __syncthreads();
    }

    if (chunked_) {
      for (int ip = tid; ip < npairs; ip += 64) {
        const int a = sep_pair_ab[ip].x;
        const int b = sep_pair_ab[ip].y;
        const double bab = B[ip];
        for (int c = 0; c <= lp - a - b; c++)
          acc[coset(a, b, c)] += bab * pz[c];
      }
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

// shared memory of integrate_separable_kernel, yrows rows per chunk and ntab
// entries in the table of the gaussian along x
static size_t integrate_separable_smem_size(const int lp_max, const int yrows,
                                            const int ntab,
                                            const bool chunked) {
  const int ls = lp_max + 1;
  const int nb = chunked ? ls * (ls + 1) / 2 : 0;
  return sizeof(double) * (ntab + ls + ncoset(lp_max) + 2 * ls * yrows + nb);
}

// largest cube edge processed by the instances with chunked_ = false: all the
// rows of the cube (yrows = mc) and for an axis aligned cell the table ex (mc
// entries) in smem_max.
static int integrate_separable_max_edge(const int lp_max,
                                        const bool axis_aligned,
                                        const size_t smem_max) {
  const size_t fixed = integrate_separable_smem_size(lp_max, 0, 0, false);
  const size_t per_edge =
      sizeof(double) * (2 * (lp_max + 1) + (axis_aligned ? 1 : 0));
  return (smem_max > fixed) ? (int)((smem_max - fixed) / per_edge) : 0;
}

// rows per chunk of the chunked_ instances, at most mc
static int integrate_separable_yrows(const int mc, const int lp_max,
                                     const size_t smem_max) {
  const size_t fixed = integrate_separable_smem_size(lp_max, 0, 0, true);
  const size_t per_row = 2 * sizeof(double) * (lp_max + 1);
  if (fixed + per_row > smem_max) {
    fprintf(stderr, "ERROR: lp %d too large for integrate_separable_kernel.\n",
            lp_max);
    abort();
  }
  return (int)std::min<size_t>(mc, (smem_max - fixed) / per_row);
}

using integrate_separable_kernel_t = void (*)(const kernel_params, const int,
                                              const int, const int,
                                              const double *);

/*******************************************************************************
 * \brief Selects the instance of integrate_separable_kernel for one value of
 *        lmax_.
 ******************************************************************************/
template <int lmax_, bool chunked_ = false>
static integrate_separable_kernel_t
select_integrate_separable(const bool distributed, const bool axis_aligned) {
  if (distributed)
    return axis_aligned
               ? integrate_separable_kernel<true, true, lmax_, chunked_>
               : integrate_separable_kernel<true, false, lmax_, chunked_>;
  return axis_aligned
             ? integrate_separable_kernel<false, true, lmax_, chunked_>
             : integrate_separable_kernel<false, false, lmax_, chunked_>;
}

/*******************************************************************************
 * \brief Selects the instance of integrate_separable_kernel with the smallest
 *        bucket of the angular momentum covering lp_max. The chunked instances
 *        are only needed for very large cubes, they are only instantiated for
 *        the largest bucket.
 ******************************************************************************/
static integrate_separable_kernel_t
select_integrate_separable(const int lp_max, const bool distributed,
                           const bool axis_aligned, const bool chunked) {
  if (chunked)
    return select_integrate_separable<GRID_SEP_MAX_L, true>(distributed,
                                                            axis_aligned);
  switch (separable_lmax_bucket(lp_max)) {
  case 4:
    return select_integrate_separable<4>(distributed, axis_aligned);
  case 8:
    return select_integrate_separable<8>(distributed, axis_aligned);
  case 12:
    return select_integrate_separable<12>(distributed, axis_aligned);
  default:
    return select_integrate_separable<GRID_SEP_MAX_L>(distributed,
                                                      axis_aligned);
  }
}

/*******************************************************************************
 * \brief Launches the Cuda kernel that integrates all tasks of one grid level.
 ******************************************************************************/
void context_info::integrate_one_grid_level(const int level, int *lp_diff) {
  if (number_of_tasks_per_level_[level] == 0)
    return;
  cp_assert(!calculate_virial || calculate_forces);

  // Compute max angular momentum.
  const ldiffs_value ldiffs =
      process_get_ldiffs(calculate_forces, calculate_virial, compute_tau);

  smem_parameters smem_params(ldiffs, lmax());

  *lp_diff = smem_params.lp_diff();
  init_constant_memory();

  kernel_params params = set_kernel_parameters(level, smem_params);

  // cartesian factorization for an axis aligned cell, lattice coordinates +
  // transformation otherwise (rotated orthogonal cells included, see
  // grid_info::check_orthogonality).
  auto &grid = grid_[level];
  const int lp_max = smem_params.lp_max();
  const int max_cube_size = std::max(grid.max_cube_size(), 1);
  const bool axis_aligned = grid.is_axis_aligned();
  // lp_max <= grid_gpu_max_l == GRID_SEP_MAX_L is checked by smem_parameters
  static_assert(grid_gpu_max_l <= GRID_SEP_MAX_L, "separable tables too small");

  init_separable_constant_memory();
  const double *lattice_to_cartesian =
      axis_aligned ? nullptr : get_lattice_to_cartesian(grid, lp_max);
  const bool distributed = grid.is_distributed();
  const int nblocks_level = number_of_tasks_per_level_[level];

  // the tasks with a cube edge of at most mc, the whole plane in one pass (see
  // sep_skip_task)
  const integrate_separable_kernel_t kernel =
      select_integrate_separable(lp_max, distributed, axis_aligned, false);
  const int mc = std::min(max_cube_size,
                          integrate_separable_max_edge(
                              lp_max, axis_aligned,
                              max_dynamic_smem(kernel, grid_gpu_smem_budget)));
  if (mc > 0) {
    const size_t smem =
        integrate_separable_smem_size(lp_max, mc, axis_aligned ? mc : 0, false);
    GRID_GPU_LAUNCH(kernel, nblocks_level, 64, smem, level_streams[level],
                    params, mc, lp_max, mc, lattice_to_cartesian);
  }

  // the larger cubes by chunks of rows, on the same stream
  if (max_cube_size > mc) {
    const integrate_separable_kernel_t chunked =
        select_integrate_separable(lp_max, distributed, axis_aligned, true);
    const int yrows = integrate_separable_yrows(
        max_cube_size, lp_max, max_dynamic_smem(chunked, grid_gpu_smem_budget));
    const size_t smem = integrate_separable_smem_size(lp_max, yrows, 0, true);
    GRID_GPU_LAUNCH(chunked, nblocks_level, 64, smem, level_streams[level],
                    params, mc, lp_max, yrows, lattice_to_cartesian);
  }
}

void context_info::compute_hab_coefficients() {

  cp_assert(!calculate_virial || calculate_forces);

  // Compute max angular momentum.
  const ldiffs_value ldiffs =
      process_get_ldiffs(calculate_forces, calculate_virial, compute_tau);

  smem_parameters smem_params(ldiffs, lmax());
  init_constant_memory();
  prepare_cab_workspace(ldiffs.la_max_diff, ldiffs.lb_max_diff);

  kernel_params params = set_kernel_parameters(-1, smem_params);

  /* WARNING if you change the block size. The number
   * of warps is hardcoded when we do the block reduction in the integrate
   * kernel. */

  const dim3 threads_per_block(4, 4, 4);

  // two step sphi transformation when its work arrays fit in the shared memory
  // next to the static shared memory of the kernel
  const int ld = ncoset(lmax());
  // (the forces kernels reuse the first 64 doubles for the block reductions)
  const size_t smem_two_steps =
      std::max(sizeof(double) * sphi_work_size(lmax(), max_nsgf_set()),
               block_reduce_64_smem);

  const auto two_steps_kernel =
      calculate_forces
          ? (compute_tau ? compute_hab_forces_two_steps<double, true>
                         : compute_hab_forces_two_steps<double, false>)
          : (compute_tau ? compute_hab_two_steps<double, true>
                         : compute_hab_two_steps<double, false>);
  if (smem_two_steps <=
      max_dynamic_smem(two_steps_kernel, grid_gpu_smem_budget)) {
    GRID_GPU_LAUNCH(two_steps_kernel, this->nblocks, threads_per_block,
                    smem_two_steps, this->main_stream, params, ld,
                    max_nsgf_set());
    return;
  }

  // fused fallback, the forces kernels need block_reduce_64_smem
  const auto fused_kernel =
      calculate_forces
          ? (compute_tau ? compute_hab_forces<double, double3, true>
                         : compute_hab_forces<double, double3, false>)
          : (compute_tau ? compute_hab<double, true>
                         : compute_hab<double, false>);
  const size_t smem_fused = calculate_forces ? block_reduce_64_smem : 0;
  GRID_GPU_LAUNCH(fused_kernel, this->nblocks, threads_per_block, smem_fused,
                  this->main_stream, params);
}
}; // namespace grid_gpu_backend
