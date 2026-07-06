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

#if defined(_OMP_H)
#error "OpenMP should not be used in .cu files to accommodate HIP."
#endif

namespace rocm_backend {

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

    // for (int i = tid / 8; i < task.nsgf_setb; i += 8) {
    //   for (int j = tid % 8; j < task.nsgf_seta; j += 8) {
    for (int loop = tid; loop < task.nsgf_seta * task.nsgf_setb;
         loop += blockDim.x * blockDim.y * blockDim.z) {
      const int j = loop % task.nsgf_seta;
      const int i = loop / task.nsgf_seta;
      T tmp = 0.0;
      for (int jco = task.first_cosetb; jco < task.ncosetb; jco++) {
        const T sphib = task.sphib[i * task.maxcob + jco];
        const auto &b = coset_inv[jco];
        for (int ico = task.first_coseta; ico < task.ncoseta; ico++) {
          const auto &a = coset_inv[ico];
          const T hab = get_hab<COMPUTE_TAU, T>(a, b, task.zeta, task.zetb,
                                                task.n1, smem_cab);
          const T sphia_times_sphib = task.sphia[j * task.maxcoa + ico] * sphib;
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
  // }
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

    // for (int i = tid / 8; i < task.nsgf_setb; i += 8) {
    //   for (int j = tid % 8; j < task.nsgf_seta; j += 8) {
    for (int loop = tid; loop < task.nsgf_seta * task.nsgf_setb;
         loop += blockDim.x * blockDim.y * blockDim.z) {
      const int j = loop % task.nsgf_seta;
      const int i = loop / task.nsgf_seta;
      T tmp = 0.0;
      T block_value = 0.0;
      if (task.block_transposed) {
        block_value = task.pab_block[j * task.nsgfb + i] * task.off_diag_twice;
      } else {
        block_value = task.pab_block[i * task.nsgfa + j] * task.off_diag_twice;
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
                (get_virial_a<COMPUTE_TAU, T>(a, b, 0, 0, task.zeta, task.zetb,
                                              task.n1, smem_cab) +
                 get_virial_b<COMPUTE_TAU, T>(a, b, 0, 0, task.zeta, task.zetb,
                                              task.rab, task.n1, smem_cab));
            virial[1] +=
                sphia_times_sphib *
                (get_virial_a<COMPUTE_TAU, T>(a, b, 0, 1, task.zeta, task.zetb,
                                              task.n1, smem_cab) +
                 get_virial_b<COMPUTE_TAU, T>(a, b, 0, 1, task.zeta, task.zetb,
                                              task.rab, task.n1, smem_cab));
            virial[2] +=
                sphia_times_sphib *
                (get_virial_a<COMPUTE_TAU, T>(a, b, 0, 2, task.zeta, task.zetb,
                                              task.n1, smem_cab) +
                 get_virial_b<COMPUTE_TAU, T>(a, b, 0, 2, task.zeta, task.zetb,
                                              task.rab, task.n1, smem_cab));
            virial[3] +=
                sphia_times_sphib *
                (get_virial_a<COMPUTE_TAU, T>(a, b, 1, 0, task.zeta, task.zetb,
                                              task.n1, smem_cab) +
                 get_virial_b<COMPUTE_TAU, T>(a, b, 1, 0, task.zeta, task.zetb,
                                              task.rab, task.n1, smem_cab));
            virial[4] +=
                sphia_times_sphib *
                (get_virial_a<COMPUTE_TAU, T>(a, b, 1, 1, task.zeta, task.zetb,
                                              task.n1, smem_cab) +
                 get_virial_b<COMPUTE_TAU, T>(a, b, 1, 1, task.zeta, task.zetb,
                                              task.rab, task.n1, smem_cab));
            virial[5] +=
                sphia_times_sphib *
                (get_virial_a<COMPUTE_TAU, T>(a, b, 1, 2, task.zeta, task.zetb,
                                              task.n1, smem_cab) +
                 get_virial_b<COMPUTE_TAU, T>(a, b, 1, 2, task.zeta, task.zetb,
                                              task.rab, task.n1, smem_cab));
            virial[6] +=
                sphia_times_sphib *
                (get_virial_a<COMPUTE_TAU, T>(a, b, 2, 0, task.zeta, task.zetb,
                                              task.n1, smem_cab) +
                 get_virial_b<COMPUTE_TAU, T>(a, b, 2, 0, task.zeta, task.zetb,
                                              task.rab, task.n1, smem_cab));
            virial[7] +=
                sphia_times_sphib *
                (get_virial_a<COMPUTE_TAU, T>(a, b, 2, 1, task.zeta, task.zetb,
                                              task.n1, smem_cab) +
                 get_virial_b<COMPUTE_TAU, T>(a, b, 2, 1, task.zeta, task.zetb,
                                              task.rab, task.n1, smem_cab));
            virial[8] +=
                sphia_times_sphib *
                (get_virial_a<COMPUTE_TAU, T>(a, b, 2, 2, task.zeta, task.zetb,
                                              task.n1, smem_cab) +
                 get_virial_b<COMPUTE_TAU, T>(a, b, 2, 2, task.zeta, task.zetb,
                                              task.rab, task.n1, smem_cab));
          }
        }
      }

      if (task.block_transposed) {
        task.hab_block[j * task.nsgfb + i] += tmp;
      } else {
        task.hab_block[i * task.nsgfa + j] += tmp;
      }
      // }
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
 * Cuda kernel for calcualting the coefficients of a potential (density,
 etc...) for a given pair of gaussian basis sets

We compute the discretized version of the following integral

$$\int _\infty ^\infty V_{ijk} P^\alpha_iP^beta_j P ^ \gamma _ k Exp(- \eta
|r_{ijk} - r_c|^2)$$

where in practice the summation is over a finite domain. The discrete form has
this shape

$$
\sum_{ijk < dmoain} V_{ijk} (x_i - x_c)^\alpha (y_j - y_c)^\beta (z_k - z_c) ^
\gamma Exp(- \eta |r_{ijk} - r_c|^2)
$$

where $0 \le \alpha + \beta + \gamma \le lmax$

It is formely the same operation than collocate (from a discrete point of view)
but the implementation differ because of technical reasons.

So most of the code is the same except the core of the routine that is
specialized to the integration.

******************************************************************************/
template <typename T, typename T3, bool distributed__, bool orthogonal_,
          int lbatch = 20>
__global__
__launch_bounds__(64) void integrate_kernel(const kernel_params dev_) {
  if (dev_.tasks[dev_.first_task + block_index()].skip_task)
    return;

  const int tid = thread_global_index();

  // __shared__ T dh_inv_[9];
  __shared__ T dh_[9];

  // for (int d = tid; d < 9; d += blockDim.x * blockDim.y * blockDim.z)
  //   dh_inv_[d] = dev_.dh_inv_[d];

  if (tid < 9)
    dh_[tid] = dev_.dh_[tid];

  __shared__ smem_task_reduced<T, T3> task;
  fill_smem_task_reduced(dev_, dev_.first_task + block_index(), task);

  if (tid == 0) {
    setup_task_cube_center<T, T3, distributed__>(dev_, task);
  }
  __syncthreads();

  __shared__ T accumulator[lbatch][64];

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

          // read the next point assuming that reading is non blocking untill
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
      __syncwarp();
    }

    if (tid < min(length - ico, lbatch)) {
      const size_t coef_offset =
          dev_.tasks[dev_.first_task + block_index()].coef_offset;
      dev_.buffers_dev.coef[coef_offset + tid + ico] = accumulator[tid][0];
    }
    __syncthreads();
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

  /* WARNING : if you change the block size please be aware of that the number
   * of warps is hardcoded when we do the block reduction in the integrate
   * kernel. The number of warps should be explicitly indicated in the
   * templating parameters or simply adjust the switch statements inside the
   * integrate kernels */

  const dim3 threads_per_block(4, 4, 4);
  if (grid_[level].is_distributed()) {
    if (grid_[level].is_orthogonal()) {
      integrate_kernel<double, double3, true, true, 20>
          <<<number_of_tasks_per_level_[level], threads_per_block, 0,
             level_streams[level]>>>(params);
    } else {
      integrate_kernel<double, double3, true, false, 20>
          <<<number_of_tasks_per_level_[level], threads_per_block, 0,
             level_streams[level]>>>(params);
    }
  } else {
    if (grid_[level].is_orthogonal()) {
      integrate_kernel<double, double3, false, true, 20>
          <<<number_of_tasks_per_level_[level], threads_per_block, 0,
             level_streams[level]>>>(params);
    } else {
      integrate_kernel<double, double3, false, false, 20>
          <<<number_of_tasks_per_level_[level], threads_per_block, 0,
             level_streams[level]>>>(params);
    }
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

  if (!compute_tau && !calculate_forces) {
    compute_hab<double, false>
        <<<this->nblocks, threads_per_block, smem_params.smem_per_block(),
           this->main_stream>>>(params);
    return;
  }

  if (compute_tau && !calculate_forces) {
    compute_hab<double, true>
        <<<this->nblocks, threads_per_block, smem_params.smem_per_block(),
           this->main_stream>>>(params);
  }

  if (!compute_tau && calculate_forces) {
    compute_hab_forces<double, double3, false>
        <<<this->nblocks, threads_per_block, smem_params.smem_per_block(),
           this->main_stream>>>(params);
    return;
  }

  if (compute_tau && calculate_forces) {
    compute_hab_forces<double, double3, true>
        <<<this->nblocks, threads_per_block, smem_params.smem_per_block(),
           this->main_stream>>>(params);
  }
}
}; // namespace rocm_backend
