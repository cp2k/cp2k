/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

/*
 * inspirartions from the gpu backend
 * Authors :
 - Dr Mathieu Taillefumier (ETH Zurich / CSCS)
 - Advanced Micro Devices, Inc.
*/

#ifdef __GRID_HIP

#include <algorithm>
#include <assert.h>
#include <hip/hip_runtime.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "grid_hip_context.h"
#include "grid_hip_internal_header.h"
#include "grid_hip_process_vab.h"

#ifdef __HIP_PLATFORM_NVIDIA__
#include <cooperative_groups.h>
#if CUDA_VERSION >= 11000
#include <cooperative_groups/reduce.h>
#endif
namespace cg = cooperative_groups;
#endif

namespace rocm_backend {
// do a warp reduction and return the final sum to thread_id = 0
template <typename T>
__device__ __inline__ T warp_reduce(T *table, const int tid) {
  // AMD GPU have warp size of 64 while nvidia GPUs have warpSize of 32 so the
  // first step is common to both platforms.
  if (tid < 32) {
    table[tid] += table[tid + 32];
  }
  __syncthreads();
  if (tid < 16) {
    table[tid] += table[tid + 16];
  }
  __syncthreads();
  if (tid < 8) {
    table[tid] += table[tid + 8];
  }
  __syncthreads();
  if (tid < 4) {
    table[tid] += table[tid + 4];
  }
  __syncthreads();
  if (tid < 2) {
    table[tid] += table[tid + 2];
  }
  __syncthreads();
  return (table[0] + table[1]);
}
// #endif

template <typename T, typename T3, bool COMPUTE_TAU, bool CALCULATE_FORCES>
__global__ __launch_bounds__(64) void compute_hab_v2(const kernel_params dev_,
                                                     const int ntasks) {
  // Copy task from global to shared memory and precompute some stuff.
  extern __shared__ T shared_memory[];
  T *coef_shared = &shared_memory[0];
  T *smem_cab = &shared_memory[dev_.smem_cab_offset];
  T *smem_alpha = &shared_memory[dev_.smem_alpha_offset];
  const int tid =
      threadIdx.x + blockDim.x * (threadIdx.y + blockDim.y * threadIdx.z);
  const int offset = dev_.sorted_blocks_offset_dev[blockIdx.x];
  const int number_of_tasks = dev_.num_tasks_per_block_dev[blockIdx.x];

  if (number_of_tasks == 0)
    return;

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

    T *coef_ = &dev_.ptr_dev[2][dev_.tasks[task_id].coef_offset];

    for (int i = tid; i < ncoset(task.lp);
         i += blockDim.x * blockDim.y * blockDim.z)
      coef_shared[i] = coef_[i];
    compute_alpha(dev_, task, smem_alpha);
    __syncthreads();
    cxyz_to_cab(dev_, task, smem_alpha, coef_shared, smem_cab);
    __syncthreads();

    for (int i = tid / 8; i < task.nsgf_setb; i += 8) {
      for (int j = tid % 8; j < task.nsgf_seta; j += 8) {
        T res = 0.0;

        T block_val1 = 0.0;

        if (CALCULATE_FORCES) {
          if (task.block_transposed) {
            block_val1 =
                task.pab_block[j * task.nsgfb + i] * task.off_diag_twice;
          } else {
            block_val1 =
                task.pab_block[i * task.nsgfa + j] * task.off_diag_twice;
          }
        }

        for (int jco = task.first_cosetb; jco < task.ncosetb; jco++) {
          const auto &b = coset_inv[jco];
          T block_val = 0.0;
          const T sphib = task.sphib[i * task.maxcob + jco];
          for (int ico = task.first_coseta; ico < task.ncoseta; ico++) {
            const auto &a = coset_inv[ico];
            const T hab = get_hab<COMPUTE_TAU, T>(a, b, task.zeta, task.zetb,
                                                  task.n1, smem_cab);
            const T sphia = task.sphia[j * task.maxcoa + ico];
            block_val += hab * sphia * sphib;

            if (CALCULATE_FORCES) {
              for (int k = 0; k < 3; k++) {
                fa[k] += block_val1 * sphia * sphib *
                         get_force_a<COMPUTE_TAU, T>(
                             a, b, k, task.zeta, task.zetb, task.n1, smem_cab);
                fb[k] +=
                    block_val1 * sphia * sphib *
                    get_force_b<COMPUTE_TAU, T>(a, b, k, task.zeta, task.zetb,
                                                task.rab, task.n1, smem_cab);
              }
              if (dev_.ptr_dev[5] != nullptr) {
                for (int k = 0; k < 3; k++) {
                  for (int l = 0; l < 3; l++) {
                    virial[3 * k + l] +=
                        (get_virial_a<COMPUTE_TAU, T>(a, b, k, l, task.zeta,
                                                      task.zetb, task.n1,
                                                      smem_cab) +
                         get_virial_b<COMPUTE_TAU, T>(a, b, k, l, task.zeta,
                                                      task.zetb, task.rab,
                                                      task.n1, smem_cab)) *
                        block_val1 * sphia * sphib;
                  }
                }
              }
            }
          }

          res += block_val;
        }

        // we can use shuffle_down if it exists for T

        if (task.block_transposed) {
          task.hab_block[j * task.nsgfb + i] += res;
        } else {
          task.hab_block[i * task.nsgfa + j] += res;
        }
      }
    }
    __syncthreads();
  }

  if (CALCULATE_FORCES) {
    const int task_id = dev_.task_sorted_by_blocks_dev[offset];
    const auto &glb_task = dev_.tasks[task_id];
    const int iatom = glb_task.iatom;
    const int jatom = glb_task.jatom;
    T *forces_a = &dev_.ptr_dev[4][3 * iatom];
    T *forces_b = &dev_.ptr_dev[4][3 * jatom];

    T *sum = (T *)shared_memory;
    if (dev_.ptr_dev[5] != nullptr) {

      for (int i = 0; i < 9; i++) {
        sum[tid] = virial[i];
        __syncthreads();

        virial[i] = warp_reduce<T>(sum, tid);
        __syncthreads();

        if (tid == 0)
          atomicAdd(dev_.ptr_dev[5] + i, virial[i]);
        __syncthreads();
      }
    }

    for (int i = 0; i < 3; i++) {
      sum[tid] = fa[i];
      __syncthreads();
      fa[i] = warp_reduce<T>(sum, tid);
      __syncthreads();
      if (tid == 0)
        atomicAdd(forces_a + i, fa[i]);
      __syncthreads();

      sum[tid] = fb[i];
      __syncthreads();
      fb[i] = warp_reduce<T>(sum, tid);
      __syncthreads();
      if (tid == 0)
        atomicAdd(forces_b + i, fb[i]);
      __syncthreads();
    }
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
template <typename T, typename T3, bool distributed__, bool orthorhombic_,
          int lbatch = 10>
__global__
__launch_bounds__(64) void integrate_kernel(const kernel_params dev_) {
  if (dev_.tasks[dev_.first_task + blockIdx.x].skip_task)
    return;

  const int tid =
      threadIdx.x + blockDim.x * (threadIdx.y + blockDim.y * threadIdx.z);

  __shared__ T dh_inv_[9];
  __shared__ T dh_[9];

  for (int d = tid; d < 9; d += blockDim.x * blockDim.y * blockDim.z)
    dh_inv_[d] = dev_.dh_inv_[d];

  for (int d = tid; d < 9; d += blockDim.x * blockDim.y * blockDim.z)
    dh_[d] = dev_.dh_[d];

  __shared__ smem_task_reduced<T, T3> task;
  fill_smem_task_reduced(dev_, dev_.first_task + blockIdx.x, task);

  if (tid == 0) {
    task.cube_center.z += task.lb_cube.z - dev_.grid_lower_corner_[0];
    task.cube_center.y += task.lb_cube.y - dev_.grid_lower_corner_[1];
    task.cube_center.x += task.lb_cube.x - dev_.grid_lower_corner_[2];

    if (distributed__) {
      if (task.apply_border_mask) {
        compute_window_size(
            dev_.grid_local_size_, dev_.grid_lower_corner_,
            dev_.grid_full_size_, /* also full size of the grid */
            dev_.tasks[dev_.first_task + blockIdx.x].border_mask,
            dev_.grid_border_width_, &task.window_size, &task.window_shift);
      }
    }
  }
  __syncthreads();

  __shared__ T accumulator[lbatch][64];
  __shared__ T sum[lbatch];

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
      int z2 = (z + task.cube_center.z) % dev_.grid_full_size_[0];

      if (z2 < 0)
        z2 += dev_.grid_full_size_[0];

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

      if (orthorhombic_ && !task.apply_border_mask) {
        ymin = (2 * (z + task.lb_cube.z) - 1) / 2;
        ymin *= ymin;
        kremain = task.discrete_radius * task.discrete_radius -
                  ymin * dh_[8] * dh_[8];
        ymin = ceil(-1.0e-8 - sqrt(fmax(0.0, kremain)) * dh_inv_[4]);
        ymax = 1 - ymin - task.lb_cube.y;
        ymin = ymin - task.lb_cube.y;
      }

      for (int y = ymin + threadIdx.y; y <= ymax; y += blockDim.y) {
        int y2 = (y + task.cube_center.y) % dev_.grid_full_size_[1];

        if (y2 < 0)
          y2 += dev_.grid_full_size_[1];

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

        if (orthorhombic_ && !task.apply_border_mask) {
          xmin = (2 * (y + task.lb_cube.y) - 1) / 2;
          xmin *= xmin;
          xmin =
              ceil(-1.0e-8 - sqrt(fmax(0.0, kremain - xmin * dh_[4] * dh_[4])) *
                                 dh_inv_[0]);
          xmax = 1 - xmin - task.lb_cube.x;
          xmin -= task.lb_cube.x;
        }

        for (int x = xmin + threadIdx.x; x <= xmax; x += blockDim.x) {
          int x2 = (x + task.cube_center.x) % dev_.grid_full_size_[2];

          if (x2 < 0)
            x2 += dev_.grid_full_size_[2];

          if (distributed__) {
            /* check if the point is within the window */
            if (task.apply_border_mask) {
              if ((x2 < task.window_shift.x) || (x2 > task.window_size.x)) {
                continue;
              }
            }
          }

          // I make no distinction between orthorhombic and non orthorhombic
          // cases
          T3 r3;

          if (orthorhombic_) {
            r3.x = (x + task.lb_cube.x + task.roffset.x) * dh_[0];
            r3.y = (y + task.lb_cube.y + task.roffset.y) * dh_[4];
            r3.z = (z + task.lb_cube.z + task.roffset.z) * dh_[8];
          } else {
            r3 = compute_coordinates(dh_, (x + task.lb_cube.x + task.roffset.x),
                                     (y + task.lb_cube.y + task.roffset.y),
                                     (z + task.lb_cube.z + task.roffset.z));
          }
          // check if the point is inside the sphere or not. Note that it does
          // not apply for the orthorhombic case when the full sphere is inside
          // the region of interest.

          if (distributed__) {
            if (((task.radius * task.radius) <=
                 (r3.x * r3.x + r3.y * r3.y + r3.z * r3.z)) &&
                (!orthorhombic_ || task.apply_border_mask))
              continue;
          } else {
            if (!orthorhombic_) {
              if ((task.radius * task.radius) <=
                  (r3.x * r3.x + r3.y * r3.y + r3.z * r3.z))
                continue;
            }
          }

          // read the next point assuming that reading is non blocking untill
          // the register is actually needed for computation. This is true on
          // NVIDIA hardware

          T grid_value =
              __ldg(&dev_.ptr_dev[1][(z2 * dev_.grid_local_size_[1] + y2) *
                                         dev_.grid_local_size_[2] +
                                     x2]);

          const T r3x2 = r3.x * r3.x;
          const T r3y2 = r3.y * r3.y;
          const T r3z2 = r3.z * r3.z;

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
          } break;
          case 1: {
            if (task.lp >= 3) {
              accumulator[0][tid] += grid_value * r3x2 * r3.x;
              accumulator[1][tid] += grid_value * r3x2 * r3.y;
              accumulator[2][tid] += grid_value * r3x2 * r3.z;
              accumulator[3][tid] += grid_value * r3.x * r3y2;
              accumulator[4][tid] += grid_value * r3xy * r3.z;
              accumulator[5][tid] += grid_value * r3.x * r3z2;
              accumulator[6][tid] += grid_value * r3y2 * r3.y;
              accumulator[7][tid] += grid_value * r3y2 * r3.z;
              accumulator[8][tid] += grid_value * r3.y * r3z2;
              accumulator[9][tid] += grid_value * r3z2 * r3.z;
            }
          } break;
          case 2: {
            if (task.lp >= 4) {
              accumulator[0][tid] += grid_value * r3x2 * r3x2;
              accumulator[1][tid] += grid_value * r3x2 * r3xy;
              accumulator[2][tid] += grid_value * r3x2 * r3xz;
              accumulator[3][tid] += grid_value * r3x2 * r3y2;
              accumulator[4][tid] += grid_value * r3x2 * r3yz;
              accumulator[5][tid] += grid_value * r3x2 * r3z2;
              accumulator[6][tid] += grid_value * r3xy * r3y2;
              accumulator[7][tid] += grid_value * r3xz * r3y2;
              accumulator[8][tid] += grid_value * r3xy * r3z2;
              accumulator[9][tid] += grid_value * r3xz * r3z2;
            }
          } break;
          case 3: {
            if (task.lp >= 4) {
              accumulator[0][tid] += grid_value * r3y2 * r3y2;
              accumulator[1][tid] += grid_value * r3y2 * r3yz;
              accumulator[2][tid] += grid_value * r3y2 * r3z2;
              accumulator[3][tid] += grid_value * r3yz * r3z2;
              accumulator[4][tid] += grid_value * r3z2 * r3z2;
            }
            if (task.lp >= 5) {
              accumulator[5][tid] += grid_value * r3x2 * r3x2 * r3.x;
              accumulator[6][tid] += grid_value * r3x2 * r3x2 * r3.y;
              accumulator[7][tid] += grid_value * r3x2 * r3x2 * r3.z;
              accumulator[8][tid] += grid_value * r3x2 * r3.x * r3y2;
              accumulator[9][tid] += grid_value * r3x2 * r3xy * r3.z;
            }
          } break;
          case 4: {
            accumulator[0][tid] += grid_value * r3x2 * r3.x * r3z2;
            accumulator[1][tid] += grid_value * r3x2 * r3y2 * r3.y;
            accumulator[2][tid] += grid_value * r3x2 * r3y2 * r3.z;
            accumulator[3][tid] += grid_value * r3x2 * r3.y * r3z2;
            accumulator[4][tid] += grid_value * r3x2 * r3z2 * r3.z;
            accumulator[5][tid] += grid_value * r3.x * r3y2 * r3y2;
            accumulator[6][tid] += grid_value * r3.x * r3y2 * r3yz;
            accumulator[7][tid] += grid_value * r3.x * r3y2 * r3z2;
            accumulator[8][tid] += grid_value * r3xy * r3z2 * r3.z;
            accumulator[9][tid] += grid_value * r3.x * r3z2 * r3z2;
          } break;
          case 5: {
            accumulator[0][tid] += grid_value * r3y2 * r3y2 * r3.y;
            accumulator[1][tid] += grid_value * r3y2 * r3y2 * r3.z;
            accumulator[2][tid] += grid_value * r3y2 * r3.y * r3z2;
            accumulator[3][tid] += grid_value * r3y2 * r3z2 * r3.z;
            accumulator[4][tid] += grid_value * r3.y * r3z2 * r3z2;
            accumulator[5][tid] += grid_value * r3z2 * r3z2 * r3.z;
            if (task.lp >= 6) {
              accumulator[6][tid] += grid_value * r3x2 * r3x2 * r3x2; // x^6
              accumulator[7][tid] += grid_value * r3x2 * r3x2 * r3xy; // x^5 y
              accumulator[8][tid] += grid_value * r3x2 * r3x2 * r3xz; // x^5 z
              accumulator[9][tid] += grid_value * r3x2 * r3x2 * r3y2; // x^4 y^2
            }
          } break;
          case 6: {
            accumulator[0][tid] += grid_value * r3x2 * r3x2 * r3yz; // x^4 y z
            accumulator[1][tid] += grid_value * r3x2 * r3x2 * r3z2; // x^4 z^2
            accumulator[2][tid] += grid_value * r3x2 * r3y2 * r3xy; // x^3 y^3
            accumulator[3][tid] += grid_value * r3x2 * r3y2 * r3xz; // x^3 y^2 z
            accumulator[4][tid] += grid_value * r3x2 * r3xy * r3z2; // x^3 y z^2
            accumulator[5][tid] += grid_value * r3x2 * r3z2 * r3xz; // x^3 z^3
            accumulator[6][tid] += grid_value * r3x2 * r3y2 * r3y2; // x^2 y^4
            accumulator[7][tid] += grid_value * r3x2 * r3y2 * r3yz; // x^3 y^2 z
            accumulator[8][tid] +=
                grid_value * r3x2 * r3y2 * r3z2; // x^2 y^2 z^2
            accumulator[9][tid] += grid_value * r3x2 * r3z2 * r3yz; // x^2 y z^3
          } break;
          case 7: {
            accumulator[0][tid] += grid_value * r3x2 * r3z2 * r3z2; // x^2 z^4
            accumulator[1][tid] += grid_value * r3y2 * r3y2 * r3xy; // x y^5
            accumulator[2][tid] += grid_value * r3y2 * r3y2 * r3xz; // x y^4 z
            accumulator[3][tid] += grid_value * r3y2 * r3xy * r3z2; // x y^3 z^2
            accumulator[4][tid] += grid_value * r3y2 * r3z2 * r3xz; // x y^2 z^3
            accumulator[5][tid] += grid_value * r3xy * r3z2 * r3z2; // x y z^4
            accumulator[6][tid] += grid_value * r3z2 * r3z2 * r3xz; // x z^5
            accumulator[7][tid] += grid_value * r3y2 * r3y2 * r3y2; // y^6
            accumulator[8][tid] += grid_value * r3y2 * r3y2 * r3yz; // y^5 z
            accumulator[9][tid] += grid_value * r3y2 * r3y2 * r3z2; // y^4 z^2
          } break;
          default:
            for (int ic = 0; (ic < lbatch) && ((ic + ico) < length); ic++) {
              auto &co = coset_inv[ic + ico];
              T tmp = 1.0;
              for (int po = 0; po < co.l[2]; po++)
                tmp *= r3.z;
              for (int po = 0; po < co.l[1]; po++)
                tmp *= r3.y;
              for (int po = 0; po < co.l[0]; po++)
                tmp *= r3.x;
              accumulator[ic][tid] += tmp * grid_value;
            }
            break;
          }
        }
      }
    }

    __syncthreads();

    // we know there is only 1 wavefront in each block
    // lbatch threads could reduce the values saved by all threads of the warp
    // and save results do a shuffle_down by hand
    for (int i = 0; i < min(length - ico, lbatch); i++) {
      if (tid < 32) {
        accumulator[i][tid] += accumulator[i][tid + 32];
      }
      __syncthreads();
      if (tid < 16) {
        accumulator[i][tid] += accumulator[i][tid + 16];
      }
      __syncthreads();
      if (tid < 8) {
        accumulator[i][tid] += accumulator[i][tid + 8];
      }
      __syncthreads();
      if (tid < 4) {
        accumulator[i][tid] += accumulator[i][tid + 4];
      }
      __syncthreads();
      if (tid < 2) {
        accumulator[i][tid] += accumulator[i][tid + 2];
      }
      __syncthreads();
      if (tid == 0)
        sum[i] = accumulator[i][0] + accumulator[i][1];
    }
    __syncthreads();

    for (int i = tid; i < min(length - ico, lbatch); i += lbatch)
      dev_.ptr_dev[2][dev_.tasks[dev_.first_task + blockIdx.x].coef_offset + i +
                      ico] = sum[i];
    __syncthreads();
    // if (tid == 0)
    //   printf("%.15lf\n", sum[0]);
  }
}

kernel_params
context_info::set_kernel_parameters(const int level,
                                    const smem_parameters &smem_params) {
  kernel_params params;
  params.smem_cab_offset = smem_params.smem_cab_offset();
  params.smem_alpha_offset = smem_params.smem_alpha_offset();
  params.first_task = 0;

  params.la_min_diff = smem_params.ldiffs().la_min_diff;
  params.lb_min_diff = smem_params.ldiffs().lb_min_diff;

  params.la_max_diff = smem_params.ldiffs().la_max_diff;
  params.lb_max_diff = smem_params.ldiffs().lb_max_diff;
  params.first_task = 0;
  params.tasks = this->tasks_dev.data();
  params.task_sorted_by_blocks_dev = task_sorted_by_blocks_dev.data();
  params.sorted_blocks_offset_dev = sorted_blocks_offset_dev.data();
  params.num_tasks_per_block_dev = this->num_tasks_per_block_dev_.data();
  params.block_offsets = this->block_offsets_dev.data();
  params.la_min_diff = smem_params.ldiffs().la_min_diff;
  params.lb_min_diff = smem_params.ldiffs().lb_min_diff;
  params.la_max_diff = smem_params.ldiffs().la_max_diff;
  params.lb_max_diff = smem_params.ldiffs().lb_max_diff;

  params.ptr_dev[0] = pab_block_.data();

  if (level >= 0) {
    params.ptr_dev[1] = grid_[level].data();
    for (int i = 0; i < 3; i++) {
      memcpy(params.dh_, grid_[level].dh(), 9 * sizeof(double));
      memcpy(params.dh_inv_, grid_[level].dh_inv(), 9 * sizeof(double));
      params.grid_full_size_[i] = grid_[level].full_size(i);
      params.grid_local_size_[i] = grid_[level].local_size(i);
      params.grid_lower_corner_[i] = grid_[level].lower_corner(i);
      params.grid_border_width_[i] = grid_[level].border_width(i);
    }
    params.first_task = first_task_per_level_[level];
  }
  params.ptr_dev[2] = this->coef_dev_.data();
  params.ptr_dev[3] = hab_block_.data();
  params.ptr_dev[4] = forces_.data();
  params.ptr_dev[5] = virial_.data();
  params.sphi_dev = this->sphi_dev.data();
  return params;
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
    if (grid_[level].is_orthorhombic()) {
      integrate_kernel<double, double3, true, true, 10>
          <<<number_of_tasks_per_level_[level], threads_per_block, 0,
             level_streams[level]>>>(params);
    } else {
      integrate_kernel<double, double3, true, false, 10>
          <<<number_of_tasks_per_level_[level], threads_per_block, 0,
             level_streams[level]>>>(params);
    }
  } else {
    if (grid_[level].is_orthorhombic()) {
      integrate_kernel<double, double3, false, true, 10>
          <<<number_of_tasks_per_level_[level], threads_per_block, 0,
             level_streams[level]>>>(params);
    } else {
      integrate_kernel<double, double3, false, false, 10>
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
    compute_hab_v2<double, double3, false, false>
        <<<this->nblocks, threads_per_block, smem_params.smem_per_block(),
           this->main_stream>>>(params, this->ntasks);
    return;
  }

  if (!compute_tau && calculate_forces) {
    compute_hab_v2<double, double3, false, true>
        <<<this->nblocks, threads_per_block, smem_params.smem_per_block(),
           this->main_stream>>>(params, this->ntasks);
    return;
  }

  if (compute_tau && calculate_forces) {
    compute_hab_v2<double, double3, true, true>
        <<<this->nblocks, threads_per_block, smem_params.smem_per_block(),
           this->main_stream>>>(params, this->ntasks);
  }

  if (compute_tau && !calculate_forces) {
    compute_hab_v2<double, double3, true, false>
        <<<this->nblocks, threads_per_block, smem_params.smem_per_block(),
           this->main_stream>>>(params, this->ntasks);
  }
}
};     // namespace rocm_backend
#endif // __GRID_ROCM
