/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2021 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

/*
 * Authors :
 - Dr Mathieu Taillefumier (ETH Zurich / CSCS)
 - Advanced Micro Devices, Inc.
*/

#ifdef __GRID_HIP

#include <algorithm>
#include <assert.h>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <hip/hip_runtime.h>

#include "grid_hip_internal_header.h"
#include "grid_hip_prepare_pab.h"

namespace rocm_backend {
/*******************************************************************************
 * \brief Decontracts the subblock, going from spherical to cartesian harmonics.
 ******************************************************************************/
template <typename T, bool IS_FUNC_AB>
__device__ __inline__ void block_to_cab(const kernel_params &params,
                                        const smem_task<T> &task, T *cab) {

  // The spherical index runs over angular momentum and then over contractions.
  // The cartesian index runs over exponents and then over angular momentum.

  // This is a T matrix product. Since the pab block can be quite large the
  // two products are fused to conserve shared memory.
  const int tid =
      threadIdx.x + blockDim.x * (threadIdx.y + blockDim.y * threadIdx.z);
  for (int i = 0; i < task.nsgf_setb; i++) {
    for (int j = 0; j < task.nsgf_seta; j++) {
      T block_val;
      if (task.block_transposed) {
        block_val = task.pab_block[j * task.nsgfb + i] * task.off_diag_twice;
      } else {
        block_val = task.pab_block[i * task.nsgfa + j] * task.off_diag_twice;
      }

      // fast path for common case
      // const int jco_start = task.first_cosetb;
      for (int jco = task.first_cosetb + tid / 8; jco < task.ncosetb;
           jco += 8) {
        const T sphib = task.sphib[i * task.maxcob + jco];
        for (int ico = task.first_coseta + (tid % 8); ico < task.ncoseta;
             ico += 8) {
          const T sphia = task.sphia[j * task.maxcoa + ico];
          const T pab_val = block_val * sphia * sphib;
          if (IS_FUNC_AB) {
            cab[jco * task.ncoseta + ico] += pab_val;
          } else {
            const auto a = coset_inv[ico];
            const auto b = coset_inv[jco];
            prepare_pab(params.func, a, b, task.zeta, task.zetb, pab_val,
                        task.n1, cab);
          }
        }
      }
      __syncthreads();
    }
  }
}

template <typename T, bool IS_FUNC_AB>
__global__ void calculate_coefficients(const kernel_params dev_) {
  __shared__ smem_task<T> task;
  if (dev_.tasks[dev_.first_task + blockIdx.x].skip_task)
    return;

  const int tid =
      threadIdx.x + blockDim.x * (threadIdx.y + blockDim.y * threadIdx.z);

  fill_smem_task_coef(dev_, dev_.first_task + blockIdx.x, task);
  extern __shared__ T shared_memory[];
  T *smem_cab = &shared_memory[dev_.smem_cab_offset];
  T *smem_alpha = &shared_memory[dev_.smem_alpha_offset];
  T *coef_ =
      &dev_.ptr_dev[2][dev_.tasks[dev_.first_task + blockIdx.x].coef_offset];
  T *coefs_ = &shared_memory[0];

  for (int z = tid; z < task.n1 * task.n2;
       z += blockDim.x * blockDim.y * blockDim.z)
    smem_cab[z] = 0.0;

  __syncthreads();
  block_to_cab<T, IS_FUNC_AB>(dev_, task, smem_cab);
  __syncthreads();
  compute_alpha(dev_, task, smem_alpha);
  __syncthreads();
  cab_to_cxyz(dev_, task, smem_alpha, smem_cab, coefs_);
  __syncthreads();

  for (int z = tid; z < ncoset(task.lp);
       z += blockDim.x * blockDim.y * blockDim.z)
    coef_[z] = coefs_[z];
}

/*
  \brief compute the real space representation of an operator expressed in the
  gaussian basis

  this kernel does the following operation

  n_{ijk} = \f[
  \sum_{p,\alpha,\beta,\gamma} C^p_{\alpha\beta\gamma} X_{\alpha,i} Y_{\beta, j}
  Z_{\gamma, k} \exp\left(- \eta (r_{ijk} - r_c) ^ 2\right)
  ]

  where $X_{\alpha,i}, Y_{\beta, j}, Z_{\gamma, k}$ are polynomials of degree
  $\alpha,\beta,\gamma$ and $r_{ijk}% a point (in cartesian coordinates) on a 3D
  grid. C^p_{\alpha\beta\gamma} are also constrained such that 0 <= \alpha +
  \beta + \gamma <= lmax. It means in practice that we need store (lmax + 1) *
  (lamax + 2) * (lmax + 3) / 6 coefficients all the other coefficients are zero

  to reduce computation, a spherical cutoff is applied such that all points
  $|r_{ijk} - r_c| > radius$ are not computed. The sum over p extends over all
  relevant pairs of gaussians (which are called task in the code).

  the kernel computes the polynomials and the gaussian then sums the result
  back to the grid.

  the coefficients $C^p_{\alpha\beta\gamma}$ are computed by
  calculate_coefficients. We only keep the non zero elements to same memory.
*/

template <typename T, typename T3, bool distributed__, bool orthorhombic_>
__global__
__launch_bounds__(64) void collocate_kernel(const kernel_params dev_) {
  // Copy task from global to shared memory and precompute some stuff.
  __shared__ smem_task_reduced<T, T3> task;

  if (dev_.tasks[dev_.first_task + blockIdx.x].skip_task)
    return;

  fill_smem_task_reduced(dev_, dev_.first_task + blockIdx.x, task);

  const int tid =
      threadIdx.x + blockDim.x * (threadIdx.y + blockDim.y * threadIdx.z);

  //  Alloc shared memory.
  extern __shared__ T coefs_[];

  T *coef_ =
      &dev_.ptr_dev[2][dev_.tasks[dev_.first_task + blockIdx.x].coef_offset];
  __shared__ T dh_[9], dh_inv_[9];

  // matrix from lattice coordinates to cartesian coordinates
  for (int i = tid; i < 9; i += blockDim.x * blockDim.y * blockDim.z) {
    dh_[i] = dev_.dh_[i];
  }

  // matrix from  cartesian coordinates to lattice coordinates.
  for (int i = tid; i < 9; i += blockDim.x * blockDim.y * blockDim.z) {
    dh_inv_[i] = dev_.dh_inv_[i];
  }

  __syncthreads();
  for (int i = tid; i < ncoset(task.lp);
       i += blockDim.x * blockDim.y * blockDim.z)
    coefs_[i] = coef_[i];

  if (tid == 0) {
    // the cube center is initialy expressed in lattice coordinates but we
    // always do something like this. x = x + lower_corner + cube_center (+
    // roffset) - grid_lower_corner so shift the cube center already
    task.cube_center.z += task.lb_cube.z - dev_.grid_lower_corner_[0];
    task.cube_center.y += task.lb_cube.y - dev_.grid_lower_corner_[1];
    task.cube_center.x += task.lb_cube.x - dev_.grid_lower_corner_[2];

    if (distributed__) {
      if (task.apply_border_mask) {
        compute_window_size(
            dev_.grid_local_size_, dev_.grid_lower_corner_,
            dev_.grid_full_size_, // also full size of the grid
            dev_.tasks[dev_.first_task + blockIdx.x].border_mask,
            dev_.grid_border_width_, &task.window_size, &task.window_shift);
      }
    }
  }
  __syncthreads();

  for (int z = threadIdx.z; z < task.cube_size.z; z += blockDim.z) {
    int z2 = (z + task.cube_center.z) % dev_.grid_full_size_[0];

    if (z2 < 0)
      z2 += dev_.grid_full_size_[0];

    if (distributed__) {
      // check if the point is within the window
      if (task.apply_border_mask) {
        // this test is only relevant when the grid is split over several mpi
        // ranks. in that case we take only the points contributing to local
        // part of the grid.
        if ((z2 < task.window_shift.z) || (z2 > task.window_size.z)) {
          continue;
        }
      }
    }

    // compute the coordinates of the point in atomic coordinates
    T kremain;
    short int ymin = 0;
    short int ymax = task.cube_size.y - 1;

    if (orthorhombic_ && !task.apply_border_mask) {
      ymin = (2 * (z + task.lb_cube.z) - 1) / 2;
      ymin *= ymin;
      kremain = task.discrete_radius * task.discrete_radius -
                ((T)ymin) * dh_[8] * dh_[8];
      ymin = ceil(-1.0e-8 - sqrt(fmax(0.0, kremain)) * dh_inv_[4]);
      ymax = 1 - ymin - task.lb_cube.y;
      ymin = ymin - task.lb_cube.y;
    }

    for (int y = ymin + threadIdx.y; y <= ymax; y += blockDim.y) {
      int y2 = (y + task.cube_center.y) % dev_.grid_full_size_[1];
      if (y2 < 0)
        y2 += dev_.grid_full_size_[1];

      if (distributed__) {
        if (task.apply_border_mask) {
          // check if the point is within the window
          if ((y2 < task.window_shift.y) || (y2 > task.window_size.y)) {
            continue;
          }
        }
      }

      short int xmin = 0;
      short int xmax = task.cube_size.x - 1;
      if (orthorhombic_ && !task.apply_border_mask) {
        xmin = (2 * (y + task.lb_cube.y) - 1) / 2;
        xmin *= xmin;
        xmin =
            ceil(-1.0e-8 - sqrt(fmax(0.0, kremain - xmin * dh_[4] * dh_[4])) *
                               dh_inv_[0]);
        xmax = 1 - xmin - task.lb_cube.x;
        xmin = xmin - task.lb_cube.x;
      }

      for (int x = xmin + threadIdx.x; x <= xmax; x += blockDim.x) {
        int x2 = (x + task.cube_center.x) % dev_.grid_full_size_[2];

        if (x2 < 0)
          x2 += dev_.grid_full_size_[2];

        if (distributed__) {
          if (task.apply_border_mask) {
            // check if the point is within the window (only true or false
            // when using mpi) otherwise MPI=1 always true
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

        if (distributed__) {
          // check if the point is inside the sphere or not. Note that it does
          // not apply for the orthorhombic case when the full sphere is inside
          // the region of interest.

          if (((task.radius * task.radius) <=
               (r3.x * r3.x + r3.y * r3.y + r3.z * r3.z)) &&
              (!orthorhombic_ || task.apply_border_mask))
            continue;
        } else {
          // we do not need to do this test for the orthorhombic case
          if ((!orthorhombic_) && ((task.radius * task.radius) <=
                                   (r3.x * r3.x + r3.y * r3.y + r3.z * r3.z)))
            continue;
        }

        // allow computation of the address in parallel to starting the
        // computations
        // T *grid_elem =
        //     dev_.ptr_dev[1] +
        //     (z2 * dev_.grid_local_size_[1] + y2) * dev_.grid_local_size_[2] +
        //     x2;

        T res = coefs_[0];

        if (task.lp >= 1) {
          res += coefs_[1] * r3.x;
          res += coefs_[2] * r3.y;
          res += coefs_[3] * r3.z;
        }
        const T r3xy = r3.x * r3.y;
        const T r3xz = r3.x * r3.z;
        const T r3yz = r3.y * r3.z;
        const T r3x2 = r3.x * r3.x;
        const T r3y2 = r3.y * r3.y;
        const T r3z2 = r3.z * r3.z;

        if (task.lp >= 2) {
          res += coefs_[4] * r3x2;
          res += coefs_[5] * r3xy;
          res += coefs_[6] * r3xz;
          res += coefs_[7] * r3y2;
          res += coefs_[8] * r3yz;
          res += coefs_[9] * r3z2;
        }

        if (task.lp >= 3) {
          res += coefs_[10] * r3x2 * r3.x;
          res += coefs_[11] * r3x2 * r3.y;
          res += coefs_[12] * r3x2 * r3.z;
          res += coefs_[13] * r3.x * r3y2;
          res += coefs_[14] * r3xy * r3.z;
          res += coefs_[15] * r3.x * r3z2;
          res += coefs_[16] * r3y2 * r3.y;
          res += coefs_[17] * r3y2 * r3.z;
          res += coefs_[18] * r3.y * r3z2;
          res += coefs_[19] * r3z2 * r3.z;
        }

        if (task.lp >= 4) {
          res += coefs_[20] * r3x2 * r3x2;
          res += coefs_[21] * r3x2 * r3xy;
          res += coefs_[22] * r3x2 * r3xz;
          res += coefs_[23] * r3x2 * r3y2;
          res += coefs_[24] * r3x2 * r3yz;
          res += coefs_[25] * r3x2 * r3z2;
          res += coefs_[26] * r3xy * r3y2;
          res += coefs_[27] * r3xz * r3y2;
          res += coefs_[28] * r3xy * r3z2;
          res += coefs_[29] * r3xz * r3z2;
          res += coefs_[30] * r3y2 * r3y2;
          res += coefs_[31] * r3y2 * r3yz;
          res += coefs_[32] * r3y2 * r3z2;
          res += coefs_[33] * r3yz * r3z2;
          res += coefs_[34] * r3z2 * r3z2;
        }

        if (task.lp > 4) {
          for (int ic = 35; ic < ncoset(task.lp); ic++) {
            auto &co = coset_inv[ic];
            T tmp = coefs_[ic];
            for (int po = 0; po < co.l[2]; po++)
              tmp *= r3.z;
            for (int po = 0; po < co.l[1]; po++)
              tmp *= r3.y;
            for (int po = 0; po < co.l[0]; po++)
              tmp *= r3.x;
            res += tmp;
          }
        }

        atomicAdd(
            dev_.ptr_dev[1] +
                (z2 * dev_.grid_local_size_[1] + y2) *
                    dev_.grid_local_size_[2] +
                x2,
            res * exp(-(r3.x * r3.x + r3.y * r3.y + r3.z * r3.z) * task.zetp));
      }
    }
  }
  __syncthreads();
}
/*******************************************************************************
 * \brief Launches the Cuda kernel that collocates all tasks of one grid level.
 ******************************************************************************/
void context_info::collocate_one_grid_level(const int level,
                                            const enum grid_func func,
                                            int *lp_diff) {

  if (number_of_tasks_per_level_[level] == 0)
    return;

  // Compute max angular momentum.
  const ldiffs_value ldiffs = prepare_get_ldiffs(func);
  smem_parameters smem_params(ldiffs, lmax());

  *lp_diff = smem_params.lp_diff();
  init_constant_memory();

  // kernel parameters
  kernel_params params = set_kernel_parameters(level, smem_params);
  params.func = func;

  // Launch !
  const dim3 threads_per_block(4, 4, 4);

  if (func == GRID_FUNC_AB) {
    calculate_coefficients<double, true>
        <<<number_of_tasks_per_level_[level], threads_per_block,
           smem_params.smem_per_block(), level_streams[level]>>>(params);
  } else {
    calculate_coefficients<double, false>
        <<<number_of_tasks_per_level_[level], threads_per_block,
           smem_params.smem_per_block(), level_streams[level]>>>(params);
  }

  if (grid_[level].is_distributed()) {
    if (grid_[level].is_orthorhombic())
      collocate_kernel<double, double3, true, true>
          <<<number_of_tasks_per_level_[level], threads_per_block,
             smem_params.cxyz_len() * sizeof(double), level_streams[level]>>>(
              params);
    else
      collocate_kernel<double, double3, true, false>
          <<<number_of_tasks_per_level_[level], threads_per_block,
             smem_params.cxyz_len() * sizeof(double), level_streams[level]>>>(
              params);
  } else {
    if (grid_[level].is_orthorhombic())
      collocate_kernel<double, double3, false, true>
          <<<number_of_tasks_per_level_[level], threads_per_block,
             smem_params.cxyz_len() * sizeof(double), level_streams[level]>>>(
              params);
    else
      collocate_kernel<double, double3, false, false>
          <<<number_of_tasks_per_level_[level], threads_per_block,
             smem_params.cxyz_len() * sizeof(double), level_streams[level]>>>(
              params);
  }
}
} // namespace rocm_backend
#endif // __GRID_ROCM
