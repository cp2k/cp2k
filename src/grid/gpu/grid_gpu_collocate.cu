/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

/*
 * Authors :
 - Mathieu Taillefumier (ETH Zurich / CSCS)
 - Advanced Micro Devices, Inc.
 - Ole Schuett
*/

#include <algorithm>
#include <assert.h>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "grid_gpu_internal_header.h"
#include "grid_gpu_prepare_pab.h"
#include "grid_gpu_separable.h"

#if defined(_OMP_H)
#error "OpenMP should not be used in .cu files to accommodate HIP."
#endif

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
  const int tid = thread_global_index();

  for (int jco = task.first_cosetb + tid / 16; jco < task.ncosetb; jco += 4) {
    for (int ico = task.first_coseta + (tid % 16); ico < task.ncoseta;
         ico += 16) {
      T pab_val = 0.0;
      for (int i = 0; i < task.nsgf_setb; i++) {
        const T sphib = task.sphib[i * task.maxcob + jco];
        T tmp_val = 0.0;
        for (int j = 0; j < task.nsgf_seta; j++) {
          T block_val;
          if (task.block_transposed) {
            block_val =
                task.pab_block[j * task.nsgfb + i] * task.off_diag_twice;
          } else {
            block_val =
                task.pab_block[i * task.nsgfa + j] * task.off_diag_twice;
          }
          tmp_val += block_val * task.sphia[j * task.maxcoa + ico];
        }
        pab_val += tmp_val * sphib;
      }

      if (IS_FUNC_AB) {
        cab[jco * task.ncoseta + ico] += pab_val;
      } else {
        const auto a = coset_inv[ico];
        const auto b = coset_inv[jco];
        // prepare_pab adds to the neighbours up/down(i, a) x up/down(i, b),
        // which other (ico, jco) pairs of the block also update: the
        // atomicAdd in prep_term is required.
        prepare_pab(params.func, a, b, task.zeta, task.zetb, pab_val, task.n1,
                    cab);
      }
    }
  }
}

template <typename T, bool IS_FUNC_AB>
__global__
__launch_bounds__(64) void calculate_coefficients(const kernel_params dev_) {
  const int number_of_tasks = dev_.num_tasks_per_block_dev[block_index()];

  if (number_of_tasks == 0)
    return;

  const int tid = thread_global_index();
  const int offset = dev_.sorted_blocks_offset_dev[block_index()];
  T *smem_cab = reinterpret_cast<double *>(
      __builtin_assume_aligned(allocate_workspace<T>(dev_), 32));

  for (int tk = 0; tk < number_of_tasks; tk++) {
    __shared__ smem_task<T> task;
    const int task_id = dev_.task_sorted_by_blocks_dev[offset + tk];
    if (dev_.tasks[task_id].skip_task)
      continue;
    fill_smem_task_coef(dev_, task_id, task);

    T *__restrict__ coef_ =
        &dev_.buffers_dev.coef[dev_.tasks[task_id].coef_offset];

    for (int z = tid; z < task.n1 * task.n2;
         z += blockDim.x * blockDim.y * blockDim.z)
      smem_cab[z] = 0.0;
    __syncthreads();

    block_to_cab<T, IS_FUNC_AB>(dev_, task, smem_cab);
    __syncthreads();

    cab_to_cxyz_factorized(task, smem_cab, coef_);
  }
}

/*******************************************************************************
 * \brief Same as calculate_coefficients, with the spherical to cartesian
 *        transformation of block_to_cab done in two steps,
 *
 *   T[i][ico]   = sum_j pab[i][j] sphia[j][ico]          (nsgfb nsgfa ncoa)
 *   P[jco][ico] = sum_i sphib[i][jco] T[i][ico]          (ncob nsgfb ncoa)
 *
 *        instead of the fused loop of cost ncoa ncob nsgfb (nsgfa + 1). The
 *        work arrays T [max_nsgf_set][ld] and P [ld][ld] live in the dynamic
 *        shared memory, see sphi_work_size.
 ******************************************************************************/
template <typename T, bool IS_FUNC_AB>
__global__ __launch_bounds__(64) void calculate_coefficients_two_steps(
    const kernel_params dev_, const int ld, const int max_nsgf_set) {
  extern __shared__ T shared_memory[];
  const int number_of_tasks = dev_.num_tasks_per_block_dev[block_index()];

  if (number_of_tasks == 0)
    return;

  T *smem_t = &shared_memory[0];          // [max_nsgf_set][ld]
  T *smem_p = smem_t + max_nsgf_set * ld; // [ld][ld], P[jco][ico]
  const int tid = thread_global_index();
  const int nthreads = blockDim.x * blockDim.y * blockDim.z;
  const int offset = dev_.sorted_blocks_offset_dev[block_index()];
  T *smem_cab = reinterpret_cast<double *>(
      __builtin_assume_aligned(allocate_workspace<T>(dev_), 32));

  for (int tk = 0; tk < number_of_tasks; tk++) {
    __shared__ smem_task<T> task;
    const int task_id = dev_.task_sorted_by_blocks_dev[offset + tk];
    if (dev_.tasks[task_id].skip_task)
      continue;
    fill_smem_task_coef(dev_, task_id, task);

    T *__restrict__ coef_ =
        &dev_.buffers_dev.coef[dev_.tasks[task_id].coef_offset];

    for (int z = tid; z < task.n1 * task.n2; z += nthreads)
      smem_cab[z] = 0.0;

    // step 1: T[i][ico] = sum_j pab[i][j] sphia[j][ico]
    const int nco_a = task.ncoseta - task.first_coseta;
    const int nco_b = task.ncosetb - task.first_cosetb;
    for (int k = tid; k < task.nsgf_setb * nco_a; k += nthreads) {
      const int i = k / nco_a;
      const int ico = task.first_coseta + (k - i * nco_a);
      T val = 0.0;
      for (int j = 0; j < task.nsgf_seta; j++) {
        const T block_val = task.block_transposed
                                ? task.pab_block[j * task.nsgfb + i]
                                : task.pab_block[i * task.nsgfa + j];
        val += block_val * task.sphia[j * task.maxcoa + ico];
      }
      smem_t[i * ld + ico] = val * task.off_diag_twice;
    }
    __syncthreads();

    // step 2: P[jco][ico] = sum_i sphib[i][jco] T[i][ico]
    for (int k = tid; k < nco_b * nco_a; k += nthreads) {
      const int jj = k / nco_a;
      const int jco = task.first_cosetb + jj;
      const int ico = task.first_coseta + (k - jj * nco_a);
      T val = 0.0;
      for (int i = 0; i < task.nsgf_setb; i++)
        val += task.sphib[i * task.maxcob + jco] * smem_t[i * ld + ico];
      if (IS_FUNC_AB) {
        smem_cab[jco * task.ncoseta + ico] += val;
      } else {
        smem_p[jco * ld + ico] = val;
      }
    }
    __syncthreads();

    if (!IS_FUNC_AB) {
      // prepare_pab adds (atomically) to several elements of cab
      for (int k = tid; k < nco_b * nco_a; k += nthreads) {
        const int jj = k / nco_a;
        const int jco = task.first_cosetb + jj;
        const int ico = task.first_coseta + (k - jj * nco_a);
        prepare_pab(dev_.func, coset_inv[ico], coset_inv[jco], task.zeta,
                    task.zetb, smem_p[jco * ld + ico], task.n1, smem_cab);
      }
      __syncthreads();
    }

    cab_to_cxyz_factorized(task, smem_cab, coef_);
    __syncthreads();
  }
}

/*******************************************************************************
 * Separable (plane streaming) collocate kernel. This is
 * integrate_separable_kernel (grid_gpu_integrate.cu) run backwards.
 *
 * Axis aligned orthorhombic cell (axis_aligned_ = true):
 *
 *   sum_{abc} C[abc] x^a y^b z^c exp(-zetp r^2)
 *     = sum_a [x^a e^{-zetp x^2}] sum_b [y^b e^{-zetp y^2}]
 *                                 sum_c [z^c e^{-zetp z^2}] C[abc]
 *
 * Other cells (axis_aligned_ = false): the
 * cartesian coefficients are first expressed in lattice monomials u^i v^j w^k
 * (u, v, w lattice coordinates of the point relative to the gaussian center),
 *
 *   C_ijk[(i,j,k)] = sum_{a+b+c=i+j+k} M[(a,b,c)][(i,j,k)] C[(a,b,c)]
 *
 * with the matrix M of build_lattice_to_cartesian_matrix (the transpose of the
 * integrate transformation). With the metric g_ij = a_i . a_j of the lattice
 * vectors, the gaussian is factored as in the dgemm backend,
 *
 *   exp(-zetp r^2) = exp(-zetp g00 u^2) exp(-zetp g11 v^2) exp(-zetp g22 w^2)
 *                    exp(-2 zetp g12 v w) exp(-2 zetp u (g01 v + g02 w))
 *
 * The first three factors go into the 1D tables, the v w term is applied once
 * per row in step 2 and the last one by a multiplicative recurrence along x,
 * so that no exponential is evaluated per grid point.
 *
 * For each z plane of the cube:
 *
 *   D[a][b] = sum_c C[coset(a, b, c)] pz[c]               a + b <= lp (step 1)
 *   E[y][a] = sum_b D[a][b] py[b][y]                                  (step 2)
 *   grid[z][y][x] += sum_a E[y][a] px[a][x]                           (step 3)
 *
 * so that a grid point costs lp + 1 FMAs (Horner in x) instead of ncoset(lp).
 * Steps 1 and 2 are O(lp^3) and O(ny lp^2) per plane. The set of grid points
 * visited is the same as for the original collocate kernel: same window masks,
 * same y / x trimming in the orthogonal case, same sphere test.
 *
 * Dynamic shared memory layout (doubles), mc = max cube size, ls = lp_max + 1:
 *   ex[mc], ey[mc], fy[mc], D[ls][ls], C[ncoset(lp_max)], E[erows][ls]
 *
 * The powers of y are computed on the fly in step 2 and the rows of a plane
 * are processed by chunks of erows rows, so that the shared memory is bounded
 * for any cube size (see collocate_separable_erows), the x coordinates are
 * computed on the fly. fy is only used when axis_aligned_ = false.
 ******************************************************************************/

template <bool distributed__, bool axis_aligned_>
__global__ __launch_bounds__(64) void collocate_separable_kernel(
    const kernel_params dev_, const int mc, const int lp_max, const int erows,
    const double *__restrict__ M) {
  // the orthogonal code paths are only used for an axis aligned cell (see
  // grid_info::check_orthogonality)
  constexpr bool orthogonal_ = axis_aligned_;

  const task_info &glb_task = dev_.tasks[dev_.first_task + block_index()];
  if (glb_task.skip_task)
    return;

  const int tid = threadIdx.x;
  const int ls = lp_max + 1;

  extern __shared__ double sep_smem[];
  double *ex = sep_smem;
  double *ey = ex + mc;
  double *fy = ey + mc;
  double *D = fy + mc;
  double *C = D + ls * ls;
  double *E = C + ncoset(lp_max); // [erows][ls]

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

  // metric of the lattice, |r|^2 = sum_ij g_ij u_i u_j
  const double *a1 = &dev_.dh_[0], *a2 = &dev_.dh_[3], *a3 = &dev_.dh_[6];
  const double g00 = a1[0] * a1[0] + a1[1] * a1[1] + a1[2] * a1[2];
  const double g11 = a2[0] * a2[0] + a2[1] * a2[1] + a2[2] * a2[2];
  const double g22 = a3[0] * a3[0] + a3[1] * a3[1] + a3[2] * a3[2];
  const double g01 = a1[0] * a2[0] + a1[1] * a2[1] + a1[2] * a2[2];
  const double g02 = a1[0] * a3[0] + a1[1] * a3[1] + a1[2] * a3[2];
  const double g12 = a2[0] * a3[0] + a2[1] * a3[1] + a2[2] * a3[2];

  const double *__restrict__ coef =
      &dev_.buffers_dev.coef[glb_task.coef_offset];

  // coefficients, cartesian or in lattice monomials
  for (int i = tid; i < ncoef; i += 64) {
    if (axis_aligned_) {
      C[i] = coef[i];
    } else {
      const orbital &co = coset_inv[i];
      const int l = co.l[0] + co.l[1] + co.l[2];
      const int base = (l > 0) ? ncoset(l - 1) : 0;
      const int nl = ncoset(l) - base;
      // column (i - base) of block l
      const double *__restrict__ Ml = M + sep_m_offset[l] + (i - base);
      double val = 0.0;
      for (int k = 0; k < nl; k++)
        val += Ml[k * nl] * coef[base + k];
      C[i] = val;
    }
  }

  // 1D tables of the diagonal part of the gaussian. The coordinates (cartesian
  // for an axis aligned cell, lattice otherwise) are computed on the fly.
  for (int x = tid; x < task.cube_size.x; x += 64) {
    const double u = x + task.lb_cube.x + task.roffset.x;
    const double r = axis_aligned_ ? u * dev_.dh_[0] : u;
    ex[x] = exp(-zetp * (axis_aligned_ ? 1.0 : g00) * r * r);
  }
  for (int y = tid; y < task.cube_size.y; y += 64) {
    const double v = y + task.lb_cube.y + task.roffset.y;
    const double r = axis_aligned_ ? v * dev_.dh_[4] : v;
    ey[y] = exp(-zetp * (axis_aligned_ ? r * r : g11 * v * v));
  }
  __syncthreads();

  const int npairs = (lp + 1) * (lp + 2) / 2;
  const int row_group = tid / GRID_SEP_COLLOCATE_ROW_LANES;
  const int lane = tid % GRID_SEP_COLLOCATE_ROW_LANES;
  constexpr int nrow_groups = 64 / GRID_SEP_COLLOCATE_ROW_LANES;

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
    const double ez = exp(-zetp * (axis_aligned_ ? rz * rz : g22 * w * w));
    const int ny = ymax - ymin + 1;

    // v w cross term of the gaussian, one value per row
    if (!axis_aligned_) {
      for (int yy = tid; yy < ny; yy += 64) {
        const double v = ymin + yy + task.lb_cube.y + task.roffset.y;
        fy[yy] = exp(-2.0 * zetp * g12 * v * w);
      }
    }

    // step 1: contract c (Horner in rz). D is only read by step 2 so the
    // barrier ending step 2 of the previous plane protects it.
    for (int ip = tid; ip < npairs; ip += 64) {
      const int a = sep_pair_ab[ip].x;
      const int b = sep_pair_ab[ip].y;
      double d = 0.0;
      for (int c = lp - a - b; c >= 0; c--)
        d = d * rz + C[coset(a, b, c)];
      D[a * ls + b] = d * ez;
    }
    __syncthreads();

    // steps 2 and 3 by chunks of erows rows, so that E fits in the shared
    // memory for any cube size. The trip count is uniform over the block.
    for (int y0 = ymin; y0 <= ymax; y0 += erows) {
      const int y1 = min(ymax, y0 + erows - 1);
      const int nyc = y1 - y0 + 1;

      // step 2: contract b, with the powers of y computed on the fly
      for (int k = tid; k < nyc * (lp + 1); k += 64) {
        const int yy = k / (lp + 1);
        const int a = k - yy * (lp + 1);
        const int y = y0 + yy;
        const double v = y + task.lb_cube.y + task.roffset.y;
        const double r = axis_aligned_ ? v * dev_.dh_[4] : v;
        double p = ey[y];
        double e = 0.0;
        for (int b = 0; b <= lp - a; b++, p *= r)
          e += D[a * ls + b] * p;
        E[yy * ls + a] = axis_aligned_ ? e : e * fy[y - ymin];
      }
      __syncthreads();

      // step 3: one group of GRID_SEP_COLLOCATE_ROW_LANES lanes per row y,
      // Horner in x.
      for (int y = y0 + row_group; y <= y1; y += nrow_groups) {
        const int y2 =
            wrap_grid_index(y + task.cube_center.y, dev_.grid_full_size_.y);
        if (distributed__ && task.apply_border_mask &&
            ((y2 < task.window_shift.y) || (y2 > task.window_size.y)))
          continue;

        int xmin = 0, xmax = task.cube_size.x - 1;
        if (trim)
          calculate_xmin_xmax_boundaries<double, double3>(task, y, kremain,
                                                          xmin, xmax);
        const double v = y + task.lb_cube.y + task.roffset.y;
        const double *__restrict__ Ey = &E[(y - y0) * ls];
        double *__restrict__ row =
            &dev_.buffers_dev.grid[(z2 * dev_.grid_local_size_.y + y2) *
                                   dev_.grid_local_size_.x];

        // u (v, w) cross term of the gaussian, exp(-2 zetp u s) with
        // s = g01 v + g02 w, by recurrence along the x of this lane.
        double t = 1.0, q = 1.0;
        if (!axis_aligned_) {
          const double s = g01 * v + g02 * w;
          const double u0 = xmin + lane + task.lb_cube.x + task.roffset.x;
          t = exp(-2.0 * zetp * s * u0);
          q = exp(-2.0 * zetp * s * GRID_SEP_COLLOCATE_ROW_LANES);
        }

        // periodic index of x, updated incrementally
        const int nx_full = dev_.grid_full_size_.x;
        int x2 = wrap_grid_index(xmin + lane + task.cube_center.x, nx_full);
        for (int x = xmin + lane; x <= xmax; x += GRID_SEP_COLLOCATE_ROW_LANES,
                 t *= q,
                 x2 = sep_advance_wrapped_index<GRID_SEP_COLLOCATE_ROW_LANES>(
                     x2, nx_full)) {
          if (distributed__ && task.apply_border_mask &&
              ((x2 < task.window_shift.x) || (x2 > task.window_size.x)))
            continue;

          const double u = x + task.lb_cube.x + task.roffset.x;
          const double r = axis_aligned_ ? u * dev_.dh_[0] : u;
          double gauss;
          if (axis_aligned_) {
            // the sphere test is never applied in this case (orthogonal_)
            gauss = ex[x];
          } else {
            const double3 r3 = compute_coordinates(dev_.dh_, r, v, w);
            const double r2 = r3.x * r3.x + r3.y * r3.y + r3.z * r3.z;
            // same sphere test as the original collocate kernel
            const bool outside = (r2max <= r2);
            if (distributed__) {
              if (outside && (!orthogonal_ || task.apply_border_mask))
                continue;
            } else if (!orthogonal_ && outside) {
              continue;
            }
            gauss = ex[x] * t;
          }
          if (axis_aligned_ && distributed__ && task.apply_border_mask) {
            const double ry = v * dev_.dh_[4];
            const double r2 = r * r + ry * ry + rz * rz;
            if (r2max <= r2)
              continue;
          }

          double val = Ey[lp];
          for (int a = lp - 1; a >= 0; a--)
            val = val * r + Ey[a];
          atomicAdd(&row[x2], val * gauss);
        }
      }
      // E is overwritten by the next chunk
      if (y1 < ymax)
        __syncthreads();
    }
  }
}

// shared memory of collocate_separable_kernel for erows rows of E per chunk
static size_t collocate_separable_smem_size(const int mc, const int lp_max,
                                            const int erows) {
  const int ls = lp_max + 1;
  return sizeof(double) * (3 * mc + ls * ls + ncoset(lp_max) + erows * ls);
}

// rows of E per chunk: all the rows of the cube when they fit in smem_max
static int collocate_separable_erows(const int mc, const int lp_max,
                                     const size_t smem_max) {
  const size_t fixed = collocate_separable_smem_size(mc, lp_max, 0);
  const size_t per_row = sizeof(double) * (lp_max + 1);
  if (fixed + per_row > smem_max) {
    fprintf(stderr,
            "ERROR: cube too large for collocate_separable_kernel (%d "
            "points, lp %d).\n",
            mc, lp_max);
    abort();
  }
  return (int)std::min<size_t>(mc, (smem_max - fixed) / per_row);
}

void context_info::calculate_all_coefficients(const enum grid_func func,
                                              int *lp_diff) {
  // Compute max angular momentum.
  const ldiffs_value ldiffs = prepare_get_ldiffs(func);
  smem_parameters smem_params(ldiffs, lmax());

  *lp_diff = smem_params.lp_diff();
  init_constant_memory();

  // kernel parameters
  kernel_params params = set_kernel_parameters(-1, smem_params);
  params.func = func;

  // Launch !
  const dim3 threads_per_block(4, 4, 4);

  // two step sphi transformation when its work arrays fit in the shared memory
  const int ld = ncoset(lmax());
  const size_t smem_two_steps =
      sizeof(double) * sphi_work_size(lmax(), max_nsgf_set());
  if (smem_two_steps <= sphi_two_step_max_smem) {
    if (func == GRID_FUNC_AB) {
      calculate_coefficients_two_steps<double, true>
          <<<nblocks, threads_per_block, smem_two_steps, main_stream>>>(
              params, ld, max_nsgf_set());
    } else {
      calculate_coefficients_two_steps<double, false>
          <<<nblocks, threads_per_block, smem_two_steps, main_stream>>>(
              params, ld, max_nsgf_set());
    }
    return;
  }

  if (func == GRID_FUNC_AB) {
    calculate_coefficients<double, true>
        <<<nblocks, threads_per_block, 0, main_stream>>>(params);
  } else {
    calculate_coefficients<double, false>
        <<<nblocks, threads_per_block, 0, main_stream>>>(params);
  }
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

  // cartesian factorization for an axis aligned cell, lattice coordinates +
  // transformation otherwise.
  auto &grid = grid_[level];
  const int lp_max = smem_params.lp_max();
  const int max_cube_size = std::max(grid.max_cube_size(), 1);
  const bool axis_aligned = grid.is_axis_aligned();
  // lp_max <= grid_gpu_max_l == GRID_SEP_MAX_L is checked by smem_parameters
  static_assert(grid_gpu_max_l <= GRID_SEP_MAX_L, "separable tables too small");

  init_separable_constant_memory();
  const double *lattice_to_cartesian =
      axis_aligned ? nullptr : get_lattice_to_cartesian(grid, lp_max);
  const int erows =
      collocate_separable_erows(max_cube_size, lp_max, separable_max_smem);
  const size_t smem =
      collocate_separable_smem_size(max_cube_size, lp_max, erows);
  const int nblocks_level = number_of_tasks_per_level_[level];
  if (grid.is_distributed()) {
    if (axis_aligned)
      collocate_separable_kernel<true, true>
          <<<nblocks_level, 64, smem, level_streams[level]>>>(
              params, max_cube_size, lp_max, erows, lattice_to_cartesian);
    else
      collocate_separable_kernel<true, false>
          <<<nblocks_level, 64, smem, level_streams[level]>>>(
              params, max_cube_size, lp_max, erows, lattice_to_cartesian);
  } else {
    if (axis_aligned)
      collocate_separable_kernel<false, true>
          <<<nblocks_level, 64, smem, level_streams[level]>>>(
              params, max_cube_size, lp_max, erows, lattice_to_cartesian);
    else
      collocate_separable_kernel<false, false>
          <<<nblocks_level, 64, smem, level_streams[level]>>>(
              params, max_cube_size, lp_max, erows, lattice_to_cartesian);
  }
}
} // namespace rocm_backend
