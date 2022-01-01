/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
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

#include "../common/grid_basis_set.h"
#include "../common/grid_common.h"
#include "grid_gpu_task_list.h"

#if (GRID_DO_COLLOCATE)
#define GRID_CONST_WHEN_COLLOCATE const
#define GRID_CONST_WHEN_INTEGRATE
#else
#define GRID_CONST_WHEN_COLLOCATE
#define GRID_CONST_WHEN_INTEGRATE const
#endif

/*******************************************************************************
 * \brief Forward declarations for inner-most loop bodies.
 *        Implementations are in grid_gpu_collocate.cu / grid_gpu_integrate.cu.
 * \author Ole Schuett
 ******************************************************************************/
#if (GRID_DO_COLLOCATE)
// collocate
typedef double cxyz_store;

__device__ static void cxyz_to_gridpoint(const double dx, const double dy,
                                         const double dz, const double zetp,
                                         const int lp, const cxyz_store *cxyz,
                                         double *gridpoint);
#else
// integrate
typedef struct {
  double *regs;
  int offset;
} cxyz_store;

__device__ static void gridpoint_to_cxyz(const double dx, const double dy,
                                         const double dz, const double zetp,
                                         const int lp, const double *gridpoint,
                                         cxyz_store *store);
#endif

/*******************************************************************************
 * \brief Atomic add for doubles that also works prior to compute capability 6.
 * \author Ole Schuett
 ******************************************************************************/
__device__ static void atomicAddDouble(double *address, double val) {
  if (val == 0.0)
    return;

#if __CUDA_ARCH__ >= 600
  atomicAdd(address, val); // part of cuda library
#else
  // https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#atomic-functions
  unsigned long long int *address_as_ull = (unsigned long long int *)address;
  unsigned long long int old = *address_as_ull, assumed;

  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,
                    __double_as_longlong(val + __longlong_as_double(assumed)));

    // Uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
  } while (assumed != old);

#endif
}

/*******************************************************************************
 * \brief Parameters of the collocate kernel.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  int smem_cab_offset;
  int smem_alpha_offset;
  int smem_cxyz_offset;
  int first_task;
  int npts_global[3];
  int npts_local[3];
  int shift_local[3];
  double dh[3][3];
  double dh_inv[3][3];
  const grid_gpu_task *tasks;
  const double *pab_blocks;
  GRID_CONST_WHEN_INTEGRATE double *grid;
  int la_min_diff;
  int lb_min_diff;
  int la_max_diff;
  int lb_max_diff;

#if (GRID_DO_COLLOCATE)
  // collocate
  enum grid_func func;
#else
  // integrate
  double *hab_blocks;
  double *forces;
  double *virial;
#endif
} kernel_params;

/*******************************************************************************
 * \brief Shared memory representation of a task.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct : grid_gpu_task_struct {
  // angular momentum range of actual collocate / integrate operation
  int la_max;
  int lb_max;
  int la_min;
  int lb_min;
  int lp;

  // size of the cab matrix
  int n1;
  int n2;

  const double *pab_block;

#if (!GRID_DO_COLLOCATE)
  // integrate
  double *hab_block;
  double *forces_a;
  double *forces_b;
#endif
} smem_task;

/*******************************************************************************
 * \brief Tabulated functions kept in constant memory to reduce register count.
 * \author Ole Schuett
 ******************************************************************************/
__constant__ orbital coset_inv[1330];
__constant__ double binomial_coef[19][19];

/*******************************************************************************
 * \brief Initializes the device's constant memory.
 * \author Ole Schuett
 ******************************************************************************/
static void init_constant_memory() {
  static bool initialized = false;
  if (initialized) {
    return; // constant memory has to be initialized only once
  }

  // Inverse coset mapping
  orbital coset_inv_host[1330];
  for (int lx = 0; lx <= 18; lx++) {
    for (int ly = 0; ly <= 18 - lx; ly++) {
      for (int lz = 0; lz <= 18 - lx - ly; lz++) {
        const int i = coset(lx, ly, lz);
        coset_inv_host[i] = {{lx, ly, lz}};
      }
    }
  }
  cudaError_t error =
      cudaMemcpyToSymbol(coset_inv, &coset_inv_host, sizeof(coset_inv_host), 0,
                         cudaMemcpyHostToDevice);
  assert(error == cudaSuccess);

  // Binomial coefficient
  double binomial_coef_host[19][19];
  memset(binomial_coef_host, 0, sizeof(binomial_coef_host));
  for (int n = 0; n <= 18; n++) {
    for (int k = 0; k <= n; k++) {
      binomial_coef_host[n][k] = fac(n) / fac(k) / fac(n - k);
    }
  }
  error =
      cudaMemcpyToSymbol(binomial_coef, &binomial_coef_host,
                         sizeof(binomial_coef_host), 0, cudaMemcpyHostToDevice);
  assert(error == cudaSuccess);

  initialized = true;
}

/*******************************************************************************
 * \brief Collocates coefficients C_xyz onto the grid for orthorhombic case.
 * \author Ole Schuett
 ******************************************************************************/
__device__ static void
ortho_cxyz_to_grid(const kernel_params *params, const smem_task *task,
                   GRID_CONST_WHEN_COLLOCATE cxyz_store *cxyz,
                   GRID_CONST_WHEN_INTEGRATE double *grid) {

  // The cube contains an even number of grid points in each direction and
  // collocation is always performed on a pair of two opposing grid points.
  // Hence, the points with index 0 and 1 are both assigned distance zero via
  // the formular distance=(2*index-1)/2.
  const int kmin = ceil(-1e-8 - task->disr_radius * params->dh_inv[2][2]);
  for (int k = threadIdx.z + kmin; k <= 1 - kmin; k += blockDim.z) {
    const int ka = task->cube_center_shifted[2] + k;
    const int kg =
        modulo(ka, params->npts_global[2]); // target location on the grid
    const int kd = (2 * k - 1) / 2; // distance from center in grid points
    const double kr = kd * params->dh[2][2]; // distance from center in a.u.
    const double kremain = task->disr_radius * task->disr_radius - kr * kr;
    const int jmin =
        ceil(-1e-8 - sqrt(fmax(0.0, kremain)) * params->dh_inv[1][1]);
    for (int j = threadIdx.y + jmin; j <= 1 - jmin; j += blockDim.y) {
      const int ja = task->cube_center_shifted[1] + j;
      const int jg =
          modulo(ja, params->npts_global[1]); // target location on the grid
      const int jd = (2 * j - 1) / 2; // distance from center in grid points
      const double jr = jd * params->dh[1][1]; // distance from center in a.u.
      const double jremain = kremain - jr * jr;
      const int imin =
          ceil(-1e-8 - sqrt(fmax(0.0, jremain)) * params->dh_inv[0][0]);

      for (int i = threadIdx.x + imin; i <= 1 - imin; i += blockDim.x) {
        const int ia = task->cube_center_shifted[0] + i;
        const int ig =
            modulo(ia, params->npts_global[0]); // target location on the grid

        // The distances above (kd, jd) do not take roffset into account,
        // ie. they always snap to the next grid point. This allowed the legacy
        // implementation to cache the loop bounds.
        // For the calculation of the grid value we'll use the true distance.
        const double dx = i * params->dh[0][0] + task->cube_offset[0];
        const double dy = j * params->dh[1][1] + task->cube_offset[1];
        const double dz = k * params->dh[2][2] + task->cube_offset[2];

        const int grid_index =
            kg * params->npts_local[1] * params->npts_local[0] +
            jg * params->npts_local[0] + ig;

#if (GRID_DO_COLLOCATE)
        // collocate
        cxyz_to_gridpoint(dx, dy, dz, task->zetp, task->lp, cxyz,
                          &grid[grid_index]);
#else
        // integrate
        gridpoint_to_cxyz(dx, dy, dz, task->zetp, task->lp, &grid[grid_index],
                          cxyz);
#endif
      }
    }
  }

  __syncthreads(); // because of concurrent writes to grid / cxyz
}

/*******************************************************************************
 * \brief Collocates coefficients C_xyz onto the grid for general case.
 * \author Ole Schuett
 ******************************************************************************/
__device__ static void
general_cxyz_to_grid(const kernel_params *params, const smem_task *task,
                     GRID_CONST_WHEN_COLLOCATE cxyz_store *cxyz,
                     GRID_CONST_WHEN_INTEGRATE double *grid) {

  // Go over the grid
  const int k_start = threadIdx.z + task->index_min[2];
  for (int k = k_start; k <= task->index_max[2]; k += blockDim.z) {
    const int kg = modulo(k - params->shift_local[2], params->npts_global[2]);
    if (kg < task->bounds_k[0] || task->bounds_k[1] < kg) {
      continue;
    }
    const int j_start = threadIdx.y + task->index_min[1];
    for (int j = j_start; j <= task->index_max[1]; j += blockDim.y) {
      const int jg = modulo(j - params->shift_local[1], params->npts_global[1]);
      if (jg < task->bounds_j[0] || task->bounds_j[1] < jg) {
        continue;
      }

      const int i_start = threadIdx.x + task->index_min[0];
      for (int i = i_start; i <= task->index_max[0]; i += blockDim.x) {
        const int ig =
            modulo(i - params->shift_local[0], params->npts_global[0]);
        if (ig < task->bounds_i[0] || task->bounds_i[1] < ig) {
          continue;
        }
        // Compute distance of current grid point from center of Gaussian.
        const double di = i - task->gp[0];
        double dx = di * params->dh[0][0];
        double dy = di * params->dh[0][1];
        double dz = di * params->dh[0][2];
        const double dj = j - task->gp[1];
        dx += dj * params->dh[1][0];
        dy += dj * params->dh[1][1];
        dz += dj * params->dh[1][2];
        const double dk = k - task->gp[2];
        dx += dk * params->dh[2][0];
        dy += dk * params->dh[2][1];
        dz += dk * params->dh[2][2];

        // Cycle if the point is not within the radius.
        const double r2 = dx * dx + dy * dy + dz * dz;
        if (r2 >= task->radius2) {
          continue;
        }

        const int stride = params->npts_local[1] * params->npts_local[0];
        const int grid_index = kg * stride + jg * params->npts_local[0] + ig;

#if (GRID_DO_COLLOCATE)
        // collocate
        cxyz_to_gridpoint(dx, dy, dz, task->zetp, task->lp, cxyz,
                          &grid[grid_index]);
#else
        // integrate
        gridpoint_to_cxyz(dx, dy, dz, task->zetp, task->lp, &grid[grid_index],
                          cxyz);
#endif
      }
    }
  }

  __syncthreads(); // because of concurrent writes to grid / cxyz
}

/*******************************************************************************
 * \brief Computes the polynomial expansion coefficients:
 *        (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
 * \author Ole Schuett
 ******************************************************************************/
__device__ static void compute_alpha(const kernel_params *params,
                                     const smem_task *task, double *alpha) {
  // strides for accessing alpha
  const int s3 = (task->lp + 1);
  const int s2 = (task->la_max + 1) * s3;
  const int s1 = (task->lb_max + 1) * s2;

  for (int idir = threadIdx.z; idir < 3; idir += blockDim.z) {
    const double drpa = task->rp[idir] - task->ra[idir];
    const double drpb = task->rp[idir] - task->rb[idir];
    for (int la = threadIdx.y; la <= task->la_max; la += blockDim.y) {
      for (int lb = threadIdx.x; lb <= task->lb_max; lb += blockDim.x) {
        for (int i = 0; i <= task->lp; i++) {
          const int base = idir * s1 + lb * s2 + la * s3;
          alpha[base + i] = 0.0;
        }
        double a = 1.0;
        for (int k = 0; k <= la; k++) {
          double b = 1.0;
          for (int l = 0; l <= lb; l++) {
            const int base = idir * s1 + lb * s2 + la * s3;
            alpha[base + la - l + lb - k] +=
                a * b * binomial_coef[la][k] * binomial_coef[lb][l];
            b *= drpb;
          }
          a *= drpa;
        }
      }
    }
  }
  __syncthreads(); // because of concurrent writes to alpha
}

/*******************************************************************************
 * \brief Transforms coefficients C_ab into C_xyz.
 * \author Ole Schuett
 ******************************************************************************/
__device__ static void cab_to_cxyz(const kernel_params *params,
                                   const smem_task *task, const double *alpha,
                                   GRID_CONST_WHEN_COLLOCATE double *cab,
                                   GRID_CONST_WHEN_INTEGRATE double *cxyz) {

  //   *** initialise the coefficient matrix, we transform the sum
  //
  // sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} *
  //         (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya
  //         (z-a_z)**lza
  //
  // into
  //
  // sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
  //
  // where p is center of the product gaussian, and lp = la_max + lb_max
  // (current implementation is l**7)

  // strides for accessing alpha
  const int s3 = (task->lp + 1);
  const int s2 = (task->la_max + 1) * s3;
  const int s1 = (task->lb_max + 1) * s2;

  // TODO: Maybe we can transpose alpha to index it directly with ico and jco.

#if (GRID_DO_COLLOCATE)
  // collocate
  for (int lzp = threadIdx.z; lzp <= task->lp; lzp += blockDim.z) {
    for (int lyp = threadIdx.y; lyp <= task->lp - lzp; lyp += blockDim.y) {
      for (int lxp = threadIdx.x; lxp <= task->lp - lzp - lyp;
           lxp += blockDim.x) {
        double reg = 0.0; // accumulate into a register
        for (int jco = 0; jco < ncoset(task->lb_max); jco++) {
          const orbital b = coset_inv[jco];
          for (int ico = 0; ico < ncoset(task->la_max); ico++) {
            const orbital a = coset_inv[ico];
#else
  // integrate
  if (threadIdx.z == 0) { // TODO: How bad is this?
    for (int jco = threadIdx.y; jco < ncoset(task->lb_max); jco += blockDim.y) {
      const orbital b = coset_inv[jco];
      for (int ico = threadIdx.x; ico < ncoset(task->la_max);
           ico += blockDim.x) {
        const orbital a = coset_inv[ico];
        double reg = 0.0; // accumulate into a register
        for (int lzp = 0; lzp <= task->lp; lzp++) {
          for (int lyp = 0; lyp <= task->lp - lzp; lyp++) {
            for (int lxp = 0; lxp <= task->lp - lzp - lyp; lxp++) {
#endif

            const double p = task->prefactor *
                             alpha[0 * s1 + b.l[0] * s2 + a.l[0] * s3 + lxp] *
                             alpha[1 * s1 + b.l[1] * s2 + a.l[1] * s3 + lyp] *
                             alpha[2 * s1 + b.l[2] * s2 + a.l[2] * s3 + lzp];

#if (GRID_DO_COLLOCATE)
            reg += p * cab[jco * task->n1 + ico]; // collocate
#else
              reg += p * cxyz[coset(lxp, lyp, lzp)]; // integrate
#endif
          }
        }

#if (GRID_DO_COLLOCATE)
        // collocate
        cxyz[coset(lxp, lyp, lzp)] = reg; // overwrite - no zeroing needed.
      }
#else
          // integrate
        }
        cab[jco * task->n1 + ico] = reg; // partial loop coverage -> zero it
      }
#endif
    }
  }
  __syncthreads(); // because of concurrent writes to cxyz / cab
}

/*******************************************************************************
 * \brief Initializes the cab matrix with zeros.
 * \author Ole Schuett
 ******************************************************************************/
__device__ static void zero_cab(const smem_task *task, double *cab) {
  if (threadIdx.z == 0) {
    for (int i = threadIdx.y; i < task->n2; i += blockDim.y) {
      for (int j = threadIdx.x; j < task->n1; j += blockDim.x) {
        cab[i * task->n1 + j] = 0.0;
      }
    }
  }
  __syncthreads(); // because of concurrent writes to cab
}

/*******************************************************************************
 * \brief Copies a task from global to shared memory and does precomputations.
 * \author Ole Schuett
 ******************************************************************************/
__device__ static void load_task(const kernel_params *params, smem_task *task) {

  // Parallel copy of base task from global to shared memory.
  if (threadIdx.y == 0 && threadIdx.z == 0) {
    const int itask = params->first_task + blockIdx.x;
    const grid_gpu_task *src_task = &params->tasks[itask];
    grid_gpu_task *dest_task = task; // Upcast to base struct.
    for (int i = threadIdx.x; i < sizeof(grid_gpu_task); i += blockDim.x) {
      ((char *)dest_task)[i] = ((const char *)src_task)[i];
    }
  }
  __syncthreads(); // because of concurrent writes to task

  if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0) {
    // angular momentum range for the actual collocate/integrate opteration.
    task->la_max = task->la_max_basis + params->la_max_diff;
    task->lb_max = task->lb_max_basis + params->lb_max_diff;
    task->la_min = imax(task->la_min_basis + params->la_min_diff, 0);
    task->lb_min = imax(task->lb_min_basis + params->lb_min_diff, 0);
    task->lp = task->la_max + task->lb_max;

    // size of the cab matrix
    task->n1 = ncoset(task->la_max);
    task->n2 = ncoset(task->lb_max);

    task->pab_block = &params->pab_blocks[task->ab_block_offset];

    // integrate
#if (!GRID_DO_COLLOCATE)
    task->hab_block = &params->hab_blocks[task->ab_block_offset];
    if (params->forces != NULL) {
      task->forces_a = &params->forces[3 * task->iatom];
      task->forces_b = &params->forces[3 * task->jatom];
    }
#endif
  }
  __syncthreads(); // because of concurrent writes to task
}

#endif // __GRID_CUDA
// EOF
