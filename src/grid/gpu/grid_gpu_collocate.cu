/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/
#include <algorithm>
#include <assert.h>
#include <cuda.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/grid_common.h"
#include "../common/grid_prepare_pab.h"
#include "grid_gpu_collocate.h"

/*******************************************************************************
 * \brief Atomic add for doubles that also works prior to compute capability 6.
 * \author Ole Schuett
 ******************************************************************************/
__device__ double atomicAddDouble(double *address, double val) {
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

  return __longlong_as_double(old);
#endif
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
 * \brief Compute value of a single grid point with distance d{xyz} from center.
 * \author Ole Schuett
 ******************************************************************************/
__device__ static double compute_point_value(const double dx, const double dy,
                                             const double dz, const double zetp,
                                             const int lp, const double *cxyz) {

  double result = 0.0; // accumulate into register
  double pow_dz = 1.0;
  for (int lzp = 0; lzp <= lp; lzp++) {
    double pow_dy = 1.0;
    for (int lyp = 0; lyp <= lp - lzp; lyp++) {
      double pow_dx = 1.0;
      for (int lxp = 0; lxp <= lp - lzp - lyp; lxp++) {
        const double pol = pow_dx * pow_dy * pow_dz;
        result += cxyz[coset(lxp, lyp, lzp)] * pol;
        pow_dx *= dx; // pow_dx = pow(dx, lxp)
      }
      pow_dy *= dy; // pow_dy = pow(dy, lyp)
    }
    pow_dz *= dz; // pow_dz = pow(dz, lzp)
  }

  // Squared distance of point from center.
  const double r2 = dx * dx + dy * dy + dz * dz;
  return result * exp(-zetp * r2);
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
  int func;
  prepare_ldiffs ldiffs;
  bool orthorhombic;
  int npts_global[3];
  int npts_local[3];
  int shift_local[3];
  int border_width[3];
  double dh[3][3];
  double dh_inv[3][3];
  // pointers
  const grid_gpu_task *tasks;
  const int *atom_kinds;
  const grid_basis_set *basis_sets;
  const int *block_offsets;
  const double *blocks_buffer;
  const double *atom_positions;
  double *grid;
} kernel_params;

/*******************************************************************************
 * \brief Shared memory representation of a task.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  int border_mask;
  bool block_transposed;
  double radius;
  double radius2;
  double ra[3];
  double rb[3];
  double rp[3];
  double rab[3];
  double gp[3];
  double rab2;
  double zeta;
  double zetb;
  double zetp;
  double prefactor;
  double dh_max;
  // angular momentum range BEFORE prepare_pab is applied
  int la_max_pab;
  int lb_max_pab;
  int la_min_pab;
  int lb_min_pab;
  // angular momentum range AFTER prepare_pab is applied
  int la_max;
  int lb_max;
  int la_min;
  int lb_min;
  int lp;
  // size of entire spherical basis
  int nsgfa;
  int nsgfb;
  // size of spherical set
  int nsgf_seta;
  int nsgf_setb;
  // size of intermediate pab matrix
  int ncoseta;
  int ncosetb;
  // size of final cab matrix
  int n1_cab;
  int n2_cab;
  // strides of the sphi transformation matrices
  int maxcoa;
  int maxcob;
  // pointers matrices
  const double *block;
  const double *sphia;
  const double *sphib;
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
__device__ static void ortho_cxyz_to_grid(const kernel_params *params,
                                          const smem_task *task,
                                          const double *cxyz, double *grid) {
  // Discretize the radius.
  const double drmin =
      fmin(params->dh[0][0], fmin(params->dh[1][1], params->dh[2][2]));
  const int imr = imax(1, (int)ceil(task->radius / drmin));
  const double disr_radius = drmin * imr;

  // *** position of the gaussian product
  //
  // this is the actual definition of the position on the grid
  // i.e. a point rp(:) gets here grid coordinates
  // MODULO(rp(:)/dr(:),npts_global(:))+1
  // hence (0.0,0.0,0.0) in real space is rsgrid%lb on the rsgrid in Fortran
  // and (1,1,1) on grid here in C.

  // cubecenter(:) = FLOOR(MATMUL(dh_inv, rp))
  int cubecenter[3];
  for (int i = 0; i < 3; i++) {
    cubecenter[i] = (int)floor(task->gp[i]);
  }

  // The cube contains an even number of grid points in each direction and
  // collocation is always performed on a pair of two opposing grid points.
  // Hence, the points with index 0 and 1 are both assigned distance zero via
  // the formular distance=(2*index-1)/2.

  const int kmin = ceil(-1e-8 - disr_radius * params->dh_inv[2][2]);
  for (int k = threadIdx.z + kmin; k <= 1 - kmin; k += blockDim.z) {
    const int ka = cubecenter[2] + k - params->shift_local[2];
    const int kg =
        modulo(ka, params->npts_global[2]); // target location on the grid
    const int kd = (2 * k - 1) / 2; // distance from center in grid points
    const double kr = kd * params->dh[2][2]; // distance from center in a.u.
    const double kremain = disr_radius * disr_radius - kr * kr;
    const int jmin =
        ceil(-1e-8 - sqrt(fmax(0.0, kremain)) * params->dh_inv[1][1]);
    for (int j = threadIdx.y + jmin; j <= 1 - jmin; j += blockDim.y) {
      const int ja = cubecenter[1] + j - params->shift_local[1];
      const int jg =
          modulo(ja, params->npts_global[1]); // target location on the grid
      const int jd = (2 * j - 1) / 2; // distance from center in grid points
      const double jr = jd * params->dh[1][1]; // distance from center in a.u.
      const double jremain = kremain - jr * jr;
      const int imin =
          ceil(-1e-8 - sqrt(fmax(0.0, jremain)) * params->dh_inv[0][0]);
      for (int i = threadIdx.x + imin; i <= 1 - imin; i += blockDim.x) {
        const int ia = cubecenter[0] + i - params->shift_local[0];
        const int ig =
            modulo(ia, params->npts_global[0]); // target location on the grid

        // The distances above (kd, jd) do not take roffset into account,
        // ie. they always snap to the next grid point. This allowed the legacy
        // implementation to cache the loop bounds.
        // For the calculation of the grid value we'll use the true distance.
        const double dx = (cubecenter[0] + i) * params->dh[0][0] - task->rp[0];
        const double dy = (cubecenter[1] + j) * params->dh[1][1] - task->rp[1];
        const double dz = (cubecenter[2] + k) * params->dh[2][2] - task->rp[2];

        // Compute grid point value and add to the global array.
        const int grid_index =
            kg * params->npts_local[1] * params->npts_local[0] +
            jg * params->npts_local[0] + ig;
        const double value =
            compute_point_value(dx, dy, dz, task->zetp, task->lp, cxyz);
        atomicAddDouble(&grid[grid_index], value);
      }
    }
  }
  __syncthreads(); // because of concurrent writes to grid
}

/*******************************************************************************
 * \brief Collocates coefficients C_xyz onto the grid for general case.
 * \author Ole Schuett
 ******************************************************************************/
__device__ static void general_cxyz_to_grid(const kernel_params *params,
                                            const smem_task *task,
                                            const double *cxyz, double *grid) {

  // get the min max indices that contain at least the cube that contains a
  // sphere around rp of radius radius if the cell is very non-orthogonal this
  // implies that many useless points are included this estimate can be
  // improved (i.e. not box but sphere should be used)
  int index_min[3] = {INT_MAX, INT_MAX, INT_MAX};
  int index_max[3] = {INT_MIN, INT_MIN, INT_MIN};
  for (int i = -1; i <= 1; i++) {
    for (int j = -1; j <= 1; j++) {
      for (int k = -1; k <= 1; k++) {
        const double x = task->rp[0] + i * task->radius;
        const double y = task->rp[1] + j * task->radius;
        const double z = task->rp[2] + k * task->radius;
        for (int idir = 0; idir < 3; idir++) {
          const double resc = params->dh_inv[0][idir] * x +
                              params->dh_inv[1][idir] * y +
                              params->dh_inv[2][idir] * z;
          index_min[idir] = imin(index_min[idir], (int)floor(resc));
          index_max[idir] = imax(index_max[idir], (int)ceil(resc));
        }
      }
    }
  }

  // Default for border_mask == 0.
  int bounds_i[2] = {0, params->npts_local[0] - 1};
  int bounds_j[2] = {0, params->npts_local[1] - 1};
  int bounds_k[2] = {0, params->npts_local[2] - 1};

  // See also rs_find_node() in task_list_methods.F.
  // If the bit is set then we need to exclude the border in that direction.
  if (task->border_mask & (1 << 0))
    bounds_i[0] += params->border_width[0];
  if (task->border_mask & (1 << 1))
    bounds_i[1] -= params->border_width[0];
  if (task->border_mask & (1 << 2))
    bounds_j[0] += params->border_width[1];
  if (task->border_mask & (1 << 3))
    bounds_j[1] -= params->border_width[1];
  if (task->border_mask & (1 << 4))
    bounds_k[0] += params->border_width[2];
  if (task->border_mask & (1 << 5))
    bounds_k[1] -= params->border_width[2];

  // Go over the grid
  const int k_start = threadIdx.z + index_min[2];
  for (int k = k_start; k <= index_max[2]; k += blockDim.z) {
    const int kg = modulo(k - params->shift_local[2], params->npts_global[2]);
    if (kg < bounds_k[0] || bounds_k[1] < kg) {
      continue;
    }
    const int j_start = threadIdx.y + index_min[1];
    for (int j = j_start; j <= index_max[1]; j += blockDim.y) {
      const int jg = modulo(j - params->shift_local[1], params->npts_global[1]);
      if (jg < bounds_j[0] || bounds_j[1] < jg) {
        continue;
      }
      const int i_start = threadIdx.x + index_min[0];
      for (int i = i_start; i <= index_max[0]; i += blockDim.x) {
        const int ig =
            modulo(i - params->shift_local[0], params->npts_global[0]);
        if (ig < bounds_i[0] || bounds_i[1] < ig) {
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

        // Compute grid point value and add to the global array.
        const double value =
            compute_point_value(dx, dy, dz, task->zetp, task->lp, cxyz);
        const int stride = params->npts_local[1] * params->npts_local[0];
        const int grid_index = kg * stride + jg * params->npts_local[0] + ig;
        atomicAddDouble(&grid[grid_index], value);
      }
    }
  }
  __syncthreads(); // because of concurrent writes to grid
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
                                   const double *cab, double *cxyz) {

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
  for (int lzp = threadIdx.z; lzp <= task->lp; lzp += blockDim.z) {
    for (int lyp = threadIdx.y; lyp <= task->lp - lzp; lyp += blockDim.y) {
      for (int lxp = threadIdx.x; lxp <= task->lp - lzp - lyp;
           lxp += blockDim.x) {
        double reg = 0.0; // accumulate into a register
        for (int jco = 0; jco < ncoset(task->lb_max); jco++) {
          const orbital b = coset_inv[jco];
          for (int ico = 0; ico < ncoset(task->la_max); ico++) {
            const orbital a = coset_inv[ico];
            reg += task->prefactor *
                   cab[jco * task->n1_cab + ico] * // [jco, ico]
                   alpha[0 * s1 + b.l[0] * s2 + a.l[0] * s3 + lxp] *
                   alpha[1 * s1 + b.l[1] * s2 + a.l[1] * s3 + lyp] *
                   alpha[2 * s1 + b.l[2] * s2 + a.l[2] * s3 + lzp];
          }
        }
        cxyz[coset(lxp, lyp, lzp)] = reg; // overwrite - no zeroing needed.
      }
    }
  }
  __syncthreads(); // because of concurrent writes to cxyz
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

  // Zero cab.
  if (threadIdx.z == 0) {
    for (int i = threadIdx.y; i < task->n2_cab; i += blockDim.y) {
      for (int j = threadIdx.x; j < task->n1_cab; j += blockDim.x) {
        cab[i * task->n1_cab + j] = 0.0;
      }
    }
  }
  __syncthreads(); // because of concurrent writes to cab

  // Decontract block, apply prepare_pab, and store in cab.
  // This is a double matrix product. Since the pab block can be quite large the
  // two products are fused to conserve shared memory.
  for (int i = threadIdx.x; i < task->nsgf_setb; i += blockDim.x) {
    for (int j = threadIdx.y; j < task->nsgf_seta; j += blockDim.y) {
      double block_val;
      if (task->block_transposed) {
        block_val = task->block[j * task->nsgfb + i];
      } else {
        block_val = task->block[i * task->nsgfa + j];
      }

      if (IS_FUNC_AB) {
        // fast path for common case
        for (int k = threadIdx.z; k < task->ncosetb; k += blockDim.z) {
          const double sphib = task->sphib[i * task->maxcob + k];
          for (int l = 0; l < task->ncoseta; l++) {
            const double sphia = task->sphia[j * task->maxcoa + l];
            const double pab_val = block_val * sphia * sphib;
            atomicAddDouble(&cab[k * task->ncoseta + l], pab_val);
          }
        }
      } else {
        // Since prepare_pab is a register hog we use it only when really needed
        for (int k = threadIdx.z; k < task->ncosetb; k += blockDim.z) {
          const orbital b = coset_inv[k];
          for (int l = 0; l < task->ncoseta; l++) {
            const orbital a = coset_inv[l];
            const double sphia = task->sphia[j * task->maxcoa + idx(a)];
            const double sphib = task->sphib[i * task->maxcob + idx(b)];
            const double pab_val = block_val * sphia * sphib;
            prepare_pab(params->func, a, b, task->zeta, task->zetb, pab_val,
                        task->n1_cab, cab);
          }
        }
      }
    }
  }
  __syncthreads(); // because of concurrent writes to cab
}

/*******************************************************************************
 * \brief Copies a task from global to shared memory and does precomputations.
 * \author Ole Schuett
 ******************************************************************************/
__device__ static void fill_smem_task(const kernel_params *params,
                                      smem_task *task) {

  if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0) {
    const int itask = params->first_task + blockIdx.x;
    const grid_gpu_task *glb_task = &params->tasks[itask];
    const int iatom = glb_task->iatom - 1;
    const int jatom = glb_task->jatom - 1;
    const int iset = glb_task->iset - 1;
    const int jset = glb_task->jset - 1;
    const int ipgf = glb_task->ipgf - 1;
    const int jpgf = glb_task->jpgf - 1;
    const int ikind = params->atom_kinds[iatom] - 1;
    const int jkind = params->atom_kinds[jatom] - 1;
    const grid_basis_set ibasis = params->basis_sets[ikind];
    const grid_basis_set jbasis = params->basis_sets[jkind];

    task->zeta = ibasis.zet[iset * ibasis.maxpgf + ipgf];
    task->zetb = jbasis.zet[jset * jbasis.maxpgf + jpgf];
    task->zetp = task->zeta + task->zetb;
    const double f = task->zetb / task->zetp;

    const double *glb_ra = &params->atom_positions[3 * iatom];
    task->rab2 = 0.0;
    for (int i = 0; i < 3; i++) {
      task->rab[i] = glb_task->rab[i];
      task->rab2 += task->rab[i] * task->rab[i];
      task->ra[i] = glb_ra[i];
      task->rb[i] = task->ra[i] + task->rab[i];
      task->rp[i] = task->ra[i] + task->rab[i] * f;
    }

    // center in grid coords, gp = MATMUL(dh_inv, rp)
    for (int i = 0; i < 3; i++) {
      task->gp[i] = params->dh_inv[0][i] * task->rp[0] +
                    params->dh_inv[1][i] * task->rp[1] +
                    params->dh_inv[2][i] * task->rp[2];
    }

    // resolution of the grid
    task->dh_max = 0.0;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        task->dh_max = fmax(task->dh_max, fabs(params->dh[i][j]));
      }
    }

    const double rscale = (iatom == jatom) ? 1.0 : 2.0;
    task->prefactor = rscale * exp(-task->zeta * f * task->rab2);
    task->border_mask = glb_task->border_mask;
    task->radius = glb_task->radius;
    task->radius2 = task->radius * task->radius;

    // angular momentum range BEFORE prepare_pab is applied
    task->la_max_pab = ibasis.lmax[iset];
    task->lb_max_pab = jbasis.lmax[jset];
    task->la_min_pab = ibasis.lmin[iset];
    task->lb_min_pab = jbasis.lmin[jset];

    // angular momentum range AFTER prepare_pab is applied
    task->la_max = ibasis.lmax[iset] + params->ldiffs.la_max_diff;
    task->lb_max = jbasis.lmax[jset] + params->ldiffs.lb_max_diff;
    task->la_min = imax(ibasis.lmin[iset] + params->ldiffs.la_min_diff, 0);
    task->lb_min = imax(jbasis.lmin[jset] + params->ldiffs.lb_min_diff, 0);
    task->lp = task->la_max + task->lb_max;

    // size of entire spherical basis
    task->nsgfa = ibasis.nsgf;
    task->nsgfb = jbasis.nsgf;

    // size of spherical set
    task->nsgf_seta = ibasis.nsgf_set[iset];
    task->nsgf_setb = jbasis.nsgf_set[jset];

    // size of intermediate pab matrix
    task->ncoseta = ncoset(task->la_max_pab);
    task->ncosetb = ncoset(task->lb_max_pab);

    // size of final cab matrix
    task->n1_cab = ncoset(task->la_max);
    task->n2_cab = ncoset(task->lb_max);

    // strides of the sphi transformation matrices
    task->maxcoa = ibasis.maxco;
    task->maxcob = jbasis.maxco;

    // start of spherical set within the basis
    const int sgfa = ibasis.first_sgf[iset] - 1;
    const int sgfb = jbasis.first_sgf[jset] - 1;

    // start of exponent within the cartesian set
    const int o1 = ipgf * task->ncoseta;
    const int o2 = jpgf * task->ncosetb;

    // transformations from contracted spherical to primitiv carthesian basis
    task->sphia = &ibasis.sphi[sgfa * task->maxcoa + o1];
    task->sphib = &jbasis.sphi[sgfb * task->maxcob + o2];

    // Locate current matrix block within the buffer.
    const int block_num = glb_task->block_num - 1;
    const int block_offset = params->block_offsets[block_num];
    task->block_transposed = (iatom > jatom);
    if (task->block_transposed) {
      task->block =
          &params->blocks_buffer[block_offset + sgfa * task->nsgfb + sgfb];
    } else {
      task->block =
          &params->blocks_buffer[block_offset + sgfb * task->nsgfa + sgfa];
    }
  }
  __syncthreads(); // because of concurrent writes to task
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

  block_to_cab<IS_FUNC_AB>(params, &task, smem_cab);

  compute_alpha(params, &task, smem_alpha);

  cab_to_cxyz(params, &task, smem_alpha, smem_cab, smem_cxyz);

  if (params->orthorhombic && task.border_mask == 0) {
    ortho_cxyz_to_grid(params, &task, smem_cxyz, params->grid);
  } else {
    general_cxyz_to_grid(params, &task, smem_cxyz, params->grid);
  }
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
    const int last_task, const bool orthorhombic, const int func,
    const int npts_global[3], const int npts_local[3], const int shift_local[3],
    const int border_width[3], const double dh[3][3], const double dh_inv[3][3],
    const cudaStream_t stream, double *grid_dev) {

  const int ntasks = last_task - first_task + 1;
  if (ntasks == 0) {
    return; // Nothing to do.
  }

  init_constant_memory();

  // Compute required shared memory.
  // TODO: Currently, cab's indicies run over 0...ncoset[lmax],
  //       however only ncoset(lmin)...ncoset(lmax) are actually needed.
  const prepare_ldiffs ldiffs = prepare_get_ldiffs(func);
  const int la_max = task_list->lmax + ldiffs.la_max_diff;
  const int lb_max = task_list->lmax + ldiffs.lb_max_diff;
  const int lp_max = la_max + lb_max;
  const int cab_len = ncoset(lb_max) * ncoset(la_max);
  const int alpha_len = 3 * (lb_max + 1) * (la_max + 1) * (lp_max + 1);
  const int cxyz_len = ncoset(lp_max);
  const int alpha_cxyz_len = alpha_len + cxyz_len;
  const size_t smem_per_block = (cab_len + alpha_cxyz_len) * sizeof(double);

  if (smem_per_block > 48 * 1024) {
    fprintf(stderr, "ERROR: Not enough shared memory.\n");
    fprintf(stderr, "alpha_len: %i, ", alpha_len);
    fprintf(stderr, "cxyz_len: %i, ", alpha_cxyz_len);
    fprintf(stderr, "cab_len: %i, ", cab_len);
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
  params.func = func;
  params.grid = grid_dev;
  params.ldiffs = ldiffs;
  params.tasks = task_list->tasks_dev;
  params.atom_kinds = task_list->atom_kinds_dev;
  params.basis_sets = task_list->basis_sets_dev;
  params.block_offsets = task_list->block_offsets_dev;
  params.atom_positions = task_list->atom_positions_dev;
  params.blocks_buffer = task_list->blocks_buffer_dev;
  memcpy(params.dh, dh, 9 * sizeof(double));
  memcpy(params.dh_inv, dh_inv, 9 * sizeof(double));
  memcpy(params.npts_global, npts_global, 3 * sizeof(int));
  memcpy(params.npts_local, npts_local, 3 * sizeof(int));
  memcpy(params.shift_local, shift_local, 3 * sizeof(int));
  memcpy(params.border_width, border_width, 3 * sizeof(int));

  // Launch !
  const int nblocks = ntasks;
  const dim3 threads_per_block(4, 8, 8);

  if (func == GRID_FUNC_AB) {
    collocate_kernel_density<<<nblocks, threads_per_block, smem_per_block,
                               stream>>>(params);
  } else {
    collocate_kernel_anyfunc<<<nblocks, threads_per_block, smem_per_block,
                               stream>>>(params);
  }
}

// EOF
