/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

/*
 * Authors :
 - Dr Mathieu Taillefumier (ETH Zurich / CSCS)
 - Advanced Micro Devices, Inc.
*/

#ifndef GRID_HIP_INTERNAL_HEADER_H
#define GRID_HIP_INTERNAL_HEADER_H

#include <algorithm>
#include <assert.h>
#include <hip/hip_runtime.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define GRID_DEVICE __device__
extern "C" {
#include "../common/grid_basis_set.h"
#include "../common/grid_constants.h"
}

#include "grid_hip_context.hpp"

namespace rocm_backend {

#if defined(__HIP_PLATFORM_NVIDIA__)
#if __CUDA_ARCH__ < 600
__device__ __inline__ double atomicAdd(double *address, double val) {
  unsigned long long int *address_as_ull = (unsigned long long int *)address;
  unsigned long long int old = *address_as_ull, assumed;

  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,
                    __double_as_longlong(val + __longlong_as_double(assumed)));

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN !=
    // NaN)
  } while (assumed != old);

  return __longlong_as_double(old);
}
#endif
#endif

/*******************************************************************************
 * \brief Orbital angular momentum.
 ******************************************************************************/
struct orbital {
  int l[3];
};

__constant__ orbital coset_inv[1330];
__constant__ int binomial_coef[19][19];

/*******************************************************************************
 * \brief Differences in angular momentum.
 ******************************************************************************/
struct ldiffs_value {
  int la_max_diff{0};
  int la_min_diff{0};
  int lb_max_diff{0};
  int lb_min_diff{0};
};

/*******************************************************************************
 * \brief data needed for calculating the coefficients, forces, and stress
 ******************************************************************************/
template <typename T> struct smem_task {
  bool block_transposed;
  T ra[3];
  T rb[3];
  T rp[3];
  T rab[3];
  T zeta;
  T zetb;
  T zetp;
  T prefactor;
  T off_diag_twice;
  char la_max;
  char lb_max;
  char la_min;
  char lb_min;
  char lp;
  // size of the cab matrix
  unsigned short int n1;
  unsigned short int n2;
  // size of entire spherical basis
  unsigned short int nsgfa;
  unsigned short int nsgfb;
  // size of spherical set
  unsigned short int nsgf_seta;
  unsigned short int nsgf_setb;
  // start of decontracted set, ie. pab and hab
  int first_coseta;
  int first_cosetb;
  // size of decontracted set, ie. pab and hab
  int ncoseta;
  int ncosetb;
  // strides of the sphi transformation matrices
  int maxcoa;
  int maxcob;
  // pointers matrices
  T *pab_block;
  T *sphia;
  T *sphib;
  // integrate
  T *hab_block;
  T *forces_a;
  T *forces_b;
};

/*******************************************************************************
 * \brief data needed for collocate and integrate kernels
 ******************************************************************************/
template <typename T, typename T3> struct smem_task_reduced {
  // bool block_transposed;
  T radius, discrete_radius;
  int3 cube_center, lb_cube, cube_size, window_size, window_shift;
  T3 roffset;
  T zetp;
  // char la_max;
  // char lb_max;
  // char la_min;
  // char lb_min;
  char lp;
  bool apply_border_mask;
};

/*******************************************************************************
 * \brief Factorial function, e.g. fac(5) = 5! = 120.
 * \author Ole Schuett
 ******************************************************************************/
__device__ __inline__ double fac(const int i) {
  static const double table[] = {
      0.10000000000000000000E+01, 0.10000000000000000000E+01,
      0.20000000000000000000E+01, 0.60000000000000000000E+01,
      0.24000000000000000000E+02, 0.12000000000000000000E+03,
      0.72000000000000000000E+03, 0.50400000000000000000E+04,
      0.40320000000000000000E+05, 0.36288000000000000000E+06,
      0.36288000000000000000E+07, 0.39916800000000000000E+08,
      0.47900160000000000000E+09, 0.62270208000000000000E+10,
      0.87178291200000000000E+11, 0.13076743680000000000E+13,
      0.20922789888000000000E+14, 0.35568742809600000000E+15,
      0.64023737057280000000E+16, 0.12164510040883200000E+18,
      0.24329020081766400000E+19, 0.51090942171709440000E+20,
      0.11240007277776076800E+22, 0.25852016738884976640E+23,
      0.62044840173323943936E+24, 0.15511210043330985984E+26,
      0.40329146112660563558E+27, 0.10888869450418352161E+29,
      0.30488834461171386050E+30, 0.88417619937397019545E+31,
      0.26525285981219105864E+33};
  return table[i];
}

/*******************************************************************************
 * \brief Number of Cartesian orbitals up to given angular momentum quantum.
 * \author Ole Schuett
 ******************************************************************************/
__host__ __device__ __inline__ int ncoset(const int l) {
  static const int table[] = {1,  // l=0
                              4,  // l=1
                              10, // l=2 ...
                              20,  35,  56,  84,  120, 165, 220,  286,
                              364, 455, 560, 680, 816, 969, 1140, 1330};
  return table[l];
}

/*******************************************************************************
 * \brief Maps three angular momentum components to a single zero based index.
 ******************************************************************************/
__host__ __device__ __inline__ int coset(int lx, int ly, int lz) {
  const int l = lx + ly + lz;
  if (l == 0) {
    return 0;
  } else {
    return ncoset(l - 1) + ((l - lx) * (l - lx + 1)) / 2 + lz;
  }
}

/*******************************************************************************
 * \brief Increase i'th component of given orbital angular momentum.
 ******************************************************************************/
__device__ __inline__ orbital up(const int i, const orbital &a) {
  orbital b = a;
  b.l[i] += 1;
  return b;
}

/*******************************************************************************
 * \brief Decrease i'th component of given orbital angular momentum.
 ******************************************************************************/
__inline__ __device__ orbital down(const int i, const orbital &a) {
  orbital b = a;
  b.l[i] = max(0, a.l[i] - 1);
  return b;
}

/*******************************************************************************
 * \brief Return coset index of given orbital angular momentum.
 ******************************************************************************/
__inline__ __device__ int idx(const orbital a) {
  return coset(a.l[0], a.l[1], a.l[2]);
}

__device__ __inline__ double power(const double x, const int expo) {
  double tmp = 1.0;
  for (int i = 1; i <= expo; i++)
    tmp *= x;
  return tmp;
}

/*******************************************************************************
 * \brief Adds given value to matrix element cab[idx(b)][idx(a)].
 ******************************************************************************/
template <typename T = double>
__device__ __inline__ void prep_term(const orbital a, const orbital b,
                                     const T value, const int n, T *cab) {
  atomicAdd(&cab[idx(b) * n + idx(a)], value);
}

/*******************************************************************************
 * \brief Initializes the device's constant memory.
 * \author Ole Schuett
 ******************************************************************************/
inline static void init_constant_memory() {
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
  hipError_t error =
      hipMemcpyToSymbol(coset_inv, &coset_inv_host, sizeof(coset_inv_host), 0,
                        hipMemcpyHostToDevice);
  assert(error == hipSuccess);

  // Binomial coefficient
  int binomial_coef_host[19][19] = {
      {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {1, 3, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {1, 4, 6, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {1, 5, 10, 10, 5, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {1, 6, 15, 20, 15, 6, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {1, 7, 21, 35, 35, 21, 7, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {1, 8, 28, 56, 70, 56, 28, 8, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {1, 9, 36, 84, 126, 126, 84, 36, 9, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {1, 10, 45, 120, 210, 252, 210, 120, 45, 10, 1, 0, 0, 0, 0, 0, 0, 0, 0},
      {1, 11, 55, 165, 330, 462, 462, 330, 165, 55, 11, 1, 0, 0, 0, 0, 0, 0, 0},
      {1, 12, 66, 220, 495, 792, 924, 792, 495, 220, 66, 12, 1, 0, 0, 0, 0, 0,
       0},
      {1, 13, 78, 286, 715, 1287, 1716, 1716, 1287, 715, 286, 78, 13, 1, 0, 0,
       0, 0, 0},
      {1, 14, 91, 364, 1001, 2002, 3003, 3432, 3003, 2002, 1001, 364, 91, 14, 1,
       0, 0, 0, 0},
      {1, 15, 105, 455, 1365, 3003, 5005, 6435, 6435, 5005, 3003, 1365, 455,
       105, 15, 1, 0, 0, 0},
      {1, 16, 120, 560, 1820, 4368, 8008, 11440, 12870, 11440, 8008, 4368, 1820,
       560, 120, 16, 1, 0, 0},
      {1, 17, 136, 680, 2380, 6188, 12376, 19448, 24310, 24310, 19448, 12376,
       6188, 2380, 680, 136, 17, 1, 0},
      {1, 18, 153, 816, 3060, 8568, 18564, 31824, 43758, 48620, 43758, 31824,
       18564, 8568, 3060, 816, 153, 18, 1}};
  error =
      hipMemcpyToSymbol(binomial_coef, &binomial_coef_host[0][0],
                        sizeof(binomial_coef_host), 0, hipMemcpyHostToDevice);
  assert(error == hipSuccess);

  initialized = true;
}

__inline__ __device__ double3
compute_coordinates(const double *__restrict__ dh_, const double x,
                    const double y, const double z) {

  double3 r3;
  // I make no distinction between orthorhombic and non orthorhombic
  // cases

  r3.x = z * dh_[6] + y * dh_[3] + x * dh_[0];
  r3.y = z * dh_[7] + y * dh_[4] + x * dh_[1];
  r3.z = z * dh_[8] + y * dh_[5] + x * dh_[2];
  return r3;
}

__inline__ __device__ float3 compute_coordinates(const float *__restrict__ dh_,
                                                 const float x, const float y,
                                                 const float z) {

  float3 r3;
  // I make no distinction between orthorhombic and non orthorhombic
  // cases

  r3.x = z * dh_[6] + y * dh_[3] + x * dh_[0];
  r3.y = z * dh_[7] + y * dh_[4] + x * dh_[1];
  r3.z = z * dh_[8] + y * dh_[5] + x * dh_[2];
  return r3;
}

/*******************************************************************************
 * \brief Computes the polynomial expansion coefficients:
 *        (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
 ******************************************************************************/
template <typename T>
__inline__ __device__ void compute_alpha(const kernel_params &params,
                                         const smem_task<T> &task,
                                         T *__restrict__ alpha) {
  // strides for accessing alpha
  const int s3 = (task.lp + 1);
  const int s2 = (task.la_max + 1) * s3;
  const int s1 = (task.lb_max + 1) * s2;
  const int tid =
      threadIdx.x + blockDim.x * (threadIdx.y + blockDim.y * threadIdx.z);
  for (int i = tid; i < 3 * s1; i += blockDim.x * blockDim.y * blockDim.z)
    alpha[i] = 0.0;

  __syncthreads();

  for (int idir = threadIdx.z; idir < 3; idir += blockDim.z) {
    const T drpa = task.rp[idir] - task.ra[idir];
    const T drpb = task.rp[idir] - task.rb[idir];
    for (int la = threadIdx.y; la <= task.la_max; la += blockDim.y) {
      for (int lb = threadIdx.x; lb <= task.lb_max; lb += blockDim.x) {
        T a = 1.0;
        for (int k = 0; k <= la; k++) {
          T b = 1.0;
          const int base = idir * s1 + lb * s2 + la * s3;
          for (int l = 0; l <= lb; l++) {
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

__host__ __inline__ __device__ void
convert_to_lattice_coordinates(const double *dh_inv_,
                               const double3 *__restrict__ const rp,
                               double3 *__restrict__ rp_c) {
  rp_c->x = dh_inv_[0] * rp->x + dh_inv_[3] * rp->y + dh_inv_[6] * rp->z;
  rp_c->y = dh_inv_[1] * rp->x + dh_inv_[4] * rp->y + dh_inv_[7] * rp->z;
  rp_c->z = dh_inv_[2] * rp->x + dh_inv_[5] * rp->y + dh_inv_[8] * rp->z;
}

__host__ __inline__ __device__ void
convert_from_lattice_coordinates_to_cartesian(
    const double *__restrict__ dh_, const double3 *__restrict__ const rp,
    double3 *__restrict__ rp_c) {
  rp_c->x = dh_[0] * rp->x + dh_[3] * rp->y + dh_[6] * rp->z;
  rp_c->y = dh_[1] * rp->x + dh_[4] * rp->y + dh_[7] * rp->z;
  rp_c->z = dh_[2] * rp->x + dh_[5] * rp->y + dh_[8] * rp->z;
}

__host__ __inline__ __device__ void
convert_to_lattice_coordinates(const float *dh_inv_,
                               const float3 *__restrict__ const rp,
                               float3 *__restrict__ rp_c) {
  rp_c->x = dh_inv_[0] * rp->x + dh_inv_[3] * rp->y + dh_inv_[6] * rp->z;
  rp_c->y = dh_inv_[1] * rp->x + dh_inv_[4] * rp->y + dh_inv_[7] * rp->z;
  rp_c->z = dh_inv_[2] * rp->x + dh_inv_[5] * rp->y + dh_inv_[8] * rp->z;
}

__host__ __inline__ __device__ void
convert_from_lattice_coordinates_to_cartesian(
    const float *__restrict__ dh_, const float3 *__restrict__ const rp,
    float3 *__restrict__ rp_c) {
  rp_c->x = dh_[0] * rp->x + dh_[3] * rp->y + dh_[6] * rp->z;
  rp_c->y = dh_[1] * rp->x + dh_[4] * rp->y + dh_[7] * rp->z;
  rp_c->z = dh_[2] * rp->x + dh_[5] * rp->y + dh_[8] * rp->z;
}

template <typename T, typename T3, bool orthorhombic_>
__inline__ T compute_cube_properties(
    const T radius, const T *const __restrict__ dh_,
    const T *const __restrict__ dh_inv_, const T3 *__restrict__ rp,
    T3 *__restrict__ roffset, int3 *__restrict__ cubecenter,
    int3 *__restrict__ lb_cube, int3 *__restrict__ cube_size) {

  /* center of the gaussian in the lattice coordinates */
  T3 rp1, rp2, rp3;

  /* it is in the lattice vector frame */
  convert_to_lattice_coordinates(dh_inv_, rp, &rp1);

  /* compute the grid point that is the closest to the sphere center. */
  cubecenter->x = std::floor(rp1.x);
  cubecenter->y = std::floor(rp1.y);
  cubecenter->z = std::floor(rp1.z);

  /* seting up the cube parameters */

  if (orthorhombic_) {

    // the cube is actually slightly bigger than the sphere of radius r. that's
    // why we need to discretize it to get the cube size "right".

    // disc_radius >= radius always. somehow despite the fact that we compile
    // things with a c++ compiler on the host side, we need to use fmin instead
    // of std::min since std::min is not allowed on the device side

    // We assume no specific form for the orthogonal matrix. (can be diaognal or
    // completely full, the only constraint is that the tree vectors are
    // orthogonal)

    T norm1, norm2, norm3;
    norm1 = dh_[0] * dh_[0] + dh_[1] * dh_[1] + dh_[2] * dh_[2];
    norm2 = dh_[3] * dh_[3] + dh_[4] * dh_[4] + dh_[5] * dh_[5];
    norm3 = dh_[6] * dh_[6] + dh_[7] * dh_[7] + dh_[8] * dh_[8];

    norm1 = std::sqrt(norm1);
    norm2 = std::sqrt(norm2);
    norm3 = std::sqrt(norm3);

    const T disr_radius =
        std::min(norm1, std::min(norm2, norm3)) *
        std::max(1,
                 (int)ceil(radius / std::min(norm1, std::min(norm2, norm3))));

    rp2.x = cubecenter->x;
    rp2.y = cubecenter->y;
    rp2.z = cubecenter->z;

    /* convert the cube center from grid points coordinates to cartesian */
    convert_from_lattice_coordinates_to_cartesian(dh_, &rp2, roffset);
    /* cube center */
    roffset->x -= rp->x;
    roffset->y -= rp->y;
    roffset->z -= rp->z;

    rp2.x = disr_radius;
    rp2.y = disr_radius;
    rp2.z = disr_radius;

    /* it is in the lattice vector frame */
    convert_to_lattice_coordinates(dh_inv_, &rp2, &rp3);
    /* lower and upper bounds */
    lb_cube->x = std::ceil(-1e-8 - rp3.x);
    lb_cube->y = std::ceil(-1e-8 - rp3.y);
    lb_cube->z = std::ceil(-1e-8 - rp3.z);

    /* it is in the lattice vector frame */
    convert_to_lattice_coordinates(dh_inv_, roffset, &rp2);

    /* Express the offset in lattice coordinates */
    roffset->x = rp2.x;
    roffset->y = rp2.y;
    roffset->z = rp2.z;

    /* compute the cube size ignoring periodicity */
    /* the interval is not symmetrical for some curious reasons. it should go
     * from [-L..L+1] so the number of points is multiple of two */
    cube_size->x = 2 - 2 * lb_cube->x;
    cube_size->y = 2 - 2 * lb_cube->y;
    cube_size->z = 2 - 2 * lb_cube->z;
    return disr_radius;
  } else {
    int3 ub_cube;

    lb_cube->x = INT_MAX;
    ub_cube.x = INT_MIN;
    lb_cube->y = INT_MAX;
    ub_cube.y = INT_MIN;
    lb_cube->z = INT_MAX;
    ub_cube.z = INT_MIN;

    for (int i = -1; i <= 1; i++) {
      for (int j = -1; j <= 1; j++) {
        for (int k = -1; k <= 1; k++) {
          T3 r;
          r.x = rp->x + ((T)i) * radius;
          r.y = rp->y + ((T)j) * radius;
          r.z = rp->z + ((T)k) * radius;
          convert_to_lattice_coordinates(dh_inv_, &r, roffset);

          lb_cube->x = std::min(lb_cube->x, (int)std::floor(roffset->x));
          ub_cube.x = std::max(ub_cube.x, (int)std::ceil(roffset->x));

          lb_cube->y = std::min(lb_cube->y, (int)std::floor(roffset->y));
          ub_cube.y = std::max(ub_cube.y, (int)std::ceil(roffset->y));

          lb_cube->z = std::min(lb_cube->z, (int)std::floor(roffset->z));
          ub_cube.z = std::max(ub_cube.z, (int)std::ceil(roffset->z));
        }
      }
    }
    /* compute the cube size ignoring periodicity */
    cube_size->x = ub_cube.x - lb_cube->x;
    cube_size->y = ub_cube.y - lb_cube->y;
    cube_size->z = ub_cube.z - lb_cube->z;

    /* compute the offset in lattice coordinates */

    roffset->x = cubecenter->x - rp1.x;
    roffset->y = cubecenter->y - rp1.y;
    roffset->z = cubecenter->z - rp1.z;

    // shift the boundaries compared to the cube center so that the
    // specialization ortho / non ortho is minimal
    lb_cube->x -= cubecenter->x;
    lb_cube->y -= cubecenter->y;
    lb_cube->z -= cubecenter->z;

    return radius;
  }
}

__inline__ __device__ void
compute_window_size(const int *const grid_size, const int *const lower_corner_,
                    const int *const period_, /* also full size of the grid */
                    const int border_mask, const int *border_width,
                    int3 *const window_size, int3 *const window_shift) {
  window_shift->x = 0;
  window_shift->y = 0;
  window_shift->z = 0;

  window_size->x = grid_size[2] - 1;
  window_size->y = grid_size[1] - 1;
  window_size->z = grid_size[0] - 1;

  if (border_mask & (1 << 0))
    window_shift->x += border_width[2];
  if (border_mask & (1 << 1))
    window_size->x -= border_width[2];
  if (border_mask & (1 << 2))
    window_shift->y += border_width[1];
  if (border_mask & (1 << 3))
    window_size->y -= border_width[1];
  if (border_mask & (1 << 4))
    window_shift->z += border_width[0];
  if (border_mask & (1 << 5))
    window_size->z -= border_width[0];
}

/*******************************************************************************
 * \brief Transforms coefficients C_ab into C_xyz.
 ******************************************************************************/
template <typename T>
__device__ __inline__ static void
cab_to_cxyz(const kernel_params &params, const smem_task<T> &task,
            const T *__restrict__ alpha, const T *__restrict__ cab,
            T *__restrict__ cxyz) {

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
  const int s3 = (task.lp + 1);
  const int s2 = (task.la_max + 1) * s3;
  const int s1 = (task.lb_max + 1) * s2;

  // TODO: Maybe we can transpose alpha to index it directly with ico and jco.
  const int tid =
      threadIdx.x + blockDim.x * (threadIdx.y + blockDim.y * threadIdx.z);

  for (int i = tid; i < ncoset(task.lp);
       i += blockDim.x * blockDim.y * blockDim.z) {
    auto &co = coset_inv[i];
    T reg = 0.0; // accumulate into a register
    for (int jco = 0; jco < ncoset(task.lb_max); jco++) {
      const auto &b = coset_inv[jco];
      for (int ico = 0; ico < ncoset(task.la_max); ico++) {
        const auto &a = coset_inv[ico];
        const T p = task.prefactor *
                    alpha[0 * s1 + b.l[0] * s2 + a.l[0] * s3 + co.l[0]] *
                    alpha[1 * s1 + b.l[1] * s2 + a.l[1] * s3 + co.l[1]] *
                    alpha[2 * s1 + b.l[2] * s2 + a.l[2] * s3 + co.l[2]];
        reg += p * cab[jco * task.n1 + ico]; // collocate
      }
    }

    cxyz[i] = reg;
  }
  __syncthreads(); // because of concurrent writes to cxyz / cab
}

/*******************************************************************************
 * \brief Transforms coefficients C_xyz into C_ab.
 ******************************************************************************/
template <typename T>
__device__ __inline__ static void
cxyz_to_cab(const kernel_params &params, const smem_task<T> &task,
            const T *__restrict__ alpha, const T *__restrict__ cxyz,
            T *__restrict__ cab) {

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
  const int s3 = (task.lp + 1);
  const int s2 = (task.la_max + 1) * s3;
  const int s1 = (task.lb_max + 1) * s2;

  const int tid =
      threadIdx.x + blockDim.x * (threadIdx.y + blockDim.y * threadIdx.z);
  for (int jco = tid / 8; jco < ncoset(task.lb_max); jco += 8) {
    const orbital b = coset_inv[jco];
    for (int ico = tid % 8; ico < ncoset(task.la_max); ico += 8) {
      const auto &a = coset_inv[ico];
      T reg = 0.0; // accumulate into a register
      for (int ic = 0; ic < ncoset(task.lp); ic++) {
        const auto &co = coset_inv[ic];
        const T p = task.prefactor *
                    alpha[b.l[0] * s2 + a.l[0] * s3 + co.l[0]] *
                    alpha[s1 + b.l[1] * s2 + a.l[1] * s3 + co.l[1]] *
                    alpha[2 * s1 + b.l[2] * s2 + a.l[2] * s3 + co.l[2]];

        reg += p * cxyz[ic]; // integrate
      }
      cab[jco * task.n1 + ico] = reg; // partial loop coverage -> zero it
    }
  }
  __syncthreads(); // because of concurrent writes to cxyz / cab
}

/*******************************************************************************
 * \brief Copies a task from global to shared memory.
 ******************************************************************************/

/* Collocate and integrate do not care about the sphi etc... computing the
 * coefficients do not care about the exponentials (the parameters are still
 * needed), so simplify the amount of information to save shared memory (it is
 * crucial on AMD GPU to max out occupancy) */
template <typename T, typename T3>
__device__ __inline__ void
fill_smem_task_reduced(const kernel_params &dev, const int task_id,
                       smem_task_reduced<T, T3> &task) {
  if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0) {
    const auto &glb_task = dev.tasks[task_id];
    task.zetp = glb_task.zeta + glb_task.zetb;
    task.radius = glb_task.radius;
    task.discrete_radius = glb_task.discrete_radius;

    // angular momentum range for the actual collocate/integrate operation.
    task.lp =
        glb_task.la_max + dev.la_max_diff + glb_task.lb_max + dev.lb_max_diff;

    task.cube_size.x = glb_task.cube_size.x;
    task.cube_size.y = glb_task.cube_size.y;
    task.cube_size.z = glb_task.cube_size.z;

    task.cube_center.x = glb_task.cube_center.x;
    task.cube_center.y = glb_task.cube_center.y;
    task.cube_center.z = glb_task.cube_center.z;

    task.roffset.x = glb_task.roffset.x;
    task.roffset.y = glb_task.roffset.y;
    task.roffset.z = glb_task.roffset.z;

    task.lb_cube.x = glb_task.lb_cube.x;
    task.lb_cube.y = glb_task.lb_cube.y;
    task.lb_cube.z = glb_task.lb_cube.z;

    task.apply_border_mask = glb_task.apply_border_mask;
  }
  __syncthreads();
}

/*******************************************************************************
 * \brief Copies a task from global to shared memory and does precomputations.
 ******************************************************************************/

/*  computing the coefficients do not care about many of the exponentials
 * parameters, etc.. so */
template <typename T>
__device__ __inline__ void fill_smem_task_coef(const kernel_params &dev,
                                               const int task_id,
                                               smem_task<T> &task) {
  if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0) {
    const auto &glb_task = dev.tasks[task_id];
    const int iatom = glb_task.iatom;
    const int jatom = glb_task.jatom;
    task.zeta = glb_task.zeta;
    task.zetb = glb_task.zetb;

    for (int i = 0; i < 3; i++) {
      task.rab[i] = glb_task.rab[i];
      task.ra[i] = glb_task.ra[i];
      task.rb[i] = task.ra[i] + task.rab[i];
      task.rp[i] =
          task.ra[i] + task.rab[i] * task.zetb / (task.zeta + task.zetb);
    }

    task.prefactor = glb_task.prefactor;

    // task.radius = glb_task.radius;
    task.off_diag_twice = glb_task.off_diag_twice;

    // angular momentum range of basis set
    const int la_max_basis = glb_task.la_max;
    const int lb_max_basis = glb_task.lb_max;
    const int la_min_basis = glb_task.la_min;
    const int lb_min_basis = glb_task.lb_min;

    // angular momentum range for the actual collocate/integrate opteration.
    task.la_max = la_max_basis + dev.la_max_diff;
    task.lb_max = lb_max_basis + dev.lb_max_diff;
    task.la_min = max(la_min_basis + (int)dev.la_min_diff, 0);
    task.lb_min = max(lb_min_basis + (int)dev.lb_min_diff, 0);
    task.lp = task.la_max + task.lb_max;

    // start of decontracted set, ie. pab and hab
    task.first_coseta = (la_min_basis > 0) ? ncoset(la_min_basis - 1) : 0;
    task.first_cosetb = (lb_min_basis > 0) ? ncoset(lb_min_basis - 1) : 0;

    // size of decontracted set, ie. pab and hab
    task.ncoseta = ncoset(la_max_basis);
    task.ncosetb = ncoset(lb_max_basis);
    // size of the cab matrix
    task.n1 = ncoset(task.la_max);
    task.n2 = ncoset(task.lb_max);

    // size of entire spherical basis
    task.nsgfa = glb_task.nsgfa; // ibasis.nsgf;
    task.nsgfb = glb_task.nsgfb;

    // size of spherical set
    task.nsgf_seta = glb_task.nsgf_seta;
    task.nsgf_setb = glb_task.nsgf_setb;

    // strides of the sphi transformation matrices
    task.maxcoa = glb_task.maxcoa;
    task.maxcob = glb_task.maxcob;

    // transformations from contracted spherical to primitive cartesian basis
    task.sphia = &dev.sphi_dev[glb_task.ikind][glb_task.sgfa * task.maxcoa +
                                               glb_task.ipgf * task.ncoseta];
    task.sphib = &dev.sphi_dev[glb_task.jkind][glb_task.sgfb * task.maxcob +
                                               glb_task.jpgf * task.ncosetb];

    // Locate current matrix block within the buffer.
    const int block_offset = dev.block_offsets[glb_task.block_num];
    task.block_transposed = glb_task.block_transposed;
    task.pab_block = dev.ptr_dev[0] + block_offset + glb_task.subblock_offset;

    if (dev.ptr_dev[3] != nullptr) {
      task.hab_block = dev.ptr_dev[3] + block_offset + glb_task.subblock_offset;
      if (dev.ptr_dev[4] != nullptr) {
        task.forces_a = &dev.ptr_dev[4][3 * iatom];
        task.forces_b = &dev.ptr_dev[4][3 * jatom];
      }
    }
  }
  __syncthreads();
}

class smem_parameters {
private:
  int la_max_{-1};
  int lb_max_{-1};
  int la_min_{-1};
  int lb_min_{-1};
  int smem_per_block_{0};
  int cxyz_len_{-1};
  int alpha_len_{-1};
  int cab_len_{-1};
  int lp_max_{-1};
  ldiffs_value ldiffs_;
  int lp_diff_{-1};

public:
  smem_parameters(const ldiffs_value ldiffs, const int lmax) {
    ldiffs_ = ldiffs;
    lp_diff_ = ldiffs.la_max_diff + ldiffs.lb_max_diff;
    la_max_ = lmax + ldiffs.la_max_diff;
    lb_max_ = lmax + ldiffs.lb_max_diff;
    lp_max_ = la_max_ + lb_max_;

    cab_len_ = ncoset(lb_max_) * ncoset(la_max_);
    alpha_len_ = 3 * (lb_max_ + 1) * (la_max_ + 1) * (lp_max_ + 1);
    cxyz_len_ = ncoset(lp_max_);
    smem_per_block_ =
        std::max(cab_len_ + alpha_len_ + cxyz_len_, 64) * sizeof(double);

    if (smem_per_block_ > 64 * 1024) {
      fprintf(stderr,
              "ERROR: Not enough shared memory in grid_gpu_collocate.\n");
      fprintf(stderr, "cab_len: %i, ", cab_len_);
      fprintf(stderr, "alpha_len: %i, ", alpha_len_);
      fprintf(stderr, "cxyz_len: %i, ", cxyz_len_);
      fprintf(stderr, "total smem_per_block: %f kb\n\n",
              smem_per_block_ / 1024.0);
      abort();
    }
  }

  ~smem_parameters(){};

  // copy and move are trivial

  inline int smem_cab_offset() const { return cxyz_len_; }

  inline int smem_alpha_offset() const { return cab_len_ + cxyz_len_; }

  inline int smem_cxyz_offset() const { return 0; }

  inline int smem_per_block() const { return smem_per_block_; }

  inline int lp_diff() const { return lp_diff_; }

  inline ldiffs_value ldiffs() const { return ldiffs_; }

  inline int lp_max() const { return lp_max_; }

  inline int cxyz_len() const { return cxyz_len_; }
};

} // namespace rocm_backend
#endif
