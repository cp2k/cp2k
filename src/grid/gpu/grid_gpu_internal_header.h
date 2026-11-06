/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

/*
 * Authors :
 - Dr Mathieu Taillefumier (ETH Zurich / CSCS)
 - Advanced Micro Devices, Inc.
*/

#ifndef GRID_GPU_INTERNAL_HEADER_H
#define GRID_GPU_INTERNAL_HEADER_H

#include <algorithm>
#include <assert.h>
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

#include "grid_gpu_context.h"

namespace rocm_backend {

#if defined(__OFFLOAD_CUDA) || defined(__HIP_PLATFORM_NVIDIA__)
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
  // signed char keeps coset_inv at 4 kB of constant memory (16 kB with int).
  // The components are promoted to int in any arithmetic expression.
  signed char l[3];
};
static_assert(sizeof(orbital) == 3, "orbital is expected to be packed");

// largest angular momentum covered by the tables of the backend (ncoset,
// coset_inv). It bounds la + lb of every operation.
constexpr int grid_gpu_max_l = 18;

__constant__ orbital coset_inv[1330];

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
  // radius: true cutoff used for the sphere-membership test.
  // discrete_radius: box-aligned cutoff (always >= radius) used only for
  // sizing the cube / trimming box boundaries.
  T radius, discrete_radius;
  T norm_lattice_vector_z_2, norm_lattice_vector_y_2;
  T norm_inverse_lattice_vector_y, norm_inverse_lattice_vector_x;
  int3 cube_center, lb_cube, cube_size, window_size, window_shift;
  T3 roffset;
  T zetp;
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
        coset_inv_host[i] = {{static_cast<signed char>(lx),
                              static_cast<signed char>(ly),
                              static_cast<signed char>(lz)}};
      }
    }
  }

  offloadMemcpyToSymbol(coset_inv, &coset_inv_host, sizeof(coset_inv_host));

  initialized = true;
}

// calculate the global index of a thread block
__inline__ __device__ unsigned int block_index() {
  return blockIdx.x + gridDim.x * (blockIdx.y + gridDim.y * blockIdx.z);
}

// Calculating the global index in the grid of any given device thread
__inline__ __device__ unsigned int thread_global_index() {
  return threadIdx.x + blockDim.x * (threadIdx.y + blockDim.y * threadIdx.z);
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
__inline__ __device__ void compute_alpha(const smem_task<T> &task,
                                         T *__restrict__ alpha) {
  // strides for accessing alpha
  const int s3 = (task.lp + 1);
  const int s2 = (task.la_max + 1) * s3;
  const int s1 = (task.lb_max + 1) * s2;
  const int tid = thread_global_index();

  for (int i = tid; i < 3 * s1; i += blockDim.x * blockDim.y * blockDim.z)
    alpha[i] = 0.0;

  __syncthreads();

  for (int idir = threadIdx.z; idir < 3; idir += blockDim.z) {
    const T drpa = task.rp[idir] - task.ra[idir];
    const T drpb = task.rp[idir] - task.rb[idir];
    for (int la = threadIdx.y; la <= task.la_max; la += blockDim.y) {
      for (int lb = threadIdx.x; lb <= task.lb_max; lb += blockDim.x) {
        // binomial coefficients by the recurrence
        // C(n, k + 1) = C(n, k) * (n - k) / (k + 1), exact in integers.
        T a = 1.0;
        int binom_a = 1;
        for (int k = 0; k <= la; k++) {
          T b = 1.0;
          int binom_b = 1;
          const int base = idir * s1 + lb * s2 + la * s3;
          for (int l = 0; l <= lb; l++) {
            alpha[base + la - l + lb - k] += a * b * binom_a * binom_b;
            b *= drpb;
            binom_b = binom_b * (lb - l) / (l + 1);
          }
          a *= drpa;
          binom_a = binom_a * (la - k) / (k + 1);
        }
      }
    }
  }
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
        (std::max(1,
                  (int)ceil(radius / std::min(norm1, std::min(norm2, norm3)))));

    rp2.x = cubecenter->x;
    rp2.y = cubecenter->y;
    rp2.z = cubecenter->z;

    /* convert the cube center from grid points coordinates to cartesian */
    convert_from_lattice_coordinates_to_cartesian(dh_, &rp2, roffset);
    /* cube center */
    roffset->x -= rp->x;
    roffset->y -= rp->y;
    roffset->z -= rp->z;

    /* This potentially makes the cube bigger than necessary because disc_radius
       is derived with the smallest lattice parameter. It is not an issue when
       the lattice is cubic but could lead to more calculations than necessary
       when the lattice is orthorhombic */
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

__inline__ __device__ void compute_window_size(const int3 grid_size,
                                               const int border_mask,
                                               const int3 &border_width,
                                               int3 &window_size,
                                               int3 &window_shift) {
  window_shift.x = 0;
  window_shift.y = 0;
  window_shift.z = 0;

  window_size.x = grid_size.x - 1;
  window_size.y = grid_size.y - 1;
  window_size.z = grid_size.z - 1;

  if (border_mask & (1 << 0))
    window_shift.x += border_width.x;
  if (border_mask & (1 << 1))
    window_size.x -= border_width.x;
  if (border_mask & (1 << 2))
    window_shift.y += border_width.y;
  if (border_mask & (1 << 3))
    window_size.y -= border_width.y;
  if (border_mask & (1 << 4))
    window_shift.z += border_width.z;
  if (border_mask & (1 << 5))
    window_size.z -= border_width.z;
}

/*******************************************************************************
 * \brief Wraps a cube-local coordinate into the periodic full grid along one
 *        axis. Single point of entry shared by collocate_kernel and
 *        integrate_kernel.
 ******************************************************************************/
__inline__ __device__ int wrap_grid_index(const int idx, const int full_size) {
  int idx2 = idx % full_size;
  if (idx2 < 0)
    idx2 += full_size;
  return idx2;
}

/*******************************************************************************
 * \brief Shifts task.cube_center into the local grid's coordinate system and,
 *        for distributed grids, derives the border window. Single point of
 *        entry shared by collocate_kernel and integrate_kernel so the two
 *        cannot drift apart on which axis maps to which.
 ******************************************************************************/
template <typename T, typename T3, bool distributed__>
__device__ __inline__ void
setup_task_cube_center(const kernel_params &dev_,
                       smem_task_reduced<T, T3> &task) {
  task.cube_center.x += task.lb_cube.x - dev_.grid_lower_corner_.x;
  task.cube_center.y += task.lb_cube.y - dev_.grid_lower_corner_.y;
  task.cube_center.z += task.lb_cube.z - dev_.grid_lower_corner_.z;

  if (distributed__) {
    if (task.apply_border_mask) {
      compute_window_size(
          dev_.grid_local_size_,
          dev_.tasks[dev_.first_task + block_index()].border_mask,
          dev_.grid_border_width_, task.window_size, task.window_shift);
    }
  }
}

/*******************************************************************************
 * \brief Transforms coefficients C_ab into C_xyz.
 ******************************************************************************/
template <typename T>
__device__ __inline__ static void
cab_to_cxyz(const smem_task<T> &task, const T *__restrict__ alpha,
            const T *__restrict__ cab, T *__restrict__ cxyz) {

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
        const T p = alpha[0 * s1 + b.l[0] * s2 + a.l[0] * s3 + co.l[0]] *
                    alpha[1 * s1 + b.l[1] * s2 + a.l[1] * s3 + co.l[1]] *
                    alpha[2 * s1 + b.l[2] * s2 + a.l[2] * s3 + co.l[2]];
        reg += p * cab[jco * task.n1 + ico]; // collocate
      }
    }

    cxyz[i] = task.prefactor * reg;
  }
  __syncthreads(); // because of concurrent writes to cxyz / cab
}

/*******************************************************************************
 * \brief Transforms coefficients C_xyz into C_ab.
 ******************************************************************************/
template <typename T>
__device__ __inline__ static void
cxyz_to_cab(const smem_task<T> &task, const T *__restrict__ alpha,
            const T *__restrict__ cxyz, T *__restrict__ cab) {

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
        const T p = alpha[b.l[0] * s2 + a.l[0] * s3 + co.l[0]] *
                    alpha[s1 + b.l[1] * s2 + a.l[1] * s3 + co.l[1]] *
                    alpha[2 * s1 + b.l[2] * s2 + a.l[2] * s3 + co.l[2]];

        reg += p * cxyz[ic]; // integrate
      }
      cab[jco * task.n1 + ico] =
          task.prefactor * reg; // partial loop coverage -> zero it
    }
  }
}

/*******************************************************************************
 * \brief Factorized cab <-> cxyz transformations.
 *
 *        In each direction alpha expands (x - a)^la (x - b)^lb around p, so
 *        cab_to_cxyz is the product of three steps: move the expansion center
 *        of the a index from a to p, of the b index from b to p, then multiply
 *        the polynomials,
 *
 *          (x - a)^n = sum_k C(n, k) (p - a)^k (x - p)^(n - k)
 *          cxyz[c]   = prefactor sum_{a' + b' = c} cab_p[b'][a']
 *
 *        The center moves are triangular, separable in x, y and z and map cab
 *        onto itself, so they are done in place, one 1D line of cab per
 *        thread. cxyz_to_cab is the transpose: a gather followed by the
 *        transposed moves. Neither needs alpha.
 ******************************************************************************/

// decodes t in [0, (l + 1) (l + 2) / 2) into (p, q) with p + q <= l
__device__ __inline__ void decode_pair(int t, int &p, int &q) {
  int s = 0;
  while (t > s) {
    t -= s + 1;
    s++;
  }
  p = t;
  q = s - t;
}

__device__ __inline__ int line_index(const int dir, const int m, const int p,
                                     const int q, const bool on_a,
                                     const int other, const int n1) {
  int l[3];
  l[dir] = m;
  l[(dir + 1) % 3] = p;
  l[(dir + 2) % 3] = q;
  const int c = coset(l[0], l[1], l[2]);
  return on_a ? other * n1 + c : c * n1 + other;
}

/*******************************************************************************
 * \brief Moves the expansion center of the a (on_a) or b cartesian index of
 *        cab along direction dir, in place:
 *
 *          forward:   new[m] = sum_{k >= 0} C(m + k, k) dr^k old[m + k]
 *          transpose: new[m] = sum_{k <= m} C(m, k)     dr^k old[m - k]
 *
 *        Each thread owns a whole line (all other indices fixed). The forward
 *        update is a Taylor shift of the polynomial sum_m old[m] x^m by dr,
 *        done in place by repeated synthetic division (Ruffini-Horner):
 *          for i = 0 .. len-1, for j = len-1 down to i: c[j] += dr c[j+1]
 *        The transpose applies the transposed elementary updates in reverse
 *        order:
 *          for i = len-1 down to 0, for j = i .. len-1: c[j+1] += dr c[j]
 *        No binomial coefficient is needed, len (len + 1) / 2 FMA per line.
 ******************************************************************************/
template <typename T, bool TRANSPOSE>
__device__ __inline__ void move_center(const int la, const int lb, const int n1,
                                       const int dir, const bool on_a,
                                       const T dr, T *cab) {
  const int tid = thread_global_index();
  const int nthreads = blockDim.x * blockDim.y * blockDim.z;
  const int l = on_a ? la : lb;
  const int nother = on_a ? ncoset(lb) : ncoset(la);
  const int ntri = (l + 1) * (l + 2) / 2;

  for (int line = tid; line < nother * ntri; line += nthreads) {
    const int other = line / ntri;
    int p, q;
    decode_pair(line - other * ntri, p, q);
    const int len = l - p - q;
    if (len == 0)
      continue; // new[0] = old[0]
    if (!TRANSPOSE) {
      for (int i = 0; i < len; i++) {
        // c[j] += dr c[j + 1] for j descending, c[j + 1] carried in a register
        T carry = cab[line_index(dir, len, p, q, on_a, other, n1)];
        for (int j = len - 1; j >= i; j--) {
          const int idx = line_index(dir, j, p, q, on_a, other, n1);
          carry = cab[idx] + dr * carry;
          cab[idx] = carry;
        }
      }
    } else {
      for (int i = len - 1; i >= 0; i--) {
        // c[j + 1] += dr c[j] for j ascending, c[j] carried in a register
        T carry = cab[line_index(dir, i, p, q, on_a, other, n1)];
        for (int j = i + 1; j <= len; j++) {
          const int idx = line_index(dir, j, p, q, on_a, other, n1);
          carry = cab[idx] + dr * carry;
          cab[idx] = carry;
        }
      }
    }
  }
}

/*******************************************************************************
 * \brief Same result as cab_to_cxyz. cab is overwritten.
 ******************************************************************************/
template <typename T>
__device__ __inline__ void cab_to_cxyz_factorized(const smem_task<T> &task,
                                                  T *__restrict__ cab,
                                                  T *__restrict__ cxyz) {
  const int la = task.la_max, lb = task.lb_max, n1 = task.n1;
  for (int dir = 0; dir < 3; dir++) {
    move_center<T, false>(la, lb, n1, dir, true, task.rp[dir] - task.ra[dir],
                          cab);
    __syncthreads();
  }
  for (int dir = 0; dir < 3; dir++) {
    move_center<T, false>(la, lb, n1, dir, false, task.rp[dir] - task.rb[dir],
                          cab);
    __syncthreads();
  }

  // cxyz[c] = prefactor sum_{a' + b' = c} cab_p[b'][a']
  const int tid = thread_global_index();
  const int nthreads = blockDim.x * blockDim.y * blockDim.z;
  for (int ic = tid; ic < ncoset(task.lp); ic += nthreads) {
    const auto &c = coset_inv[ic];
    T val = 0.0;
    for (int ax = 0; ax <= min(c.l[0], la); ax++) {
      for (int ay = 0; ay <= min(c.l[1], la - ax); ay++) {
        for (int az = 0; az <= min(c.l[2], la - ax - ay); az++) {
          const int bx = c.l[0] - ax, by = c.l[1] - ay, bz = c.l[2] - az;
          if (bx + by + bz <= lb)
            val += cab[coset(bx, by, bz) * n1 + coset(ax, ay, az)];
        }
      }
    }
    cxyz[ic] = task.prefactor * val;
  }
  __syncthreads();
}

/*******************************************************************************
 * \brief Same result as cxyz_to_cab.
 ******************************************************************************/
template <typename T>
__device__ __inline__ void cxyz_to_cab_factorized(const smem_task<T> &task,
                                                  const T *__restrict__ cxyz,
                                                  T *__restrict__ cab) {
  const int la = task.la_max, lb = task.lb_max, n1 = task.n1;
  const int tid = thread_global_index();
  const int nthreads = blockDim.x * blockDim.y * blockDim.z;

  // cab_p[b'][a'] = prefactor cxyz[a' + b']
  for (int k = tid; k < n1 * task.n2; k += nthreads) {
    const int jco = k / n1;
    const int ico = k - jco * n1;
    const auto &a = coset_inv[ico];
    const auto &b = coset_inv[jco];
    cab[k] = task.prefactor *
             cxyz[coset(a.l[0] + b.l[0], a.l[1] + b.l[1], a.l[2] + b.l[2])];
  }
  __syncthreads();

  for (int dir = 0; dir < 3; dir++) {
    move_center<T, true>(la, lb, n1, dir, true, task.rp[dir] - task.ra[dir],
                         cab);
    __syncthreads();
  }
  for (int dir = 0; dir < 3; dir++) {
    move_center<T, true>(la, lb, n1, dir, false, task.rp[dir] - task.rb[dir],
                         cab);
    __syncthreads();
  }
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

    task.norm_lattice_vector_z_2 = dev.dh_[6] * dev.dh_[6] +
                                   dev.dh_[7] * dev.dh_[7] +
                                   dev.dh_[8] * dev.dh_[8];
    task.norm_lattice_vector_y_2 = dev.dh_[3] * dev.dh_[3] +
                                   dev.dh_[4] * dev.dh_[4] +
                                   dev.dh_[5] * dev.dh_[5];

    task.norm_inverse_lattice_vector_y =
        sqrt(dev.dh_inv_[3] * dev.dh_inv_[3] + dev.dh_inv_[4] * dev.dh_inv_[4] +
             dev.dh_inv_[5] * dev.dh_inv_[5]);
    task.norm_inverse_lattice_vector_x =
        sqrt(dev.dh_inv_[0] * dev.dh_inv_[0] + dev.dh_inv_[1] * dev.dh_inv_[1] +
             dev.dh_inv_[2] * dev.dh_inv_[2]);

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
    task.pab_block =
        dev.buffers_dev.pab_block + block_offset + glb_task.subblock_offset;

    if (dev.buffers_dev.hab_block != nullptr) {
      task.hab_block =
          dev.buffers_dev.hab_block + block_offset + glb_task.subblock_offset;
      if (dev.buffers_dev.forces != nullptr) {
        task.forces_a = &dev.buffers_dev.forces[3 * iatom];
        task.forces_b = &dev.buffers_dev.forces[3 * jatom];
      }
    }
  }
  __syncthreads();
}

/*******************************************************************************
 * \brief Size (in doubles) of the work arrays of the two step spherical <->
 *        cartesian transformations of calculate_coefficients_two_step and
 *        compute_hab_two_step: [max_nsgf_set][ld] and [ld][ld] with
 *        ld = ncoset(lmax).
 ******************************************************************************/
inline int sphi_work_size(const int lmax, const int max_nsgf_set) {
  const int ld = ncoset(lmax);
  return max_nsgf_set * ld + ld * ld;
}

// shared memory budget of the two step kernels (default limit without opting
// in on NVIDIA hardware). Above, the fused kernels are used.
constexpr size_t sphi_two_step_max_smem = 48 * 1024;

class smem_parameters {
private:
  int la_max_{-1};
  int lb_max_{-1};
  int smem_per_block_{0};
  int alpha_size_{-1};
  int cab_size_{-1};
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

    if (lp_max_ > grid_gpu_max_l) {
      fprintf(stderr,
              "ERROR: the GPU grid backend supports la + lb <= %d, this "
              "operation needs %d (lmax %d, la_max_diff %d, lb_max_diff %d).\n",
              grid_gpu_max_l, lp_max_, lmax, ldiffs.la_max_diff,
              ldiffs.lb_max_diff);
      abort();
    }

    // NB: cab is allocated in global memory not shared memory. Each block has
    // its own cab space

    cab_size_ = (rocm_backend::ncoset(la_max_) * rocm_backend::ncoset(lb_max_));
    alpha_size_ = 3 * (lb_max_ + 1) * (la_max_ + 1) * (lp_max_ + 1);
    smem_per_block_ = std::max(alpha_size_, 64) * sizeof(double);

    if (smem_per_block_ > 64 * 1024) {
      fprintf(stderr,
              "ERROR: Not enough shared memory in grid_gpu_collocate.\n");
      fprintf(stderr, "alpha_len: %i, ", alpha_size_);
      fprintf(stderr, "total smem_per_block: %f kB\n\n",
              smem_per_block_ / 1024.0);
      abort();
    }
  }

  ~smem_parameters() {};

  // copy and move are trivial

  inline int alpha_size() const { return alpha_size_; }

  inline int cab_size() const { return cab_size_; }

  inline int smem_per_block() const { return smem_per_block_; }

  inline int lp_diff() const { return lp_diff_; }

  inline ldiffs_value ldiffs() const { return ldiffs_; }

  inline int lp_max() const { return lp_max_; }

  inline int cxyz_size() const { return ncoset(lp_max_); }
};
template <typename T>
__inline__ __device__ T *allocate_workspace(const kernel_params &dev_) {
  unsigned int offset = dev_.cab_block_offset_dev[block_index()];
  return (T *)(dev_.buffers_dev.cab + offset);
}

template <typename T, typename T3>
__device__ __inline__ T
calculate_ymix_ymax_boundaries(smem_task_reduced<T, T3> &task, const int z,
                               int &ymin, int &ymax) {
  T kremain = 0.0;
  ymin = (2 * (z + task.lb_cube.z) - 1) / 2;
  ymin *= ymin;
  kremain = task.discrete_radius * task.discrete_radius -
            ((T)ymin) * task.norm_lattice_vector_z_2;
  ymin = ceil(-1.0e-8 -
              sqrt(fmax(0.0, kremain)) * task.norm_inverse_lattice_vector_y);
  ymax = 1 - ymin - task.lb_cube.y;
  ymin = ymin - task.lb_cube.y;
  return kremain;
}

template <typename T, typename T3>
__device__ __inline__ void
calculate_xmin_xmax_boundaries(smem_task_reduced<T, T3> &task, const int y,
                               const T kremain, int &xmin, int &xmax) {
  xmin = (2 * (y + task.lb_cube.y) - 1) / 2;
  xmin *= xmin;
  xmin = ceil(-1.0e-8 -
              sqrt(fmax(0.0, kremain - xmin * task.norm_lattice_vector_y_2)) *
                  task.norm_inverse_lattice_vector_x);
  xmax = 1 - xmin - task.lb_cube.x;
  xmin -= task.lb_cube.x;
}

} // namespace rocm_backend
#endif
