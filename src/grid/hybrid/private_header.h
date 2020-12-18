/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#ifndef PRIVATE_HEADER_H
#define PRIVATE_HEADER_H

#include <cuda.h>

class Matrix {
  int row;
  int col;
  int ld;
  double *data;

public:
  __device__ __host__ Matrix(const int row_, const int col_, void *data) {
    row = row_;
    col = col_;
    ld = col_;
    data = (double *)data;
  }

  __device__ __host__ Matrix(const int row_, const int col_, const int ld_,
                             void *data_) {
    row = row_;
    col = col_;
    ld = ld_;
    data = (double *)data_;
  }

  __device__ __host__ ~Matrix() {}

  __device__ __host__ inline double &operator()(int const r, int const c) {
    return data[r * ld + c];
  }

  __device__ __host__ inline double operator()(int const r, int const c) const {
    if ((r >= row) || (c >= col))
      return 0.0;
    return data[r * ld + c];
  }

  __device__ __host__ inline double *at(int const r, int const c) const {
    return &data[r * ld + c];
  }

  __device__ __host__ inline double *at(int const r, int const c) {
    return &data[r * ld + c];
  }
};

class tensor3 {
  int l1, l2, l3;
  int ld;
  double *data;

public:
  __device__ __host__ tensor3(int const l1_, int const l2_, int const l3_,
                              int const ld_, void *data_) {
    l1 = l1_;
    l2 = l2_;
    l3 = l3_;
    ld = ld_;

    data = (double *)data_;
  }

  __device__ __host__ ~tensor3() {}

  __device__ __host__ inline double operator()(int const i, int const j,
                                               int const k) {
    if ((i >= l1) || (j >= l2) || (k >= l3))
      return 0.0;
    return data[(i * l2 + j) * ld + k];
  }

  __device__ __host__ inline double operator()(int const i, int const j,
                                               int const k) const {
    if ((i >= l1) || (j >= l2) || (k >= l3))
      return 0;
    return data[(i * l2 + j) * ld + k];
  }

  __device__ __host__ inline double *at(int const i, int const j,
                                        int const k) const {
    return &data[(i * l2 + j) * ld + k];
  }

  __device__ __host__ inline double *at(int const i, int const j, int const k) {
    return &data[(i * l2 + j) * ld + k];
  }
};

struct tensor4 {
  int l1, l2, l3, l4;
  int ld;
  double *data;

  __device__ __host__ tensor4(int const l1_, int const l2_, int const l3_,
                              int const l4_, void *data_) {
    l1 = l1_;
    l2 = l2_;
    l3 = l3_;
    l4 = l4_;
    ld = l4_;

    data = (double *)data_;
  }

  __device__ __host__ tensor4(int const l1_, int const l2_, int const l3_,
                              int const l4_, int const ld_, void *data_) {
    l1 = l1_;
    l2 = l2_;
    l3 = l3_;
    l4 = l4_;
    ld = ld_;

    data = (double *)data_;
  }

  __device__ __host__ ~tensor4() {}

  __device__ __host__ int size() const { return l1 * l2 * l3 * ld; }
  __device__ __host__ inline double &operator()(int const i, int const j,
                                                int const k, int const l) {
    return data[((i * l2 + j) * l3 + k) * ld + l];
  }

  __device__ __host__ inline double operator()(int const i, int const j,
                                               int const k, int const l) const {
    if ((i >= l1) || (j >= l2) || (k >= l3) || (l >= l4))
      return 0.0;
    return data[((i * l2 + j) * l3 + k) * ld + l];
  }

  __device__ __host__ inline double *at(int const i, int const j, int const k,
                                        int const l) const {
    return &data[((i * l2 + j) * l3 + k) * ld + l];
  }

  __device__ __host__ inline double *at(int const i, int const j, int const k,
                                        int const l) {
    return &data[((i * l2 + j) * l3 + k) * ld + l];
  }

  __device__ void zero(cooperative_groups::thread_block &block) {
    for (int i = block.thread_rank(); i < block.size(); i++)
      data[i] = 0.0;
  }
};

#if __CUDA_ARCH__ < 600
__device__ double atomicAdd1(double *address, double val) {
  unsigned long long int *address_as_ull = (unsigned long long int *)address;
  unsigned long long int old = *address_as_ull, assumed;
  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,
                    __double_as_longlong(val + __longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
}
#endif

__inline__ __device__ void
return_cube_position(const int3 *__restrict__ const cube_center,
                     const int3 *__restrict__ const lower_boundaries_cube,
                     const int3 period, int3 *const position) {
  position->x = (cube_center->x + lower_boundaries_cube->x) % period.x;
  position->y = (cube_center->y + lower_boundaries_cube->y) % period.y;
  position->z = (cube_center->z + lower_boundaries_cube->z) % period.z;

  if (position->x < 0)
    position->x += period.x;
  if (position->y < 0)
    position->y += period.y;
  if (position->x < 0)
    position->y += period.z;
}

__inline__ __device__ void
convert_to_lattice_coordinates(const double3 *__restrict__ const rp,
                               double3 *__restrict__ rp_c) {
  rp_c->x = dh_inv_[0] * rp->x + dh_inv_[3] * rp->y + dh_inv_[6] * rp->z;
  rp_c->y = dh_inv_[1] * rp->x + dh_inv_[4] * rp->y + dh_inv_[7] * rp->z;
  rp_c->z = dh_inv_[2] * rp->x + dh_inv_[5] * rp->y + dh_inv_[8] * rp->z;
}

__inline__ __device__ void convert_from_lattice_coordinates_to_cartesian(
    const double3 *__restrict__ const rp, double3 *__restrict__ rp_c) {
  rp_c->x = dh_[0] * rp->x + dh_[3] * rp->y + dh_[6] * rp->z;
  rp_c->y = dh_[1] * rp->x + dh_[4] * rp->y + dh_[7] * rp->z;
  rp_c->z = dh_[2] * rp->x + dh_[5] * rp->y + dh_[8] * rp->z;
}

__device__ void compute_cube_properties(const double radius,
                                        const double3 *__restrict__ rp,
                                        double3 *__restrict__ roffset,
                                        int3 *__restrict__ cubecenter,
                                        int3 *__restrict__ lb_cube,
                                        int3 *__restrict__ cube_size) {
  int3 ub_cube;

  /* center of the gaussian in the lattice coordinates */
  double3 rp1;

  /* it is in the lattice vector frame */
  convert_to_lattice_coordinates(rp, &rp1);

  cubecenter->x = floor(rp1.x);
  cubecenter->y = floor(rp1.y);
  cubecenter->z = floor(rp1.z);

  if (is_orthorhombic_[0]) {
    /* seting up the cube parameters */
    const double3 dr = {.x = dh_[0], .y = dh_[4], .z = dh_[8]};
    const double3 dr_inv = {.x = dh_inv_[0], .y = dh_inv_[4], .z = dh_inv_[8]};
    /* cube center */

    /* lower and upper bounds */

    // Historically, the radius gets discretized.
    const double drmin = min(dr.x, min(dr.y, dr.z));
    const double disr_radius = drmin * max(1.0, ceil(radius / drmin));

    roffset->x = rp->x - cubecenter->x * dr.x;
    roffset->y = rp->y - cubecenter->y * dr.y;
    roffset->z = rp->z - cubecenter->z * dr.z;

    lb_cube->x = ceil(-1e-8 - disr_radius * dr_inv.x);
    lb_cube->y = ceil(-1e-8 - disr_radius * dr_inv.y);
    lb_cube->z = ceil(-1e-8 - disr_radius * dr_inv.z);

    roffset->x *= dr_inv.x;
    roffset->y *= dr_inv.y;
    roffset->z *= dr_inv.z;

    // Symetric interval
    ub_cube.x = 1 - lb_cube->x;
    ub_cube.y = 1 - lb_cube->y;
    ub_cube.z = 1 - lb_cube->z;

  } else {

    lb_cube->x = INT_MAX;
    ub_cube.x = INT_MIN;
    lb_cube->y = INT_MAX;
    ub_cube.y = INT_MIN;
    lb_cube->z = INT_MAX;
    ub_cube.z = INT_MIN;

    for (int i = -1; i <= 1; i++) {
      for (int j = -1; j <= 1; j++) {
        for (int k = -1; k <= 1; k++) {
          double3 r = make_double3(((double)i) * radius, ((double)j) * radius,
                                   ((double)k) * radius);
          convert_to_lattice_coordinates(&r, roffset);

          lb_cube->x = min(lb_cube->x, (int)floor(roffset->x));
          ub_cube.x = max(ub_cube.x, (int)ceil(roffset->x));

          lb_cube->y = min(lb_cube->y, (int)floor(roffset->y));
          ub_cube.y = max(ub_cube.y, (int)ceil(roffset->y));

          lb_cube->z = min(lb_cube->z, (int)floor(roffset->z));
          ub_cube.z = max(ub_cube.z, (int)ceil(roffset->z));
        }
      }
    }

    /* compute the offset in lattice coordinates */

    roffset->x = rp1.x - cubecenter->x;
    roffset->y = rp1.y - cubecenter->y;
    roffset->z = rp1.z - cubecenter->z;
  }

  /* compute the cube size ignoring periodicity */
  cube_size->x = ub_cube.x - lb_cube->x + 1;
  cube_size->y = ub_cube.y - lb_cube->y + 1;
  cube_size->z = ub_cube.z - lb_cube->z + 1;
}

__inline__ __device__ void
compute_window_size(const int3 *const grid_size,
                    const int3 *const lower_corner_,
                    const int3 *const period_, /* also full size of the grid */
                    const int border_mask, const int3 *border_width,
                    int3 *const window_size, int3 *const window_shift) {
  window_shift->x = 0;
  window_shift->y = 0;
  window_shift->z = 0;

  window_size->x = grid_size->x;
  window_size->y = grid_size->y;
  window_size->z = grid_size->z;

  if (grid_size->x != period_->x)
    window_size->x--;

  if (grid_size->y != period_->y)
    window_size->y--;

  if (grid_size->z != period_->z)
    window_size->z--;

  if ((grid_size->x != period_->x) || (grid_size->y != period_->y) ||
      (grid_size->z != period_->z)) {
    if (border_mask & (1 << 0))
      window_shift->x += border_width->x;
    if (border_mask & (1 << 1))
      window_size->x -= border_width->x;
    if (border_mask & (1 << 2))
      window_shift->y += border_width->y;
    if (border_mask & (1 << 3))
      window_size->y -= border_width->y;
    if (border_mask & (1 << 4))
      window_shift->z += border_width->z;
    if (border_mask & (1 << 5))
      window_size->z -= border_width->z;
  }
}

#endif
