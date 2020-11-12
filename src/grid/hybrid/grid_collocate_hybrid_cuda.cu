/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#ifdef __GRID_CUDA

#include <assert.h>
#include <cooperative_groups.h>
#include <cuda.h>

#include "../cpu/collocation_integration.h"

extern "C" void reset_list_gpu(pgf_list_gpu *lst);
extern "C" void return_dh(void *const ptr, const int level, double *const dh);
extern "C" void return_dh_inv(void *const ptr, const int level,
                              double *const dh);
extern "C" int return_num_devs(void *const ptr);
extern "C" int return_device_id(void *const ptr, const int device_id);
extern "C" int is_grid_orthorhombic(void *const ptr);

namespace cg = cooperative_groups;

__constant__ __device__ double dh_[9];
__constant__ __device__ double dh_inv_[9];
__constant__ __device__ int is_orthorhombic_[1];
/* lower corner of the grid block when the grid is divided over multiple mpi
 * ranks */
__constant__ __device__ int grid_lower_boundaries_[3];

extern __shared__ double array[];

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

template <typename T> __device__ T reduce(cg::thread_group g, T *temp, T val) {
  int lane = g.thread_rank();

  // Each iteration halves the number of active threads
  // Each thread adds its partial sum[i] to sum[lane+i]
  for (int i = g.size() / 2; i > 0; i /= 2) {
    temp[lane] = val;
    g.sync(); // wait for all threads to store
    if (lane < i)
      val += temp[lane + i];
    g.sync(); // wait for all threads to load
  }
  return val; // note: only thread 0 will return full sum
}

__device__ void compute_cube_properties(const double radius,
                                        const double3 *__restrict__ rp,
                                        double *__restrict__ disr_radius,
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
    *disr_radius = drmin * max(1.0, ceil(radius / drmin));

    roffset->x = rp->x - cubecenter->x * dr.x;
    roffset->y = rp->y - cubecenter->y * dr.y;
    roffset->z = rp->z - cubecenter->z * dr.z;

    lb_cube->x = ceil(-1e-8 - *disr_radius * dr_inv.x);
    lb_cube->y = ceil(-1e-8 - *disr_radius * dr_inv.y);
    lb_cube->z = ceil(-1e-8 - *disr_radius * dr_inv.z);

    roffset->x /= dr.x;
    roffset->y /= dr.y;
    roffset->z /= dr.z;

    // Symetric interval
    ub_cube.x = -lb_cube->x;
    ub_cube.y = -lb_cube->y;
    ub_cube.z = -lb_cube->z;

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
          // roffset.x = dh_inv[0][0] * x[0] + dh_inv[1][idir] * x[1] +
          // dh_inv[2][idir] * x[2];
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

__global__ void compute_collocation_gpu_(
    const int apply_cutoff, const int3 grid_size_,
    const int3 grid_lower_corner_pos_, const int3 period_,
    const int *__restrict__ border_mask, const int3 border_width,
    const int *__restrict__ lmax_gpu_, const double *__restrict__ zeta_gpu,
    const double3 *__restrict__ rp, const double *__restrict__ radius_gpu_,
    const int *__restrict__ coef_offset_gpu_,
    const double *__restrict__ coef_gpu_, double *__restrict__ grid_gpu_) {
  /* the period is sotred in constant memory */
  /* the displacement vectors as well */

  int lmax = lmax_gpu_[blockIdx.x];

  int3 cube_size, cube_center, lb_cube, window_size, window_shift;

  double3 roffset;
  double disr_radius = 0;
  const double radius = radius_gpu_[blockIdx.x];

  compute_cube_properties(radius, rp + blockIdx.x, &disr_radius, &roffset,
                          &cube_center, &lb_cube, &cube_size);

  compute_window_size(&grid_size_, &grid_lower_corner_pos_,
                      &period_, /* also full size of the grid */
                      border_mask[blockIdx.x], &border_width, &window_size,
                      &window_shift);

  const double *__restrict__ coef = coef_gpu_ + coef_offset_gpu_[blockIdx.x];

  const double zeta = zeta_gpu[blockIdx.x];

  double *coefs_ = (double *)array;

  int id = (threadIdx.z * blockDim.y + threadIdx.y) * blockDim.x + threadIdx.x;

  for (int i = id; i < ((lmax + 1) * (lmax + 1) * (lmax + 1));
       i += (blockDim.x * blockDim.y * blockDim.z))
    coefs_[i] = coef[i];
  __syncthreads();

  for (int z = threadIdx.z; z <= cube_size.z; z += blockDim.z) {
    const double z1 = z + lb_cube.z - roffset.z;
    const int z2 = (z + cube_center.z + lb_cube.z - grid_lower_corner_pos_.z +
                    32 * period_.z) %
                   period_.z;

    /* check if the point is within the window */
    if ((z2 < window_shift.z) || (z2 > window_size.z)) {
      continue;
    }

    for (int y = threadIdx.y; y <= cube_size.y; y += blockDim.y) {
      double y1 = y + lb_cube.y - roffset.y;
      const int y2 = (y + lb_cube.y + cube_center.y - grid_lower_corner_pos_.y +
                      32 * period_.y) %
                     period_.y;

      /* check if the point is within the window */
      if ((y2 < window_shift.y) || (y2 > window_size.y)) {
        continue;
      }

      for (int x = threadIdx.x; x <= cube_size.x; x += blockDim.x) {
        const double x1 = x + lb_cube.x - roffset.x;
        const int x2 = (x + lb_cube.x + cube_center.x -
                        grid_lower_corner_pos_.x + 32 * period_.x) %
                       period_.x;

        /* check if the point is within the window */
        if ((x2 < window_shift.x) || (x2 > window_size.x)) {
          continue;
        }

        double3 r3;
        r3.x = z1 * dh_[6] + y1 * dh_[3] + x1 * dh_[0];
        r3.y = z1 * dh_[7] + y1 * dh_[4] + x1 * dh_[1];
        r3.z = z1 * dh_[8] + y1 * dh_[5] + x1 * dh_[2];

        if (apply_cutoff &&
            ((radius * radius) < (r3.x * r3.x + r3.y * r3.y + r3.z * r3.z)))
          continue;

        /* compute the coordinates of the point in atomic coordinates */
        double exp_factor =
            exp(-(r3.x * r3.x + r3.y * r3.y + r3.z * r3.z) * zeta);

        double res = 0.0;
        double dx = 1;

        /* NOTE: the coefficients are stored as lx,lz,ly */
        /* It is suboptimal right now because i do more operations than needed
         * (a lot of coefs_ are zero). Moreover, it is a dgemm underneath and
         * could be treated with tensor cores */

        for (int alpha = 0; alpha <= lmax; alpha++) {
          double dz = 1;
          for (int gamma = 0; gamma <= (lmax - alpha); gamma++) {
            double dy = dx * dz;
            const int off = (alpha * (lmax + 1) + gamma) * (lmax + 1);
            for (int beta = 0; beta <= (lmax - alpha - gamma); beta++) {
              res += coefs_[off + beta] * dy;
              dy *= r3.y;
            }
            dz *= r3.z;
          }
          dx *= r3.x;
        }

        res *= exp_factor;

#if __CUDA_ARCH__ < 600
        atomicAdd1(&grid_gpu_[(z2 * grid_size_.y + y2) * grid_size_.x + x2],
                   res);
#else
        atomicAdd(&grid_gpu_[(z2 * grid_size_.y + y2) * grid_size_.x + x2],
                  res);
#endif
      }
    }
  }
}

__global__ void collocate_add_cube_to_grid(
    const int zeroing_grid, const int list_length, /* size of the batch */
    const int number_of_replicas,                  /* number of grids */
    const int3 grid_size_, const int3 grid_lower_corner_pos_,
    const int3 period_, const int3 window_shift, const int3 window_size,
    const int cmax,
    const double3 *__restrict__ rp,         // center of the cubes
    const double *__restrict__ radius_gpu_, // radius of the spheres
    const size_t *const __restrict__ cube_offset_,
    const double *const __restrict__ cube_gpu_,
    double *const __restrict__ grid_gpu_) {
  if ((list_length * number_of_replicas) >= blockDim.z * gridDim.z)
    return;

  int grid_index = blockIdx.z * blockDim.z / grid_size_.z;

  double *const __restrict__ grid =
      grid_gpu_ + grid_index * grid_size_.z * grid_size_.y * grid_size_.x;
  int list_length_local = list_length / number_of_replicas;

  const int lt_min = grid_index * list_length_local;
  const int lt_max = min((grid_index + 1) * list_length_local, list_length);
  const int3 grid_idx =
      make_int3((threadIdx.z + blockIdx.z * blockDim.z) % grid_size_.z,
                threadIdx.y + blockIdx.y * blockDim.y,
                threadIdx.x + blockIdx.x * blockDim.x);

  const int id =
      (grid_idx.z * grid_size_.y + grid_idx.y) * grid_size_.x + grid_idx.x;
  double grid_point = 0.0;

  if (!zeroing_grid)
    grid_point = grid[id];

  int *__restrict__ px = (int *)array;
  int *__restrict__ py = px + cmax;
  int *__restrict__ pz = py + cmax;

  for (int lt = lt_min; lt < lt_max; lt++) {
    for (int i = 0; i < cmax; i += (blockDim.x * blockDim.y * blockDim.z)) {
      px[i] = 0xFFFFFFFF;
      py[i] = 0xFFFFFFFF;
      pz[i] = 0xFFFFFFFF;
    }

    int3 position;

    int3 cube_size, cube_center, lb_cube;

    double3 roffset;
    double disr_radius = 0;
    const double radius = radius_gpu_[blockIdx.x];

    compute_cube_properties(radius, rp + blockIdx.x, &disr_radius, &roffset,
                            &cube_center, &lb_cube, &cube_size);

    return_cube_position(&cube_center, &lb_cube, period_, &position);

    // TODO we should exclude the cube that do not contribute this block right
    // away

    int dont_ignore_cube = 0;
    for (int i = 0; i <= cube_size.z;
         i += (blockDim.x * blockDim.y * blockDim.z)) {
      int res = (i + lb_cube.z + grid_lower_corner_pos_.z) % period_.z;
      res = (res) < 0 ? (res + period_.z) : (res);
      res -= grid_lower_corner_pos_.z;
      pz[i] = res;
      dont_ignore_cube = (res == grid_idx.z) || dont_ignore_cube;
    }

    __syncthreads();

    /* check if the intersection between [z_min..z_max] of the cube
     * intersect [grid^z_min..grid^z_max]
     */

    int sum = 0;
#if __CUDA_ARCH__ < 700
    sum = __any_sync(0xffffffff, dont_ignore_cube);
#else
    sum = __match_any_sync(0xffffffff, dont_ignore_cube);
#endif

    // all threads should have the same value. If equal to zero then take
    // the next cube on the list
    if (sum == 0) {
      continue;
    }

    dont_ignore_cube = 0;
    for (int i = 0; i <= cube_size.y;
         i += (blockDim.x * blockDim.y * blockDim.z)) {
      int res = (i + lb_cube.y + grid_lower_corner_pos_.y) % period_.y;
      res = (res) < 0 ? (res + period_.y) : (res);
      res -= grid_lower_corner_pos_.y;
      py[i] = res;
      dont_ignore_cube = (res == grid_idx.y) || dont_ignore_cube;
    }

#if __CUDA_ARCH__ < 700
    sum = __any_sync(0xffffffff, dont_ignore_cube);
#else
    sum = __match_any_sync(0xffffffff, dont_ignore_cube);
#endif

    if (sum == 0) {
      continue;
    }

    for (int i = 0; i <= cube_size.x;
         i += (blockDim.x * blockDim.y * blockDim.z)) {
      int res = (i + lb_cube.x + grid_lower_corner_pos_.x) % period_.x;
      res = (res) < 0 ? (res + period_.x) : (res);
      res -= grid_lower_corner_pos_.x;
      px[i] = res;
      dont_ignore_cube = (res == grid_idx.x) || dont_ignore_cube;
    }

#if __CUDA_ARCH__ < 700
    sum = __any_sync(0xffffffff, dont_ignore_cube);
#else
    sum = __match_any_sync(0xffffffff, dont_ignore_cube);
#endif

    if (sum == 0) {
      continue;
    }

    const double *const __restrict__ cube = cube_gpu_ + cube_offset_[lt];

    for (int k = 0; k <= cube_size.z; k++) {
      const int z = pz[k + threadIdx.z];

      /* make the code very divergent */
      if (z != grid_idx.z)
        continue;

      for (int j = 0; j <= cube_size.y; j++) {
        const int y = py[j + threadIdx.y];

        /* make the code very divergent */
        if (y != grid_idx.y)
          continue;

        for (int i = 0; i < cube_size.x; i++) {
          const int x = px[i + threadIdx.x];

          if (x == grid_idx.x)
            grid_point +=
                cube[((k + threadIdx.z) * cube_size.y + (j + threadIdx.y)) *
                         cube_size.x +
                     i + threadIdx.x];
        }
      }
    }
  }
  grid[id] = grid_point;
}

__global__ void
compute_collocation_cube_gpu_(const int apply_cutoff, const int3 period_,
                              const int *const __restrict__ lmax_gpu_,
                              const double *const __restrict__ zeta_gpu,
                              const double3 *const __restrict__ rp,
                              const double *const __restrict__ radius_gpu_,
                              const int *const __restrict__ coef_offset_gpu_,
                              const double *const __restrict__ coef_gpu_,
                              const size_t *const __restrict__ cube_offset_,
                              double *const __restrict__ cube_gpu_) {
  /* the period is sotred in constant memory */
  /* the displacement vectors as well */

  int lmax = lmax_gpu_[blockIdx.x];

  int3 position;

  int3 cube_size, cube_center, lb_cube;

  double3 roffset;
  double disr_radius = 0;
  const double radius = radius_gpu_[blockIdx.x];

  compute_cube_properties(radius, rp + blockIdx.x, &disr_radius, &roffset,
                          &cube_center, &lb_cube, &cube_size);

  return_cube_position(&cube_center, &lb_cube, period_, &position);

  const double *__restrict__ coef = coef_gpu_ + coef_offset_gpu_[blockIdx.x];

  const double zeta = zeta_gpu[blockIdx.x];

  double *coefs_ = (double *)array;

  int id = (threadIdx.z * blockDim.y + threadIdx.y) * blockDim.x + threadIdx.x;

  for (int i = id; i < ((lmax + 1) * (lmax + 1) * (lmax + 1));
       i += (blockDim.x * blockDim.y * blockDim.z))
    coefs_[i] = coef[i];
  __syncthreads();

  double *__restrict__ const cube_ = cube_gpu_ + cube_offset_[blockIdx.x];

  for (int z = threadIdx.z; z <= cube_size.z; z += blockDim.z) {
    const double z1 = z - (cube_size.z / 2) - roffset.z;

    for (int y = threadIdx.y; y <= cube_size.y; y += blockDim.y) {
      double y1 = y - (cube_size.y / 2) - roffset.y;

      for (int x = threadIdx.x; x <= cube_size.x; x += blockDim.x) {
        const double x1 = x - (cube_size.x / 2) - roffset.x;

        /* compute the coordinates of the point in atomic coordinates */
        double3 r3;
        r3.x = z1 * dh_[6] + y1 * dh_[3] + x1 * dh_[0];
        r3.y = z1 * dh_[7] + y1 * dh_[4] + x1 * dh_[1];
        r3.z = z1 * dh_[8] + y1 * dh_[5] + x1 * dh_[2];

        if ((apply_cutoff) &&
            ((radius * radius) < (r3.x * r3.x + r3.y * r3.y + r3.z * r3.z))) {
          continue;
        }

        double exp_factor =
            exp(-(r3.x * r3.x + r3.y * r3.y + r3.z * r3.z) * zeta);

        double res = 0.0;
        double dx = 1;

        /* NOTE: the coefficients are stored as lx,lz,ly */
        /* It is suboptimal right now because i do more operations than needed
         * (a lot of coefs_ are zero). Moreover, it is a dgemm underneath and
         * could be treated with tensor cores */

        for (int alpha = 0; alpha <= lmax; alpha++) {
          double dz = 1;
          for (int gamma = 0; gamma <= (lmax - alpha); gamma++) {
            double dy = dx * dz;
            const int off = (alpha * (lmax + 1) + gamma) * (lmax + 1);
            for (int beta = 0; beta <= (lmax - alpha - gamma); beta++) {
              res += coefs_[off + beta] * dy;
              dy *= r3.y;
            }
            dz *= r3.z;
          }
          dx *= r3.x;
        }

        res *= exp_factor;

        cube_[(z * cube_size.z + y) * cube_size.y + x] += res;
      }
    }
  }
}

__global__ void compute_integration_gpu_(
    const int3 grid_size, const int3 grid_lower_corner_pos, const int3 period,
    const int *__restrict__ lmax_gpu_, const double *__restrict__ zeta_gpu,
    const int *cube_size_gpu_, const int *cube_position_,
    const double *__restrict__ roffset_gpu_,
    const int *__restrict__ coef_offset_gpu_,
    const double *__restrict__ grid_gpu_, double *__restrict__ coef_gpu_) {
  /* the period is stored in constant memory */
  /* the displacement vectors as well */

  int lmax = lmax_gpu_[blockIdx.x];

  const int position[3] = {cube_position_[3 * blockIdx.x],
                           cube_position_[3 * blockIdx.x + 1],
                           cube_position_[3 * blockIdx.x + 2]};

  const int cube_size[3] = {cube_size_gpu_[3 * blockIdx.x],
                            cube_size_gpu_[3 * blockIdx.x + 1],
                            cube_size_gpu_[3 * blockIdx.x + 2]};

  const double roffset[3] = {roffset_gpu_[3 * blockIdx.x],
                             roffset_gpu_[3 * blockIdx.x + 1],
                             roffset_gpu_[3 * blockIdx.x + 2]};

  const double zeta = zeta_gpu[blockIdx.x];

  double *coefs_shared = (double *)array;

  int id = (threadIdx.z * blockDim.y + threadIdx.y) * blockDim.x + threadIdx.x;

  for (int i = id; i < ((lmax + 1) * (lmax + 1) * (lmax + 1));
       i += (blockDim.x * blockDim.y * blockDim.z))
    coefs_shared[i] = 0.0;
  __syncthreads();

  for (int z = threadIdx.z; z < cube_size[0]; z += blockDim.z) {
    const double z1 = z - (cube_size[0] / 2) - roffset[0];
    const int z2 = (z + position[0] + 32 * period.z) % period.z;
    for (int y = threadIdx.y; y < cube_size[1]; y += blockDim.y) {
      double y1 = y - (cube_size[1] / 2) - roffset[1];
      const int y2 = (y + position[1] + 32 * period.y) % period.y;
      for (int x = threadIdx.x; x < cube_size[2]; x += blockDim.x) {
        const double x1 = x - (cube_size[2] / 2) - roffset[2];
        const int x2 = (x + position[2] + 32 * period.x) % period.x;
        double3 r3;

        /* compute the coordinates of the point in atomic coordinates */
        r3.x = z1 * dh_[6] + y1 * dh_[3] + x1 * dh_[0];
        r3.y = z1 * dh_[7] + y1 * dh_[4] + x1 * dh_[1];
        r3.z = z1 * dh_[8] + y1 * dh_[5] + x1 * dh_[2];
        double exp_factor =
            exp(-(r3.x * r3.x + r3.y * r3.y + r3.z * r3.z) * zeta);
        double dx = 1;
        const double grid_value =
            grid_gpu_[(z2 * grid_size.y + y2) * grid_size.x + x2] * exp_factor;
        /* NOTE: the coefficients are stored as lx,lz,ly */

        for (int alpha = 0; alpha <= lmax; alpha++) {
          double dz = 1;
          for (int gamma = 0; gamma <= lmax; gamma++) {
            double dy = dx * dz;
            const int off = (alpha * (lmax + 1) + gamma) * (lmax + 1);
            for (int beta = 0; beta <= lmax; beta++) {
#if __CUDA_ARCH__ < 600
              atomicAdd1(&coefs_shared[off + beta], grid_value * dy);
#else
              atomicAdd(&coefs_shared[off + beta], grid_value * dy);
#endif
              dy *= r3.y;
            }
            dz *= r3.z;
          }
          dx *= r3.x;
        }
      }
    }
  }

  double *__restrict__ coef = coef_gpu_ + coef_offset_gpu_[blockIdx.x];
  for (int i = id; i < ((lmax + 1) * (lmax + 1) * (lmax + 1));
       i += (blockDim.x * blockDim.y * blockDim.z))
    coef[i] = coefs_shared[i];
}

extern "C" void compute_collocation_gpu(pgf_list_gpu *handler) {
  cudaSetDevice(handler->device_id);
  cudaStreamSynchronize(handler->stream);

  if (handler->durty) {
    cudaFree(handler->coef_gpu_);
    cudaMalloc(&handler->coef_gpu_,
               sizeof(double) * handler->coef_alloc_size_gpu_);
    handler->durty = false;
  }

  cudaMemcpyAsync(handler->rp_gpu_, handler->rp_cpu_,
                  sizeof(double3) * handler->list_length,
                  cudaMemcpyHostToDevice, handler->stream);
  cudaMemcpyAsync(handler->radius_gpu_, handler->radius_cpu_,
                  sizeof(double) * handler->list_length, cudaMemcpyHostToDevice,
                  handler->stream);
  cudaMemcpyAsync(handler->zeta_gpu_, handler->zeta_cpu_,
                  sizeof(double) * handler->list_length, cudaMemcpyHostToDevice,
                  handler->stream);
  cudaMemcpyAsync(handler->lmax_gpu_, handler->lmax_cpu_,
                  sizeof(int) * handler->list_length, cudaMemcpyHostToDevice,
                  handler->stream);
  cudaMemcpyAsync(handler->coef_offset_gpu_, handler->coef_offset_cpu_,
                  sizeof(int) * handler->list_length, cudaMemcpyHostToDevice,
                  handler->stream);

  cudaMemcpyAsync(handler->border_mask_gpu_, handler->border_mask_cpu_,
                  sizeof(int) * handler->list_length, cudaMemcpyHostToDevice,
                  handler->stream);

  cudaMemcpyAsync(handler->coef_gpu_, handler->coef_cpu_,
                  sizeof(double) * handler->coef_dynamic_alloc_size_gpu_,
                  cudaMemcpyHostToDevice, handler->stream);

  // cudaEventRecord(handler->event, handler->stream);
  dim3 gridSize, blockSize;

  gridSize.x = handler->list_length;

  blockSize.x = 4;
  blockSize.y = 4;
  blockSize.z = 4;

  compute_collocation_gpu_<<<gridSize, blockSize,
                             (handler->lmax + 1) * (handler->lmax + 1) *
                                 (handler->lmax + 1) * sizeof(double),
                             handler->stream>>>(
      handler->apply_cutoff, handler->grid_size,
      handler->grid_lower_corner_position, handler->grid_full_size,
      handler->border_mask_gpu_, handler->border_width, handler->lmax_gpu_,
      handler->zeta_gpu_, handler->rp_gpu_, handler->radius_gpu_,
      handler->coef_offset_gpu_, handler->coef_gpu_, handler->data_gpu_);
}

extern "C" void compute_collocation_gpu_variant(pgf_list_gpu *handler) {
  cudaSetDevice(handler->device_id);
  cudaStreamSynchronize(handler->stream);

  if (handler->durty) {
    cudaFree(handler->coef_gpu_);
    cudaMalloc(&handler->coef_gpu_,
               sizeof(double) * handler->coef_alloc_size_gpu_);
    handler->durty = false;
  }

  if (handler->cube_data_gpu_ == NULL) {
    cudaMalloc(&handler->cube_data_gpu_,
               sizeof(double) * handler->cube_alloc_size_);
  } else {
    if (handler->max_cube_alloc_size_ < handler->cube_alloc_size_) {
      cudaFree(handler->cube_data_gpu_);
      cudaMalloc(&handler->cube_data_gpu_,
                 sizeof(double) * handler->cube_alloc_size_);
      handler->max_cube_alloc_size_ = handler->cube_alloc_size_;
    }
  }

  cudaMemcpyAsync(handler->rp_gpu_, handler->rp_cpu_,
                  sizeof(double3) * handler->list_length,
                  cudaMemcpyHostToDevice, handler->stream);

  cudaMemcpyAsync(handler->radius_gpu_, handler->radius_cpu_,
                  sizeof(double) * handler->list_length, cudaMemcpyHostToDevice,
                  handler->stream);

  cudaMemcpyAsync(handler->zeta_gpu_, handler->zeta_cpu_,
                  sizeof(double) * handler->list_length, cudaMemcpyHostToDevice,
                  handler->stream);

  cudaMemcpyAsync(handler->lmax_gpu_, handler->lmax_cpu_,
                  sizeof(int) * handler->list_length, cudaMemcpyHostToDevice,
                  handler->stream);

  cudaMemcpyAsync(handler->coef_offset_gpu_, handler->coef_offset_cpu_,
                  sizeof(int) * handler->list_length, cudaMemcpyHostToDevice,
                  handler->stream);

  cudaMemcpyAsync(handler->coef_gpu_, handler->coef_cpu_,
                  sizeof(double) * handler->coef_dynamic_alloc_size_gpu_,
                  cudaMemcpyHostToDevice, handler->stream);

  cudaMemcpyAsync(handler->cube_offset_gpu_, handler->cube_offset_cpu_,
                  sizeof(size_t) * handler->list_length, cudaMemcpyHostToDevice,
                  handler->stream);
  dim3 gridSize, blockSize;

  blockSize.x = 4;
  blockSize.y = 4;
  blockSize.z = 4;

  gridSize.x = handler->list_length;

  compute_collocation_cube_gpu_<<<gridSize, blockSize,
                                  (handler->lmax + 1) * (handler->lmax + 1) *
                                      (handler->lmax + 1) * sizeof(double),
                                  handler->stream>>>(
      handler->apply_cutoff, handler->grid_full_size, handler->lmax_gpu_,
      handler->zeta_gpu_, handler->rp_gpu_, handler->radius_gpu_,
      handler->coef_offset_gpu_, handler->coef_gpu_, handler->cube_offset_gpu_,
      handler->cube_data_gpu_);

  gridSize.x = handler->grid_size.x;
  gridSize.y = handler->grid_size.y;
  gridSize.z = handler->grid_size.z * handler->number_of_replicas;

  /* block size can be rearranged to something maybe better. right now 4x4x4  */
  int shared_mem_size = 3 * handler->cmax * sizeof(int);
  collocate_add_cube_to_grid<<<gridSize, blockSize, shared_mem_size,
                               handler->stream>>>(
      handler->zeroing_grid, handler->list_length, /* size of the batch */
      handler->number_of_replicas,                 /* number of grids */
      handler->grid_size, handler->grid_lower_corner_position,
      handler->grid_full_size, handler->window_shift, handler->window_size,
      handler->cmax,
      handler->rp_gpu_,     // center of the cubes
      handler->radius_gpu_, // radius of the spheres
      handler->cube_offset_gpu_, handler->cube_data_gpu_, handler->data_gpu_);

  handler->zeroing_grid = false;
}

extern "C" void initialize_grid_parameters_on_gpu_step1(void *const ctx,
                                                        const int level) {
  double dh[9], dh_inv[9];

  return_dh(ctx, level, dh);
  return_dh_inv(ctx, level, dh_inv);

  int orthorhombic = is_grid_orthorhombic(ctx);
  for (int device = 0; device < return_num_devs(ctx); device++) {
    cudaSetDevice(return_device_id(ctx, device));
    cudaMemcpyToSymbol(dh_, dh, sizeof(double) * 9);
    cudaMemcpyToSymbol(dh_inv_, dh_inv, sizeof(double) * 9);
    cudaMemcpyToSymbol(is_orthorhombic_, &orthorhombic, sizeof(int));
  }
}

#endif
