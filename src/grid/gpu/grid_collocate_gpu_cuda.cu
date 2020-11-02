#ifdef __COLLOCATE_GPU

#include <assert.h>
#include <cooperative_groups.h>
#include <cuda.h>

#include "../cpu/collocation_integration.h"

extern "C" void reset_list_gpu(pgf_list_gpu *lst);
extern "C" void return_dh(void *const ptr, const int level, double *const dh);
extern "C" void return_dh_inv(void *const ptr, const int level, double *const dh);
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
    ub_cube.x =  - lb_cube->x;
    ub_cube.y =  - lb_cube->y;
    ub_cube.z =  - lb_cube->z;

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

__global__ void compute_collocation_gpu_spherical_cutoff_(
    const int3 grid_size_, const int3 grid_lower_corner_pos_,
    const int3 period_, const int3 window_shift, const int3 window_size,
    const int *__restrict__ lmax_gpu_, const double *__restrict__ zeta_gpu,
    const double3 *__restrict__ rp, const double *__restrict__ radius_gpu_,
    const int *__restrict__ coef_offset_gpu_,
    const double *__restrict__ coef_gpu_, double *__restrict__ grid_gpu_) {
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

  for (int z = threadIdx.z; z <= cube_size.z; z += blockDim.z) {
    const double z1 = z - (cube_size.z / 2) - roffset.z;
    const int z2 = (z + position.z + 32 * period_.z) % period_.z -
                   grid_lower_corner_pos_.z;

    /* check if the point is within the window */
    // if ((z2 < (window_shift.z - grid_lower_corner_pos_.z)) || (z2 >=
    // window_size.z)) {
    //     continue;
    // }

    for (int y = threadIdx.y; y <= cube_size.y; y += blockDim.y) {
      double y1 = y - (cube_size.y / 2) - roffset.y;
      const int y2 = (y + position.y + 32 * period_.y) % period_.y -
                     grid_lower_corner_pos_.y;

      /* check if the point is within the window */
      // if ((y2 < window_shift.y - grid_lower_corner_pos_.y) || (y2 >=
      // window_size.y)) {
      //     continue;
      // }

      for (int x = threadIdx.x; x <= cube_size.x; x += blockDim.x) {
        const double x1 = x - (cube_size.x / 2) - roffset.x;
        const int x2 = (x + position.x + 32 * period_.x) % period_.x -
                       grid_lower_corner_pos_.x;

        /* check if the point is within the window */
        // if ((x2 < (window_shift.x - grid_lower_corner_pos_.x)) || (x2 >=
        // window_size.x)) {
        //      continue;
        // }

        double3 r3;
        r3.x = z1 * dh_[6] + y1 * dh_[3] + x1 * dh_[0];
        r3.y = z1 * dh_[7] + y1 * dh_[4] + x1 * dh_[1];
        r3.z = z1 * dh_[8] + y1 * dh_[5] + x1 * dh_[2];
        /* compute the coordinates of the point in atomic coordinates */
        double exp_factor =
            exp(-(r3.x * r3.x + r3.y * r3.y + r3.z * r3.z) * zeta);

        double res = 0.0;
        double dx = 1.0;

        /* NOTE: the coefficients are stored as lx,lz,ly */

        for (int alpha = 0; alpha <= lmax; alpha++) {
          double dz = 1;
          for (int gamma = 0; gamma <= lmax; gamma++) {
            double dy = dx * dz;
            const int off = (alpha * (lmax + 1) + gamma) * (lmax + 1);
            for (int beta = 0; beta <= lmax; beta++) {
              res += coefs_[off + beta] * dy;
              dy *= r3.y;
            }
            dz *= r3.z;
          }
          dx *= r3.x;
        }

        if (radius * radius >= (r3.x * r3.x + r3.y * r3.y + r3.z * r3.z)) {
            res *= exp_factor;

            /* this operationm is the limiting factor of this code. AtomicAdd
             * act at the thread level while actually we do not need to act at
             * that level but at the grid level. Update of each block should be
             * sequential but not within the block */
#if __CUDA_ARCH__ < 600
            /* kepler does not have hardware atomic operations on double */
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
}

__global__ void compute_collocation_gpu_(
    const int3 grid_size_, const int3 grid_lower_corner_pos_,
    const int3 period_, const int3 window_shift, const int3 window_size,
    const int *__restrict__ lmax_gpu_, const double *__restrict__ zeta_gpu,
    const double3 *__restrict__ rp, const double *__restrict__ radius_gpu_,
    const int *__restrict__ coef_offset_gpu_,
    const double *__restrict__ coef_gpu_, double *__restrict__ grid_gpu_) {
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

  for (int z = threadIdx.z; z <= cube_size.z; z += blockDim.z) {
    const double z1 = z - (cube_size.z / 2) - roffset.z;
    const int z2 = (z + position.z + 32 * period_.z) % period_.z -
                   grid_lower_corner_pos_.z;

    /* check if the point is within the window */
    // if ((z2 < (window_shift.z - grid_lower_corner_pos_.z)) || (z2 >=
    // window_size.z)) {
    //     continue;
    // }

    for (int y = threadIdx.y; y <= cube_size.y; y += blockDim.y) {
      double y1 = y - (cube_size.y / 2) - roffset.y;
      const int y2 = (y + position.y + 32 * period_.y) % period_.y -
                     grid_lower_corner_pos_.y;

      /* check if the point is within the window */
      // if ((y2 < window_shift.y - grid_lower_corner_pos_.y) || (y2 >=
      // window_size.y)) {
      //     continue;
      // }

      for (int x = threadIdx.x; x <= cube_size.x; x += blockDim.x) {
        const double x1 = x - (cube_size.x / 2) - roffset.x;
        const int x2 = (x + position.x + 32 * period_.x) % period_.x -
                       grid_lower_corner_pos_.x;

        /* check if the point is within the window */
        // if ((x2 < (window_shift.x - grid_lower_corner_pos_.x)) || (x2 >=
        // window_size.x)) {
        //      continue;
        // }

        double3 r3;
        r3.x = z1 * dh_[6] + y1 * dh_[3] + x1 * dh_[0];
        r3.y = z1 * dh_[7] + y1 * dh_[4] + x1 * dh_[1];
        r3.z = z1 * dh_[8] + y1 * dh_[5] + x1 * dh_[2];
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
          for (int gamma = 0; gamma <= lmax; gamma++) {
            double dy = dx * dz;
            const int off = (alpha * (lmax + 1) + gamma) * (lmax + 1);
            for (int beta = 0; beta <= lmax; beta++) {
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

__global__ void
compute_integration_gpu_(const int3 grid_size, const int3 grid_lower_corner_pos,
                         const int3 period, const int *__restrict__ lmax_gpu_,
                         const double *__restrict__ zeta_gpu,
                         const int *cube_size_gpu_, const int *cube_position_,
                         const double *__restrict__ roffset_gpu_,
                         const int *__restrict__ coef_offset_gpu_,
                         const double *__restrict__ grid_gpu_,
                         const double *__restrict__ coef_gpu_) {
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
  const double *__restrict__ coef = coef_gpu_ + coef_offset_gpu_[blockIdx.x];

  const double zeta = zeta_gpu[blockIdx.x];

  double *coefs_ = (double *)array;

  int id = (threadIdx.z * blockDim.y + threadIdx.y) * blockDim.x + threadIdx.x;

  for (int i = id; i < ((lmax + 1) * (lmax + 1) * (lmax + 1));
       i += (blockDim.x * blockDim.y * blockDim.z))
    coefs_[i] = 0.0;
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
        double res = 0.0;
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
              coefs_[off + beta] += grid_value * dy;
              dy *= r3.y;
            }
            dz *= r3.z;
          }
          dx *= r3.x;
        }
      }
    }
  }
}

extern "C" void compute_collocation_gpu(pgf_list_gpu *handler) {
  cudaSetDevice(handler->device_id);
  cudaStreamSynchronize(handler->stream);
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

  if (handler->durty) {
    cudaFree(handler->coef_gpu_);
    cudaMalloc(&handler->coef_gpu_,
               sizeof(double) * handler->coef_alloc_size_gpu_);
    handler->durty = false;
  }

  cudaMemcpyAsync(handler->coef_gpu_, handler->coef_cpu_,
                  sizeof(double) * handler->coef_dynamic_alloc_size_gpu_,
                  cudaMemcpyHostToDevice, handler->stream);

  // cudaEventRecord(handler->event, handler->stream);
  dim3 block, thread;

  block.x = handler->list_length;

  thread.x = 4;
  thread.y = 4;
  thread.z = 4;

  if (handler->apply_cutoff) {
      compute_collocation_gpu_spherical_cutoff_<<<block, thread,
          (handler->lmax + 1) * (handler->lmax + 1) *
          (handler->lmax + 1) * sizeof(double),
          handler->stream>>>(
              handler->grid_size, handler->grid_lower_corner_position,
              handler->grid_full_size, handler->window_shift, handler->window_size,
              handler->lmax_gpu_, handler->zeta_gpu_, handler->rp_gpu_,
              handler->radius_gpu_, handler->coef_offset_gpu_, handler->coef_gpu_,
              handler->data_gpu_);
  }  else {
      compute_collocation_gpu_<<<block, thread,
        (handler->lmax + 1) * (handler->lmax + 1) *
        (handler->lmax + 1) * sizeof(double),
        handler->stream>>>(
            handler->grid_size, handler->grid_lower_corner_position,
            handler->grid_full_size, handler->window_shift, handler->window_size,
            handler->lmax_gpu_, handler->zeta_gpu_, handler->rp_gpu_,
            handler->radius_gpu_, handler->coef_offset_gpu_, handler->coef_gpu_,
            handler->data_gpu_);
}
}

extern "C" void
initialize_grid_parameters_on_gpu(collocation_integration *const handler) {
    for (int worker = 0; worker < handler->worker_list_size; worker++) {
      assert(handler->worker_list[worker].device_id >= 0);
      cudaSetDevice(handler->worker_list[worker].device_id);

      handler->worker_list[worker].grid_size.x = handler->grid.size[2];
      handler->worker_list[worker].grid_size.y = handler->grid.size[1];
      handler->worker_list[worker].grid_size.z = handler->grid.size[0];

      handler->worker_list[worker].grid_full_size.x = handler->grid.full_size[2];
      handler->worker_list[worker].grid_full_size.y = handler->grid.full_size[1];
      handler->worker_list[worker].grid_full_size.z = handler->grid.full_size[0];

      handler->worker_list[worker].window_size.x = handler->grid.window_size[2];
      handler->worker_list[worker].window_size.y = handler->grid.window_size[1];
      handler->worker_list[worker].window_size.z = handler->grid.window_size[0];

      handler->worker_list[worker].window_shift.x = handler->grid.window_shift[2];
      handler->worker_list[worker].window_shift.y = handler->grid.window_shift[1];
      handler->worker_list[worker].window_shift.z = handler->grid.window_shift[0];

      handler->worker_list[worker].grid_lower_corner_position.x =
          handler->grid.lower_corner[2];
      handler->worker_list[worker].grid_lower_corner_position.y =
          handler->grid.lower_corner[1];
      handler->worker_list[worker].grid_lower_corner_position.z =
          handler->grid.lower_corner[0];

      if (handler->worker_list[worker].data_gpu_ == NULL) {
          cudaMalloc(&handler->worker_list[worker].data_gpu_,
                     sizeof(double) * handler->grid.alloc_size_);
          handler->worker_list[worker].data_gpu_old_size_ = handler->grid.alloc_size_;
      } else {
          if (handler->worker_list[worker].data_gpu_old_size_ < handler->grid.alloc_size_) {
              cudaFree(handler->worker_list[worker].data_gpu_);
              cudaMalloc(&handler->worker_list[worker].data_gpu_,
                         sizeof(double) * handler->grid.alloc_size_);
              handler->worker_list[worker].data_gpu_old_size_ = handler->grid.alloc_size_;
          }
      }

      cudaMemset(handler->worker_list[worker].data_gpu_, 0,
                 sizeof(double) * handler->grid.alloc_size_);
      reset_list_gpu(handler->worker_list + worker);
  }
}

extern "C" void initialize_grid_parameters_on_gpu_step1(void *const  ctx, const int level) {
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
