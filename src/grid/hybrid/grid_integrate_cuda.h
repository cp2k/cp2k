/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

__device__ void
grid_compute_pab_gpu(cg::thread_block &block, const int *lmin, const int *lmax,
                     const int lp, const double prefactor,
                     const tensor4 &alpha, // [3][lb_max+1][la_max+1][lp+1]
                     const tensor3 &coef_xyz,
                     Matrix &pab) //[lp+1][lp+1][lp+1]
{
  for (int lzb = 0; lzb <= lmax[1]; lzb++) {
    for (int lyb = 0; lyb <= lmax[1] - lzb; lyb++) {
      const int lxb_min = imax(lmin[1] - lzb - lyb, 0);
      for (int lxb = lxb_min; lxb <= lmax[1] - lzb - lyb; lxb++) {
        const int jco = coset(lxb, lyb, lzb);
        for (int lza = 0; lza <= lmax[0]; lza++) {
          for (int lya = 0; lya <= lmax[0] - lza; lya++) {
            const int lxa_min = imax(lmin[0] - lza - lya, 0);
            for (int lxa = lxa_min; lxa <= lmax[1] - lza - lya; lxa++) {
              const int ico = coset(lxa, lya, lza);
              double pab_ = 0.0;
              for (int lyp = 0; lyp <= lya + lyb; lyp++) {
                double p = prefactor * alpha(1, lyb, lya, lyp);
                for (int lxp = 0; lxp <= lp - lya - lyb; lxp++) {
                  double p1 = p * alpha(0, lxb, lxa, lxp);
                  const double *__restrict__ const src1 =
                      coef_xyz.at(lyp, lxp, 0);
                  const double *__restrict__ const src2 =
                      alpha.at(0, lzb, lza, 0);
                  for (int lzp = block.thread_rank();
                       lzp <= lp - lya - lyb - lxp; lyp += block.size()) {
                    pab_ += src1[lyp] * p1 * src2[lyp]; // collocate
                  }
                }
              }

              // need reduction over the group thread

              if (block.thread_rank() == 0) {
                pab(jco, ico) += pab_;
              }
            }
          }
        }
      }
    }
  }
}

__global__ void compute_integration_gpu_(
    const int apply_cutoff, const int3 grid_size_,
    const int3 grid_lower_corner_pos_, const int3 period_,
    const int *__restrict__ border_mask, const int3 border_width,
    const int *__restrict__ lmax_gpu_, const double *__restrict__ zeta_gpu,
    const double3 *__restrict__ rp, const double *__restrict__ radius_gpu_,
    const int *__restrict__ coef_offset_gpu_,
    const double *__restrict__ grid_gpu_, double *__restrict__ coef_gpu_) {
  /* the period is sotred in constant memory */
  /* the displacement vectors as well */

  int lmax = lmax_gpu_[blockIdx.x];

  int3 cube_size, cube_center, lb_cube, window_size, window_shift;

  double3 roffset;
  const double radius = radius_gpu_[blockIdx.x];

  compute_cube_properties(radius, rp + blockIdx.x, &roffset, &cube_center,
                          &lb_cube, &cube_size);

  compute_window_size(&grid_size_, &grid_lower_corner_pos_,
                      &period_, /* also full size of the grid */
                      border_mask[blockIdx.x], &border_width, &window_size,
                      &window_shift);

  const double zeta = zeta_gpu[blockIdx.x];

  double *coefs_ = (double *)array;

  int id = (threadIdx.z * blockDim.y + threadIdx.y) * blockDim.x + threadIdx.x;

  for (int i = id; i < ((lmax + 1) * (lmax + 2) * (lmax + 3)) / 6;
       i += (blockDim.x * blockDim.y * blockDim.z))
    coefs_[i] = 0;
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
            grid_gpu_[(z2 * grid_size_.y + y2) * grid_size_.x + x2] *
            exp(-(r3.x * r3.x + r3.y * r3.y + r3.z * r3.z) * zeta);

        double dx = 1;

        /* NOTE: the coefficients are stored as lx,lz,ly */
        /* It is suboptimal right now because i do more operations than needed
         * (a lot of coefs_ are zero). Moreover, it is a dgemm underneath and
         * could be treated with tensor cores */

        int off = 0;
        for (int alpha = 0; alpha <= lmax; alpha++) {
          double dz = 1;
          for (int gamma = 0; gamma <= (lmax - alpha); gamma++) {
            double dy = dx * dz;
            for (int beta = 0; beta <= (lmax - alpha - gamma); beta++) {
#if __CUDA_ARCH__ < 600
              atomicAdd1(&coefs_[off], exp_factor * dy);
#else
              atomicAdd_block(&coefs_[off], exp_factor * dy);
#endif
              dy *= r3.y;
              off++;
            }
            dz *= r3.z;
          }
          dx *= r3.x;
        }
      }
    }
  }
  __syncthreads();
  double *__restrict__ coef = coef_gpu_ + coef_offset_gpu_[blockIdx.x];

  for (int i = id; i < ((lmax + 1) * (lmax + 2) * (lmax + 3)) / 6;
       i += (blockDim.x * blockDim.y * blockDim.z))
    coef[i] = coefs_[i];
  __syncthreads();
}
