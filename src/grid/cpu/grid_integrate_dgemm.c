/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <omp.h>

#ifdef __LIBXSMM
#include <libxsmm.h>
#endif

#include "../common/grid_common.h"
#include "coefficients.h"
#include "collocation_integration.h"
#include "grid_collocate_dgemm.h"
#include "non_orthorombic_corrections.h"
#include "private_header.h"
#include "tensor_local.h"
#include "utils.h"

/* It is a sub-optimal version of the mapping in case of a cubic cutoff. But is
 * very general and does not depend on the orthorombic nature of the grid. for
 * orthorombic cases, it is faster to apply PCB directly on the polynomials. */

void extract_cube(struct collocation_integration_ *handler,
                  const int *lower_boundaries_cube, const int *cube_center) {
  int position[3];
  /* return the cube position in global coordinates */

  return_cube_position(handler->grid.size, handler->grid.window_shift,
                       cube_center, lower_boundaries_cube,
                       handler->grid.full_size, position);

  int z1 = position[0];
  int z_offset = 0;

  int lower_corner[3];
  int upper_corner[3];
  int diff[3];

  memset(handler->cube.data, 0, sizeof(double) * handler->cube.alloc_size_);

  for (int z = 0; (z < handler->cube.size[0]); z++, z1++) {
    /* lower boundary is within the window */
    lower_corner[0] = z1;
    /* now compute the upper corner */
    /* needs to be as large as possible but still within the region of interest
     */
    upper_corner[0] = compute_next_boundaries(
        &z1, z, handler->grid.full_size[0], handler->cube.size[0]);

    /* now check if the intersection between this interval and the window is
     * empty or not */
    if ((upper_corner[0] < handler->grid.window_shift[0]) ||
        (lower_corner[0] >=
         (handler->grid.window_shift[0] + handler->grid.window_size[0]))) {
      update_loop_index(lower_corner[0], upper_corner[0],
                        handler->grid.full_size[0], &z_offset, &z, &z1);
      continue;
    }

    diff[0] = imax(handler->grid.window_shift[0] - lower_corner[0], 0);
    lower_corner[0] = imax(lower_corner[0], handler->grid.window_shift[0]);
    upper_corner[0] = imin(upper_corner[0], handler->grid.window_shift[0] +
                                                handler->grid.window_size[0]);

    /* // We have a full plane. */
    if ((upper_corner[0] - lower_corner[0]) <= handler->grid.window_size[0]) {
      int y1 = position[1];
      int y_offset = 0;
      for (int y = 0; y < handler->cube.size[1]; y1++, y++) {
        /* lower boundary is within the window */
        lower_corner[1] = y1;
        /* now compute the upper corner */
        /* needs to be as large as possible but still within the region of
         * interest */
        upper_corner[1] = compute_next_boundaries(
            &y1, y, handler->grid.full_size[1], handler->cube.size[1]);

        /* now check if the intersection between this interval and the window is
         * empty or not */
        if ((upper_corner[1] < handler->grid.window_shift[1]) ||
            (lower_corner[1] >=
             (handler->grid.window_shift[1] + handler->grid.window_size[1]))) {
          update_loop_index(lower_corner[1], upper_corner[1],
                            handler->grid.full_size[1], &y_offset, &y, &y1);
          continue;
        }
        diff[1] = imax(handler->grid.window_shift[1] - lower_corner[1], 0);
        lower_corner[1] = imax(lower_corner[1], handler->grid.window_shift[1]);
        upper_corner[1] =
            imin(upper_corner[1],
                 handler->grid.window_shift[1] + handler->grid.window_size[1]);

        if (upper_corner[1] - lower_corner[1]) {
          int x1 = position[2];
          int x_offset = 0;
          for (int x = 0; x < handler->cube.size[2]; x1++, x++) {
            /* lower boundary is within the window */
            lower_corner[2] = x1;
            /* now compute the upper corner */
            /* needs to be as large as possible but still within the region of
             * interest */
            upper_corner[2] = compute_next_boundaries(
                &x1, x, handler->grid.full_size[2], handler->cube.size[2]);

            /* now check if the intersection between this interval and the
             * window is empty or not */
            if ((upper_corner[2] < handler->grid.window_shift[2]) ||
                (lower_corner[2] >= (handler->grid.window_shift[2] +
                                     handler->grid.window_size[2]))) {
              update_loop_index(lower_corner[2], upper_corner[2],
                                handler->grid.full_size[2], &x_offset, &x, &x1);
              continue;
            }

            diff[2] = imax(handler->grid.window_shift[2] - lower_corner[2], 0);
            lower_corner[2] =
                imax(lower_corner[2], handler->grid.window_shift[2]);
            upper_corner[2] =
                imin(upper_corner[2], handler->grid.window_shift[2] +
                                          handler->grid.window_size[2]);

            if (upper_corner[2] - lower_corner[2]) {
              const int position1[3] = {z + diff[0], y + diff[1], x + diff[2]};

              /* the function will internally take care of the local vx global
               * grid */
              extract_sub_grid(
                  lower_corner, // lower corner of the portion of cube (in the
                                // full grid)
                  upper_corner, // upper corner of the portion of cube (in the
                                // full grid)
                  position1,    // starting position subblock inside the cube
                  &handler->grid,
                  &handler->cube); // the grid to add data from

              update_loop_index(lower_corner[2], upper_corner[2],
                                handler->grid.full_size[2], &x_offset, &x, &x1);
            }
          }
          update_loop_index(lower_corner[1], upper_corner[1],
                            handler->grid.full_size[1], &y_offset, &y, &y1);
        }
      }
      update_loop_index(lower_corner[0], upper_corner[0],
                        handler->grid.full_size[0], &z_offset, &z, &z1);
    }
  }
}

double integrate_l0_z(double *scratch, const int z0,
                      const struct tensor_ *const p_alpha_beta_reduced_,
                      struct tensor_ *const grid) {
  const double *__restrict const pz =
      &idx3(p_alpha_beta_reduced_[0], 0, 0, z0); /* k indice */
  double *__restrict__ dst = &idx3(grid[0], 0, 0, 0);

  for (int s = 0; s < grid->size[1] * grid->ld_; s++) {
    dst[s] *= pz[0];
  }

  for (int z1 = 1; z1 < grid->size[0]; z1++) {
    double *__restrict__ src = &idx3(grid[0], z1, 0, 0);
    const double pzz = pz[z1];
    cblas_daxpy(grid->size[1] * grid->ld_, pzz, src, 1, dst, 1);
  }

  return cblas_ddot(grid->size[1] * grid->ld_, &idx3(grid[0], 0, 0, 0), 1,
                    scratch, 1);
}

void grid_integrate(collocation_integration *const handler,
                    const bool use_ortho, const double zetp, const double rp[3],
                    const double radius) {
  if (handler == NULL)
    abort();

  const int lp = handler->coef.size[0] - 1;
  // *** position of the gaussian product
  //
  // this is the actual definition of the position on the grid
  // i.e. a point rp(:) gets here grid coordinates
  // MODULO(rp(:)/dr(:),npts(:))+1
  // hence (0.0,0.0,0.0) in real space is rsgrid%lb on the rsgrid ((1,1,1) on
  // grid)

  // cubecenter(:) = FLOOR(MATMUL(dh_inv, rp))
  int cubecenter[3];
  int cube_size[3];
  int lb_cube[3], ub_cube[3];
  double roffset[3];
  double disr_radius;

  /* cube : grid comtaining pointlike product between polynomials
   *
   * pol : grid  containing the polynomials in all three directions
   *
   * pol_folded : grid containing the polynomials after folding for periodic
   * boundaries conditions
   */

  /* seting up the cube parameters */
  int cmax = compute_cube_properties(use_ortho, radius, handler->dh,
                                     handler->dh_inv, rp, &disr_radius, roffset,
                                     cubecenter, lb_cube, ub_cube, cube_size);

  /* initialize the multidimensional array containing the polynomials */
  if (lp != 0) {
    initialize_tensor_3(&handler->pol, 3, cmax, handler->coef.size[0]);
  } else {
    initialize_tensor_3(&handler->pol, 3, handler->coef.size[0], cmax);
  }
  handler->pol_alloc_size = realloc_tensor(&handler->pol);

  /* allocate memory for the polynomial and the cube */

  initialize_tensor_3(&handler->cube, cube_size[0], cube_size[1], cube_size[2]);

  handler->cube_alloc_size = realloc_tensor(&handler->cube);

  initialize_W_and_T(handler, &handler->coef, &handler->cube);

  /* compute the polynomials */

  // WARNING : do not reverse the order in pol otherwise you will have to
  // reverse the order in collocate_dgemm as well.

  /* the tensor contraction is done for a given order so I have to be careful
   * how the tensors X, Y, Z are stored. In collocate, they are stored
   * normally 0 (Z), (1) Y, (2) X in the table pol. but in integrate (which
   * uses the same tensor reduction), I have a special treatment for l = 0. In
   * that case the order *should* be the same than for collocate. For l > 0,
   * we need a different storage order which is (X) 2, (Y) 0, and (Z) 1. */

  int perm[3] = {2, 0, 1};

  if (lp == 0) {
    /* I need to restore the same order than for collocate */
    perm[0] = 0;
    perm[1] = 1;
    perm[2] = 2;
  }

  bool use_ortho_forced = handler->orthogonal[0] && handler->orthogonal[1] &&
                          handler->orthogonal[2];
  if (use_ortho) {
    grid_fill_pol_dgemm((lp != 0), handler->dh[0][0], roffset[2], 0, lb_cube[2],
                        ub_cube[2], lp, cmax, zetp,
                        &idx3(handler->pol, perm[2], 0, 0)); /* i indice */
    grid_fill_pol_dgemm((lp != 0), handler->dh[1][1], roffset[1], 0, lb_cube[1],
                        ub_cube[1], lp, cmax, zetp,
                        &idx3(handler->pol, perm[1], 0, 0)); /* j indice */
    grid_fill_pol_dgemm((lp != 0), handler->dh[2][2], roffset[0], 0, lb_cube[0],
                        ub_cube[0], lp, cmax, zetp,
                        &idx3(handler->pol, perm[0], 0, 0)); /* k indice */
  } else {
    double dx[3];

    dx[2] = handler->dh[0][0] * handler->dh[0][0] +
            handler->dh[0][1] * handler->dh[0][1] +
            handler->dh[0][2] * handler->dh[0][2];

    dx[1] = handler->dh[1][0] * handler->dh[1][0] +
            handler->dh[1][1] * handler->dh[1][1] +
            handler->dh[1][2] * handler->dh[1][2];

    dx[0] = handler->dh[2][0] * handler->dh[2][0] +
            handler->dh[2][1] * handler->dh[2][1] +
            handler->dh[2][2] * handler->dh[2][2];

    grid_fill_pol_dgemm((lp != 0), 1.0, roffset[2], 0, lb_cube[2], ub_cube[2],
                        lp, cmax, zetp * dx[2],
                        &idx3(handler->pol, perm[2], 0, 0)); /* i indice */
    grid_fill_pol_dgemm((lp != 0), 1.0, roffset[1], 0, lb_cube[1], ub_cube[1],
                        lp, cmax, zetp * dx[1],
                        &idx3(handler->pol, perm[1], 0, 0)); /* j indice */
    grid_fill_pol_dgemm((lp != 0), 1.0, roffset[0], 0, lb_cube[0], ub_cube[0],
                        lp, cmax, zetp * dx[0],
                        &idx3(handler->pol, perm[0], 0, 0)); /* k indice */

    /* the three remaining tensors are initialized in the function */
    calculate_non_orthorombic_corrections_tensor(
        zetp, roffset, handler->dh, lb_cube, ub_cube, handler->orthogonal,
        &handler->Exp);
  }

  extract_cube(handler, lb_cube, cubecenter);

  if (!use_ortho && !use_ortho_forced)
    apply_non_orthorombic_corrections(handler->orthogonal, &handler->Exp,
                                      &handler->cube);

  if (lp != 0) {
    tensor_reduction_for_collocate_integrate(handler->scratch,
                                             // pointer to scratch memory
                                             1.0, handler->orthogonal, NULL,
                                             &handler->cube, &handler->pol,
                                             &handler->coef);
  } else {
    /* it is very specific to integrate because we might end up with a
     * single element after the tensor product/contraction. In that case, I
     * compute the cube and then do a scalar product between the two. */

    /* we could also do this with 2 matrix-vector multiplications and a scalar
     * product
     *
     * H_{jk} = C_{ijk} . P_i (along x) C_{ijk} is *stored C[k][j][i]* !!!!!!
     * L_{k} = H_{jk} . P_j (along y)
     * v_{ab} = L_k . P_k (along z)
     */

    tensor cube_tmp;
    initialize_tensor_2(&cube_tmp, cube_size[0], cube_size[1]);
    alloc_tensor(&cube_tmp);

    /* first along x */
    cblas_dgemv(CblasRowMajor, CblasNoTrans,
                handler->cube.size[0] * handler->cube.size[1],
                handler->cube.size[2], 1.0, handler->cube.data,
                handler->cube.ld_, &idx3(handler->pol, 2, 0, 0), 1, 0.0,
                cube_tmp.data, 1);

    /* second along j */
    cblas_dgemv(CblasRowMajor, CblasNoTrans, handler->cube.size[0],
                handler->cube.size[1], 1.0, cube_tmp.data, cube_tmp.ld_,
                &idx3(handler->pol, 1, 0, 0), 1, 0.0, handler->scratch, 1);

    /* finally along k, it is a scalar product.... */
    handler->coef.data[0] = cblas_ddot(handler->cube.size[0], handler->scratch,
                                       1, &idx3(handler->pol, 0, 0, 0), 1);

    free(cube_tmp.data);
  }

  /* go from ijk -> xyz */
  if (!use_ortho)
    grid_transform_coef_jik_to_yxz(handler->dh, &handler->coef);
}

// *****************************************************************************
void grid_integrate_pgf_product_dgemm(
    void *const handle, const bool use_ortho, const int lp, const double zeta,
    const double zetb, const double dh[3][3], const double dh_inv[3][3],
    const double ra[3], const double rab[3], const int grid_global_size[3],
    const int grid_local_size[3], const int shift_local[3], const double radius,
    double *const grid_, double *const coef) {
  if (!handle) {
    abort();
  }

  collocation_integration *handler = (collocation_integration *)handle;
  const double zetp = zeta + zetb;
  const double f = zetb / zetp;
  double rp[3];

  initialize_basis_vectors(handler, dh, dh_inv);
  verify_orthogonality(dh, handler->orthogonal);

  initialize_tensor_3(&(handler->grid), grid_local_size[2], grid_local_size[1],
                      grid_local_size[0]);

  handler->grid.ld_ = grid_local_size[0];
  handler->grid.data = grid_;
  handler->blocked_grid.blocked_decomposition = false;

  setup_global_grid_size(&handler->grid, (const int *const)grid_global_size);

  initialize_tensor_3(&handler->grid, grid_local_size[2], grid_local_size[1],
                      grid_local_size[0]);
  handler->grid.ld_ = grid_local_size[0];

  setup_grid_window(&handler->grid, shift_local, NULL, 0);

  for (int i = 0; i < 3; i++) {
    rp[i] = ra[i] + f * rab[i];
  }

  initialize_tensor_3(&(handler->coef), lp + 1, lp + 1, lp + 1);
  alloc_tensor(&(handler->coef));

  grid_integrate(handler, use_ortho, zetp, rp, radius);

  /* I need to transpose the coefficients because they are computed as ly, lx,
   * lz while we want them in the format lz, ly, lx. Fortunately it is a
   * single transpose. So either I include it in the next tranformation or I
   * do it separately. */

  transform_yxz_to_triangular(&handler->coef, coef);
  handler->grid.data = NULL;
  /* Return the result to cp2k for now */
}
