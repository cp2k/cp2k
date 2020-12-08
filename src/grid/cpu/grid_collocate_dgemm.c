/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __LIBXSMM
#include <libxsmm.h>
#endif

#include "../common/grid_basis_set.h"
#include "../common/grid_common.h"
#include "../common/grid_constants.h"
#include "coefficients.h"
#include "collocation_integration.h"
#include "grid_collocate_dgemm.h"
#include "grid_cpu_task_list.h"
#include "grid_prepare_pab_dgemm.h"
#include "non_orthorombic_corrections.h"
#include "private_header.h"
#include "tensor_local.h"

void collocate_l0(double *scratch, const double alpha, const bool orthogonal,
                  const struct tensor_ *exp_xy,
                  const struct tensor_ *p_alpha_beta_reduced_,
                  struct tensor_ *cube);

void rotate_to_cartesian_harmonics(const grid_basis_set *ibasis,
                                   const grid_basis_set *jbasis,
                                   const int iatom, const int jatom,
                                   const int iset, const int jset,
                                   double *const block, tensor *work,
                                   tensor *pab) {
  dgemm_params m1, m2;
  memset(&m1, 0, sizeof(dgemm_params));
  memset(&m2, 0, sizeof(dgemm_params));

  // Define some more convenient aliases.
  const int nsgf_seta = ibasis->nsgf_set[iset]; // size of spherical set
  const int nsgf_setb = jbasis->nsgf_set[jset];
  const int nsgfa = ibasis->nsgf; // size of entire spherical basis
  const int nsgfb = jbasis->nsgf;
  const int sgfa = ibasis->first_sgf[iset] - 1; // start of spherical set
  const int sgfb = jbasis->first_sgf[jset] - 1;
  const int maxcoa = ibasis->maxco;
  const int maxcob = jbasis->maxco;
  const int ncoseta = ncoset(ibasis->lmax[iset]);
  const int ncosetb = ncoset(jbasis->lmax[jset]);
  const int ncoa = ibasis->npgf[iset] * ncoseta; // size of carthesian set
  const int ncob = jbasis->npgf[jset] * ncosetb;

  initialize_tensor_2(work, nsgf_setb, ncoa);
  realloc_tensor(work);

  initialize_tensor_2(pab, ncob, ncoa);
  realloc_tensor(pab);

  // the rotations happen here.
  if (iatom <= jatom) {
    m1.op1 = 'N';
    m1.op2 = 'N';
    m1.m = work->size[0];
    m1.n = work->size[1];
    m1.k = nsgf_seta;
    m1.alpha = 1.0;
    m1.beta = 0.0;
    m1.a = block + sgfb * nsgfa + sgfa;
    m1.lda = nsgfa;
    m1.b = &ibasis->sphi[sgfa * maxcoa];
    m1.ldb = maxcoa;
    m1.c = work->data;
    m1.ldc = work->ld_;
  } else {
    m1.op1 = 'T';
    m1.op2 = 'N';
    m1.m = work->size[0];
    m1.n = work->size[1];
    m1.k = nsgf_seta;
    m1.alpha = 1.0;
    m1.beta = 0.0;
    m1.a = block + sgfa * nsgfb + sgfb;
    m1.lda = nsgfb;
    m1.b = &ibasis->sphi[sgfa * maxcoa];
    m1.ldb = maxcoa;
    m1.c = work->data;
    m1.ldc = work->ld_;
  }
  m1.use_libxsmm = true;
  dgemm_simplified(&m1);

  m2.op1 = 'T';
  m2.op2 = 'N';
  m2.m = pab->size[0];
  m2.n = pab->size[1];
  m2.k = work->size[0];
  m2.alpha = 1.0;
  m2.beta = 0.0;
  m2.a = &jbasis->sphi[sgfb * maxcob];
  m2.lda = maxcob;
  m2.b = work->data;
  m2.ldb = work->ld_;
  m2.c = pab->data;
  m2.ldc = pab->ld_;
  m2.use_libxsmm = true;
  dgemm_simplified(&m2);
}

void tensor_reduction_for_collocate_integrate(
    double *scratch, const double alpha, const bool *const orthogonal,
    const struct tensor_ *Exp, const struct tensor_ *co,
    const struct tensor_ *p_alpha_beta_reduced_, struct tensor_ *cube);

/* compute the functions (x - x_i)^l exp (-eta (x - x_i)^2) for l = 0..lp using
 * a recursive relation to avoid computing the exponential on each grid point. I
 * think it is not really necessary anymore since it is *not* the dominating
 * contribution to computation of collocate and integrate */

void grid_fill_pol_dgemm(const bool transpose, const double dr,
                         const double roffset, const int pol_offset,
                         const int xmin, const int xmax, const int lp,
                         const int cmax, const double zetp, double *pol_) {
  tensor pol;
  const double t_exp_1 = exp(-zetp * dr * dr);
  const double t_exp_2 = t_exp_1 * t_exp_1;

  double t_exp_min_1 = exp(-zetp * (dr - roffset) * (dr - roffset));
  double t_exp_min_2 = exp(2.0 * zetp * (dr - roffset) * dr);

  double t_exp_plus_1 = exp(-zetp * roffset * roffset);
  double t_exp_plus_2 = exp(2.0 * zetp * roffset * dr);

  if (transpose) {
    initialize_tensor_2(&pol, cmax, lp + 1);
    pol.data = pol_;
    /* It is original Ole code. I need to transpose the polynomials for the
     * integration routine and Ole code already does it. */
    for (int ig = 0; ig >= xmin; ig--) {
      const double rpg = ig * dr - roffset;
      t_exp_min_1 *= t_exp_min_2 * t_exp_1;
      t_exp_min_2 *= t_exp_2;
      double pg = t_exp_min_1;
      for (int icoef = 0; icoef <= lp; icoef++) {
        idx2(pol, pol_offset + ig - xmin, icoef) = pg;
        pg *= rpg;
      }
    }

    double t_exp_plus_1 = exp(-zetp * roffset * roffset);
    double t_exp_plus_2 = exp(2 * zetp * roffset * dr);
    for (int ig = 1; ig <= xmax; ig++) {
      const double rpg = ig * dr - roffset;
      t_exp_plus_1 *= t_exp_plus_2 * t_exp_1;
      t_exp_plus_2 *= t_exp_2;
      double pg = t_exp_plus_1;
      for (int icoef = 0; icoef <= lp; icoef++) {
        idx2(pol, pol_offset + ig - xmin, icoef) = pg;
        pg *= rpg;
      }
    }

  } else {
    initialize_tensor_2(&pol, lp + 1, cmax);
    pol.data = pol_;
    /* memset(pol.data, 0, sizeof(double) * pol.alloc_size_); */
    /*
     *   compute the values of all (x-xp)**lp*exp(..)
     *
     *  still requires the old trick:
     *  new trick to avoid to many exps (reuse the result from the previous
     * gridpoint): exp( -a*(x+d)**2)=exp(-a*x**2)*exp(-2*a*x*d)*exp(-a*d**2)
     *  exp(-2*a*(x+d)*d)=exp(-2*a*x*d)*exp(-2*a*d**2)
     */

    /* compute the exponential recursively and store the polynomial prefactors
     * as well */
    for (int ig = 0; ig >= xmin; ig--) {
      const double rpg = ig * dr - roffset;
      t_exp_min_1 *= t_exp_min_2 * t_exp_1;
      t_exp_min_2 *= t_exp_2;
      double pg = t_exp_min_1;
      idx2(pol, 0, pol_offset + ig - xmin) = pg;
      if (lp > 0)
        idx2(pol, 1, pol_offset + ig - xmin) = rpg;
    }

    for (int ig = 1; ig <= xmax; ig++) {
      const double rpg = ig * dr - roffset;
      t_exp_plus_1 *= t_exp_plus_2 * t_exp_1;
      t_exp_plus_2 *= t_exp_2;
      double pg = t_exp_plus_1;
      idx2(pol, 0, pol_offset + ig - xmin) = pg;
      if (lp > 0)
        idx2(pol, 1, pol_offset + ig - xmin) = rpg;
    }

    /* compute the remaining powers using previously computed stuff */
    if (lp >= 2) {
      double *__restrict__ poly = &idx2(pol, 1, 0);
      double *__restrict__ src1 = &idx2(pol, 0, 0);
      double *__restrict__ dst = &idx2(pol, 2, 0);
//#pragma omp simd
#pragma GCC ivdep
      for (int ig = 0; ig < (xmax - xmin + 1 + pol_offset); ig++)
        dst[ig] = src1[ig] * poly[ig] * poly[ig];
    }

    for (int icoef = 3; icoef <= lp; icoef++) {
      const double *__restrict__ poly = &idx2(pol, 1, 0);
      const double *__restrict__ src1 = &idx2(pol, icoef - 1, 0);
      double *__restrict__ dst = &idx2(pol, icoef, 0);
//#pragma omp simd
#pragma GCC ivdep
      for (int ig = 0; ig < (xmax - xmin + 1 + pol_offset); ig++) {
        dst[ig] = src1[ig] * poly[ig];
      }
    }

    //
    if (lp > 0) {
      double *__restrict__ dst = &idx2(pol, 1, 0);
      const double *__restrict__ src = &idx2(pol, 0, 0);
#pragma GCC ivdep
      for (int ig = 0; ig < (xmax - xmin + 1 + pol_offset); ig++) {
        dst[ig] *= src[ig];
      }
    }
  }
}

void apply_sphere_cutoff_ortho(struct collocation_integration_ *const handler,
                               const double disr_radius, const int cmax,
                               const int *const lb_cube,
                               const int *const cube_center) {
  // a mapping so that the ig corresponds to the right grid point
  int **map = handler->map;
  map[1] = map[0] + 2 * cmax + 1;
  map[2] = map[1] + 2 * cmax + 1;
  memset(map[0], 0xff, sizeof(int) * 3 * (2 * cmax + 1));

  for (int i = 0; i < 3; i++) {
    for (int ig = 0; ig < handler->cube.size[i]; ig++) {
      map[i][ig] = modulo(cube_center[i] + lb_cube[i] + ig -
                              handler->grid.lower_corner[i],
                          handler->grid.full_size[i]);
    }
  }

  const Interval zwindow = {.xmin = handler->grid.window_shift[0],
                            .xmax = handler->grid.window_size[0]};
  const Interval ywindow = {.xmin = handler->grid.window_shift[1],
                            .xmax = handler->grid.window_size[1]};
  const Interval xwindow = {.xmin = handler->grid.window_shift[2],
                            .xmax = handler->grid.window_size[2]};

  for (int kg = 0; kg < handler->cube.size[0]; kg++) {
    const int k = map[0][kg];
    const int kd =
        (2 * (kg + lb_cube[0]) - 1) / 2; // distance from center in grid points
    const double kr = kd * handler->dh[2][2]; // distance from center in a.u.
    const double kremain = disr_radius * disr_radius - kr * kr;

    if ((kremain >= 0.0) && is_point_in_interval(k, zwindow)) {

      const int jgmin = ceil(-1e-8 - sqrt(kremain) * handler->dh_inv[1][1]);
      for (int jg = jgmin; jg <= (1 - jgmin); jg++) {
        const int j = map[1][jg - lb_cube[1]];
        const double jr = ((2 * jg - 1) >> 1) *
                          handler->dh[1][1]; // distance from center in a.u.
        const double jremain = kremain - jr * jr;

        if ((jremain >= 0.0) && is_point_in_interval(j, ywindow)) {
          const int xmin = ceil(-1e-8 - sqrt(jremain) * handler->dh_inv[0][0]);
          const int xmax = 1 - xmin;

          // printf("xmin %d, xmax %d\n", xmin, xmax);
          for (int x = xmin - lb_cube[2];
               x < imin((xmax - lb_cube[2]), handler->cube.size[2]); x++) {
            const int x1 = map[2][x];

            if (!is_point_in_interval(x1, xwindow))
              continue;

            int lower_corner[3] = {k, j, x1};
            int upper_corner[3] = {k + 1, j + 1, x1 + 1};

            compute_interval(map[2], handler->grid.full_size[2],
                             handler->grid.size[2], handler->cube.size[2], x1,
                             &x, lower_corner + 2, upper_corner + 2, xwindow);

            if (upper_corner[2] - lower_corner[2]) {
              const int position1[3] = {kg, jg - lb_cube[1], x};

              double *restrict dst = &idx3(handler->grid, lower_corner[0],
                                           lower_corner[1], lower_corner[2]);
              double *restrict src = &idx3(handler->cube, position1[0],
                                           position1[1], position1[2]);

              const int sizex = upper_corner[2] - lower_corner[2];
#pragma GCC ivdep
              for (int x = 0; x < sizex; x++) {
                dst[x] += src[x];
              }
            }

            if (handler->grid.size[2] == handler->grid.full_size[2])
              update_loop_index(handler->grid.full_size[2], x1, &x);
            else
              x += upper_corner[2] - lower_corner[2] - 1;
          }
        }
      }
    }
  }
}

void apply_spherical_cutoff_generic(
    struct collocation_integration_ *const handler, const double disr_radius,
    const int cmax, const int *const lb_cube, const int *const ub_cube,
    const double *const roffset, const int *const cube_center) {

  const double a = handler->dh[0][0] * handler->dh[0][0] +
                   handler->dh[0][1] * handler->dh[0][1] +
                   handler->dh[0][2] * handler->dh[0][2];
  const double a_inv = 1.0 / a;
  // a mapping so that the ig corresponds to the right grid point
  int **map = handler->map;
  map[1] = map[0] + 2 * cmax + 1;
  map[2] = map[1] + 2 * cmax + 1;
  memset(map[0], 0xff, sizeof(int) * 3 * (2 * cmax + 1));

  for (int i = 0; i < 3; i++) {
    for (int ig = 0; ig < handler->cube.size[i]; ig++) {
      map[i][ig] = modulo(cube_center[i] + ig + lb_cube[i] -
                              handler->grid.lower_corner[i],
                          handler->grid.full_size[i]);
    }
  }

  const Interval zwindow = {.xmin = handler->grid.window_shift[0],
                            .xmax = handler->grid.window_size[0]};
  const Interval ywindow = {.xmin = handler->grid.window_shift[1],
                            .xmax = handler->grid.window_size[1]};
  const Interval xwindow = {.xmin = handler->grid.window_shift[2],
                            .xmax = handler->grid.window_size[2]};

  for (int k = 0; k < handler->cube.size[0]; k++) {
    const int iz = map[0][k];

    if (!is_point_in_interval(iz, zwindow))
      continue;

    const double tz = (k + lb_cube[0] - roffset[0]);
    const double z[3] = {tz * handler->dh[2][0], tz * handler->dh[2][1],
                         tz * handler->dh[2][2]};

    for (int j = 0; j < handler->cube.size[1]; j++) {
      const int iy = map[1][j];

      if (!is_point_in_interval(iy, ywindow))
        continue;

      const double ty = (j - roffset[1] + lb_cube[1]);
      const double y[3] = {z[0] + ty * handler->dh[1][0],
                           z[1] + ty * handler->dh[1][1],
                           z[2] + ty * handler->dh[1][2]};

      /* Sqrt[(-2 x1 \[Alpha] - 2 y1 \[Beta] - 2 z1 \[Gamma])^2 - */
      /*                                            4 (x1^2 + y1^2 + z1^2)
       * (\[Alpha]^2 + \[Beta]^2 + \[Gamma]^2)] */

      const double b =
          -2.0 * (handler->dh[0][0] * (roffset[2] * handler->dh[0][0] - y[0]) +
                  handler->dh[0][1] * (roffset[2] * handler->dh[0][1] - y[1]) +
                  handler->dh[0][2] * (roffset[2] * handler->dh[0][2] - y[2]));

      const double c = (roffset[2] * handler->dh[0][0] - y[0]) *
                           (roffset[2] * handler->dh[0][0] - y[0]) +
                       (roffset[2] * handler->dh[0][1] - y[1]) *
                           (roffset[2] * handler->dh[0][1] - y[1]) +
                       (roffset[2] * handler->dh[0][2] - y[2]) *
                           (roffset[2] * handler->dh[0][2] - y[2]) -
                       disr_radius * disr_radius;

      double delta = b * b - 4.0 * a * c;

      if (delta < 0.0)
        continue;

      delta = sqrt(delta);

      const int xmin = imax(ceil((-b - delta) * 0.5 * a_inv), lb_cube[2]);
      const int xmax = imin(floor((-b + delta) * 0.5 * a_inv), ub_cube[2]);

      int lower_corner[3] = {iz, iy, xmin};
      int upper_corner[3] = {iz + 1, iy + 1, xmin};

      for (int x = xmin - lb_cube[2];
           x < imin((xmax - lb_cube[2]), handler->cube.size[2]); x++) {
        const int x1 = map[2][x];

        if (!is_point_in_interval(x1, xwindow))
          continue;

        compute_interval(map[2], handler->grid.full_size[2],
                         handler->grid.size[2], handler->cube.size[2], x1, &x,
                         lower_corner + 2, upper_corner + 2, xwindow);

        if (upper_corner[2] - lower_corner[2]) {
          const int position1[3] = {k, j, x};

          /* the function will internally take care of the local vs global
           * grid */

          double *__restrict__ dst = &idx3(handler->grid, lower_corner[0],
                                           lower_corner[1], lower_corner[2]);
          double *__restrict__ src =
              &idx3(handler->cube, position1[0], position1[1], position1[2]);

          const int sizex = upper_corner[2] - lower_corner[2];
#pragma GCC ivdep
          for (int x = 0; x < sizex; x++) {
            dst[x] += src[x];
          }

          if (handler->grid.size[0] == handler->grid.full_size[0])
            update_loop_index(handler->grid.full_size[2], x1, &x);
          else
            x += upper_corner[2] - lower_corner[2] - 1;
        }
      }
    }
  }
}

void collocate_l0(double *scratch, const double alpha, const bool orthogonal_xy,
                  const struct tensor_ *exp_xy,
                  const struct tensor_ *p_alpha_beta_reduced_,
                  struct tensor_ *cube) {
  const double *__restrict pz =
      &idx3(p_alpha_beta_reduced_[0], 0, 0, 0); /* k indice */
  const double *__restrict py =
      &idx3(p_alpha_beta_reduced_[0], 1, 0, 0); /* j indice */
  const double *__restrict px =
      &idx3(p_alpha_beta_reduced_[0], 2, 0, 0); /* i indice */

  memset(&idx3(cube[0], 0, 0, 0), 0, sizeof(double) * cube->alloc_size_);
  memset(scratch, 0, sizeof(double) * cube->size[1] * cube->ld_);

  cblas_dger(CblasRowMajor, cube->size[1], cube->size[2], alpha, py, 1, px, 1,
             scratch, cube->ld_);

  if (exp_xy && !orthogonal_xy) {
    for (int y = 0; y < cube->size[1]; y++) {
      const double *__restrict src = &idx2(exp_xy[0], y, 0);
      double *__restrict dst = &scratch[y * cube->ld_];
#pragma GCC ivdep
      for (int x = 0; x < cube->size[2]; x++) {
        dst[x] *= src[x];
      }
    }
  }

  cblas_dger(CblasRowMajor, cube->size[0], cube->size[1] * cube->ld_, 1.0, pz,
             1, scratch, 1, &idx3(cube[0], 0, 0, 0), cube->size[1] * cube->ld_);
}

/* compute the following operation (variant) it is a tensor contraction

   V_{kji} = \sum_{\alpha\beta\gamma} C_{\alpha\gamma\beta} T_{2,\alpha,i}
   T_{1,\beta,j} T_{0,\gamma,k}

*/
void tensor_reduction_for_collocate_integrate(
    double *scratch, const double alpha, const bool *const orthogonal,
    const struct tensor_ *Exp, const struct tensor_ *co,
    const struct tensor_ *p_alpha_beta_reduced_, struct tensor_ *cube) {
  if (co->size[0] > 1) {
    dgemm_params m1, m2, m3;

    memset(&m1, 0, sizeof(dgemm_params));
    memset(&m2, 0, sizeof(dgemm_params));
    memset(&m3, 0, sizeof(dgemm_params));
    tensor T;
    tensor W;

    double *__restrict const pz =
        &idx3(p_alpha_beta_reduced_[0], 0, 0, 0); /* k indice */
    double *__restrict const py =
        &idx3(p_alpha_beta_reduced_[0], 1, 0, 0); /* j indice */
    double *__restrict const px =
        &idx3(p_alpha_beta_reduced_[0], 2, 0, 0); /* i indice */

    initialize_tensor_3(&T, co->size[0] /* alpha */, co->size[1] /* gamma */,
                        cube->size[1] /* j */);
    initialize_tensor_3(&W, co->size[1] /* gamma */, cube->size[1] /* j */,
                        cube->size[2] /* i */);

    T.data = scratch;
    W.data = scratch + T.alloc_size_;
    /* WARNING we are in row major layout. cblas allows it and it is more
     * natural to read left to right than top to bottom
     *
     * we do first T_{\alpha,\gamma,j} = \sum_beta C_{alpha\gamma\beta}
     * Y_{\beta, j}
     *
     * keep in mind that Y_{\beta, j} = p_alpha_beta_reduced(1, \beta, j)
     * and the order of indices is also important. the last indice is the
     * fastest one. it can be done with one dgemm.
     */

    m1.op1 = 'N';
    m1.op2 = 'N';
    m1.alpha = alpha;
    m1.beta = 0.0;
    m1.m = co->size[0] * co->size[1]; /* alpha gamma */
    m1.n = cube->size[1];             /* j */
    m1.k = co->size[2];               /* beta */
    m1.a = co->data;                  // Coef_{alpha,gamma,beta} Coef_xzy
    m1.lda = co->ld_;
    m1.b = py; // Y_{beta, j} = p_alpha_beta_reduced(1, beta, j)
    m1.ldb = p_alpha_beta_reduced_->ld_;
    m1.c = T.data; // T_{\alpha, \gamma, j} = T(alpha, gamma, j)
    m1.ldc = T.ld_;
    m1.use_libxsmm = true;
    /*
     * the next step is a reduction along the alpha index.
     *
     * We compute then
     *
     * W_{gamma, j, i} = sum_{\alpha} T_{\gamma, j, alpha} X_{\alpha, i}
     *
     * which means we need to transpose T_{\alpha, \gamma, j} to get
     * T_{\gamma, j, \alpha}. Fortunately we can do it while doing the
     * matrix - matrix multiplication
     */

    m2.op1 = 'T';
    m2.op2 = 'N';
    m2.alpha = 1.0;
    m2.beta = 0.0;
    m2.m = cube->size[1] * co->size[1]; // (\gamma j) direction
    m2.n = cube->size[2];               // i
    m2.k = co->size[0];                 // alpha
    m2.a = T.data;                      // T_{\alpha, \gamma, j}
    m2.lda = T.ld_ * co->size[1];
    m2.b = px; // X_{alpha, i}  = p_alpha_beta_reduced(0, alpha, i)
    m2.ldb = p_alpha_beta_reduced_->ld_;
    m2.c = W.data; // W_{\gamma, j, i}
    m2.ldc = W.ld_;
    m2.use_libxsmm = true;
    /* the final step is again a reduction along the gamma indice. It can
     * again be done with one dgemm. The operation is simply
     *
     * Cube_{k, j, i} = \sum_{alpha} Z_{k, \gamma} W_{\gamma, j, i}
     *
     * which means we need to transpose Z_{\gamma, k}.
     */

    m3.op1 = 'T';
    m3.op2 = 'N';
    m3.alpha = alpha;
    m3.beta = 0.0;
    m3.m = cube->size[0];                 // Z_{k \gamma}
    m3.n = cube->size[1] * cube->size[2]; // (ji) direction
    m3.k = co->size[1];                   // \gamma
    m3.a = pz;                            // p_alpha_beta_reduced(0, gamma, i)
    m3.lda = p_alpha_beta_reduced_->ld_;
    m3.b = &idx3(W, 0, 0, 0); // W_{\gamma, j, i}
    m3.ldb = W.size[1] * W.ld_;
    m3.c = &idx3(cube[0], 0, 0, 0); // cube_{kji}
    m3.ldc = cube->ld_ * cube->size[1];
    m3.use_libxsmm = true;
    dgemm_simplified(&m1);
    dgemm_simplified(&m2);

    // apply the non orthorombic corrections in the xy plane
    if (Exp && !orthogonal[2]) {
      tensor exp_xy;
      initialize_tensor_2(&exp_xy, Exp->size[1], Exp->size[2]);
      exp_xy.data = &idx3(Exp[0], 2, 0, 0);
      apply_non_orthorombic_corrections_xy_blocked(&exp_xy, &W);
    }

    dgemm_simplified(&m3);
  } else {
    if (Exp && !orthogonal[2]) {
      tensor exp_xy;
      initialize_tensor_2(&exp_xy, Exp->size[1], Exp->size[2]);

      exp_xy.data = &idx3(Exp[0], 2, 0, 0);
      collocate_l0(scratch, co->data[0] * alpha, orthogonal[2], &exp_xy,
                   p_alpha_beta_reduced_, cube);
    } else {
      collocate_l0(scratch, co->data[0] * alpha, true, NULL,
                   p_alpha_beta_reduced_, cube);
    }
  }

  if (Exp && (!orthogonal[0] || !orthogonal[1])) {
    tensor exp_xy;
    initialize_tensor_2(&exp_xy, Exp->size[1], Exp->size[2]);
    if (!orthogonal[0]) {
      exp_xy.data = &idx3(Exp[0], 0, 0, 0);
      apply_non_orthorombic_corrections_xz_blocked(&exp_xy, cube);
    }

    if (!orthogonal[1]) {
      exp_xy.data = &idx3(Exp[0], 1, 0, 0);
      apply_non_orthorombic_corrections_yz_blocked(&exp_xy, cube);
    }
  }

  return;
}

/* It is a sub-optimal version of the mapping in case of a cubic cutoff. But is
 * very general and does not depend on the orthorombic nature of the grid. for
 * orthorombic cases, it is faster to apply PCB directly on the polynomials. */

void apply_mapping_cubic(struct collocation_integration_ *handler,
                         const int cmax, const int *const lower_boundaries_cube,
                         const int *const cube_center) {

  // a mapping so that the ig corresponds to the right grid point
  int **map = handler->map;
  map[1] = map[0] + 2 * cmax + 1;
  map[2] = map[1] + 2 * cmax + 1;
  memset(map[0], 0xff, sizeof(int) * 3 * (2 * cmax + 1));
  for (int i = 0; i < 3; i++) {
    for (int ig = 0; ig < handler->cube.size[i]; ig++) {
      map[i][ig] = modulo(cube_center[i] + ig + lower_boundaries_cube[i] -
                              handler->grid.lower_corner[i],
                          handler->grid.full_size[i]);
    }
  }

  int lower_corner[3];
  int upper_corner[3];

  const Interval zwindow = {.xmin = handler->grid.window_shift[0],
                            .xmax = handler->grid.window_size[0]};
  const Interval ywindow = {.xmin = handler->grid.window_shift[1],
                            .xmax = handler->grid.window_size[1]};
  const Interval xwindow = {.xmin = handler->grid.window_shift[2],
                            .xmax = handler->grid.window_size[2]};

  /* this code makes a decomposition of the cube such that we can add block of
   * datas in a vectorized way. */

  /* the decomposition depends unfortunately on the way the grid is split over
   * mpi ranks. If each mpi rank has the full grid then it is simple. A 1D
   * example of the decomposition will explain it better. We have an interval
   * [x1, x1 + cube_size - 1] (and a second index x [0, cube_size -1]) and a
   * grid that goes from [0.. grid_size - 1].
   *
   * We start from x1 and compute the largest interval [x1.. x1 + diff] that fit
   * to [0.. grid_size - 1]. Computing the difference diff is simply
   * min(grid_size - x1, cube_size - x). then we add the result in a vectorized
   * fashion. we itterate the processus by reducing the interval [x1, x1 +
   * cube_size - 1] until it is empty. */

  for (int z = 0; (z < handler->cube.size[0]); z++) {
    const int z1 = map[0][z];

    if (!is_point_in_interval(z1, zwindow))
      continue;

    compute_interval(map[0], handler->grid.full_size[0], handler->grid.size[0],
                     handler->cube.size[0], z1, &z, lower_corner, upper_corner,
                     zwindow);

    /* // We have a full plane. */
    if (upper_corner[0] - lower_corner[0]) {
      for (int y = 0; y < handler->cube.size[1]; y++) {
        const int y1 = map[1][y];

        // this check is completely irrelevant when running without MPI.
        if (!is_point_in_interval(y1, ywindow))
          continue;

        compute_interval(map[1], handler->grid.full_size[1],
                         handler->grid.size[1], handler->cube.size[1], y1, &y,
                         lower_corner + 1, upper_corner + 1, ywindow);

        if (upper_corner[1] - lower_corner[1]) {
          for (int x = 0; x < handler->cube.size[2]; x++) {
            const int x1 = map[2][x];

            if (!is_point_in_interval(x1, xwindow))
              continue;

            compute_interval(map[2], handler->grid.full_size[2],
                             handler->grid.size[2], handler->cube.size[2], x1,
                             &x, lower_corner + 2, upper_corner + 2, xwindow);

            if (upper_corner[2] - lower_corner[2]) {
              const int position1[3] = {z, y, x};

              /* the function will internally take care of the local vx global
               * grid */
              add_sub_grid(
                  lower_corner, // lower corner of the portion of cube (in the
                  // full grid)
                  upper_corner, // upper corner of the portion of cube (in the
                  // full grid)
                  position1,       // starting position subblock inside the cube
                  &handler->cube,  // the cube to extract data from
                  &handler->grid); // the grid to add data from
              if (handler->grid.size[2] == handler->grid.full_size[2])
                update_loop_index(handler->grid.full_size[2], x1, &x);
              else
                x += upper_corner[2] - lower_corner[2] - 1;
            }
          }
          if (handler->grid.size[1] == handler->grid.full_size[1])
            update_loop_index(handler->grid.full_size[1], y1, &y);
          else
            y += upper_corner[1] - lower_corner[1] - 1;
        }
      }
      if (handler->grid.size[0] == handler->grid.full_size[0])
        update_loop_index(handler->grid.full_size[0], z1, &z);
      else
        z += upper_corner[0] - lower_corner[0] - 1;
    }
  }
}

// *****************************************************************************
void grid_collocate(collocation_integration *const handler,
                    const bool use_ortho, const double zetp, const double rp[3],
                    const double radius) {
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
  int pol_offset[3] = {0, 0, 0};
  double roffset[3];
  double disr_radius;
  /* cube : grid containing pointlike product between polynomials
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
  initialize_tensor_3(&handler->pol, 3, handler->coef.size[0], cmax);
  realloc_tensor(&handler->pol);
  memset(handler->pol.data, 0, sizeof(double) * handler->pol.alloc_size_);

  /* compute the polynomials */

  // WARNING : do not reverse the order in pol otherwise you will have to
  // reverse the order in collocate_dgemm as well.

  if (use_ortho) {
    grid_fill_pol_dgemm(false, handler->dh[0][0], roffset[2], pol_offset[2],
                        lb_cube[2], ub_cube[2], handler->coef.size[2] - 1, cmax,
                        zetp, &idx3(handler->pol, 2, 0, 0)); /* i indice */
    grid_fill_pol_dgemm(false, handler->dh[1][1], roffset[1], pol_offset[1],
                        lb_cube[1], ub_cube[1], handler->coef.size[1] - 1, cmax,
                        zetp, &idx3(handler->pol, 1, 0, 0)); /* j indice */
    grid_fill_pol_dgemm(false, handler->dh[2][2], roffset[0], pol_offset[0],
                        lb_cube[0], ub_cube[0], handler->coef.size[0] - 1, cmax,
                        zetp, &idx3(handler->pol, 0, 0, 0)); /* k indice */
  } else {
    grid_fill_pol_dgemm(false, 1.0, roffset[0], pol_offset[2], lb_cube[0],
                        ub_cube[0], handler->coef.size[0] - 1, cmax,
                        zetp * handler->dx[0],
                        &idx3(handler->pol, 0, 0, 0)); /* k indice */
    grid_fill_pol_dgemm(false, 1.0, roffset[1], pol_offset[1], lb_cube[1],
                        ub_cube[1], handler->coef.size[1] - 1, cmax,
                        zetp * handler->dx[1],
                        &idx3(handler->pol, 1, 0, 0)); /* j indice */
    grid_fill_pol_dgemm(false, 1.0, roffset[2], pol_offset[0], lb_cube[2],
                        ub_cube[2], handler->coef.size[2] - 1, cmax,
                        zetp * handler->dx[2],
                        &idx3(handler->pol, 2, 0, 0)); /* i indice */

    calculate_non_orthorombic_corrections_tensor(
        zetp, roffset, handler->dh, lb_cube, ub_cube, handler->orthogonal,
        &handler->Exp);

    /* Use a slightly modified version of Ole code */
    grid_transform_coef_xzy_to_ikj(handler->dh, &handler->coef);
  }

  /* allocate memory for the polynomial and the cube */

  initialize_tensor_3(&handler->cube, cube_size[0], cube_size[1], cube_size[2]);

  realloc_tensor(&handler->cube);
  initialize_W_and_T(handler, &handler->cube, &handler->coef);

  tensor_reduction_for_collocate_integrate(
      handler->scratch, // pointer to scratch memory
      1.0, handler->orthogonal, &handler->Exp, &handler->coef, &handler->pol,
      &handler->cube);

  if (handler->apply_cutoff) {
    if (use_ortho) {
      apply_sphere_cutoff_ortho(handler, disr_radius, cmax, lb_cube,
                                cubecenter);
    } else {
      apply_spherical_cutoff_generic(handler, radius, cmax, lb_cube, ub_cube,
                                     roffset, cubecenter);
    }
    return;
  }
  apply_mapping_cubic(handler, cmax, lb_cube, cubecenter);
}

//******************************************************************************
void grid_collocate_pgf_product_cpu_dgemm(
    const bool use_ortho, const int border_mask, const int func,
    const int la_max, const int la_min, const int lb_max, const int lb_min,
    const double zeta, const double zetb, const double rscale,
    const double dh[3][3], const double dh_inv[3][3], const double ra[3],
    const double rab[3], const int grid_global_size[3],
    const int grid_local_size[3], const int shift_local[3],
    const int border_width[3], const double radius, const int o1, const int o2,
    const int n1, const int n2, const double pab_[n2][n1],
    double *const grid_) {
  collocation_integration *handler = collocate_create_handle();

  tensor pab;
  initialize_tensor_2(&pab, n2, n1);
  alloc_tensor(&pab);
  memcpy(pab.data, pab_, sizeof(double) * n1 * n2);
  // Uncomment this to dump all tasks to file.
  // #define __GRID_DUMP_TASKS
  int offset[2] = {o1, o2};

  int lmax[2] = {la_max, lb_max};
  int lmin[2] = {la_min, lb_min};

  const double zetp = zeta + zetb;
  const double f = zetb / zetp;
  const double rab2 = rab[0] * rab[0] + rab[1] * rab[1] + rab[2] * rab[2];
  const double prefactor = rscale * exp(-zeta * f * rab2);
  const double zeta_pair[2] = {zeta, zetb};
  initialize_basis_vectors(handler, dh, dh_inv);
  verify_orthogonality(dh, handler->orthogonal);

  initialize_tensor_3(&(handler->grid), grid_local_size[2], grid_local_size[1],
                      grid_local_size[0]);

  handler->grid.ld_ = grid_local_size[0];
  handler->grid.data = grid_;

  setup_global_grid_size(&handler->grid, (const int *const)grid_global_size);

  initialize_tensor_3(&handler->grid, grid_local_size[2], grid_local_size[1],
                      grid_local_size[0]);
  handler->grid.ld_ = grid_local_size[0];

  setup_grid_window(&handler->grid, shift_local, border_width, border_mask);

  double rp[3], rb[3];
  for (int i = 0; i < 3; i++) {
    rp[i] = ra[i] + f * rab[i];
    rb[i] = ra[i] + rab[i];
  }

  int lmin_diff[2], lmax_diff[2];
  grid_prepare_get_ldiffs_dgemm(func, lmin_diff, lmax_diff);

  int lmin_prep[2];
  int lmax_prep[2];

  lmin_prep[0] = imax(lmin[0] + lmin_diff[0], 0);
  lmin_prep[1] = imax(lmin[1] + lmin_diff[1], 0);

  lmax_prep[0] = lmax[0] + lmax_diff[0];
  lmax_prep[1] = lmax[1] + lmax_diff[1];

  const int n1_prep = ncoset(lmax_prep[0]);
  const int n2_prep = ncoset(lmax_prep[1]);

  /* I really do not like this. This will disappear */
  tensor pab_prep;
  initialize_tensor_2(&pab_prep, n2_prep, n1_prep);
  alloc_tensor(&pab_prep);
  memset(pab_prep.data, 0, pab_prep.alloc_size_ * sizeof(double));

  grid_prepare_pab_dgemm(func, offset, lmax, lmin, &zeta_pair[0], &pab,
                         &pab_prep);

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
  //

  /* precautionary tail since I probably intitialize data to NULL when I
   * initialize a new tensor. I want to keep the memory (I put a ridiculous
   * amount already) */

  initialize_tensor_4(&(handler->alpha), 3, lmax_prep[1] + 1, lmax_prep[0] + 1,
                      lmax_prep[0] + lmax_prep[1] + 1);
  realloc_tensor(&(handler->alpha));

  const int lp = lmax_prep[0] + lmax_prep[1];

  initialize_tensor_3(&(handler->coef), lp + 1, lp + 1, lp + 1);
  realloc_tensor(&(handler->coef));

  // initialy cp2k stores coef_xyz as coef[z][y][x]. this is fine but I
  // need them to be stored as

  grid_prepare_alpha_dgemm(ra, rb, rp, lmax_prep, &(handler->alpha));

  //
  //   compute P_{lxp,lyp,lzp} given P_{lxa,lya,lza,lxb,lyb,lzb} and
  //   alpha(ls,lxa,lxb,1) use a three step procedure we don't store zeros, so
  //   counting is done using lxyz,lxy in order to have contiguous memory access
  //   in collocate_fast.F
  //

  // coef[x][z][y]
  grid_compute_coefficients_dgemm(lmin_prep, lmax_prep, lp, prefactor,
                                  &(handler->alpha), &pab_prep,
                                  &(handler->coef));

  grid_collocate(handler, use_ortho, zetp, rp, radius);

  collocate_destroy_handle(handler);
  free(pab_prep.data);
}

void extract_blocks(grid_context *const ctx, const _task *const task,
                    const grid_buffer *pab_blocks, tensor *const work,
                    tensor *const pab) {
  const int iatom = task->iatom;
  const int jatom = task->jatom;
  const int iset = task->iset;
  const int jset = task->jset;
  const int ikind = ctx->atom_kinds[iatom];
  const int jkind = ctx->atom_kinds[jatom];
  const grid_basis_set *ibasis = ctx->basis_sets[ikind];
  const grid_basis_set *jbasis = ctx->basis_sets[jkind];

  const int block_num = task->block_num;

  // Locate current matrix block within the buffer. This block
  // contains the weights of the gaussian pairs in the spherical
  // harmonic basis, but we do computation in the cartesian
  // harmonic basis so we have to rotate the coefficients. It is nothing
  // else than a basis change and it done with two dgemm.

  const int block_offset = ctx->block_offsets[block_num]; // zero based
  double *const block = &pab_blocks->host_buffer[block_offset];

  rotate_to_cartesian_harmonics(ibasis, jbasis, iatom, jatom, iset, jset, block,
                                work, pab);
}

void compute_coefficients(grid_context *const ctx,
                          struct collocation_integration_ *handler,
                          const _task *const previous_task,
                          const _task *const task,
                          const grid_buffer *pab_blocks, tensor *const pab,
                          tensor *const work, tensor *const pab_prep) {
  // Load subblock from buffer and decontract into Cartesian sublock pab.
  // The previous pab can be reused when only ipgf or jpgf has changed.
  if (task->update_block_ || (previous_task == NULL)) {
    extract_blocks(ctx, task, pab_blocks, work, pab);
  }

  int lmin_prep[2];
  int lmax_prep[2];

  lmin_prep[0] = imax(task->lmin[0] + handler->lmin_diff[0], 0);
  lmin_prep[1] = imax(task->lmin[1] + handler->lmin_diff[1], 0);

  lmax_prep[0] = task->lmax[0] + handler->lmax_diff[0];
  lmax_prep[1] = task->lmax[1] + handler->lmax_diff[1];

  const int n1_prep = ncoset(lmax_prep[0]);
  const int n2_prep = ncoset(lmax_prep[1]);

  /* we do not reallocate memory. We initialized the structure with the
   * maximum lmax of the all list already.
   */
  initialize_tensor_2(pab_prep, n2_prep, n1_prep);
  realloc_tensor(pab_prep);

  grid_prepare_pab_dgemm(handler->func, task->offset, task->lmin, task->lmax,
                         &task->zeta[0], pab, pab_prep);

  //   *** initialise the coefficient matrix, we transform the sum
  //
  // sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} *
  //         (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb
  //         (y-a_y)**lya (z-a_z)**lza
  //
  // into
  //
  // sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp
  // (z-p_z)**lzp
  //
  // where p is center of the product gaussian, and lp = la_max + lb_max
  // (current implementation is l**7)
  //

  /* precautionary tail since I probably intitialize data to NULL when I
   * initialize a new tensor. I want to keep the memory (I put a ridiculous
   * amount already) */

  initialize_tensor_4(&handler->alpha, 3, lmax_prep[1] + 1, lmax_prep[0] + 1,
                      lmax_prep[0] + lmax_prep[1] + 1);
  realloc_tensor(&handler->alpha);

  const int lp = lmax_prep[0] + lmax_prep[1];

  initialize_tensor_3(&handler->coef, lp + 1, lp + 1, lp + 1);
  realloc_tensor(&handler->coef);

  // these two functions can be done with dgemm again....

  grid_prepare_alpha_dgemm(task->ra, task->rb, task->rp, lmax_prep,
                           &handler->alpha);

  // compute the coefficients after applying the function of interest
  // coef[x][z][y]
  grid_compute_coefficients_dgemm(
      lmin_prep, lmax_prep, lp,
      task->prefactor * ((task->iatom == task->jatom) ? 1.0 : 2.0),
      &handler->alpha, pab_prep, &handler->coef);
}

void collocate_one_grid_level_dgemm(grid_context *const ctx,
                                    const int *const border_width,
                                    const int *const shift_local,
                                    const int func, const int level,
                                    const grid_buffer *pab_blocks) {
  tensor *const grid = &ctx->grid[level];
  // Using default(shared) because with GCC 9 the behavior around const changed:
  // https://www.gnu.org/software/gcc/gcc-9/porting_to.html
#pragma omp parallel default(shared)
  {
    const int num_threads = omp_get_num_threads();
    const int thread_id = omp_get_thread_num();

    tensor work, pab, pab_prep;

    struct collocation_integration_ *handler = ctx->handler[thread_id];

    handler->func = func;
    grid_prepare_get_ldiffs_dgemm(handler->func, handler->lmin_diff,
                                  handler->lmax_diff);

    handler->apply_cutoff = ctx->apply_cutoff;

    // Allocate pab matrix for re-use across tasks.
    initialize_tensor_2(&pab, ctx->maxco, ctx->maxco);
    alloc_tensor(&pab);

    initialize_tensor_2(&work, ctx->maxco, ctx->maxco);
    alloc_tensor(&work);

    initialize_tensor_2(&pab_prep, ctx->maxco, ctx->maxco);
    alloc_tensor(&pab_prep);

    initialize_basis_vectors(handler, grid->dh, grid->dh_inv);

    /* setup the grid parameters, window parameters (if the grid is split), etc
     */

    tensor_copy(&handler->grid, grid);

    for (int d = 0; d < 3; d++)
      handler->orthogonal[d] = handler->grid.orthogonal[d];

    if ((thread_id == 0) || (num_threads == 1)) {
      // thread id 0 directly store the results in the final storage space
      handler->grid.data = ctx->grid[level].data;
    }

    if ((num_threads > 1) && (thread_id > 0)) {
      handler->grid.data = ((double *)ctx->scratch) +
                           (thread_id - 1) * handler->grid.alloc_size_;
      memset(handler->grid.data, 0, sizeof(double) * grid->alloc_size_);
    }

    /* it is only useful when we split the list over multiple threads. The first
     * iteration should load the block whatever status the task->block_update_
     * has */
    const _task *prev_task = NULL;
#pragma omp for schedule(static)
    for (int itask = 0; itask < ctx->tasks_per_level[level]; itask++) {
      // Define some convenient aliases.
      const _task *task = &ctx->tasks[level][itask];

      if (task->level != level) {
        printf("level %d, %d\n", task->level, level);
        abort();
      }
      /* the grid is divided over several ranks or not periodic */
      if ((handler->grid.size[0] != handler->grid.full_size[0]) ||
          (handler->grid.size[1] != handler->grid.full_size[1]) ||
          (handler->grid.size[2] != handler->grid.full_size[2])) {
        /* unfortunately the window where the gaussian should be added depends
         * on the bonds. So I have to adjust the window all the time. */

        setup_grid_window(&handler->grid, shift_local, border_width,
                          task->border_mask);
      }

      /* this is a three steps procedure. pab_blocks contains the coefficients
       * of the operator in the spherical harmonic basis while we do computation
       * in the cartesian harmonic basis.
       *
       * step 1 : rotate the coefficients from the harmonic to the cartesian
       * basis

       * step 2 : extract the subblock and apply additional transformations
       * corresponding the spatial derivatives of the operator (the last is not
       * always done)

       * step 3 : change from (x - x_1)^\alpha (x - x_2)^\beta to (x -
       * x_{12})^k. It is a transformation which involves binomial
       * coefficients.
       *
       * \f[ (x - x_1) ^\alpha (x - x_2) ^ beta = \sum_{k_{1} k_{2}} ^
       *     {\alpha\beta} \text{Binomial}(\alpha,k_1)
       *     \text{Binomial}(\beta,k_2) (x - x_{12})^{k_1 + k_2} (x_12 - x_1)
       *     ^{\alpha - k_1} (x_12 - x_2) ^{\beta - k_2} ]
       *
       * step 1 is done only when necessary, the two remaining steps are done
       * for each bond.
       */

      compute_coefficients(ctx, handler, prev_task, task, pab_blocks, &pab,
                           &work, &pab_prep);

      grid_collocate(handler, ctx->orthorhombic, task->zetp, task->rp,
                     task->radius);
      prev_task = task;
    }

    // Merge thread local grids into shared grid. Could be improved though....

    // thread 0 does nothing since the data are already placed in the final
    // destination
    if (num_threads > 1) {
      if ((grid->alloc_size_ / (num_threads - 1)) >= 2) {
        const int block_size = grid->alloc_size_ / (num_threads - 1) +
                               (grid->alloc_size_ % (num_threads - 1));

        for (int bk = 0; bk < num_threads; bk++) {
          if (thread_id > 0) {
            int bk_id = (bk + thread_id - 1) % num_threads;
            int begin = bk_id * block_size;
            int end = imin((bk_id + 1) * block_size, grid->alloc_size_);
            cblas_daxpy(end - begin, 1.0, handler->grid.data + begin, 1,
                        grid->data + begin, 1);
          }
#pragma omp barrier
        }
      }
    } else {
#pragma omp critical
      if (thread_id > 0)
        cblas_daxpy(handler->grid.alloc_size_, 1.0,
                    &idx3(handler->grid, 0, 0, 0), 1, &idx3(grid[0], 0, 0, 0),
                    1);
    }
    handler->grid.data = NULL;
    free(pab.data);
    free(pab_prep.data);
    free(work.data);
  }
}

/*******************************************************************************
 * \brief Collocate all tasks of a given list onto given grids.
 *        See grid_task_list.h for details.
 ******************************************************************************/
void grid_cpu_collocate_task_list(
    grid_cpu_task_list *const ptr, const bool orthorhombic,
    const enum grid_func func, const int nlevels,
    const int npts_global[nlevels][3], const int npts_local[nlevels][3],
    const int shift_local[nlevels][3], const int border_width[nlevels][3],
    const double dh[nlevels][3][3], const double dh_inv[nlevels][3][3],
    const grid_buffer *pab_blocks, double *grid[nlevels]) {

  grid_context *const ctx = (grid_context *const)ptr;

  assert(ctx->checksum == ctx_checksum);

  ctx->orthorhombic = orthorhombic;

  const int max_threads = omp_get_max_threads();

  assert(ctx->nlevels == nlevels);

  //#pragma omp parallel for
  for (int level = 0; level < ctx->nlevels; level++) {
    set_grid_parameters(&ctx->grid[level], orthorhombic, npts_global[level],
                        npts_local[level], shift_local[level],
                        border_width[level], dh[level], dh_inv[level],
                        grid[level]);
    memset(ctx->grid[level].data, 0,
           sizeof(double) * ctx->grid[level].alloc_size_);
  }

  if (ctx->scratch == NULL) {
    int max_size = ctx->grid[0].alloc_size_;

    /* compute the size of the largest grid. It is used afterwards to allocate
     * scratch memory for the grid on each omp thread */
    for (int x = 1; x < nlevels; x++) {
      max_size = imax(ctx->grid[x].alloc_size_, max_size);
    }

    max_size = ((max_size / 4096) + (max_size % 4096 != 0)) * 4096;

    /* scratch is a void pointer !!!!! */
    ctx->scratch =
        grid_allocate_scratch(max_size * max_threads * sizeof(double));
  }

  for (int level = 0; level < ctx->nlevels; level++) {
    collocate_one_grid_level_dgemm(ctx, border_width[level], shift_local[level],
                                   func, level, pab_blocks);
  }

  grid_free_scratch(ctx->scratch);
  ctx->scratch = NULL;
}
