/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2021 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include "non_orthorombic_corrections.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/grid_common.h"
#include "utils.h"

double exp_recursive(const double c_exp, const double c_exp_minus_1,
                     const int index) {
  if (index == -1)
    return c_exp_minus_1;

  if (index == 1)
    return c_exp;

  double res = 1.0;

  if (index < 0) {
    for (int i = 0; i < -index; i++) {
      res *= c_exp_minus_1;
    }
    return res;
  }

  if (index > 0) {
    for (int i = 0; i < index; i++) {
      res *= c_exp;
    }
    return res;
  }

  return 1.0;
}

void exp_i(const double alpha, const int imin, const int imax,
           double *restrict const res) {
  const double c_exp_co = exp(alpha);
  /* const double c_exp_minus_1 = 1/ c_exp; */
  res[0] = exp(imin * alpha);
  for (int i = 1; i < (imax - imin); i++) {
    res[i] = res[i - 1] *
             c_exp_co; // exp_recursive(c_exp_co, 1.0 / c_exp_co, i + imin);
  }
}

void exp_ij(const double alpha, const int offset_i, const int imin,
            const int imax, const int offset_j, const int jmin, const int jmax,
            tensor *exp_ij_) {
  double c_exp = exp(alpha * imin);
  const double c_exp_co = exp(alpha);

  for (int i = 0; i < (imax - imin); i++) {
    double *restrict dst = &idx2(exp_ij_[0], i + offset_i, offset_j);
    double ctmp = exp_recursive(c_exp, 1.0 / c_exp, jmin);

    for (int j = 0; j < (jmax - jmin); j++) {
      dst[j] *= ctmp;
      ctmp *= c_exp;
    }
    c_exp *= c_exp_co;
  }
}

void calculate_non_orthorombic_corrections_tensor(
    const double mu_mean, const double *r_ab, const double basis[3][3],
    const int *const xmin, const int *const xmax, bool *plane,
    tensor *const Exp) {
  // zx, zy, yx
  const int n[3][2] = {{0, 2}, {0, 1}, {1, 2}};

  // need to review this
  const double c[3] = {
      /* alpha gamma */
      -2.0 * mu_mean *
          (basis[0][0] * basis[2][0] + basis[0][1] * basis[2][1] +
           basis[0][2] * basis[2][2]),
      /* beta gamma */
      -2.0 * mu_mean *
          (basis[1][0] * basis[2][0] + basis[1][1] * basis[2][1] +
           basis[1][2] * basis[2][2]),
      /* alpha beta */
      -2.0 * mu_mean *
          (basis[0][0] * basis[1][0] + basis[0][1] * basis[1][1] +
           basis[0][2] * basis[1][2])};

  /* a naive implementation of the computation of exp(-2 (v_i . v_j) (i
   * - r_i) (j _ r_j)) requires n m exponentials but we can do much much
   * better with only 7 exponentials
   *
   * first expand it. we get exp(2 (v_i . v_j) i j) exp(2 (v_i . v_j) i r_j)
   * exp(2 (v_i . v_j) j r_i) exp(2 (v_i . v_j) r_i r_j). we can use the fact
   * that the sum of two terms in an exponential is equal to the product of
   * the two exponentials.
   *
   * this means that exp (a i) with i integer can be computed recursively with
   * one exponential only
   */

  /* we have a orthorombic case */
  if (plane[0] && plane[1] && plane[2])
    return;

  tensor exp_tmp;
  double *x1, *x2;
  /* printf("%d %d %d\n", plane[0], plane[1], plane[2]); */
  const int max_elem =
      imax(imax(xmax[0] - xmin[0], xmax[1] - xmin[1]), xmax[2] - xmin[2]) + 1;
  initialize_tensor_3(Exp, 3, max_elem, max_elem);
  realloc_tensor(Exp);

  x1 = grid_allocate_scratch(sizeof(double) * max_elem);
  x2 = grid_allocate_scratch(sizeof(double) * max_elem);
  initialize_tensor_2(&exp_tmp, Exp->size[1], Exp->size[2]);

  memset(&idx3(Exp[0], 0, 0, 0), 0, sizeof(double) * Exp->alloc_size_);
  for (int dir = 0; dir < 3; dir++) {
    int d1 = n[dir][0];
    int d2 = n[dir][1];

    if (!plane[dir]) {
      const double c_exp_const = exp(c[dir] * r_ab[d1] * r_ab[d2]);

      exp_i(-r_ab[d2] * c[dir], xmin[d1], xmax[d1] + 1, x1);
      exp_i(-r_ab[d1] * c[dir], xmin[d2], xmax[d2] + 1, x2);

      exp_tmp.data = &idx3(Exp[0], dir, 0, 0);
      cblas_dger(CblasRowMajor, xmax[d1] - xmin[d1] + 1,
                 xmax[d2] - xmin[d2] + 1, c_exp_const, x1, 1, x2, 1,
                 &idx2(exp_tmp, 0, 0), exp_tmp.ld_);
      exp_ij(c[dir], 0, xmin[d1], xmax[d1] + 1, 0, xmin[d2], xmax[d2] + 1,
             &exp_tmp);
    }
  }
  grid_free_scratch(x1);
  grid_free_scratch(x2);
}

void calculate_non_orthorombic_corrections_tensor_blocked(
    const double mu_mean, const double *r_ab, const double basis[3][3],
    const int *const lower_corner, const int *const upper_corner,
    const int *const block_size, const int *const offset, const int *const xmin,
    const int *const xmax, bool *plane, tensor *const Exp) {
  // zx, zy, yx
  const int n[3][2] = {{0, 2}, {0, 1}, {1, 2}};

  // need to review this
  const double c[3] = {
      /* alpha gamma */
      -2.0 * mu_mean *
          (basis[0][0] * basis[2][0] + basis[0][1] * basis[2][1] +
           basis[0][2] * basis[2][2]),
      /* beta gamma */
      -2.0 * mu_mean *
          (basis[1][0] * basis[2][0] + basis[1][1] * basis[2][1] +
           basis[1][2] * basis[2][2]),
      /* alpha beta */
      -2.0 * mu_mean *
          (basis[0][0] * basis[1][0] + basis[0][1] * basis[1][1] +
           basis[0][2] * basis[1][2])};

  /* a naive implementation of the computation of exp(-2 (v_i . v_j) (i
   * - r_i) (j _ r_j)) requires n m exponentials but we can do much much
   * better with only 7 exponentials
   *
   * first expand it. we get exp(2 (v_i . v_j) i j) exp(2 (v_i . v_j) i r_j)
   * exp(2 (v_i . v_j) j r_i) exp(2 (v_i . v_j) r_i r_j). we can use the fact
   * that the sum of two terms in an exponential is equal to the product of
   * the two exponentials.
   *
   * this means that exp (a i) with i integer can be computed recursively with
   * one exponential only
   */

  /* we have a orthorombic case */
  if (plane[0] && plane[1] && plane[2])
    return;

  tensor exp_blocked;
  double *x1, *x2;
  /* printf("%d %d %d\n", plane[0], plane[1], plane[2]); */
  initialize_tensor_2(&exp_blocked, imax(block_size[0], block_size[1]),
                      imax(block_size[1], block_size[2]));

  const int cube_size[3] = {(upper_corner[0] - lower_corner[0]) * block_size[0],
                            (upper_corner[1] - lower_corner[1]) * block_size[1],
                            (upper_corner[2] - lower_corner[2]) *
                                block_size[2]};

  const int max_elem = imax(imax(cube_size[0], cube_size[1]), cube_size[2]);
  x1 = grid_allocate_scratch(sizeof(double) * max_elem);
  x2 = grid_allocate_scratch(sizeof(double) * max_elem);

  initialize_tensor_4(Exp, 3,
                      imax(upper_corner[0] - lower_corner[0],
                           upper_corner[1] - lower_corner[1]),
                      imax(upper_corner[2] - lower_corner[2],
                           upper_corner[1] - lower_corner[1]),
                      exp_blocked.alloc_size_);

  realloc_tensor(Exp);
  /* memset(Exp->data, 0, sizeof(double) * Exp->alloc_size_); */

  for (int dir = 0; dir < 3; dir++) {
    int d1 = n[dir][0];
    int d2 = n[dir][1];

    if (!plane[dir]) {
      const double c_exp_const = exp(c[dir] * r_ab[d1] * r_ab[d2]);
      memset(x1, 0, sizeof(double) * max_elem);
      memset(x2, 0, sizeof(double) * max_elem);
      /* memset(exp_tmp.data, 0, sizeof(double) * exp_tmp.alloc_size_); */
      exp_i(-r_ab[d2] * c[dir], xmin[d1], xmax[d1] + 1, x1 + offset[d1]);
      exp_i(-r_ab[d1] * c[dir], xmin[d2], xmax[d2] + 1, x2 + offset[d2]);

      for (int y = 0; y < (upper_corner[d1] - lower_corner[d1]); y++) {
        const int y_1 = y * block_size[d1];
        for (int x = 0; x < (upper_corner[d2] - lower_corner[d2]); x++) {
          const int x_1 = x * block_size[d2];
          exp_blocked.data = &idx4(Exp[0], dir, y, x, 0);

          for (int y2 = 0; y2 < block_size[d1]; y2++) {
            double *restrict dst = &idx2(exp_blocked, y2, 0);
            const double scal = x1[y_1 + y2] * c_exp_const;
            const double *restrict src = &x2[x_1];
            //#pragma omp simd linear(dst, src) simdlen(8)
            GRID_PRAGMA_SIMD((dst, src), 8)
            for (int x3 = 0; x3 < block_size[d2]; x3++) {
              dst[x3] = scal * src[x3];
            }
          }

          exp_ij(c[dir], 0, xmin[d1] - offset[d1] + y_1,
                 xmin[d1] - offset[d1] + y_1 + block_size[d1], 0,
                 xmin[d2] - offset[d2] + x_1,
                 xmin[d2] - offset[d2] + x_1 + block_size[d2], &exp_blocked);
        }
      }
    }
  }

  grid_free_scratch(x1);
  grid_free_scratch(x2);
  /* free(exp_tmp.data); */
}

void apply_non_orthorombic_corrections(const bool *restrict plane,
                                       const tensor *const Exp,
                                       tensor *const cube) {
  // Well we should never call non orthorombic corrections if everything is
  // orthorombic
  if (plane[0] && plane[1] && plane[2])
    return;

  /*k and i are orthogonal, k and j as well */
  if (plane[0] && plane[1]) {
    for (int z = 0; z < cube->size[0]; z++) {
      for (int y = 0; y < cube->size[1]; y++) {
        const double *restrict yx = &idx3(Exp[0], 2, y, 0);
        double *restrict dst = &idx3(cube[0], z, y, 0);

        //#pragma omp simd linear(dst, yx) simdlen(8)
        GRID_PRAGMA_SIMD((dst, yx), 8)
        for (int x = 0; x < cube->size[2]; x++) {
          dst[x] *= yx[x];
        }
      }
    }
    return;
  }

  /* k and i are orhogonal, i and j as well */
  if (plane[0] && plane[2]) {
    for (int z = 0; z < cube->size[0]; z++) {
      for (int y = 0; y < cube->size[1]; y++) {
        const double zy = idx3(Exp[0], 1, z, y);
        double *restrict dst = &idx3(cube[0], z, y, 0);

        //#pragma omp simd linear(dst) simdlen(8)
        GRID_PRAGMA_SIMD((dst), 8)
        for (int x = 0; x < cube->size[2]; x++) {
          dst[x] *= zy;
        }
      }
    }
    return;
  }

  /* j, k are orthognal, i and j are orthognal */
  if (plane[1] && plane[2]) {
    for (int z = 0; z < cube->size[0]; z++) {
      double *restrict zx = &idx3(Exp[0], 0, z, 0);
      for (int y = 0; y < cube->size[1]; y++) {
        double *restrict dst = &idx3(cube[0], z, y, 0);
        //#pragma omp simd linear(dst, zx) simdlen(8)
        GRID_PRAGMA_SIMD((dst, zx), 8)
        for (int x = 0; x < cube->size[2]; x++) {
          dst[x] *= zx[x];
        }
      }
    }
    return;
  }

  if (plane[0]) {
    // z perpendicular to x. but y non perpendicular to any
    for (int z = 0; z < cube->size[0]; z++) {
      for (int y = 0; y < cube->size[1]; y++) {
        const double zy = idx3(Exp[0], 1, z, y);
        const double *restrict yx = &idx3(Exp[0], 2, y, 0);
        double *restrict dst = &idx3(cube[0], z, y, 0);

        //#pragma omp simd linear(dst, yx) simdlen(8)
        GRID_PRAGMA_SIMD((dst, yx), 8)
        for (int x = 0; x < cube->size[2]; x++) {
          dst[x] *= zy * yx[x];
        }
      }
    }
    return;
  }

  if (plane[1]) {
    // z perpendicular to y, but x and z are not and y and x neither
    for (int z = 0; z < cube->size[0]; z++) {
      double *restrict zx = &idx3(Exp[0], 0, z, 0);
      for (int y = 0; y < cube->size[1]; y++) {
        const double *restrict yx = &idx3(Exp[0], 2, y, 0);
        double *restrict dst = &idx3(cube[0], z, y, 0);
        //#pragma omp simd linear(dst, yx) simdlen(8)
        GRID_PRAGMA_SIMD((dst, yx), 8)
        for (int x = 0; x < cube->size[2]; x++) {
          dst[x] *= zx[x] * yx[x];
        }
      }
    }
    return;
  }

  if (plane[2]) {
    // x perpendicular to y, but x and z are not and y and z neither
    for (int z = 0; z < cube->size[0]; z++) {
      double *restrict zx = &idx3(Exp[0], 0, z, 0);
      for (int y = 0; y < cube->size[1]; y++) {
        const double zy = idx3(Exp[0], 1, z, y);
        double *restrict dst = &idx3(cube[0], z, y, 0);

        //#pragma omp simd linear(dst) simdlen(8)
        GRID_PRAGMA_SIMD((dst), 8)
        for (int x = 0; x < cube->size[2]; x++) {
          dst[x] *= zx[x] * zy;
        }
      }
    }
    return;
  }

  /* generic  case */

  for (int z = 0; z < cube->size[0]; z++) {
    double *restrict zx = &idx3(Exp[0], 0, z, 0);
    for (int y = 0; y < cube->size[1]; y++) {
      const double zy = idx3(Exp[0], 1, z, y);
      const double *restrict yx = &idx3(Exp[0], 2, y, 0);
      double *restrict dst = &idx3(cube[0], z, y, 0);

      //#pragma omp simd linear(dst, zx, yx) simdlen(8)
      GRID_PRAGMA_SIMD((dst, zx), 8)
      for (int x = 0; x < cube->size[2]; x++) {
        dst[x] *= zx[x] * zy * yx[x];
      }
    }
  }
  return;
}

void apply_non_orthorombic_corrections_xy_blocked(
    const struct tensor_ *const Exp, struct tensor_ *const m) {
  for (int gamma = 0; gamma < m->size[0]; gamma++) {
    for (int y1 = 0; y1 < m->size[1]; y1++) {
      double *restrict dst = &idx3(m[0], gamma, y1, 0);
      const double *restrict src = &idx2(Exp[0], y1, 0);

      //#pragma omp simd linear(dst, src) simdlen(8)
      GRID_PRAGMA_SIMD((dst, src), 8)
      for (int x1 = 0; x1 < m->size[2]; x1++) {
        dst[x1] *= src[x1];
      }
    }
  }
}

void apply_non_orthorombic_corrections_xz_blocked(
    const struct tensor_ *const Exp, struct tensor_ *const m) {
  for (int z1 = 0; z1 < m->size[0]; z1++) {
    const double *restrict src = &idx2(Exp[0], z1, 0);
    for (int y1 = 0; y1 < m->size[1]; y1++) {
      double *restrict dst = &idx3(m[0], z1, y1, 0);
      //#pragma omp simd linear(dst, src) simdlen(8)
      GRID_PRAGMA_SIMD((dst, src), 8)
      for (int x1 = 0; x1 < m->size[2]; x1++) {
        dst[x1] *= src[x1];
      }
    }
  }
}

void apply_non_orthorombic_corrections_yz_blocked(
    const struct tensor_ *const Exp, struct tensor_ *const m) {
  for (int z1 = 0; z1 < m->size[0]; z1++) {
    for (int y1 = 0; y1 < m->size[1]; y1++) {
      const double src = idx2(Exp[0], z1, y1);
      double *restrict dst = &idx3(m[0], z1, y1, 0);
      //#pragma omp simd linear(dst) simdlen(8)
      GRID_PRAGMA_SIMD((dst), 8)
      for (int x1 = 0; x1 < m->size[2]; x1++) {
        dst[x1] *= src;
      }
    }
  }
}

void apply_non_orthorombic_corrections_xz_yz_blocked(
    const struct tensor_ *const Exp_xz, const struct tensor_ *const Exp_yz,
    struct tensor_ *const m) {
  for (int z1 = 0; z1 < m->size[0]; z1++) {
    const double *restrict src_xz = &idx2(Exp_xz[0], z1, 0);
    for (int y1 = 0; y1 < m->size[1]; y1++) {
      const double src = idx2(Exp_yz[0], z1, y1);
      double *restrict dst = &idx3(m[0], z1, y1, 0);
      //#pragma omp simd linear(dst) simdlen(8)
      GRID_PRAGMA_SIMD((dst), 8)
      for (int x1 = 0; x1 < m->size[2]; x1++) {
        dst[x1] *= src * src_xz[x1];
      }
    }
  }
}
