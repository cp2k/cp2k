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
#include <string.h>

#ifdef __LIBXSMM
#include <libxsmm.h>
#endif

#include "../common/grid_common.h"
#include "tensor_local.h"
#include "utils.h"

void convert_to_lattice_coordinates(const double dh_inv_[3][3],
                                    const double *restrict const rp,
                                    double *restrict rp_c) {
  rp_c[0] =
      dh_inv_[0][0] * rp[0] + dh_inv_[1][0] * rp[1] + dh_inv_[0][0] * rp[2];
  rp_c[1] =
      dh_inv_[0][1] * rp[0] + dh_inv_[1][1] * rp[1] + dh_inv_[1][1] * rp[2];
  rp_c[2] =
      dh_inv_[0][2] * rp[0] + dh_inv_[1][2] * rp[1] + dh_inv_[2][2] * rp[2];
}

/* DO NOT CHANGE THIS. */

void dgemm_simplified(dgemm_params *const m, const bool use_libxsmm) {
  if (m == NULL)
    abort();

#if defined(__LIBXSMM)
  if (use_libxsmm && m->op2 == 'N') {
    /* we are in row major but xsmm is in column major */
    m->prefetch = LIBXSMM_PREFETCH_AUTO;
    /* in the future, more flags can be or'd together (like NONE | TRANS_B,
     * etc.)*/
    m->flags =
        (m->op1 != 'T' ? LIBXSMM_GEMM_FLAG_NONE : LIBXSMM_GEMM_FLAG_TRANS_B);

    if (m->kernel == NULL) {
      m->kernel =
          libxsmm_dmmdispatch(m->n, m->m, m->k, &m->ldb, &m->lda, &m->ldc,
                              &m->alpha, &m->beta, &m->flags, &m->prefetch);
    }

    if (m->kernel) {
      m->kernel(m->b, m->a, m->c, m->b, m->a, m->c);
      return;
    }
  }
#endif

#if defined(__MKL)
  // fall back to mkl
  if ((m->op1 == 'N') && (m->op2 == 'N'))
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m->m, m->n, m->k,
                m->alpha, m->a, m->lda, m->b, m->ldb, m->beta, m->c, m->ldc);

  if ((m->op1 == 'T') && (m->op2 == 'N'))
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, m->m, m->n, m->k,
                m->alpha, m->a, m->lda, m->b, m->ldb, m->beta, m->c, m->ldc);

  if ((m->op1 == 'N') && (m->op2 == 'T'))
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m->m, m->n, m->k,
                m->alpha, m->a, m->lda, m->b, m->ldb, m->beta, m->c, m->ldc);

  if ((m->op1 == 'T') && (m->op2 == 'T'))
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, m->m, m->n, m->k,
                m->alpha, m->a, m->lda, m->b, m->ldb, m->beta, m->c, m->ldc);

#else

  if ((m->op1 == 'N') && (m->op2 == 'N'))
    dgemm_("N", "N", &m->n, &m->m, &m->k, &m->alpha, m->b, &m->ldb, m->a,
           &m->lda, &m->beta, m->c, &m->ldc);

  if ((m->op1 == 'T') && (m->op2 == 'N'))
    dgemm_("N", "T", &m->n, &m->m, &m->k, &m->alpha, m->b, &m->ldb, m->a,
           &m->lda, &m->beta, m->c, &m->ldc);

  if ((m->op1 == 'T') && (m->op2 == 'T'))
    dgemm_("T", "T", &m->n, &m->m, &m->k, &m->alpha, m->b, &m->ldb, m->a,
           &m->lda, &m->beta, m->c, &m->ldc);

  if ((m->op1 == 'N') && (m->op2 == 'T'))
    dgemm_("T", "N", &m->n, &m->m, &m->k, &m->alpha, m->b, &m->ldb, m->a,
           &m->lda, &m->beta, m->c, &m->ldc);

#endif
}

void batched_dgemm_simplified(dgemm_params *const m, const int batch_size,
                              const bool use_libxsmm) {
  assert(m != NULL);
  assert(batch_size > 0);

#if defined(__LIBXSMM)
  if (use_libxsmm && m->op2 == 'N') {
    /* we are in row major but xsmm is in column major */
    m->prefetch = LIBXSMM_PREFETCH_AUTO;
    /* in the future, more flags can be or'd together (like NONE | TRANS_B,
     * etc.)*/
    m->flags =
        (m->op1 != 'T' ? LIBXSMM_GEMM_FLAG_NONE : LIBXSMM_GEMM_FLAG_TRANS_B);

    if (m->kernel == NULL) {
      m->kernel =
          libxsmm_dmmdispatch(m->n, m->m, m->k, &m->ldb, &m->lda, &m->ldc,
                              &m->alpha, &m->beta, &m->flags, &m->prefetch);
    }

    if (m->kernel) {
      for (int s = 0; s < batch_size - 1; s++) {
        m->kernel(m[s].b, m[s].a, m[s].c, m[s + 1].b, m[s + 1].a, m[s + 1].c);
      }
      m->kernel(m[batch_size - 1].b, m[batch_size - 1].a, m[batch_size - 1].c,
                m[batch_size - 1].b, m[batch_size - 1].a, m[batch_size - 1].c);
      return;
    }
  }
#endif

#if defined(__MKL)
  // fall back to mkl
  for (int s = 0; s < batch_size; s++) {
    if ((m[s].op1 == 'N') && (m[s].op2 == 'N'))
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m[s].m, m[s].n,
                  m[s].k, m[s].alpha, m[s].a, m[s].lda, m[s].b, m[s].ldb,
                  m[s].beta, m[s].c, m[s].ldc);

    if ((m[s].op1 == 'T') && (m[s].op2 == 'N'))
      cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, m[s].m, m[s].n,
                  m[s].k, m[s].alpha, m[s].a, m[s].lda, m[s].b, m[s].ldb,
                  m[s].beta, m[s].c, m[s].ldc);

    if ((m[s].op1 == 'N') && (m[s].op2 == 'T'))
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m[s].m, m[s].n,
                  m[s].k, m[s].alpha, m[s].a, m[s].lda, m[s].b, m[s].ldb,
                  m[s].beta, m[s].c, m[s].ldc);

    if ((m[s].op1 == 'T') && (m[s].op2 == 'T'))
      cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, m[s].m, m[s].n, m[s].k,
                  m[s].alpha, m[s].a, m[s].lda, m[s].b, m[s].ldb, m[s].beta,
                  m[s].c, m[s].ldc);
  }
#else
  for (int s = 0; s < batch_size; s++) {
    if ((m[s].op1 == 'N') && (m[s].op2 == 'N'))
      dgemm_("N", "N", &m[s].n, &m[s].m, &m[s].k, &m[s].alpha, m[s].b,
             &m[s].ldb, m[s].a, &m[s].lda, &m[s].beta, m[s].c, &m[s].ldc);

    if ((m[s].op1 == 'T') && (m[s].op2 == 'N'))
      dgemm_("N", "T", &m[s].n, &m[s].m, &m[s].k, &m[s].alpha, m[s].b,
             &m[s].ldb, m[s].a, &m[s].lda, &m[s].beta, m[s].c, &m[s].ldc);

    if ((m[s].op1 == 'T') && (m[s].op2 == 'T'))
      dgemm_("T", "T", &m[s].n, &m[s].m, &m[s].k, &m[s].alpha, m[s].b,
             &m[s].ldb, m[s].a, &m[s].lda, &m[s].beta, m[s].c, &m[s].ldc);

    if ((m[s].op1 == 'N') && (m[s].op2 == 'T'))
      dgemm_("T", "N", &m[s].n, &m[s].m, &m[s].k, &m[s].alpha, m[s].b,
             &m[s].ldb, m[s].a, &m[s].lda, &m[s].beta, m[s].c, &m[s].ldc);
  }
#endif
}

void extract_sub_grid(const int *lower_corner, const int *upper_corner,
                      const int *position, const tensor *const grid,
                      tensor *const subgrid) {
  int position1[3] = {0, 0, 0};

  if (position) {
    position1[0] = position[0];
    position1[1] = position[1];
    position1[2] = position[2];
  }

  const int sizex = upper_corner[2] - lower_corner[2];
  const int sizey = upper_corner[1] - lower_corner[1];
  const int sizez = upper_corner[0] - lower_corner[0];

  for (int z = 0; z < sizez; z++) {
    /* maybe use matcopy from libxsmm if possible */
    for (int y = 0; y < sizey; y++) {
      double *restrict src =
          &idx3(grid[0], lower_corner[0] + z - grid->window_shift[0],
                lower_corner[1] + y - grid->window_shift[1],
                lower_corner[2] - grid->window_shift[2]);
      double *restrict dst =
          &idx3(subgrid[0], position1[0] + z, position1[1] + y, position1[2]);
#ifdef __LIBXSMM
      LIBXSMM_PRAGMA_SIMD
#else
#pragma GCC ivdep
#endif
      for (int x = 0; x < sizex; x++) {
        dst[x] = src[x];
      }
    }
  }

  return;
}

void add_sub_grid(const int *lower_corner, const int *upper_corner,
                  const int *position, const tensor *subgrid, tensor *grid) {
  int position1[3] = {0, 0, 0};

  if (position) {
    position1[0] = position[0];
    position1[1] = position[1];
    position1[2] = position[2];
  }
  for (int d = 0; d < 3; d++) {
    if ((lower_corner[d] < grid->window_shift[d]) || (lower_corner[d] < 0) ||
        (lower_corner[d] >= upper_corner[d]) || (upper_corner[d] <= 0) ||
        (upper_corner[d] - lower_corner[d] > subgrid->size[d]) ||
        (grid == NULL) || (subgrid == NULL)) {
      printf("Error : invalid parameters. Values of the given parameters along "
             "the first wrong dimension\n");
      printf("      : lorner corner  [%d] = %d\n", d, lower_corner[d]);
      printf("      : upper  corner  [%d] = %d\n", d, upper_corner[d]);
      printf("      : diff           [%d] = %d\n", d,
             upper_corner[d] - lower_corner[d]);
      printf("      : src grid size  [%d] = %d\n", d, subgrid->size[d]);
      printf("      : dst grid size  [%d] = %d\n", d, grid->size[d]);
      printf("      : window dst grid size  [%d] = %d\n", d,
             grid->window_size[d]);
      printf("      : window dst shift  [%d] = %d\n", d, grid->window_shift[d]);
      abort();
    }
  }

  const int sizex = upper_corner[2] - lower_corner[2];
  const int sizey = upper_corner[1] - lower_corner[1];
  const int sizez = upper_corner[0] - lower_corner[0];

  for (int z = 0; z < sizez; z++) {
    double *restrict dst =
        &idx3(grid[0], lower_corner[0] + z, lower_corner[1], lower_corner[2]);
    double *restrict src =
        &idx3(subgrid[0], position1[0] + z, position1[1], position1[2]);
    for (int y = 0; y < sizey - 1; y++) {
      //__builtin_prefetch(dst + grid->ld_);
#ifdef __LIBXSMM
      LIBXSMM_PRAGMA_SIMD
#else
#pragma GCC ivdep
#endif
      for (int x = 0; x < sizex; x++) {
        dst[x] += src[x];
      }

      dst += grid->ld_;
      src += subgrid->ld_;
    }

#pragma GCC ivdep
    for (int x = 0; x < sizex; x++) {
      dst[x] += src[x];
    }
  }
  return;
}

int compute_cube_properties(const bool ortho, const double radius,
                            const double dh[3][3], const double dh_inv[3][3],
                            const double *rp, double *disr_radius,
                            double *roffset, int *cubecenter, int *lb_cube,
                            int *ub_cube, int *cube_size) {
  int cmax = 0;

  /* center of the gaussian in the lattice coordinates */
  double rp1[3];

  /* it is in the lattice vector frame */
  for (int i = 0; i < 3; i++) {
    double dh_inv_rp = 0.0;
    for (int j = 0; j < 3; j++) {
      dh_inv_rp += dh_inv[j][i] * rp[j];
    }
    rp1[2 - i] = dh_inv_rp;
    cubecenter[2 - i] = floor(dh_inv_rp);
  }

  if (ortho) {
    /* seting up the cube parameters */
    const double dx[3] = {dh[2][2], dh[1][1], dh[0][0]};
    const double dx_inv[3] = {dh_inv[2][2], dh_inv[1][1], dh_inv[0][0]};
    /* cube center */

    /* lower and upper bounds */

    // Historically, the radius gets discretized.
    const double drmin = fmin(dh[0][0], fmin(dh[1][1], dh[2][2]));
    *disr_radius = drmin * fmax(1.0, ceil(radius / drmin));

    for (int i = 0; i < 3; i++) {
      roffset[i] = rp[2 - i] - ((double)cubecenter[i]) * dx[i];
    }

    for (int i = 0; i < 3; i++) {
      lb_cube[i] = ceil(-1e-8 - *disr_radius * dx_inv[i]);
    }

    // Symmetric interval
    for (int i = 0; i < 3; i++) {
      ub_cube[i] = 1 - lb_cube[i];
    }

  } else {
    for (int idir = 0; idir < 3; idir++) {
      lb_cube[idir] = INT_MAX;
      ub_cube[idir] = INT_MIN;
    }

    /* compute the size of the box. It is a fairly trivial way to compute
     * the box and it may have far more point than needed */
    for (int i = -1; i <= 1; i++) {
      for (int j = -1; j <= 1; j++) {
        for (int k = -1; k <= 1; k++) {
          double x[3] = {/* rp[0] + */ ((double)i) * radius,
                         /* rp[1] + */ ((double)j) * radius,
                         /* rp[2] + */ ((double)k) * radius};
          /* convert_to_lattice_coordinates(dh_inv, x, y); */
          for (int idir = 0; idir < 3; idir++) {
            const double resc = dh_inv[0][idir] * x[0] +
                                dh_inv[1][idir] * x[1] + dh_inv[2][idir] * x[2];
            lb_cube[2 - idir] = imin(lb_cube[2 - idir], floor(resc));
            ub_cube[2 - idir] = imax(ub_cube[2 - idir], ceil(resc));
          }
        }
      }
    }

    /* compute the offset in lattice coordinates */

    for (int i = 0; i < 3; i++) {
      roffset[i] = rp1[i] - cubecenter[i];
    }

    *disr_radius = radius;
  }

  /* compute the cube size ignoring periodicity */

  /* the +1 is normal here */
  cube_size[0] = ub_cube[0] - lb_cube[0] + 1;
  cube_size[1] = ub_cube[1] - lb_cube[1] + 1;
  cube_size[2] = ub_cube[2] - lb_cube[2] + 1;

  for (int i = 0; i < 3; i++) {
    cmax = imax(cmax, cube_size[i]);
  }

  return cmax;
}

void return_cube_position(const int *const grid_size, const int *const lb_grid,
                          const int *const cube_center,
                          const int *const lower_boundaries_cube,
                          const int *const period, int *const position) {
  for (int i = 0; i < 3; i++)
    position[i] = modulo(cube_center[i] - lb_grid[i] + lower_boundaries_cube[i],
                         period[i]);

  assert(!(position[0] >= grid_size[0]) || (position[1] >= grid_size[1]) ||
         (position[2] >= grid_size[2]));
}

void verify_orthogonality(const double dh[3][3], bool orthogonal[3]) {
  double norm1, norm2, norm3;

  norm1 = dh[0][0] * dh[0][0] + dh[0][1] * dh[0][1] + dh[0][2] * dh[0][2];
  norm2 = dh[1][0] * dh[1][0] + dh[1][1] * dh[1][1] + dh[1][2] * dh[1][2];
  norm3 = dh[2][0] * dh[2][0] + dh[2][1] * dh[2][1] + dh[2][2] * dh[2][2];

  norm1 = 1.0 / sqrt(norm1);
  norm2 = 1.0 / sqrt(norm2);
  norm3 = 1.0 / sqrt(norm3);

  /* x z */
  orthogonal[0] =
      ((fabs(dh[0][0] * dh[2][0] + dh[0][1] * dh[2][1] + dh[0][2] * dh[2][2]) *
        norm1 * norm3) < 1e-12);
  /* y z */
  orthogonal[1] =
      ((fabs(dh[1][0] * dh[2][0] + dh[1][1] * dh[2][1] + dh[1][2] * dh[2][2]) *
        norm2 * norm3) < 1e-12);
  /* x y */
  orthogonal[2] =
      ((fabs(dh[0][0] * dh[1][0] + dh[0][1] * dh[1][1] + dh[0][2] * dh[1][2]) *
        norm1 * norm2) < 1e-12);
}

#ifndef __MKL
void cblas_daxpy(const int N, const double alpha, const double *X,
                 const int incX, double *Y, const int incY) {
  if ((incX == 1) && (incY == 1)) {
#pragma GCC ivdep
    for (int i = 0; i < N; i++)
      Y[i] += alpha * X[i];
    return;
  }

  if (incX == 1) {
#pragma GCC ivdep
    for (int i = 0; i < N; i++)
      Y[i + incY] += alpha * X[i];
    return;
  }

  if (incY == 1) {
#pragma GCC ivdep
    for (int i = 0; i < N; i++)
      Y[i] += alpha * X[i + incX];
    return;
  }

#pragma GCC ivdep
  for (int i = 0; i < N; i++)
    Y[i + incY] += alpha * X[i + incX];
  return;
}

double cblas_ddot(const int N, const double *X, const int incX, const double *Y,
                  const int incY) {
  if ((incX == incY) && (incY == 1)) {
    double res = 0.0;

#pragma GCC ivdep
    for (int i = 0; i < N; i++) {
      res += X[i] * Y[i];
    }
    return res;
  }

  if (incX == 1) {
    double res = 0.0;

#pragma GCC ivdep
    for (int i = 0; i < N; i++) {
      res += X[i] * Y[i + incY];
    }
    return res;
  }

  if (incY == 1) {
    double res = 0.0;

#pragma GCC ivdep
    for (int i = 0; i < N; i++) {
      res += X[i + incX] * Y[i];
    }
    return res;
  }

  double res = 0.0;

#pragma GCC ivdep
  for (int i = 0; i < N; i++) {
    res += X[i + incX] * Y[i + incY];
  }
  return res;
}

void cblas_dger(const CBLAS_LAYOUT Layout, const int M, const int N,
                const double alpha, const double *X, const int incX,
                const double *Y, const int incY, double *A, const int lda) {
  if (Layout == CblasRowMajor) {
    for (int i = 0; i < M; i++) {
      const double x = alpha * X[i + incX];
      if (incY == 1) {
        double *restrict dst = &A[i * lda];
        const double *restrict y = Y;
        for (int k = 0; k < N; k++) {
          dst[k] += y[k] * x;
        }
      } else {
        double *restrict dst = &A[i * lda];
        const double *restrict y = Y;
        for (int k = 0; k < N; k++) {
          dst[k] += y[k + incY] * x;
        }
      }
    }
  } else {
    for (int j = 0; j < N; j++) {
      const double y = alpha * Y[j + incY];
      if (incX == 1) {
        double *restrict dst = &A[j * lda];
        const double *restrict x = X;
        for (int k = 0; k < M; k++) {
          dst[k] += x[k] * y;
        }
      } else {
        double *restrict dst = &A[j * lda];
        const double *restrict x = X;
        for (int k = 0; k < M; k++) {
          dst[k] += x[k + incX] * y;
        }
      }
    }
  }
}

/* code taken from gsl_cblas. We really need to use a proper cblas interface and
 * build system.... */
void cblas_dgemv(const CBLAS_LAYOUT order, const CBLAS_TRANSPOSE TransA,
                 const int M, const int N, const double alpha, const double *A,
                 const int lda, const double *X, const int incX,
                 const double beta, double *Y, const int incY) {
  int i, j;
  int lenX, lenY;
#define OFFSET(N, incX) ((incX) > 0 ? 0 : ((N)-1) * (-(incX)))

  const int Trans = (TransA != CblasConjTrans) ? TransA : CblasTrans;

  if (M == 0 || N == 0)
    return;

  if (alpha == 0.0 && beta == 1.0)
    return;

  if (Trans == CblasNoTrans) {
    lenX = N;
    lenY = M;
  } else {
    lenX = M;
    lenY = N;
  }

  /* form  y := beta*y */
  if (beta == 0.0) {
    int iy = OFFSET(lenY, incY);
    for (i = 0; i < lenY; i++) {
      Y[iy] = 0.0;
      iy += incY;
    }
  } else if (beta != 1.0) {
    int iy = OFFSET(lenY, incY);
    for (i = 0; i < lenY; i++) {
      Y[iy] *= beta;
      iy += incY;
    }
  }

  if (alpha == 0.0)
    return;

  if ((order == CblasRowMajor && Trans == CblasNoTrans) ||
      (order == CblasColMajor && Trans == CblasTrans)) {
    /* form  y := alpha*A*x + y */
    int iy = OFFSET(lenY, incY);
    for (i = 0; i < lenY; i++) {
      double temp = 0.0;
      int ix = OFFSET(lenX, incX);
      for (j = 0; j < lenX; j++) {
        temp += X[ix] * A[lda * i + j];
        ix += incX;
      }
      Y[iy] += alpha * temp;
      iy += incY;
    }
  } else if ((order == CblasRowMajor && Trans == CblasTrans) ||
             (order == CblasColMajor && Trans == CblasNoTrans)) {
    /* form  y := alpha*A'*x + y */
    int ix = OFFSET(lenX, incX);
    for (j = 0; j < lenX; j++) {
      const double temp = alpha * X[ix];
      if (temp != 0.0) {
        int iy = OFFSET(lenY, incY);
        for (i = 0; i < lenY; i++) {
          Y[iy] += temp * A[lda * j + i];
          iy += incY;
        }
      }
      ix += incX;
    }
  }
#undef OFFSET
}
#endif

int compute_interval(const int *const map, const int full_size, const int size,
                     const int cube_size, const int x1, int *x,
                     int *const lower_corner, int *const upper_corner,
                     const Interval window) {
  if (size == full_size) {
    /* we have the full grid in that direction */
    /* lower boundary is within the window */
    *lower_corner = x1;
    /* now compute the upper corner */
    /* needs to be as large as possible. basically I take [x1..
     * min(grid.full_size, cube_size - x)] */

    *upper_corner = compute_next_boundaries(x1, *x, full_size, cube_size);

    {
      Interval tz = create_interval(*lower_corner, *upper_corner);
      Interval res = intersection_interval(tz, window);
      *lower_corner = res.xmin;
      *upper_corner = res.xmax;
    }
  } else {
    *lower_corner = x1;
    *upper_corner = x1 + 1;

    // the map is always increasing by 1 except when we cross the boundaries of
    // the grid and pbc are applied. Since we are only interested in by a
    // subwindow of the full table we check that the next point is inside the
    // window of interest and is also equal to the previous point + 1. The last
    // check is pointless in practice.

    for (int i = *x + 2; (i < size) && ((*upper_corner + 1) == map[i]) &&
                         is_point_in_interval(map[i], window);
         i++) {
      (*upper_corner)++;
    }
    return (*upper_corner) - (*lower_corner);
  }

  return imax(*lower_corner - *upper_corner, 0);
}
