/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#ifndef UTILS_H
#define UTILS_H

#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#if defined(__MKL)
#include <mkl.h>
#include <mkl_cblas.h>
#endif

#if defined(__LIBXSMM)
#include <libxsmm.h>
#endif

#include "../common/grid_common.h"
#include "private_header.h"
#include "tensor_local.h"

/* inverse of the factorials */
static const double inv_fac[] = {1.0,
                                 1.0,
                                 0.5,
                                 0.166666666666666666666666666667,
                                 0.0416666666666666666666666666667,
                                 0.00833333333333333333333333333333,
                                 0.00138888888888888888888888888889,
                                 0.000198412698412698412698412698413,
                                 0.0000248015873015873015873015873016,
                                 2.7557319223985890652557319224e-6,
                                 2.7557319223985890652557319224e-7,
                                 2.50521083854417187750521083854e-8,
                                 2.08767569878680989792100903212e-9,
                                 1.60590438368216145993923771702e-10,
                                 1.14707455977297247138516979787e-11,
                                 7.64716373181981647590113198579e-13,
                                 4.77947733238738529743820749112e-14,
                                 2.81145725434552076319894558301e-15,
                                 1.56192069685862264622163643501e-16,
                                 8.22063524662432971695598123687e-18,
                                 4.11031762331216485847799061844e-19,
                                 1.95729410633912612308475743735e-20,
                                 8.8967913924505732867488974425e-22,
                                 3.86817017063068403771691193152e-23,
                                 1.6117375710961183490487133048e-24,
                                 6.4469502843844733961948532192e-26,
                                 2.47959626322479746007494354585e-27,
                                 9.18368986379554614842571683647e-29,
                                 3.27988923706983791015204172731e-30,
                                 1.13099628864477169315587645769e-31,
                                 3.76998762881590564385292152565e-33};

inline int coset_without_offset(int lx, int ly, int lz) {
  const int l = lx + ly + lz;
  if (l == 0) {
    return 0;
  } else {
    return ((l - lx) * (l - lx + 1)) / 2 + lz;
  }
}

typedef struct dgemm_params_ {
  char storage;
  char op1;
  char op2;
  double alpha;
  double beta;
  double *a, *b, *c;
  int m, n, k, lda, ldb, ldc;
  int x, y, z;
  int x1, y1, z1;
  bool use_libxsmm;
#if defined(__LIBXSMM)
  libxsmm_dmmfunction kernel;
  int prefetch;
  int flags;
#endif
} dgemm_params;

extern void dgemm_simplified(dgemm_params *const m);
extern void batched_dgemm_simplified(dgemm_params *const m,
                                     const int batch_size);

/*******************************************************************************
 * \brief Prototype for BLAS dgemm.
 * \author Ole Schuett
 ******************************************************************************/
void dgemm_(const char *transa, const char *transb, const int *m, const int *n,
            const int *k, const double *alpha, const double *a, const int *lda,
            const double *b, const int *ldb, const double *beta, double *c,
            const int *ldc);

extern void extract_sub_grid(const int *lower_corner, const int *upper_corner,
                             const int *position, const tensor *const grid,
                             tensor *const subgrid);
extern void add_sub_grid(const int *lower_corner, const int *upper_corner,
                         const int *position, const tensor *subgrid,
                         tensor *grid);
extern void return_cube_position(const int *grid_size, const int *lb_grid,
                                 const int *cube_center,
                                 const int *lower_boundaries_cube,
                                 const int *period, int *const position);

extern void verify_orthogonality(const double dh[3][3], bool *orthogonal);

extern int compute_cube_properties(const bool ortho, const double radius,
                                   const double dh[3][3],
                                   const double dh_inv[3][3], const double *rp,
                                   double *disr_radius, double *roffset,
                                   int *cubecenter, int *lb_cube, int *ub_cube,
                                   int *cube_size);

inline int return_offset_l(const int l) {
  static const int offset_[] = {1,   4,   7,   11,  16,  22,  29,
                                37,  46,  56,  67,  79,  92,  106,
                                121, 137, 154, 172, 191, 211, 232};
  return offset_[l];
}

inline int return_linear_index_from_exponents(const int alpha, const int beta,
                                              const int gamma) {
  const int l = alpha + beta + gamma;
  return return_offset_l(l) + (l - alpha) * (l - alpha + 1) / 2 + gamma;
}

static inline void *grid_allocate_scratch(size_t size) {
#ifdef __LIBXSMM
  return libxsmm_aligned_scratch(size, 0 /*auto-alignment*/);
#else
  return malloc(size);
#endif
}

static inline void grid_free_scratch(void *ptr) {
#ifdef __LIBXSMM
  libxsmm_free(ptr);
#else
  free(ptr);
#endif
}

/* even openblas has cblas versions of lapack and blas. */
#ifndef __MKL
enum CBLAS_LAYOUT { CblasRowMajor = 101, CblasColMajor = 102 };
enum CBLAS_TRANSPOSE {
  CblasNoTrans = 111,
  CblasTrans = 112,
  CblasConjTrans = 113
};
enum CBLAS_UPLO { CblasUpper = 121, CblasLower = 122 };
enum CBLAS_DIAG { CblasNonUnit = 131, CblasUnit = 132 };
enum CBLAS_SIDE { CblasLeft = 141, CblasRight = 142 };

typedef enum CBLAS_LAYOUT CBLAS_LAYOUT;
typedef enum CBLAS_TRANSPOSE CBLAS_TRANSPOSE;
typedef enum CBLAS_UPLO CBLAS_UPLO;
typedef enum CBLAS_DIAG CBLAS_DIAG;

double cblas_ddot(const int N, const double *X, const int incX, const double *Y,
                  const int incY);

void cblas_dger(const CBLAS_LAYOUT Layout, const int M, const int N,
                const double alpha, const double *X, const int incX,
                const double *Y, const int incY, double *A, const int lda);

void cblas_daxpy(const int N, const double alpha, const double *X,
                 const int incX, double *Y, const int incY);

void cblas_dgemv(const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE TransA,
                 const int M, const int N, const double alpha, const double *A,
                 const int lda, const double *X, const int incX,
                 const double beta, double *Y, const int incY);

#endif

extern void compute_interval(const int *const map, const int full_size,
                             const int size, const int cube_size, const int x1,
                             int *x, int *const lower_corner,
                             int *const upper_corner, Interval window);
#endif
