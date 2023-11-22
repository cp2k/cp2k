/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2023 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#if defined(__DLAF)

#include <dlaf_c/eigensolver/eigensolver.h>
#include <dlaf_c/factorization/cholesky.h>
#include <dlaf_c/init.h>

void dlaf_init() {
  const char *pika_argv[] = {"cp2k"};
  const char *dlaf_argv[] = {"cp2k"};
  dlaf_initialize(1, pika_argv, 1, dlaf_argv);
}

// Eigensolver (double)
void dlaf_pdsyevd_wrapper(int n, double *a, int ia, int ja, int desca[9],
                          double *w, double *z, int iz, int jz, int descz[9],
                          int *info) {
  char uplo = 'L'; // Only lower triangular matrices are supported
  dlaf_pdsyevd(uplo, n, a, ia, ja, desca, w, z, iz, jz, descz, info);
}

// Cholesky decomposition (double)
// Wrapper without uplo parameter, avoids passing characters
void dlaf_pdpotrf_wrapper(int n, double *a, int ia, int ja, int desca[9],
                          int *info) {
  dlaf_pdpotrf('U', n, a, ia, ja, desca, info);
}

// Cholesky decomposition (float)
// Wrapper without uplo parameter, avoids passing characters
void dlaf_pspotrf_wrapper(int n, float *a, int ia, int ja, int desca[9],
                          int *info) {
  dlaf_pspotrf('U', n, a, ia, ja, desca, info);
}

#endif // __DLAF
