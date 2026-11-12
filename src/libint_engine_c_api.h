/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#ifndef LIBINT_ENGINE_C_API_H
#define LIBINT_ENGINE_C_API_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Stable bridge identifiers, independent of Libint's C++ enum values.
 * parameters[0] is omega for erf/erfc/erfx, the radius for truncated Coulomb,
 * or alpha in exp(-alpha*r12^2) for Gaussian. For erfx, parameters[1:2] are
 * the erf and erfc coefficients. All three entries must be finite; the
 * length/attenuation parameter must be nonnegative. Other entries are ignored.
 */
enum cp2k_libint_engine_operator {
  CP2K_LIBINT_ENGINE_COULOMB = 1,
  CP2K_LIBINT_ENGINE_ERF_COULOMB = 2,
  CP2K_LIBINT_ENGINE_ERFC_COULOMB = 3,
  CP2K_LIBINT_ENGINE_TRUNCATED_COULOMB = 4,
  CP2K_LIBINT_ENGINE_GAUSSIAN = 5,
  CP2K_LIBINT_ENGINE_DELTA = 6,
  CP2K_LIBINT_ENGINE_ERFX_COULOMB = 7
};

/* Status functions return 0 on success and 1 on error. C++ exceptions are
 * caught here. A non-NULL error buffer is NUL-terminated if error_capacity > 0.
 * A NULL error buffer is allowed. Initialization/finalization belongs to the
 * process lifecycle, outside threaded regions. Destroy all handles before
 * finalize.
 */
int cp2k_libint_engine_initialize(char *error, size_t error_capacity);
int cp2k_libint_engine_finalize(char *error, size_t error_capacity);

/* Creates one owning handle, initially NULL, for use by one thread at a time.
 * ncenters = 2, 3, or 4 selects (A|B), (A|BC), or (AB|CD).
 * derivative_order = 0 or 1. Capacities max_nprim > 0 and max_l >= 0 are fixed
 * for the handle's lifetime. precision >= 0 is Engine's absolute primitive
 * screening precision; zero disables this screening.
 */
int cp2k_libint_engine_create(int operator_id, int ncenters,
                              int derivative_order, int max_nprim, int max_l,
                              const double parameters[3], double precision,
                              void **handle, char *error,
                              size_t error_capacity);

/* Sets one 1-based physical center slot. Each shell has one Cartesian
 * contraction. Exponents are strictly non-negative. Coefficients multiply
 * normalization-free Cartesian primitives; no primitive or contraction
 * normalization is applied. Input arrays are copied into reusable storage.
 * All physical slots must be set before compute. Dummy centers are internal.
 * Libint must supply Cartesian kernels for the requested rank and derivative
 * order.
 */
int cp2k_libint_engine_set_shell(void *handle, int slot, int l, int nprim,
                                 const double *exponents,
                                 const double *coefficients,
                                 const double center[3], char *error,
                                 size_t error_capacity);

/* Returns borrowed, read-only component pointers in Engine's native layout:
 * row-major, last physical shell index fastest; within each shell, the
 * Cartesian ordering configured in Libint. ncart is the product of shell
 * sizes. ncomponents is 1 for values or 3*ncenters for center derivatives
 * ordered Ax,Ay,Az,Bx,..., without dummy-center components.
 * Successful screened-zero shell sets have *targets == NULL and valid sizes.
 * A NULL component pointer also denotes zeros. On failure all outputs are zero.
 * Every set_shell/compute/destroy call on this handle invalidates earlier
 * borrowed results, including calls which fail. Copy results before reuse.
 */
int cp2k_libint_engine_compute(void *handle, const double *const **targets,
                               size_t *ncart, int *ncomponents, char *error,
                               size_t error_capacity);

/* Releases an owning handle and sets it to NULL. A NULL handle is allowed. */
void cp2k_libint_engine_destroy(void **handle);

#ifdef __cplusplus
}
#endif

#endif
