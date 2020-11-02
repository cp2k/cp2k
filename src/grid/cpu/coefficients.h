#ifndef COEFFICIENTS_H
#define COEFFICIENTS_H
#include <stdbool.h>
#include <string.h>

#include "../common/grid_common.h"
#include "utils.h"

// *****************************************************************************
extern void grid_prepare_alpha_dgemm(const double ra[3], const double rb[3],
                                     const double rp[3], const int *lmax,
                                     tensor *alpha);

extern void grid_prepare_coef_dgemm(const int *const lmin,
                                    const int *const lmax, const int lp,
                                    const double prefactor,
                                    const tensor *const alpha,
                                    const tensor *const pab, tensor *coef_xyz);

extern void compute_compact_polynomial_coefficients(
    const tensor *coef, const int *coef_offset_, const int *lmin,
    const int *lmax, const double *ra, const double *rb, const double *rab,
    const double prefactor, tensor *co);

extern void grid_transform_coef_xyz_to_ijk(const double dh[3][3],
                                           const tensor *coef_xyz);

extern void grid_transform_coef_jik_to_yxz(const double dh[3][3],
                                           const tensor *coef_xyz);
extern void transform_triangular_to_xyz(const double *const coef_xyz,
                                        tensor *const coef);
extern void transform_xyz_to_triangular(const tensor *const coef,
                                        double *const coef_xyz);
extern void transform_yxz_to_triangular(const tensor *const coef,
                                        double *const coef_xyz);
extern void grid_transform_coef_ijk_to_xyz_cp2k(
    const int lp, const double dh[3][3], const double *__restrict coef_ijk,
    double *__restrict coef_xyz);

#endif
