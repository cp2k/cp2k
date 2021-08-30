/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2021 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(__AVX2__) && defined(__FMA__)
#include <immintrin.h>
#endif

#include "../common/grid_common.h"
#include "../common/grid_library.h"
#include "../common/grid_sphere_cache.h"

#define GRID_MAX_LP_OPTIMIZED 9

#if (GRID_DO_COLLOCATE)
#define GRID_CONST_WHEN_COLLOCATE const
#define GRID_CONST_WHEN_INTEGRATE
#else
#define GRID_CONST_WHEN_COLLOCATE
#define GRID_CONST_WHEN_INTEGRATE const
#endif

/*******************************************************************************
 * \brief Simple loop body for ortho_cx_to_grid using plain C.
 * \author Ole Schuett
 ******************************************************************************/
static inline void __attribute__((always_inline))
ortho_cx_to_grid_scalar(const int lp, const int cmax, const int i,
                        const double pol[3][lp + 1][2 * cmax + 1],
                        GRID_CONST_WHEN_COLLOCATE double *cx,
                        GRID_CONST_WHEN_INTEGRATE double *grid_0,
                        GRID_CONST_WHEN_INTEGRATE double *grid_1,
                        GRID_CONST_WHEN_INTEGRATE double *grid_2,
                        GRID_CONST_WHEN_INTEGRATE double *grid_3) {

#if (GRID_DO_COLLOCATE)
  // collocate
  double reg[4] = {0.0, 0.0, 0.0, 0.0};
  for (int lxp = 0; lxp <= lp; lxp++) {
    const double p = pol[0][lxp][i + cmax];
    reg[0] += cx[lxp * 4 + 0] * p;
    reg[1] += cx[lxp * 4 + 1] * p;
    reg[2] += cx[lxp * 4 + 2] * p;
    reg[3] += cx[lxp * 4 + 3] * p;
  }
  *grid_0 += reg[0];
  *grid_1 += reg[1];
  *grid_2 += reg[2];
  *grid_3 += reg[3];

#else
  // integrate
  const double reg[4] = {*grid_0, *grid_1, *grid_2, *grid_3};
  for (int lxp = 0; lxp <= lp; lxp++) {
    const double p = pol[0][lxp][i + cmax];
    cx[lxp * 4 + 0] += reg[0] * p;
    cx[lxp * 4 + 1] += reg[1] * p;
    cx[lxp * 4 + 2] += reg[2] * p;
    cx[lxp * 4 + 3] += reg[3] * p;
  }
#endif
}

/*******************************************************************************
 * \brief Optimized loop body for ortho_cx_to_grid using AVX2 Intel Intrinsics.
 *        This routine always processes four consecutive grid elements at once.
 * \author Ole Schuett
 ******************************************************************************/
#if defined(__AVX2__) && defined(__FMA__)
static inline void __attribute__((always_inline))
ortho_cx_to_grid_avx2(const int lp, const int cmax, const int i,
                      const double pol[3][lp + 1][2 * cmax + 1],
                      GRID_CONST_WHEN_COLLOCATE double *cx,
                      GRID_CONST_WHEN_INTEGRATE double *grid_0,
                      GRID_CONST_WHEN_INTEGRATE double *grid_1,
                      GRID_CONST_WHEN_INTEGRATE double *grid_2,
                      GRID_CONST_WHEN_INTEGRATE double *grid_3) {

  const int icmax = i + cmax;

#if (GRID_DO_COLLOCATE)
  // collocate
  // First iteration for lxp == 0 does not need add instructions.
  __m256d p_vec = _mm256_loadu_pd(&pol[0][0][icmax]);
  __m256d r_vec_0 = _mm256_mul_pd(p_vec, _mm256_set1_pd(cx[0]));
  __m256d r_vec_1 = _mm256_mul_pd(p_vec, _mm256_set1_pd(cx[1]));
  __m256d r_vec_2 = _mm256_mul_pd(p_vec, _mm256_set1_pd(cx[2]));
  __m256d r_vec_3 = _mm256_mul_pd(p_vec, _mm256_set1_pd(cx[3]));

  // Remaining iterations for lxp > 0 use fused multiply adds.
  GRID_PRAGMA_UNROLL_UP_TO(GRID_MAX_LP_OPTIMIZED)
  for (int lxp = 1; lxp <= lp; lxp++) {
    const double *cx_base = &cx[lxp * 4];
    p_vec = _mm256_loadu_pd(&pol[0][lxp][icmax]);
    r_vec_0 = _mm256_fmadd_pd(p_vec, _mm256_set1_pd(cx_base[0]), r_vec_0);
    r_vec_1 = _mm256_fmadd_pd(p_vec, _mm256_set1_pd(cx_base[1]), r_vec_1);
    r_vec_2 = _mm256_fmadd_pd(p_vec, _mm256_set1_pd(cx_base[2]), r_vec_2);
    r_vec_3 = _mm256_fmadd_pd(p_vec, _mm256_set1_pd(cx_base[3]), r_vec_3);
  }

  // Add vectors to grid one at a time, because they can aliase when cube wraps.
  _mm256_storeu_pd(grid_0, _mm256_add_pd(_mm256_loadu_pd(grid_0), r_vec_0));
  _mm256_storeu_pd(grid_1, _mm256_add_pd(_mm256_loadu_pd(grid_1), r_vec_1));
  _mm256_storeu_pd(grid_2, _mm256_add_pd(_mm256_loadu_pd(grid_2), r_vec_2));
  _mm256_storeu_pd(grid_3, _mm256_add_pd(_mm256_loadu_pd(grid_3), r_vec_3));

#else
  // integrate
  __m256d grid_vec_0 = _mm256_loadu_pd(grid_0);
  __m256d grid_vec_1 = _mm256_loadu_pd(grid_1);
  __m256d grid_vec_2 = _mm256_loadu_pd(grid_2);
  __m256d grid_vec_3 = _mm256_loadu_pd(grid_3);

  GRID_PRAGMA_UNROLL_UP_TO(GRID_MAX_LP_OPTIMIZED + 1)
  for (int lxp = 0; lxp <= lp; lxp++) {
    __m256d p_vec = _mm256_loadu_pd(&pol[0][lxp][icmax]);

    // Do 4 dot products at once. https://stackoverflow.com/a/10454420
    __m256d xy0 = _mm256_mul_pd(p_vec, grid_vec_0);
    __m256d xy1 = _mm256_mul_pd(p_vec, grid_vec_1);
    __m256d xy2 = _mm256_mul_pd(p_vec, grid_vec_2);
    __m256d xy3 = _mm256_mul_pd(p_vec, grid_vec_3);

    // low to high: xy00+xy01 xy10+xy11 xy02+xy03 xy12+xy13
    __m256d temp01 = _mm256_hadd_pd(xy0, xy1);

    // low to high: xy20+xy21 xy30+xy31 xy22+xy23 xy32+xy33
    __m256d temp23 = _mm256_hadd_pd(xy2, xy3);

    // low to high: xy02+xy03 xy12+xy13 xy20+xy21 xy30+xy31
    __m256d swapped = _mm256_permute2f128_pd(temp01, temp23, 0x21);

    // low to high: xy00+xy01 xy10+xy11 xy22+xy23 xy32+xy33
    __m256d blended = _mm256_blend_pd(temp01, temp23, 0b1100);

    __m256d r_vec = _mm256_add_pd(swapped, blended);

    // cx += r_vec
    double *cx_base = &cx[lxp * 4];
    _mm256_storeu_pd(cx_base, _mm256_add_pd(r_vec, _mm256_loadu_pd(cx_base)));
  }
#endif
}
#endif // __AVX2__ && __FMA__

/*******************************************************************************
 * \brief Collocates coefficients C_x onto the grid for orthorhombic case.
 * \author Ole Schuett
 ******************************************************************************/
static inline void __attribute__((always_inline))
ortho_cx_to_grid(const int lp, const int kg1, const int kg2, const int jg1,
                 const int jg2, const int cmax,
                 const double pol[3][lp + 1][2 * cmax + 1],
                 const int map[3][2 * cmax + 1],
                 const int sections[3][2 * cmax + 1], const int npts_local[3],
                 int **sphere_bounds_iter, GRID_CONST_WHEN_COLLOCATE double *cx,
                 GRID_CONST_WHEN_INTEGRATE double *grid) {

  // Lower and upper sphere bounds relative to center, ie. in cube coordinates.
  const int lb = *((*sphere_bounds_iter)++);
  const int ub = 1 - lb;

  // AVX instructions can only load/store from evenly spaced memory locations.
  // Since the sphere bounds can wrap around due to the grid's periodicity,
  // the inner loop runs over sections with homogeneous cube to grid mapping.
  for (int istart = lb; istart <= ub; istart++) {
    const int istop = imin(ub, istart + sections[0][istart + cmax]);
    const int cube2grid = map[0][istart + cmax] - istart;

    const int stride = npts_local[1] * npts_local[0];
    const int grid_index_0 = kg1 * stride + jg1 * npts_local[0];
    const int grid_index_1 = kg2 * stride + jg1 * npts_local[0];
    const int grid_index_2 = kg1 * stride + jg2 * npts_local[0];
    const int grid_index_3 = kg2 * stride + jg2 * npts_local[0];
    GRID_CONST_WHEN_INTEGRATE double *grid_base_0 = &grid[grid_index_0];
    GRID_CONST_WHEN_INTEGRATE double *grid_base_1 = &grid[grid_index_1];
    GRID_CONST_WHEN_INTEGRATE double *grid_base_2 = &grid[grid_index_2];
    GRID_CONST_WHEN_INTEGRATE double *grid_base_3 = &grid[grid_index_3];

    // Use AVX2 to process grid points in chunks of four, ie. 256 bit vectors.
#if defined(__AVX2__) && defined(__FMA__)
    const int istop_vec = istart + 4 * ((istop - istart + 1) / 4) - 1;
    for (int i = istart; i <= istop_vec; i += 4) {
      const int ig = i + cube2grid;
      ortho_cx_to_grid_avx2(lp, cmax, i, pol, cx, &grid_base_0[ig],
                            &grid_base_1[ig], &grid_base_2[ig],
                            &grid_base_3[ig]);
    }
    istart = istop_vec + 1;
#endif

    // Process up to 3 remaining points - or everything if AVX2 isn't available.
    for (int i = istart; i <= istop; i++) {
      const int ig = i + cube2grid;
      ortho_cx_to_grid_scalar(lp, cmax, i, pol, cx, &grid_base_0[ig],
                              &grid_base_1[ig], &grid_base_2[ig],
                              &grid_base_3[ig]);
    }
    istart = istop;
  }
}

/*******************************************************************************
 * \brief Transforms coefficients C_xy into C_x by fixing grid index j.
 * \author Ole Schuett
 ******************************************************************************/
static inline void __attribute__((always_inline))
ortho_cxy_to_cx(const int lp, const int j1, const int j2, const int cmax,
                const double pol[3][lp + 1][2 * cmax + 1],
                GRID_CONST_WHEN_COLLOCATE double *cxy,
                GRID_CONST_WHEN_INTEGRATE double *cx) {

  for (int lyp = 0; lyp <= lp; lyp++) {
    for (int lxp = 0; lxp <= lp - lyp; lxp++) {
      const double p1 = pol[1][lyp][j1 + cmax];
      const double p2 = pol[1][lyp][j2 + cmax];
      const int cxy_index = lyp * (lp + 1) * 2 + lxp * 2; // [lyp, lxp, 0]

#if (GRID_DO_COLLOCATE)
      // collocate
      cx[lxp * 4 + 0] += cxy[cxy_index + 0] * p1;
      cx[lxp * 4 + 1] += cxy[cxy_index + 1] * p1;
      cx[lxp * 4 + 2] += cxy[cxy_index + 0] * p2;
      cx[lxp * 4 + 3] += cxy[cxy_index + 1] * p2;
#else
      // integrate
      cxy[cxy_index + 0] += cx[lxp * 4 + 0] * p1;
      cxy[cxy_index + 1] += cx[lxp * 4 + 1] * p1;
      cxy[cxy_index + 0] += cx[lxp * 4 + 2] * p2;
      cxy[cxy_index + 1] += cx[lxp * 4 + 3] * p2;
#endif
    }
  }
}

/*******************************************************************************
 * \brief Loop body of ortho_cxy_to_grid to be inlined for low values of lp.
 * \author Ole Schuett
 ******************************************************************************/
static inline void __attribute__((always_inline))
ortho_cxy_to_grid_low(const int lp, const int j1, const int j2, const int kg1,
                      const int kg2, const int jg1, const int jg2,
                      const int cmax, const double pol[3][lp + 1][2 * cmax + 1],
                      const int map[3][2 * cmax + 1],
                      const int sections[3][2 * cmax + 1],
                      const int npts_local[3], int **sphere_bounds_iter,
                      double *cx, GRID_CONST_WHEN_COLLOCATE double *cxy,
                      GRID_CONST_WHEN_INTEGRATE double *grid) {

#if (GRID_DO_COLLOCATE)
  // collocate
  ortho_cxy_to_cx(lp, j1, j2, cmax, pol, cxy, cx);
  ortho_cx_to_grid(lp, kg1, kg2, jg1, jg2, cmax, pol, map, sections, npts_local,
                   sphere_bounds_iter, cx, grid);
#else
  // integrate
  ortho_cx_to_grid(lp, kg1, kg2, jg1, jg2, cmax, pol, map, sections, npts_local,
                   sphere_bounds_iter, cx, grid);
  ortho_cxy_to_cx(lp, j1, j2, cmax, pol, cxy, cx);
#endif
}

/*******************************************************************************
 * \brief Collocates coefficients C_xy onto the grid for orthorhombic case.
 * \author Ole Schuett
 ******************************************************************************/
static inline void ortho_cxy_to_grid(
    const int lp, const int kg1, const int kg2, const int cmax,
    const double pol[3][lp + 1][2 * cmax + 1], const int map[3][2 * cmax + 1],
    const int sections[3][2 * cmax + 1], const int npts_local[3],
    int **sphere_bounds_iter, GRID_CONST_WHEN_COLLOCATE double *cxy,
    GRID_CONST_WHEN_INTEGRATE double *grid) {

  // The cube contains an even number of grid points in each direction and
  // collocation is always performed on a pair of two opposing grid points.
  // Hence, the points with index 0 and 1 are both assigned distance zero via
  // the formular distance=(2*index-1)/2.

  const int jstart = *((*sphere_bounds_iter)++);
  const size_t cx_size = (lp + 1) * 4;
  double cx[cx_size];
  for (int j1 = jstart; j1 <= 0; j1++) {
    const int j2 = 1 - j1;
    const int jg1 = map[1][j1 + cmax];
    const int jg2 = map[1][j2 + cmax];

    memset(cx, 0, cx_size * sizeof(double));

    // Generate separate branches for low values of lp gives up to 30% speedup.
    if (lp <= GRID_MAX_LP_OPTIMIZED) {
      GRID_PRAGMA_UNROLL(GRID_MAX_LP_OPTIMIZED + 1)
      for (int ilp = 0; ilp <= GRID_MAX_LP_OPTIMIZED; ilp++) {
        if (lp == ilp) {
          ortho_cxy_to_grid_low(ilp, j1, j2, kg1, kg2, jg1, jg2, cmax, pol, map,
                                sections, npts_local, sphere_bounds_iter, cx,
                                cxy, grid);
        }
      }
    } else {
      ortho_cxy_to_grid_low(lp, j1, j2, kg1, kg2, jg1, jg2, cmax, pol, map,
                            sections, npts_local, sphere_bounds_iter, cx, cxy,
                            grid);
    }
  }
}

/*******************************************************************************
 * \brief Transforms coefficients C_xyz into C_xz by fixing grid index k.
 * \author Ole Schuett
 ******************************************************************************/
static inline void ortho_cxyz_to_cxy(const int lp, const int k1, const int k2,
                                     const int cmax,
                                     const double pol[3][lp + 1][2 * cmax + 1],
                                     GRID_CONST_WHEN_COLLOCATE double *cxyz,
                                     GRID_CONST_WHEN_INTEGRATE double *cxy) {

  for (int lzp = 0; lzp <= lp; lzp++) {
    for (int lyp = 0; lyp <= lp - lzp; lyp++) {
      for (int lxp = 0; lxp <= lp - lzp - lyp; lxp++) {
        const double p1 = pol[2][lzp][k1 + cmax];
        const double p2 = pol[2][lzp][k2 + cmax];
        const int cxyz_index =
            lzp * (lp + 1) * (lp + 1) + lyp * (lp + 1) + lxp; // [lzp, lyp, lxp]
        const int cxy_index = lyp * (lp + 1) * 2 + lxp * 2;   // [lyp, lxp, 0]

#if (GRID_DO_COLLOCATE)
        // collocate
        cxy[cxy_index + 0] += cxyz[cxyz_index] * p1;
        cxy[cxy_index + 1] += cxyz[cxyz_index] * p2;
#else
        // integrate
        cxyz[cxyz_index] += cxy[cxy_index + 0] * p1;
        cxyz[cxyz_index] += cxy[cxy_index + 1] * p2;
#endif
      }
    }
  }
}

/*******************************************************************************
 * \brief Collocates coefficients C_xyz onto the grid for orthorhombic case.
 * \author Ole Schuett
 ******************************************************************************/
static inline void
ortho_cxyz_to_grid(const int lp, const double zetp, const double dh[3][3],
                   const double dh_inv[3][3], const double rp[3],
                   const int npts_global[3], const int npts_local[3],
                   const int shift_local[3], const double radius,
                   GRID_CONST_WHEN_COLLOCATE double *cxyz,
                   GRID_CONST_WHEN_INTEGRATE double *grid) {

  // *** position of the gaussian product
  //
  // this is the actual definition of the position on the grid
  // i.e. a point rp(:) gets here grid coordinates
  // MODULO(rp(:)/dr(:),npts_global(:))+1
  // hence (0.0,0.0,0.0) in real space is rsgrid%lb on the rsgrid in Fortran
  // and (1,1,1) on grid here in C.

  // cubecenter(:) = FLOOR(MATMUL(dh_inv, rp))
  int cubecenter[3];
  for (int i = 0; i < 3; i++) {
    double dh_inv_rp = 0.0;
    for (int j = 0; j < 3; j++) {
      dh_inv_rp += dh_inv[j][i] * rp[j];
    }
    cubecenter[i] = (int)floor(dh_inv_rp);
  }

  double roffset[3];
  for (int i = 0; i < 3; i++) {
    roffset[i] = rp[i] - ((double)cubecenter[i]) * dh[i][i];
  }

  // Lookup loop bounds for spherical cutoff.
  int *sphere_bounds;
  double disr_radius;
  grid_sphere_cache_lookup(radius, dh, dh_inv, &sphere_bounds, &disr_radius);
  int **sphere_bounds_iter = &sphere_bounds;

  // Cube bounds.
  int lb_cube[3], ub_cube[3];
  for (int i = 0; i < 3; i++) {
    lb_cube[i] = (int)ceil(-1e-8 - disr_radius * dh_inv[i][i]);
    ub_cube[i] = 1 - lb_cube[i];
    // If grid is not period check that cube fits without wrapping.
    if (npts_global[i] != npts_local[i]) {
      const int offset =
          modulo(cubecenter[i] + lb_cube[i] - shift_local[i], npts_global[i]) -
          lb_cube[i];
      assert(offset + ub_cube[i] < npts_local[i]);
      assert(offset + lb_cube[i] >= 0);
    }
  }

  // cmax = MAXVAL(ub_cube)
  const int cmax = imax(imax(ub_cube[0], ub_cube[1]), ub_cube[2]);

  // Precompute (x-xp)**lp*exp(..) for each direction.
  double pol_mutable[3][lp + 1][2 * cmax + 1];
  for (int idir = 0; idir < 3; idir++) {
    const double dr = dh[idir][idir];
    const double ro = roffset[idir];
    //  Reuse the result from the previous gridpoint to avoid to many exps:
    //  exp( -a*(x+d)**2) = exp(-a*x**2)*exp(-2*a*x*d)*exp(-a*d**2)
    //  exp(-2*a*(x+d)*d) = exp(-2*a*x*d)*exp(-2*a*d**2)
    const double t_exp_1 = exp(-zetp * pow(dr, 2));
    const double t_exp_2 = pow(t_exp_1, 2);
    double t_exp_min_1 = exp(-zetp * pow(+dr - ro, 2));
    double t_exp_min_2 = exp(-2 * zetp * (+dr - ro) * (-dr));
    for (int ig = 0; ig >= lb_cube[idir]; ig--) {
      const double rpg = ig * dr - ro;
      t_exp_min_1 *= t_exp_min_2 * t_exp_1;
      t_exp_min_2 *= t_exp_2;
      double pg = t_exp_min_1;
      for (int icoef = 0; icoef <= lp; icoef++) {
        pol_mutable[idir][icoef][ig + cmax] = pg; // exp(-zetp*rpg**2)
        pg *= rpg;
      }
    }
    double t_exp_plus_1 = exp(-zetp * pow(-ro, 2));
    double t_exp_plus_2 = exp(-2 * zetp * (-ro) * (+dr));
    for (int ig = 0; ig >= lb_cube[idir]; ig--) {
      const double rpg = (1 - ig) * dr - ro;
      t_exp_plus_1 *= t_exp_plus_2 * t_exp_1;
      t_exp_plus_2 *= t_exp_2;
      double pg = t_exp_plus_1;
      for (int icoef = 0; icoef <= lp; icoef++) {
        pol_mutable[idir][icoef][1 - ig + cmax] = pg; // exp(-zetp*rpg**2)
        pg *= rpg;
      }
    }
  }
  const double(*pol)[lp + 1][2 * cmax + 1] =
      (const double(*)[lp + 1][2 * cmax + 1]) pol_mutable;

  // Precompute mapping from cube to grid indices for each direction
  int map_mutable[3][2 * cmax + 1];
  for (int i = 0; i < 3; i++) {
    for (int k = -cmax; k <= +cmax; k++) {
      map_mutable[i][k + cmax] =
          modulo(cubecenter[i] + k - shift_local[i], npts_global[i]);
    }
  }
  const int(*map)[2 * cmax + 1] = (const int(*)[2 * cmax + 1]) map_mutable;

  // Precompute lenght of sections with homogeneous cube to grid mapping.
  int sections_mutable[3][2 * cmax + 1];
  for (int i = 0; i < 3; i++) {
    for (int kg = 2 * cmax; kg >= 0; kg--) {
      if (kg == 2 * cmax || map[i][kg] != map[i][kg + 1] - 1) {
        sections_mutable[i][kg] = 0;
      } else {
        sections_mutable[i][kg] = sections_mutable[i][kg + 1] + 1;
      }
    }
  }
  const int(*sections)[2 * cmax + 1] =
      (const int(*)[2 * cmax + 1]) sections_mutable;

  // Loop over k dimension of the cube.
  const int kstart = *((*sphere_bounds_iter)++);
  const size_t cxy_size = (lp + 1) * (lp + 1) * 2;
  double cxy[cxy_size];
  for (int k1 = kstart; k1 <= 0; k1++) {
    const int k2 = 1 - k1;
    const int kg1 = map[2][k1 + cmax];
    const int kg2 = map[2][k2 + cmax];

    memset(cxy, 0, cxy_size * sizeof(double));

#if (GRID_DO_COLLOCATE)
    // collocate
    ortho_cxyz_to_cxy(lp, k1, k2, cmax, pol, cxyz, cxy);
    ortho_cxy_to_grid(lp, kg1, kg2, cmax, pol, map, sections, npts_local,
                      sphere_bounds_iter, cxy, grid);
#else
    // integrate
    ortho_cxy_to_grid(lp, kg1, kg2, cmax, pol, map, sections, npts_local,
                      sphere_bounds_iter, cxy, grid);
    ortho_cxyz_to_cxy(lp, k1, k2, cmax, pol, cxyz, cxy);
#endif
  }
}

/*******************************************************************************
 * \brief Collocates coefficients C_i onto the grid for general case.
 * \author Ole Schuett
 ******************************************************************************/
static inline void __attribute__((always_inline))
general_ci_to_grid(const int lp, const int jg, const int kg, const int ismin,
                   const int ismax, const int npts_local[3],
                   const int index_min[3], const int index_max[3],
                   const int map_i[], const int sections_i[],
                   const double gp[3], const int k, const int j,
                   const double exp_ij[], const double exp_jk[],
                   const double exp_ki[], GRID_CONST_WHEN_COLLOCATE double *ci,
                   GRID_CONST_WHEN_INTEGRATE double *grid) {

  const int base = kg * npts_local[1] * npts_local[0] + jg * npts_local[0];

  // AVX instructions can only load/store from evenly spaced memory locations.
  // Since the cube can wrap around due to the grid's periodicity,
  // the inner loop runs over sections with homogeneous cube to grid mapping.
  for (int istart = ismin; istart <= ismax; istart++) {
    const int istop = imin(ismax, istart + sections_i[istart - index_min[0]]);
    if (map_i[istart - index_min[0]] < 0) {
      istart = istop; // skip over out-of-bounds indicies
      continue;
    }

    const int cube2grid = map_i[istart - index_min[0]] - istart;
    for (int i = istart; i <= istop; i++) {
      const int ig = i + cube2grid;
      const double di = i - gp[0];

      const int stride_i = index_max[0] - index_min[0] + 1;
      const int stride_j = index_max[1] - index_min[1] + 1;
      const int stride_k = index_max[2] - index_min[2] + 1;
      const int idx_ij = (j - index_min[1]) * stride_i + i - index_min[0];
      const int idx_jk = (k - index_min[2]) * stride_j + j - index_min[1];
      const int idx_ki = (i - index_min[0]) * stride_k + k - index_min[2];

      // Mathieu's trick: Calculate 3D Gaussian from three precomputed 2D tables
      //
      // r   =  (i-gp[0])*dh[0,:] + (j-gp[1])*dh[1,:] + (k-gp[2])*dh[2,:]
      //     =  a                 + b                 + c
      //
      // r**2  =  (a + b + c)**2  =  a**2 + b**2 + c**2 + 2ab + 2bc + 2ca
      //
      // exp(-r**2)  =  exp(-a(a+2b)) * exp(-b*(b+2c)) * exp(-c*(c+2a))
      //
      const double gaussian = exp_ij[idx_ij] * exp_jk[idx_jk] * exp_ki[idx_ki];

      const int grid_index = base + ig; // [kg, jg, ig]
      double dip = gaussian;

#if (GRID_DO_COLLOCATE)
      // collocate
      double reg = 0.0;
      for (int il = 0; il <= lp; il++) {
        reg += ci[il] * dip;
        dip *= di;
      }
      grid[grid_index] += reg;
#else
      // integrate
      const double reg = grid[grid_index];
      for (int il = 0; il <= lp; il++) {
        ci[il] += reg * dip;
        dip *= di;
      }
#endif
    }
    istart = istop;
  }
}

/*******************************************************************************
 * \brief Transforms coefficients C_ij into C_i by fixing grid index j.
 * \author Ole Schuett
 ******************************************************************************/
static inline void __attribute__((always_inline))
general_cij_to_ci(const int lp, const double dj,
                  GRID_CONST_WHEN_COLLOCATE double *cij,
                  GRID_CONST_WHEN_INTEGRATE double *ci) {
  double djp = 1.0;
  for (int jl = 0; jl <= lp; jl++) {
    for (int il = 0; il <= lp - jl; il++) {
      const int cij_index = jl * (lp + 1) + il; // [jl, il]
#if (GRID_DO_COLLOCATE)
      ci[il] += cij[cij_index] * djp; // collocate
#else
      cij[cij_index] += ci[il] * djp; // integrate
#endif
    }
    djp *= dj;
  }
}

/*******************************************************************************
 * \brief Loop body of general_cij_to_grid to be inlined for low values of lp.
 * \author Ole Schuett
 ******************************************************************************/
static inline void __attribute__((always_inline)) general_cij_to_grid_low(
    const int lp, const int jg, const int kg, const int ismin, const int ismax,
    const int npts_local[3], const int index_min[3], const int index_max[3],
    const int map_i[], const int sections_i[], const double gp[3], const int k,
    const int j, const double exp_ij[], const double exp_jk[],
    const double exp_ki[], const double dj, double *ci,
    GRID_CONST_WHEN_COLLOCATE double *cij,
    GRID_CONST_WHEN_INTEGRATE double *grid) {

#if (GRID_DO_COLLOCATE)
  // collocate
  general_cij_to_ci(lp, dj, cij, ci);
  general_ci_to_grid(lp, jg, kg, ismin, ismax, npts_local, index_min, index_max,
                     map_i, sections_i, gp, k, j, exp_ij, exp_jk, exp_ki, ci,
                     grid);
#else
  // integrate
  general_ci_to_grid(lp, jg, kg, ismin, ismax, npts_local, index_min, index_max,
                     map_i, sections_i, gp, k, j, exp_ij, exp_jk, exp_ki, ci,
                     grid);
  general_cij_to_ci(lp, dj, cij, ci);
#endif
}

/*******************************************************************************
 * \brief Collocates coefficients C_ij onto the grid for general case.
 * \author Ole Schuett
 ******************************************************************************/
static inline void general_cij_to_grid(
    const int lp, const int k, const int kg, const int npts_local[3],
    const int index_min[3], const int index_max[3], const int map_i[],
    const int map_j[], const int sections_i[], const int sections_j[],
    const double dh[3][3], const double gp[3], const double radius,
    const double exp_ij[], const double exp_jk[], const double exp_ki[],
    GRID_CONST_WHEN_COLLOCATE double *cij,
    GRID_CONST_WHEN_INTEGRATE double *grid) {

  for (int j = index_min[1]; j <= index_max[1]; j++) {
    const int jg = map_j[j - index_min[1]];
    if (jg < 0) {
      j += sections_j[j - index_min[1]]; // skip over out-of-bounds indicies
      continue;
    }

    //--------------------------------------------------------------------
    // Find bounds for the inner loop based on a quadratic equation in i.
    //
    // The real-space vector from the center of the gaussian to the
    // grid point i,j,k is given by:
    //   r = (i-gp[0])*dh[0,:] + (j-gp[1])*dh[1,:] + (k-gp[2])*dh[2,:]
    //
    // Separating the term that depends on i:
    //   r = i*dh[0,:] - gp[0]*dh[0,:] + (j-gp[1])*dh[1,:] + (k-gp[2])*dh[2,:]
    //     = i*dh[0,:] + v
    //
    // The squared distance works out to:
    //   r**2 = dh[0,:]**2 * i**2  +  2 * v * dh[0,:] * i  +  v**2
    //        = a * i**2           +  b * i                +  c
    //
    // Solving r**2==radius**2 for i yields:
    //    d =  b**2  -  4 * a * (c - radius**2)
    //    i = (-b \pm sqrt(d)) / (2*a)
    //
    double a = 0.0, b = 0.0, c = 0.0;
    for (int i = 0; i < 3; i++) {
      const double v = (0 - gp[0]) * dh[0][i] + (j - gp[1]) * dh[1][i] +
                       (k - gp[2]) * dh[2][i];
      a += dh[0][i] * dh[0][i];
      b += 2.0 * v * dh[0][i];
      c += v * v;
    }
    const double d = b * b - 4.0 * a * (c - radius * radius);

    if (0.0 < d) {
      const double sqrt_d = sqrt(d);
      const double inv_2a = 1.0 / (2.0 * a);
      const int ismin = (int)ceil((-b - sqrt_d) * inv_2a);
      const int ismax = (int)floor((-b + sqrt_d) * inv_2a);
      const double dj = j - gp[1];

      double ci[lp + 1];
      memset(ci, 0, sizeof(ci));

      // Generate separate branches for low values of lp.
      if (lp <= GRID_MAX_LP_OPTIMIZED) {
        GRID_PRAGMA_UNROLL(GRID_MAX_LP_OPTIMIZED + 1)
        for (int ilp = 0; ilp <= GRID_MAX_LP_OPTIMIZED; ilp++) {
          if (lp == ilp) {
            general_cij_to_grid_low(ilp, jg, kg, ismin, ismax, npts_local,
                                    index_min, index_max, map_i, sections_i, gp,
                                    k, j, exp_ij, exp_jk, exp_ki, dj, ci, cij,
                                    grid);
          }
        }
      } else {
        general_cij_to_grid_low(lp, jg, kg, ismin, ismax, npts_local, index_min,
                                index_max, map_i, sections_i, gp, k, j, exp_ij,
                                exp_jk, exp_ki, dj, ci, cij, grid);
      }
    }
  }
}

/*******************************************************************************
 * \brief Transforms coefficients C_ijk into C_ij by fixing grid index k.
 * \author Ole Schuett
 ******************************************************************************/
static inline void general_cijk_to_cij(const int lp, const double dk,
                                       GRID_CONST_WHEN_COLLOCATE double *cijk,
                                       GRID_CONST_WHEN_INTEGRATE double *cij) {
  double dkp = 1.0;
  for (int kl = 0; kl <= lp; kl++) {
    for (int jl = 0; jl <= lp - kl; jl++) {
      for (int il = 0; il <= lp - kl - jl; il++) {
        const int cij_index = jl * (lp + 1) + il; // [jl, il]
        const int cijk_index =
            kl * (lp + 1) * (lp + 1) + jl * (lp + 1) + il; // [kl, jl, il]
#if (GRID_DO_COLLOCATE)
        cij[cij_index] += cijk[cijk_index] * dkp; // collocate
#else
        cijk[cijk_index] += cij[cij_index] * dkp; // integrate
#endif
      }
    }
    dkp *= dk;
  }
}

/*******************************************************************************
 * \brief Precompute mapping of grid indices and its homogeneous sections.
 * \author Ole Schuett
 ******************************************************************************/
static inline void precompute_mapping(const int index_min, const int index_max,
                                      const int shift_local,
                                      const int npts_global,
                                      const int bounds[2], int map[],
                                      int sections[]) {

  // Precompute mapping from continous grid indices to pbc wraped.
  for (int k = index_min; k <= index_max; k++) {
    const int kg = modulo(k - shift_local, npts_global);
    if (bounds[0] <= kg && kg <= bounds[1]) {
      map[k - index_min] = kg;
    } else {
      map[k - index_min] = INT_MIN; // out of bounds - not mapped
    }
  }

  // Precompute lenght of sections with homogeneous cube to grid mapping.
  const int range = index_max - index_min + 1;
  for (int kg = range - 1; kg >= 0; kg--) {
    if (kg == range - 1 || map[kg] != map[kg + 1] - 1) {
      sections[kg] = 0;
    } else {
      sections[kg] = sections[kg + 1] + 1;
    }
  }
}

/*******************************************************************************
 * \brief Fill one of the 2D tables that is used to assemble the 3D Gaussian.
 * \author Ole Schuett
 ******************************************************************************/
static inline void
precompute_exp_table(const int idir, const int jdir, const int index_min[3],
                     const int index_max[3], const double zetp,
                     const double dh[3][3], const double gp[3],
                     double exp_table[]) {

  const int stride_i = index_max[idir] - index_min[idir] + 1;
  const double h_ii = dh[idir][0] * dh[idir][0] + dh[idir][1] * dh[idir][1] +
                      dh[idir][2] * dh[idir][2];
  const double h_ij = dh[idir][0] * dh[jdir][0] + dh[idir][1] * dh[jdir][1] +
                      dh[idir][2] * dh[jdir][2];

  for (int i = index_min[idir]; i <= index_max[idir]; i++) {
    const double di = i - gp[idir];
    const double rii = di * di * h_ii;
    const double rij_unit = di * h_ij;
    const double exp_ij_unit = exp(-zetp * 2.0 * rij_unit);

    // compute exponentials symmetrically around cube center
    const int j_center = (int)gp[jdir];
    const double dj_center = j_center - gp[jdir];
    const double rij_center = dj_center * rij_unit;
    const double exp_ij_center = exp(-zetp * (rii + 2.0 * rij_center));

    // above center
    double exp_ij = exp_ij_center;
    for (int j = j_center; j <= index_max[jdir]; j++) {
      const int idx = (j - index_min[jdir]) * stride_i + i - index_min[idir];
      exp_table[idx] = exp_ij; // exp(-zetp * (di*di*h_ii + 2*di*dj*h_ij));
      exp_ij *= exp_ij_unit;
    }

    // below center
    const double exp_ij_unit_inv = 1.0 / exp_ij_unit;
    exp_ij = exp_ij_center * exp_ij_unit_inv;
    for (int j = j_center - 1; j >= index_min[jdir]; j--) {
      const int idx = (j - index_min[jdir]) * stride_i + i - index_min[idir];
      exp_table[idx] = exp_ij; // exp(-zetp * (di*di*h_ii + 2*di*dj*h_ij));
      exp_ij *= exp_ij_unit_inv;
    }
  }
}

/*******************************************************************************
 * \brief Collocates coefficients C_ijk onto the grid for general case.
 * \author Ole Schuett
 ******************************************************************************/
static inline void
general_cijk_to_grid(const int border_mask, const int lp, const double zetp,
                     const double dh[3][3], const double dh_inv[3][3],
                     const double rp[3], const int npts_global[3],
                     const int npts_local[3], const int shift_local[3],
                     const int border_width[3], const double radius,
                     GRID_CONST_WHEN_COLLOCATE double *cijk,
                     GRID_CONST_WHEN_INTEGRATE double *grid) {

  // Default for border_mask == 0.
  int bounds_i[2] = {0, npts_local[0] - 1};
  int bounds_j[2] = {0, npts_local[1] - 1};
  int bounds_k[2] = {0, npts_local[2] - 1};

  // See also rs_find_node() in task_list_methods.F.
  // If the bit is set then we need to exclude the border in that direction.
  if (border_mask & (1 << 0))
    bounds_i[0] += border_width[0];
  if (border_mask & (1 << 1))
    bounds_i[1] -= border_width[0];
  if (border_mask & (1 << 2))
    bounds_j[0] += border_width[1];
  if (border_mask & (1 << 3))
    bounds_j[1] -= border_width[1];
  if (border_mask & (1 << 4))
    bounds_k[0] += border_width[2];
  if (border_mask & (1 << 5))
    bounds_k[1] -= border_width[2];

  // center in grid coords
  // gp = MATMUL(dh_inv, rp)
  double gp[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      gp[i] += dh_inv[j][i] * rp[j];
    }
  }

  // Get the min max indices that contain at least the cube that contains a
  // sphere around rp of radius radius if the cell is very non-orthogonal this
  // implies that many useless points are included this estimate can be improved
  // (i.e. not box but sphere should be used)
  int index_min[3] = {INT_MAX, INT_MAX, INT_MAX};
  int index_max[3] = {INT_MIN, INT_MIN, INT_MIN};
  for (int i = -1; i <= 1; i++) {
    for (int j = -1; j <= 1; j++) {
      for (int k = -1; k <= 1; k++) {
        const double x = rp[0] + i * radius;
        const double y = rp[1] + j * radius;
        const double z = rp[2] + k * radius;
        for (int idir = 0; idir < 3; idir++) {
          const double resc =
              dh_inv[0][idir] * x + dh_inv[1][idir] * y + dh_inv[2][idir] * z;
          index_min[idir] = imin(index_min[idir], (int)floor(resc));
          index_max[idir] = imax(index_max[idir], (int)ceil(resc));
        }
      }
    }
  }

  // Precompute mappings
  const int range_i = index_max[0] - index_min[0] + 1;
  int map_i[range_i], sections_i[range_i];
  precompute_mapping(index_min[0], index_max[0], shift_local[0], npts_global[0],
                     bounds_i, map_i, sections_i);
  const int range_j = index_max[1] - index_min[1] + 1;
  int map_j[range_j], sections_j[range_j];
  precompute_mapping(index_min[1], index_max[1], shift_local[1], npts_global[1],
                     bounds_j, map_j, sections_j);
  const int range_k = index_max[2] - index_min[2] + 1;
  int map_k[range_k], sections_k[range_k];
  precompute_mapping(index_min[2], index_max[2], shift_local[2], npts_global[2],
                     bounds_k, map_k, sections_k);

  // Precompute exponentials
  double exp_ij[range_i * range_j];
  precompute_exp_table(0, 1, index_min, index_max, zetp, dh, gp, exp_ij);
  double exp_jk[range_j * range_k];
  precompute_exp_table(1, 2, index_min, index_max, zetp, dh, gp, exp_jk);
  double exp_ki[range_k * range_i];
  precompute_exp_table(2, 0, index_min, index_max, zetp, dh, gp, exp_ki);

  // go over the grid, but cycle if the point is not within the radius
  const int cij_size = (lp + 1) * (lp + 1);
  double cij[cij_size];
  for (int k = index_min[2]; k <= index_max[2]; k++) {
    const int kg = map_k[k - index_min[2]];
    if (kg < 0) {
      k += sections_k[k - index_min[2]]; // skip over out-of-bounds indicies
      continue;
    }

    // zero coef_xyt
    memset(cij, 0, cij_size * sizeof(double));

#if (GRID_DO_COLLOCATE)
    // collocate
    general_cijk_to_cij(lp, (double)k - gp[2], cijk, cij);
    general_cij_to_grid(lp, k, kg, npts_local, index_min, index_max, map_i,
                        map_j, sections_i, sections_j, dh, gp, radius, exp_ij,
                        exp_jk, exp_ki, cij, grid);
#else
    // integrate
    general_cij_to_grid(lp, k, kg, npts_local, index_min, index_max, map_i,
                        map_j, sections_i, sections_j, dh, gp, radius, exp_ij,
                        exp_jk, exp_ki, cij, grid);
    general_cijk_to_cij(lp, (double)k - gp[2], cijk, cij);
#endif
  }
}

/*******************************************************************************
 * \brief Transforms coefficients C_xyz into C_ijk.
 * \author Ole Schuett
 ******************************************************************************/
static inline void
general_cxyz_to_cijk(const int lp, const double dh[3][3],
                     GRID_CONST_WHEN_COLLOCATE double *cxyz,
                     GRID_CONST_WHEN_INTEGRATE double *cijk) {

  // transform P_{lxp,lyp,lzp} into a P_{lip,ljp,lkp} such that
  // sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-x_p)**lxp (y-y_p)**lyp (z-z_p)**lzp =
  // sum_{lip,ljp,lkp} P_{lip,ljp,lkp} (i-i_p)**lip (j-j_p)**ljp (k-k_p)**lkp

  // transform using multinomials
  double hmatgridp[lp + 1][3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      hmatgridp[0][j][i] = 1.0;
      for (int k = 1; k <= lp; k++) {
        hmatgridp[k][j][i] = hmatgridp[k - 1][j][i] * dh[j][i];
      }
    }
  }

  const int lpx = lp;
  for (int klx = 0; klx <= lpx; klx++) {
    for (int jlx = 0; jlx <= lpx - klx; jlx++) {
      for (int ilx = 0; ilx <= lpx - klx - jlx; ilx++) {
        const int lx = ilx + jlx + klx;
        const int lpy = lp - lx;
        for (int kly = 0; kly <= lpy; kly++) {
          for (int jly = 0; jly <= lpy - kly; jly++) {
            for (int ily = 0; ily <= lpy - kly - jly; ily++) {
              const int ly = ily + jly + kly;
              const int lpz = lp - lx - ly;
              for (int klz = 0; klz <= lpz; klz++) {
                for (int jlz = 0; jlz <= lpz - klz; jlz++) {
                  for (int ilz = 0; ilz <= lpz - klz - jlz; ilz++) {
                    const int lz = ilz + jlz + klz;
                    const int il = ilx + ily + ilz;
                    const int jl = jlx + jly + jlz;
                    const int kl = klx + kly + klz;
                    const int lp1 = lp + 1;
                    const int cijk_index =
                        kl * lp1 * lp1 + jl * lp1 + il; // [kl,jl,il]
                    const int cxyz_index =
                        lz * lp1 * lp1 + ly * lp1 + lx; // [lz,ly,lx]
                    const double p =
                        hmatgridp[ilx][0][0] * hmatgridp[jlx][1][0] *
                        hmatgridp[klx][2][0] * hmatgridp[ily][0][1] *
                        hmatgridp[jly][1][1] * hmatgridp[kly][2][1] *
                        hmatgridp[ilz][0][2] * hmatgridp[jlz][1][2] *
                        hmatgridp[klz][2][2] * fac(lx) * fac(ly) * fac(lz) /
                        (fac(ilx) * fac(ily) * fac(ilz) * fac(jlx) * fac(jly) *
                         fac(jlz) * fac(klx) * fac(kly) * fac(klz));
#if (GRID_DO_COLLOCATE)
                    cijk[cijk_index] += cxyz[cxyz_index] * p; // collocate
#else
                    cxyz[cxyz_index] += cijk[cijk_index] * p; // integrate
#endif
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

/*******************************************************************************
 * \brief Collocates coefficients C_xyz onto the grid for general case.
 * \author Ole Schuett
 ******************************************************************************/
static inline void
general_cxyz_to_grid(const int border_mask, const int lp, const double zetp,
                     const double dh[3][3], const double dh_inv[3][3],
                     const double rp[3], const int npts_global[3],
                     const int npts_local[3], const int shift_local[3],
                     const int border_width[3], const double radius,
                     GRID_CONST_WHEN_COLLOCATE double *cxyz,
                     GRID_CONST_WHEN_INTEGRATE double *grid) {

  const size_t cijk_size = (lp + 1) * (lp + 1) * (lp + 1);
  double cijk[cijk_size];
  memset(cijk, 0, cijk_size * sizeof(double));

#if (GRID_DO_COLLOCATE)
  // collocate
  general_cxyz_to_cijk(lp, dh, cxyz, cijk);
  general_cijk_to_grid(border_mask, lp, zetp, dh, dh_inv, rp, npts_global,
                       npts_local, shift_local, border_width, radius, cijk,
                       grid);
#else
  // integrate
  general_cijk_to_grid(border_mask, lp, zetp, dh, dh_inv, rp, npts_global,
                       npts_local, shift_local, border_width, radius, cijk,
                       grid);
  general_cxyz_to_cijk(lp, dh, cxyz, cijk);
#endif
}

/*******************************************************************************
 * \brief Collocates coefficients C_xyz onto the grid.
 * \author Ole Schuett
 ******************************************************************************/
static inline void
cxyz_to_grid(const bool orthorhombic, const int border_mask, const int lp,
             const double zetp, const double dh[3][3],
             const double dh_inv[3][3], const double rp[3],
             const int npts_global[3], const int npts_local[3],
             const int shift_local[3], const int border_width[3],
             const double radius, GRID_CONST_WHEN_COLLOCATE double *cxyz,
             GRID_CONST_WHEN_INTEGRATE double *grid) {

  enum grid_library_kernel k;
  if (orthorhombic && border_mask == 0) {
    k = (GRID_DO_COLLOCATE) ? GRID_COLLOCATE_ORTHO : GRID_INTEGRATE_ORTHO;
    ortho_cxyz_to_grid(lp, zetp, dh, dh_inv, rp, npts_global, npts_local,
                       shift_local, radius, cxyz, grid);
  } else {
    k = (GRID_DO_COLLOCATE) ? GRID_COLLOCATE_GENERAL : GRID_INTEGRATE_GENERAL;
    general_cxyz_to_grid(border_mask, lp, zetp, dh, dh_inv, rp, npts_global,
                         npts_local, shift_local, border_width, radius, cxyz,
                         grid);
  }
  grid_library_counter_add(lp, GRID_BACKEND_REF, k, 1);
}

/*******************************************************************************
 * \brief Transforms coefficients C_ab into C_xyz.
 * \author Ole Schuett
 ******************************************************************************/
static inline void cab_to_cxyz(const int la_max, const int la_min,
                               const int lb_max, const int lb_min,
                               const double prefactor, const double ra[3],
                               const double rb[3], const double rp[3],
                               GRID_CONST_WHEN_COLLOCATE double *cab,
                               GRID_CONST_WHEN_INTEGRATE double *cxyz) {

  // Computes the polynomial expansion coefficients:
  //     (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
  const int lp = la_max + lb_max;
  double alpha[3][lb_max + 1][la_max + 1][lp + 1];
  memset(alpha, 0, 3 * (lb_max + 1) * (la_max + 1) * (lp + 1) * sizeof(double));

  for (int i = 0; i < 3; i++) {
    const double drpa = rp[i] - ra[i];
    const double drpb = rp[i] - rb[i];
    for (int lxa = 0; lxa <= la_max; lxa++) {
      for (int lxb = 0; lxb <= lb_max; lxb++) {
        double binomial_k_lxa = 1.0;
        double a = 1.0;
        for (int k = 0; k <= lxa; k++) {
          double binomial_l_lxb = 1.0;
          double b = 1.0;
          for (int l = 0; l <= lxb; l++) {
            alpha[i][lxb][lxa][lxa - l + lxb - k] +=
                binomial_k_lxa * binomial_l_lxb * a * b;
            binomial_l_lxb *= ((double)(lxb - l)) / ((double)(l + 1));
            b *= drpb;
          }
          binomial_k_lxa *= ((double)(lxa - k)) / ((double)(k + 1));
          a *= drpa;
        }
      }
    }
  }

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

  for (int lzb = 0; lzb <= lb_max; lzb++) {
    for (int lza = 0; lza <= la_max; lza++) {
      for (int lyb = 0; lyb <= lb_max - lzb; lyb++) {
        for (int lya = 0; lya <= la_max - lza; lya++) {
          const int lxb_min = imax(lb_min - lzb - lyb, 0);
          const int lxa_min = imax(la_min - lza - lya, 0);
          for (int lxb = lxb_min; lxb <= lb_max - lzb - lyb; lxb++) {
            for (int lxa = lxa_min; lxa <= la_max - lza - lya; lxa++) {
              const int ico = coset(lxa, lya, lza);
              const int jco = coset(lxb, lyb, lzb);
              const int cab_index = jco * ncoset(la_max) + ico; // [jco, ico]
              for (int lzp = 0; lzp <= lza + lzb; lzp++) {
                for (int lyp = 0; lyp <= lp - lza - lzb; lyp++) {
                  for (int lxp = 0; lxp <= lp - lza - lzb - lyp; lxp++) {
                    const double p = alpha[0][lxb][lxa][lxp] *
                                     alpha[1][lyb][lya][lyp] *
                                     alpha[2][lzb][lza][lzp] * prefactor;
                    const int lp1 = lp + 1;
                    const int cxyz_index =
                        lzp * lp1 * lp1 + lyp * lp1 + lxp; // [lzp, lyp, lxp]
#if (GRID_DO_COLLOCATE)
                    cxyz[cxyz_index] += cab[cab_index] * p; // collocate
#else
                    cab[cab_index] += cxyz[cxyz_index] * p; // integrate
#endif
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

/*******************************************************************************
 * \brief Collocates coefficients C_ab onto the grid.
 * \author Ole Schuett
 ******************************************************************************/
static inline void
cab_to_grid(const bool orthorhombic, const int border_mask, const int la_max,
            const int la_min, const int lb_max, const int lb_min,
            const double zeta, const double zetb, const double rscale,
            const double dh[3][3], const double dh_inv[3][3],
            const double ra[3], const double rab[3], const int npts_global[3],
            const int npts_local[3], const int shift_local[3],
            const int border_width[3], const double radius,
            GRID_CONST_WHEN_COLLOCATE double *cab,
            GRID_CONST_WHEN_INTEGRATE double *grid) {

  // Check if radius is too small to be mapped onto grid of given resolution.
  double dh_max = 0.0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      dh_max = fmax(dh_max, fabs(dh[i][j]));
    }
  }
  if (2.0 * radius < dh_max) {
    return;
  }

  const double zetp = zeta + zetb;
  const double f = zetb / zetp;
  const double rab2 = rab[0] * rab[0] + rab[1] * rab[1] + rab[2] * rab[2];
  const double prefactor = rscale * exp(-zeta * f * rab2);
  double rp[3], rb[3];
  for (int i = 0; i < 3; i++) {
    rp[i] = ra[i] + f * rab[i];
    rb[i] = ra[i] + rab[i];
  }

  const int lp = la_max + lb_max;
  const size_t cxyz_size = (lp + 1) * (lp + 1) * (lp + 1);
  double cxyz[cxyz_size];
  memset(cxyz, 0, cxyz_size * sizeof(double));

#if (GRID_DO_COLLOCATE)
  // collocate
  cab_to_cxyz(la_max, la_min, lb_max, lb_min, prefactor, ra, rb, rp, cab, cxyz);
  cxyz_to_grid(orthorhombic, border_mask, lp, zetp, dh, dh_inv, rp, npts_global,
               npts_local, shift_local, border_width, radius, cxyz, grid);
#else
  // integrate
  cxyz_to_grid(orthorhombic, border_mask, lp, zetp, dh, dh_inv, rp, npts_global,
               npts_local, shift_local, border_width, radius, cxyz, grid);
  cab_to_cxyz(la_max, la_min, lb_max, lb_min, prefactor, ra, rb, rp, cab, cxyz);
#endif
}

// EOF
