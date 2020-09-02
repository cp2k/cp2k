/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/grid_common.h"
#include "../common/grid_library.h"
#include "grid_ref_collocate.h"
#include "grid_ref_prepare_pab.h"

//******************************************************************************
// \brief Computes the polynomial expansion coefficients:
//        (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
//        Results are passed to grid_prepare_coef.
// \author Ole Schuett
//******************************************************************************
static void
prepare_alpha(const double ra[3], const double rb[3], const double rp[3],
              const int la_max, const int lb_max,
              double alpha[3][lb_max + 1][la_max + 1][la_max + lb_max + 1]) {

  // Initialize with zeros.
  for (int iaxis = 0; iaxis < 3; iaxis++) {
    for (int lxb = 0; lxb <= lb_max; lxb++) {
      for (int lxa = 0; lxa <= la_max; lxa++) {
        for (int lxp = 0; lxp <= la_max + lb_max; lxp++) {
          alpha[iaxis][lxb][lxa][lxp] = 0.0;
        }
      }
    }
  }

  for (int iaxis = 0; iaxis < 3; iaxis++) {
    const double drpa = rp[iaxis] - ra[iaxis];
    const double drpb = rp[iaxis] - rb[iaxis];
    for (int lxa = 0; lxa <= la_max; lxa++) {
      for (int lxb = 0; lxb <= lb_max; lxb++) {
        double binomial_k_lxa = 1.0;
        double a = 1.0;
        for (int k = 0; k <= lxa; k++) {
          double binomial_l_lxb = 1.0;
          double b = 1.0;
          for (int l = 0; l <= lxb; l++) {
            alpha[iaxis][lxb][lxa][lxa - l + lxb - k] +=
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
}

//******************************************************************************
// \brief Compute coefficients for all combinations of angular momentum.
//        Results are passed to collocate_ortho and collocate_general.
// \author Ole Schuett
//******************************************************************************
static void prepare_coef(const int la_max, const int la_min, const int lb_max,
                         const int lb_min, const int lp, const double prefactor,
                         const double alpha[3][lb_max + 1][la_max + 1][lp + 1],
                         const double pab[ncoset[lb_max]][ncoset[la_max]],
                         double coef_xyz[lp + 1][lp + 1][lp + 1]) {

  memset(coef_xyz, 0, (lp + 1) * (lp + 1) * (lp + 1) * sizeof(double));

  double coef_xyt[lp + 1][lp + 1];
  double coef_xtt[lp + 1];

  for (int lzb = 0; lzb <= lb_max; lzb++) {
    for (int lza = 0; lza <= la_max; lza++) {
      for (int lyp = 0; lyp <= lp - lza - lzb; lyp++) {
        for (int lxp = 0; lxp <= lp - lza - lzb - lyp; lxp++) {
          coef_xyt[lyp][lxp] = 0.0;
        }
      }
      for (int lyb = 0; lyb <= lb_max - lzb; lyb++) {
        for (int lya = 0; lya <= la_max - lza; lya++) {
          const int lxpm = (lb_max - lzb - lyb) + (la_max - lza - lya);
          for (int i = 0; i <= lxpm; i++) {
            coef_xtt[i] = 0.0;
          }
          for (int lxb = LIBGRID_MAX(lb_min - lzb - lyb, 0); lxb <= lb_max - lzb - lyb;
               lxb++) {
            for (int lxa = LIBGRID_MAX(la_min - lza - lya, 0);
                 lxa <= la_max - lza - lya; lxa++) {
              const int ico = coset(lxa, lya, lza);
              const int jco = coset(lxb, lyb, lzb);
              const double p_ele = prefactor * pab[jco][ico];
              for (int lxp = 0; lxp <= lxa + lxb; lxp++) {
                coef_xtt[lxp] += p_ele * alpha[0][lxb][lxa][lxp];
              }
            }
          }
          for (int lyp = 0; lyp <= lya + lyb; lyp++) {
            for (int lxp = 0; lxp <= lp - lza - lzb - lya - lyb; lxp++) {
              coef_xyt[lyp][lxp] += alpha[1][lyb][lya][lyp] * coef_xtt[lxp];
            }
          }
        }
      }
      for (int lzp = 0; lzp <= lza + lzb; lzp++) {
        for (int lyp = 0; lyp <= lp - lza - lzb; lyp++) {
          for (int lxp = 0; lxp <= lp - lza - lzb - lyp; lxp++) {
            coef_xyz[lzp][lyp][lxp] +=
                alpha[2][lzb][lza][lzp] * coef_xyt[lyp][lxp];
          }
        }
      }
    }
  }
}

//******************************************************************************
// \brief Computes mapping from cube to grid indices.
//        Used only in the orthorhombic case.
// \author Ole Schuett
//******************************************************************************
static void fill_map(const int lb_cube, const int ub_cube, const int cubecenter,
                     const int npts_global, const int shift_local,
                     const int cmax, int map[2 * cmax + 1]) {

  for (int i = 0; i < 2 * cmax + 1; i++) {
    map[i] = INT_MAX; // Safety net, will trigger out of bounds.
  }

  for (int ig = lb_cube; ig <= ub_cube; ig++) {
    map[ig + cmax] = LIBGRID_MOD(cubecenter + ig - shift_local, npts_global);
  }
}

//******************************************************************************
// \brief Computes (x-xp)**lp*exp(..) for all cube points in one dimension.
//        Used only in the orthorhombic case.
// \author Ole Schuett
//******************************************************************************
static void fill_pol(const double dr, const double roffset, const int lb_cube,
                     const int lp, const int cmax, const double zetp,
                     double pol[lp + 1][2 * cmax + 1]) {

  //  Reuse the result from the previous gridpoint to avoid to many exps:
  //  exp( -a*(x+d)**2) = exp(-a*x**2)*exp(-2*a*x*d)*exp(-a*d**2)
  //  exp(-2*a*(x+d)*d) = exp(-2*a*x*d)*exp(-2*a*d**2)

  const double t_exp_1 = exp(-zetp * pow(dr, 2));
  const double t_exp_2 = pow(t_exp_1, 2);

  double t_exp_min_1 = exp(-zetp * pow(+dr - roffset, 2));
  double t_exp_min_2 = exp(-2 * zetp * (+dr - roffset) * (-dr));
  for (int ig = 0; ig >= lb_cube; ig--) {
    const double rpg = ig * dr - roffset;
    t_exp_min_1 *= t_exp_min_2 * t_exp_1;
    t_exp_min_2 *= t_exp_2;
    double pg = t_exp_min_1;
    // pg  = EXP(-zetp*rpg**2)
    for (int icoef = 0; icoef <= lp; icoef++) {
      pol[icoef][ig - lb_cube] = pg;
      pg *= rpg;
    }
  }

  double t_exp_plus_1 = exp(-zetp * pow(-roffset, 2));
  double t_exp_plus_2 = exp(-2 * zetp * (-roffset) * (+dr));
  for (int ig = 0; ig >= lb_cube; ig--) {
    const double rpg = (1 - ig) * dr - roffset;
    t_exp_plus_1 *= t_exp_plus_2 * t_exp_1;
    t_exp_plus_2 *= t_exp_2;
    double pg = t_exp_plus_1;
    // pg  = EXP(-zetp*rpg**2)
    for (int icoef = 0; icoef <= lp; icoef++) {
      pol[icoef][1 - ig - lb_cube] = pg;
      pg *= rpg;
    }
  }
}

//******************************************************************************
// \brief  A much simpler but also slower implementation of collocate_core.
// \author Ole Schuett
//******************************************************************************
static void collocate_core_simple(
    const int lp, const int cmax, const double coef_xyz[lp + 1][lp + 1][lp + 1],
    const double pol[3][lp + 1][2 * cmax + 1], const int lb_cube[3],
    const int ub_cube[3], const double dh[3][3], const double dh_inv[3][3],
    const double disr_radius, const int cubecenter[3], const int npts_global[3],
    const int npts_local[3], const int shift_local[3], double *grid) {

  // Create the full cube, ignoring periodicity for now.
  const int nz = ub_cube[2] - lb_cube[2] + 1;
  const int ny = ub_cube[1] - lb_cube[1] + 1;
  const int nx = ub_cube[0] - lb_cube[0] + 1;

  double *cube = malloc(nz * ny * nx * sizeof(double));
  memset(cube, 0, nz * ny * nx * sizeof(double));

  for (int lzp = 0; lzp <= lp; lzp++) {
    for (int lyp = 0; lyp <= lp - lzp; lyp++) {
      for (int lxp = 0; lxp <= lp - lzp - lyp; lxp++) {
        for (int k = 0; k < nz; k++) {
          for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
              cube[k * ny * nx + j * nx + i] += coef_xyz[lzp][lyp][lxp] *
                                                pol[2][lzp][k] *
                                                pol[1][lyp][j] * pol[0][lxp][i];
            }
          }
        }
      }
    }
  }

  //
  // Write cube back to large grid taking periodicity and radius into account.
  //

  // The cube contains an even number of grid points in each direction and
  // collocation is always performed on a pair of two opposing grid points.
  // Hence, the points with index 0 and 1 are both assigned distance zero via
  // the formular distance=(2*index-1)/2.

  const int kgmin = ceil(-1e-8 - disr_radius * dh_inv[2][2]);
  for (int kg = kgmin; kg <= 1 - kgmin; kg++) {
    const int ka = cubecenter[2] + kg - shift_local[2];
    const int k = LIBGRID_MOD(ka, npts_global[2]); // target location on the grid
    const int kd = (2 * kg - 1) / 2; // distance from center in grid points
    const double kr = kd * dh[2][2]; // distance from center in a.u.
    const double kremain = disr_radius * disr_radius - kr * kr;
    const int jgmin = ceil(-1e-8 - sqrt(fmax(0.0, kremain)) * dh_inv[1][1]);
    for (int jg = jgmin; jg <= 1 - jgmin; jg++) {
      const int ja = cubecenter[1] + jg - shift_local[1];
      const int j = LIBGRID_MOD(ja, npts_global[1]); // target location on the grid
      const int jd = (2 * jg - 1) / 2; // distance from center in grid points
      const double jr = jd * dh[1][1]; // distance from center in a.u.
      const double jremain = kremain - jr * jr;
      const int igmin = ceil(-1e-8 - sqrt(fmax(0.0, jremain)) * dh_inv[0][0]);
      for (int ig = igmin; ig <= 1 - igmin; ig++) {
        const int ia = cubecenter[0] + ig - shift_local[0];
        const int i = LIBGRID_MOD(ia, npts_global[0]); // target location on the grid
        const int grid_index =
            k * npts_local[1] * npts_local[0] + j * npts_local[0] + i;
        const int cube_index = (kg - lb_cube[2]) * ny * nx +
                               (jg - lb_cube[1]) * nx + (ig - lb_cube[0]);
        grid[grid_index] += cube[cube_index];
      }
    }
  }

  free(cube);
}

//******************************************************************************
// \brief Fills the 3D cube by taking the outer product of the 1D pol arrays.
//        The majority of cpu cycles are spend in this routine.
//        Used only in the orthorhombic case.
// \author Ole Schuett
//******************************************************************************
static void collocate_core(const int lp, const int cmax,
                           const double coef_xyz[lp + 1][lp + 1][lp + 1],
                           const double pol[3][lp + 1][2 * cmax + 1],
                           const int map[3][2 * cmax + 1], const int lb_cube[3],
                           const double dh[3][3], const double dh_inv[3][3],
                           const double disr_radius, const int npts_local[3],
                           double *grid) {

  const int kgmin = ceil(-1e-8 - disr_radius * dh_inv[2][2]);
  for (int kg = kgmin; kg <= 0; kg++) {
    const int kg2 = 1 - kg;
    const int k = map[2][kg + cmax];
    const int k2 = map[2][kg2 + cmax];

    // initialize coef_xy
    const int ncoef_xy = (lp + 1) * (lp + 2) / 2;
    double coef_xy[ncoef_xy][2];
    for (int i = 0; i < ncoef_xy; i++) {
      coef_xy[i][0] = 0.0;
      coef_xy[i][1] = 0.0;
    }

    for (int lzp = 0; lzp <= lp; lzp++) {
      int lxy = 0;
      for (int lyp = 0; lyp <= lp - lzp; lyp++) {
        for (int lxp = 0; lxp <= lp - lzp - lyp; lxp++) {
          coef_xy[lxy][0] +=
              coef_xyz[lzp][lyp][lxp] * pol[2][lzp][kg - lb_cube[2]];
          coef_xy[lxy][1] +=
              coef_xyz[lzp][lyp][lxp] * pol[2][lzp][kg2 - lb_cube[2]];
          lxy++;
        }
        lxy += lzp;
      }
    }

    const int kd = (2 * kg - 1) / 2; // distance from center in grid points
    const double kr = kd * dh[2][2]; // distance from center in a.u.
    const double kremain = disr_radius * disr_radius - kr * kr;
    const int jgmin = ceil(-1e-8 - sqrt(fmax(0.0, kremain)) * dh_inv[1][1]);
    for (int jg = jgmin; jg <= 0; jg++) {
      const int jg2 = 1 - jg;
      const int j = map[1][jg + cmax];
      const int j2 = map[1][jg2 + cmax];

      // initialize coef_x
      double coef_x[lp + 1][4];
      for (int i = 0; i < lp + 1; i++) {
        for (int j = 0; j < 4; j++) {
          coef_x[i][j] = 0.0;
        }
      }

      int lxy = 0;
      for (int lyp = 0; lyp <= lp; lyp++) {
        for (int lxp = 0; lxp <= lp - lyp; lxp++) {
          coef_x[lxp][0] += coef_xy[lxy][0] * pol[1][lyp][jg - lb_cube[1]];
          coef_x[lxp][1] += coef_xy[lxy][1] * pol[1][lyp][jg - lb_cube[1]];
          coef_x[lxp][2] += coef_xy[lxy][0] * pol[1][lyp][jg2 - lb_cube[1]];
          coef_x[lxp][3] += coef_xy[lxy][1] * pol[1][lyp][jg2 - lb_cube[1]];
          lxy++;
        }
      }

      const int jd = (2 * jg - 1) / 2; // distance from center in grid points
      const double jr = jd * dh[1][1]; // distance from center in a.u.
      const double jremain = kremain - jr * jr;
      const int igmin = ceil(-1e-8 - sqrt(fmax(0.0, jremain)) * dh_inv[0][0]);
      for (int ig = igmin; ig <= 0; ig++) {
        const int ig2 = 1 - ig;
        const int i = map[0][ig + cmax];
        const int i2 = map[0][ig2 + cmax];

        double s01 = 0.0;
        double s02 = 0.0;
        double s03 = 0.0;
        double s04 = 0.0;
        double s05 = 0.0;
        double s06 = 0.0;
        double s07 = 0.0;
        double s08 = 0.0;

        for (int lxp = 0; lxp <= lp; lxp++) {
          s01 += coef_x[lxp][0] * pol[0][lxp][ig - lb_cube[0]];
          s02 += coef_x[lxp][1] * pol[0][lxp][ig - lb_cube[0]];
          s03 += coef_x[lxp][2] * pol[0][lxp][ig - lb_cube[0]];
          s04 += coef_x[lxp][3] * pol[0][lxp][ig - lb_cube[0]];
          s05 += coef_x[lxp][0] * pol[0][lxp][ig2 - lb_cube[0]];
          s06 += coef_x[lxp][1] * pol[0][lxp][ig2 - lb_cube[0]];
          s07 += coef_x[lxp][2] * pol[0][lxp][ig2 - lb_cube[0]];
          s08 += coef_x[lxp][3] * pol[0][lxp][ig2 - lb_cube[0]];
        }

        const int stride = npts_local[1] * npts_local[0];
        grid[k * stride + j * npts_local[0] + i] += s01;
        grid[k2 * stride + j * npts_local[0] + i] += s02;
        grid[k * stride + j2 * npts_local[0] + i] += s03;
        grid[k2 * stride + j2 * npts_local[0] + i] += s04;
        grid[k * stride + j * npts_local[0] + i2] += s05;
        grid[k2 * stride + j * npts_local[0] + i2] += s06;
        grid[k * stride + j2 * npts_local[0] + i2] += s07;
        grid[k2 * stride + j2 * npts_local[0] + i2] += s08;
      }
    }
  }
}

//******************************************************************************
// \brief Collocate kernel for the orthorhombic case.
// \author Ole Schuett
//******************************************************************************
static void collocate_ortho(const int lp, const double zetp,
                            const double coef_xyz[lp + 1][lp + 1][lp + 1],
                            const double dh[3][3], const double dh_inv[3][3],
                            const double rp[3], const int npts_global[3],
                            const int npts_local[3], const int shift_local[3],
                            const double radius, double *grid) {

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
    cubecenter[i] = floor(dh_inv_rp);
  }

  double roffset[3];
  for (int i = 0; i < 3; i++) {
    roffset[i] = rp[i] - ((double)cubecenter[i]) * dh[i][i];
  }

  // Historically, the radius gets discretized.
  const double drmin = fmin(dh[0][0], fmin(dh[1][1], dh[2][2]));
  const double disr_radius = drmin * fmax(1.0, ceil(radius / drmin));

  int lb_cube[3], ub_cube[3];
  for (int i = 0; i < 3; i++) {
    lb_cube[i] = ceil(-1e-8 - disr_radius * dh_inv[i][i]);
    ub_cube[i] = 1 - lb_cube[i];
    // If grid is not period check that cube fits without wrapping.
    if (npts_global[i] != npts_local[i]) {
      const int offset =
          LIBGRID_MOD(cubecenter[i] + lb_cube[i] - shift_local[i], npts_global[i]) -
          lb_cube[i];
      assert(offset + ub_cube[i] < npts_local[i]);
      assert(offset + lb_cube[i] >= 0);
    }
  }

  // cmax = MAXVAL(ub_cube)
  int cmax = INT_MIN;
  for (int i = 0; i < 3; i++) {
    cmax = LIBGRID_MAX(cmax, ub_cube[i]);
  }

  double pol_mutable[3][lp + 1][2 * cmax + 1];
  for (int i = 0; i < 3; i++) {
    fill_pol(dh[i][i], roffset[i], lb_cube[i], lp, cmax, zetp, pol_mutable[i]);
  }
  const double(*pol)[lp + 1][2 * cmax + 1] =
      (const double(*)[lp + 1][2 * cmax + 1]) pol_mutable;

  // Enable to run a much simpler, but also slower implementation.
  if (false) {
    collocate_core_simple(lp, cmax, coef_xyz, pol, lb_cube, ub_cube, dh, dh_inv,
                          disr_radius, cubecenter, npts_global, npts_local,
                          shift_local, grid);
  } else {
    // a mapping so that the ig corresponds to the right grid point
    int map_mutable[3][2 * cmax + 1];
    for (int i = 0; i < 3; i++) {
      fill_map(lb_cube[i], ub_cube[i], cubecenter[i], npts_global[i],
               shift_local[i], cmax, map_mutable[i]);
    }
    const int(*map)[2 * cmax + 1] = (const int(*)[2 * cmax + 1]) map_mutable;

    collocate_core(lp, cmax, coef_xyz, pol, map, lb_cube, dh, dh_inv,
                   disr_radius, npts_local, grid);
  }
}

//******************************************************************************
// \brief Collocate kernel for general case, ie. non-ortho or with subpatches.
// \author Ole Schuett
//******************************************************************************
static void collocate_general(const int border_mask, const int lp,
                              const double zetp,
                              const double coef_xyz[lp + 1][lp + 1][lp + 1],
                              const double dh[3][3], const double dh_inv[3][3],
                              const double rp[3], const int npts_global[3],
                              const int npts_local[3], const int shift_local[3],
                              const int border_width[3], const double radius,
                              double *grid) {

  int bounds[3][2] = {{0, npts_local[0] - 1}, // Default for border_mask == 0.
                      {0, npts_local[1] - 1},
                      {0, npts_local[2] - 1}};

  // See also rs_find_node() in task_list_methods.F.
  // If the bit is set then we need to exclude the border in that direction.
  if (border_mask & (1 << 0))
    bounds[0][0] += border_width[0];
  if (border_mask & (1 << 1))
    bounds[0][1] -= border_width[0];
  if (border_mask & (1 << 2))
    bounds[1][0] += border_width[1];
  if (border_mask & (1 << 3))
    bounds[1][1] -= border_width[1];
  if (border_mask & (1 << 4))
    bounds[2][0] += border_width[2];
  if (border_mask & (1 << 5))
    bounds[2][1] -= border_width[2];

  // Translated from collocate_general_opt()
  //
  // transform P_{lxp,lyp,lzp} into a P_{lip,ljp,lkp} such that
  // sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-x_p)**lxp (y-y_p)**lyp (z-z_p)**lzp =
  // sum_{lip,ljp,lkp} P_{lip,ljp,lkp} (i-i_p)**lip (j-j_p)**ljp (k-k_p)**lkp
  //

  // aux mapping array to simplify life
  int coef_map[lp + 1][lp + 1][lp + 1];

  // Safety net, will trigger out of bounds.
  for (int lzp = 0; lzp <= lp; lzp++) {
    for (int lyp = 0; lyp <= lp; lyp++) {
      for (int lxp = 0; lxp <= lp; lxp++) {
        coef_map[lzp][lyp][lxp] = INT_MAX;
      }
    }
  }

  int lxyz = 0;
  for (int lzp = 0; lzp <= lp; lzp++) {
    for (int lyp = 0; lyp <= lp - lzp; lyp++) {
      for (int lxp = 0; lxp <= lp - lzp - lyp; lxp++) {
        coef_map[lzp][lyp][lxp] = ++lxyz;
      }
    }
  }

  // center in grid coords
  // gp = MATMUL(dh_inv, rp)
  double gp[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      gp[i] += dh_inv[j][i] * rp[j];
    }
  }

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

  // zero coef_ijk
  const int ncoef_ijk = ((lp + 1) * (lp + 2) * (lp + 3)) / 6;
  double coef_ijk[ncoef_ijk];
  for (int i = 0; i < ncoef_ijk; i++) {
    coef_ijk[i] = 0.0;
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
                    const int lijk = coef_map[kl][jl][il];
                    coef_ijk[lijk - 1] +=
                        coef_xyz[lz][ly][lx] * hmatgridp[ilx][0][0] *
                        hmatgridp[jlx][1][0] * hmatgridp[klx][2][0] *
                        hmatgridp[ily][0][1] * hmatgridp[jly][1][1] *
                        hmatgridp[kly][2][1] * hmatgridp[ilz][0][2] *
                        hmatgridp[jlz][1][2] * hmatgridp[klz][2][2] * fac[lx] *
                        fac[ly] * fac[lz] /
                        (fac[ilx] * fac[ily] * fac[ilz] * fac[jlx] * fac[jly] *
                         fac[jlz] * fac[klx] * fac[kly] * fac[klz]);
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // CALL return_cube_nonortho(cube_info, radius, index_min, index_max, rp)
  //
  // get the min max indices that contain at least the cube that contains a
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
          index_min[idir] = LIBGRID_MIN(index_min[idir], floor(resc));
          index_max[idir] = LIBGRID_MAX(index_max[idir], ceil(resc));
        }
      }
    }
  }

  // go over the grid, but cycle if the point is not within the radius
  for (int k = index_min[2]; k <= index_max[2]; k++) {
    const int k_index = LIBGRID_MOD(k - shift_local[2], npts_global[2]);
    if (k_index < bounds[2][0] || bounds[2][1] < k_index) {
      continue;
    }

    // zero coef_xyt
    const int ncoef_xyt = ((lp + 1) * (lp + 2)) / 2;
    double coef_xyt[ncoef_xyt];
    for (int i = 0; i < ncoef_xyt; i++) {
      coef_xyt[i] = 0.0;
    }

    int lxyz = 0;
    double dkp = 1.0;
    const double dk = k - gp[2];
    for (int kl = 0; kl <= lp; kl++) {
      int lxy = 0;
      for (int jl = 0; jl <= lp - kl; jl++) {
        for (int il = 0; il <= lp - kl - jl; il++) {
          coef_xyt[lxy++] += coef_ijk[lxyz++] * dkp;
        }
        lxy += kl;
      }
      dkp *= dk;
    }

    for (int j = index_min[1]; j <= index_max[1]; j++) {
      const int j_index = LIBGRID_MOD(j - shift_local[1], npts_global[1]);
      if (j_index < bounds[1][0] || bounds[1][1] < j_index) {
        continue;
      }

      double coef_xtt[lp + 1];
      for (int i = 0; i <= lp; i++) {
        coef_xtt[i] = 0.0;
      }
      int lxy = 0;
      double djp = 1.0;
      const double dj = j - gp[1];
      for (int jl = 0; jl <= lp; jl++) {
        for (int il = 0; il <= lp - jl; il++) {
          coef_xtt[il] += coef_xyt[lxy++] * djp;
        }
        djp *= dj;
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
      if (d <= 0.0) {
        continue;
      }
      const double sqrt_d = sqrt(d);
      const int ismin = ceil((-b - sqrt_d) / (2.0 * a));
      const int ismax = floor((-b + sqrt_d) / (2.0 * a));

      for (int i = ismin; i <= ismax; i++) {
        const int i_index = LIBGRID_MOD(i - shift_local[0], npts_global[0]);
        if (i_index < bounds[0][0] || bounds[0][1] < i_index) {
          continue;
        }

        // polynomial terms
        double res = 0.0;
        double dip = 1.0;
        const double di = i - gp[0];
        for (int il = 0; il <= lp; il++) {
          res += coef_xtt[il] * dip;
          dip *= di;
        }

        // Could be done with a recursion, but would obfuscate code a lot.
        res *= exp(-zetp * ((a * i + b) * i + c));

        const int grid_index = k_index * npts_local[1] * npts_local[0] +
                               j_index * npts_local[0] + i_index;
        grid[grid_index] += res;
      }
    }
  }
}

//******************************************************************************
// \brief Collocates a single product of primitiv Gaussians.
//        See grid_collocate.h for details.
// \author Ole Schuett
//******************************************************************************
void grid_ref_collocate_pgf_product(
    const bool orthorhombic, const int border_mask, const int func,
    const int la_max, const int la_min, const int lb_max, const int lb_min,
    const double zeta, const double zetb, const double rscale,
    const double dh[3][3], const double dh_inv[3][3], const double ra[3],
    const double rab[3], const int npts_global[3], const int npts_local[3],
    const int shift_local[3], const int border_width[3], const double radius,
    const int o1, const int o2, const int n1, const int n2,
    const double pab[n2][n1], double *grid) {

  // Check if radius is too small to be mapped onto grid of given resolution.
  double dh_max = 0.0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      dh_max = fmax(dh_max, fabs(dh[i][j]));
  if (2.0 * radius < dh_max)
    return;

  const double zetp = zeta + zetb;
  const double f = zetb / zetp;
  const double rab2 = rab[0] * rab[0] + rab[1] * rab[1] + rab[2] * rab[2];
  const double prefactor = rscale * exp(-zeta * f * rab2);
  double rp[3], rb[3];
  for (int i = 0; i < 3; i++) {
    rp[i] = ra[i] + f * rab[i];
    rb[i] = ra[i] + rab[i];
  }

  int la_min_diff, la_max_diff, lb_min_diff, lb_max_diff;
  grid_ref_prepare_get_ldiffs(func, &la_min_diff, &la_max_diff, &lb_min_diff,
                              &lb_max_diff);

  const int la_min_prep = LIBGRID_MAX(la_min + la_min_diff, 0);
  const int lb_min_prep = LIBGRID_MAX(lb_min + lb_min_diff, 0);
  const int la_max_prep = la_max + la_max_diff;
  const int lb_max_prep = lb_max + lb_max_diff;
  const int lp = la_max_prep + lb_max_prep;

  const int n1_prep = ncoset[la_max_prep];
  const int n2_prep = ncoset[lb_max_prep];
  double pab_prep_mutable[n2_prep][n1_prep];
  memset(pab_prep_mutable, 0, n2_prep * n1_prep * sizeof(double));
  grid_ref_prepare_pab(func, o1, o2, la_max, la_min, lb_max, lb_min, zeta, zetb,
                       n1, n2, pab, n1_prep, n2_prep, pab_prep_mutable);
  const double(*pab_prep)[n1_prep] = (const double(*)[n1_prep])pab_prep_mutable;

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

  double alpha_mutable[3][lb_max_prep + 1][la_max_prep + 1][lp + 1];
  prepare_alpha(ra, rb, rp, la_max_prep, lb_max_prep, alpha_mutable);
  const double(*alpha)[lb_max_prep + 1][la_max_prep + 1][lp + 1] =
      (const double(*)[lb_max_prep + 1][la_max_prep + 1][lp + 1]) alpha_mutable;

  //
  //   compute P_{lxp,lyp,lzp} given P_{lxa,lya,lza,lxb,lyb,lzb} and
  //   alpha(ls,lxa,lxb,1) use a three step procedure we don't store zeros, so
  //   counting is done using lxyz,lxy in order to have contiguous memory access
  //   in collocate_fast.F
  //

  double coef_xyz_mutable[lp + 1][lp + 1][lp + 1];
  prepare_coef(la_max_prep, la_min_prep, lb_max_prep, lb_min_prep, lp,
               prefactor, alpha, pab_prep, coef_xyz_mutable);
  const double(*coef_xyz)[lp + 1][lp + 1] =
      (const double(*)[lp + 1][lp + 1]) coef_xyz_mutable;

  if (orthorhombic && border_mask == 0) {
    // Here we ignore bounds_owned and always collocate the entire cube,
    // thereby assuming that the cube fits into the local grid.
    grid_library_gather_stats((grid_library_stats){.ref_collocate_ortho = 1});
    collocate_ortho(lp, zetp, coef_xyz, dh, dh_inv, rp, npts_global, npts_local,
                    shift_local, radius, grid);
  } else {
    grid_library_gather_stats((grid_library_stats){.ref_collocate_general = 1});
    collocate_general(border_mask, lp, zetp, coef_xyz, dh, dh_inv, rp,
                      npts_global, npts_local, shift_local, border_width,
                      radius, grid);
  }
}

// EOF
