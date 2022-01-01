/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define GRID_DO_COLLOCATE 0
#include "../common/grid_common.h"
#include "../common/grid_process_vab.h"
#include "grid_ref_collint.h"
#include "grid_ref_integrate.h"

/*******************************************************************************
 * \brief Integrates a single task. See grid_ref_integrate.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_ref_integrate_pgf_product(
    const bool orthorhombic, const bool compute_tau, const int border_mask,
    const int la_max, const int la_min, const int lb_max, const int lb_min,
    const double zeta, const double zetb, const double dh[3][3],
    const double dh_inv[3][3], const double ra[3], const double rab[3],
    const int npts_global[3], const int npts_local[3], const int shift_local[3],
    const int border_width[3], const double radius, const int o1, const int o2,
    const int n1, const int n2, const double *grid, double hab[n2][n1],
    const double pab[n2][n1], double forces[2][3], double virials[2][3][3],
    double hdab[n2][n1][3], double hadb[n2][n1][3],
    double a_hdab[n2][n1][3][3]) {

  const bool calculate_forces = (forces != NULL || hdab != NULL);
  const bool calculate_virial = (virials != NULL || a_hdab != NULL);
  assert(!calculate_virial || calculate_forces);
  const process_ldiffs ldiffs =
      process_get_ldiffs(calculate_forces, calculate_virial, compute_tau);
  int la_max_local = la_max + ldiffs.la_max_diff;
  int lb_max_local = lb_max + ldiffs.lb_max_diff;
  int la_min_local = imax(0, la_min + ldiffs.la_min_diff);
  int lb_min_local = imax(0, lb_min + ldiffs.lb_min_diff);

  const int m1 = ncoset(la_max_local);
  const int m2 = ncoset(lb_max_local);
  double cab[m2 * m1];
  memset(cab, 0, m2 * m1 * sizeof(double));

  const double rscale = 1.0; // TODO: remove rscale from cab_to_grid
  cab_to_grid(orthorhombic, border_mask, la_max_local, la_min_local,
              lb_max_local, lb_min_local, zeta, zetb, rscale, dh, dh_inv, ra,
              rab, npts_global, npts_local, shift_local, border_width, radius,
              cab, grid);

  //  cab contains all the information needed to find the elements of hab
  //  and optionally of derivatives of these elements
  for (int la = la_min; la <= la_max; la++) {
    for (int ax = 0; ax <= la; ax++) {
      for (int ay = 0; ay <= la - ax; ay++) {
        const int az = la - ax - ay;
        const orbital a = {{ax, ay, az}};
        for (int lb = lb_min; lb <= lb_max; lb++) {
          for (int bx = 0; bx <= lb; bx++) {
            for (int by = 0; by <= lb - bx; by++) {
              const int bz = lb - bx - by;
              const orbital b = {{bx, by, bz}};

              // Update hab block.
              hab[o2 + idx(b)][o1 + idx(a)] +=
                  get_hab(a, b, zeta, zetb, m1, cab, compute_tau);

              // Update forces.
              if (forces != NULL) {
                const double pabval = pab[o2 + idx(b)][o1 + idx(a)];
                for (int i = 0; i < 3; i++) {
                  forces[0][i] += pabval * get_force_a(a, b, i, zeta, zetb, m1,
                                                       cab, compute_tau);
                  forces[1][i] += pabval * get_force_b(a, b, i, zeta, zetb, rab,
                                                       m1, cab, compute_tau);
                }
              }

              // Update virials.
              if (virials != NULL) {
                const double pabval = pab[o2 + idx(b)][o1 + idx(a)];
                for (int i = 0; i < 3; i++) {
                  for (int j = 0; j < 3; j++) {
                    virials[0][i][j] +=
                        pabval * get_virial_a(a, b, i, j, zeta, zetb, m1, cab,
                                              compute_tau);
                    virials[1][i][j] +=
                        pabval * get_virial_b(a, b, i, j, zeta, zetb, rab, m1,
                                              cab, compute_tau);
                  }
                }
              }

              // Update hdab, hadb, and a_hdab (not used in batch mode).
              if (hdab != NULL) {
                assert(!compute_tau);
                for (int i = 0; i < 3; i++) {
                  hdab[o2 + idx(b)][o1 + idx(a)][i] +=
                      get_force_a(a, b, i, zeta, zetb, m1, cab, false);
                }
              }
              if (hadb != NULL) {
                assert(!compute_tau);
                for (int i = 0; i < 3; i++) {
                  hadb[o2 + idx(b)][o1 + idx(a)][i] +=
                      get_force_b(a, b, i, zeta, zetb, rab, m1, cab, false);
                }
              }
              if (a_hdab != NULL) {
                assert(!compute_tau);
                for (int i = 0; i < 3; i++) {
                  for (int j = 0; j < 3; j++) {
                    a_hdab[o2 + idx(b)][o1 + idx(a)][i][j] +=
                        get_virial_a(a, b, i, j, zeta, zetb, m1, cab, false);
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

// EOF
