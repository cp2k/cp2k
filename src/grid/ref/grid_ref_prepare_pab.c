/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "../common/grid_common.h"
#include "../common/grid_prepare_pab.h"
#include "grid_ref_prepare_pab.h"

/*******************************************************************************
 * \brief Adds given value to matrix element cab[idx(b)][idx(a)].
 * \author Ole Schuett
 ******************************************************************************/
static inline void prep_term(const orbital a, const orbital b,
                             const double value, const int n, double *cab) {
  cab[idx(b) * n + idx(a)] += value;
}

/*******************************************************************************
 * \brief Returns block size changes due to transformation grid_prepare_pab.
 * \author Ole Schuett
 ******************************************************************************/
void grid_ref_prepare_get_ldiffs(const int func, int *la_min_diff,
                                 int *la_max_diff, int *lb_min_diff,
                                 int *lb_max_diff) {
  const prepare_ldiffs ldiffs = prepare_get_ldiffs(func);
  *la_min_diff = ldiffs.la_min_diff;
  *la_max_diff = ldiffs.la_max_diff;
  *lb_min_diff = ldiffs.lb_min_diff;
  *lb_max_diff = ldiffs.lb_max_diff;
}

/*******************************************************************************
 * \brief Selects and transforms a sub-block of the given density matrix block.
 *        See grid_ref_prepare_pab.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_ref_prepare_pab(const int func, const int o1, const int o2,
                          const int la_max, const int la_min, const int lb_max,
                          const int lb_min, const double zeta,
                          const double zetb, const int n1, const int n2,
                          const double pab[n2][n1], const int n1_prep,
                          const int n2_prep,
                          double pab_prep[n2_prep][n1_prep]) {

  for (int lxa = 0; lxa <= la_max; lxa++) {
    for (int lxb = 0; lxb <= lb_max; lxb++) {
      for (int lya = 0; lya <= la_max - lxa; lya++) {
        for (int lyb = 0; lyb <= lb_max - lxb; lyb++) {
          const int lza_start = imax(la_min - lxa - lya, 0);
          for (int lza = lza_start; lza <= la_max - lxa - lya; lza++) {
            const int lzb_start = imax(lb_min - lxb - lyb, 0);
            for (int lzb = lzb_start; lzb <= lb_max - lxb - lyb; lzb++) {
              const orbital a = {{lxa, lya, lza}};
              const orbital b = {{lxb, lyb, lzb}};
              const double pab_val = pab[o2 + idx(b)][o1 + idx(a)];
              prepare_pab(func, a, b, zeta, zetb, pab_val, n1_prep,
                          (double *)pab_prep);
            }
          }
        }
      }
    }
  }
}

// EOF
