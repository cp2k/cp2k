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

#define GRID_DO_COLLOCATE 0
#include "../common/grid_common.h"
#include "grid_ref_collint.h"
#include "grid_ref_integrate.h"

/*******************************************************************************
 * \brief Orbital angular momentum.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  int l[3];
} orbital;

/*******************************************************************************
 * \brief Increase i'th component of given orbital angular momentum.
 * \author Ole Schuett
 ******************************************************************************/
static inline orbital up(const int i, const orbital a) {
  orbital b = a;
  b.l[i] += 1;
  return b;
}

/*******************************************************************************
 * \brief Decrease i'th component of given orbital angular momentum.
 * \author Ole Schuett
 ******************************************************************************/
static inline orbital down(const int i, const orbital a) {
  orbital b = a;
  b.l[i] = imax(0, a.l[i] - 1);
  return b;
}

/*******************************************************************************
 * \brief Return coset index of given orbital angular momentum.
 * \author Ole Schuett
 ******************************************************************************/
static inline int index(const orbital a) {
  return coset(a.l[0], a.l[1], a.l[2]);
}

/*******************************************************************************
 * \brief Contracts given matrix elements to obtain the forces for atom a.
 * \author Ole Schuett
 ******************************************************************************/
static inline void update_force_a(const orbital a, const orbital b,
                                  const double pab, const double ftza,
                                  const int m1, const int m2,
                                  const double vab[m2][m1], double force_a[3]) {

  for (int i = 0; i < 3; i++) {
    const double aip1 = vab[index(b)][index(up(i, a))];
    const double aim1 = vab[index(b)][index(down(i, a))];
    force_a[i] += pab * (ftza * aip1 - a.l[i] * aim1);
  }
}

/*******************************************************************************
 * \brief Contracts given matrix elements to obtain the forces for atom b.
 * \author Ole Schuett
 ******************************************************************************/
static inline void update_force_b(const orbital a, const orbital b,
                                  const double pab, const double ftzb,
                                  const double rab[3], const int m1,
                                  const int m2, const double vab[m2][m1],
                                  double force_b[3]) {

  const double axpm0 = vab[index(b)][index(a)];
  for (int i = 0; i < 3; i++) {
    const double aip1 = vab[index(b)][index(up(i, a))];
    const double bim1 = vab[index(down(i, b))][index(a)];
    force_b[i] += pab * (ftzb * (aip1 - rab[i] * axpm0) - b.l[i] * bim1);
  }
}

/*******************************************************************************
 * \brief Contracts given matrix elements to obtain the virial for atom a.
 * \author Ole Schuett
 ******************************************************************************/
static inline void update_virial_a(const orbital a, const orbital b,
                                   const double pab, const double ftza,
                                   const int m1, const int m2,
                                   const double vab[m2][m1],
                                   double virial_a[3][3]) {

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      virial_a[i][j] += pab * ftza * vab[index(b)][index(up(i, up(j, a)))] -
                        pab * a.l[j] * vab[index(b)][index(up(i, down(j, a)))];
    }
  }
}

/*******************************************************************************
 * \brief Contracts given matrix elements to obtain the virial for atom b.
 * \author Ole Schuett
 ******************************************************************************/
static inline void update_virial_b(const orbital a, const orbital b,
                                   const double pab, const double ftzb,
                                   const double rab[3], const int m1,
                                   const int m2, const double vab[m2][m1],
                                   double virial_b[3][3]) {

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      virial_b[i][j] += pab * ftzb *
                            (vab[index(b)][index(up(i, up(j, a)))] -
                             vab[index(b)][index(up(i, a))] * rab[j] -
                             vab[index(b)][index(up(j, a))] * rab[i] +
                             vab[index(b)][index(a)] * rab[j] * rab[i]) -
                        pab * b.l[j] * vab[index(up(i, down(j, b)))][index(a)];
    }
  }
}

/*******************************************************************************
 * \brief Contracts given matrix elements to obtain forces and virials.
 * \author Ole Schuett
 ******************************************************************************/
static void update_all(const orbital a, const orbital b, const double f,
                       const double ftza, const double ftzb,
                       const double rab[3], const int m1, const int m2,
                       const double vab[m2][m1], const double pab, double *hab,
                       double forces[2][3], double virials[2][3][3]) {

  *hab += f * vab[index(b)][index(a)];

  if (forces != NULL) {
    update_force_a(a, b, f * pab, ftza, m1, m2, vab, forces[0]);
    update_force_b(a, b, f * pab, ftzb, rab, m1, m2, vab, forces[1]);
  }

  if (virials != NULL) {
    update_virial_a(a, b, f * pab, ftza, m1, m2, vab, virials[0]);
    update_virial_b(a, b, f * pab, ftzb, rab, m1, m2, vab, virials[1]);
  }
}

/*******************************************************************************
 * \brief Contracts given matrix elements to obtain forces and virials for tau.
 * \author Ole Schuett
 ******************************************************************************/
static void update_tau(const orbital a, const orbital b, const double ftza,
                       const double ftzb, const double rab[3], const int m1,
                       const int m2, const double vab[m2][m1], const double pab,
                       double *hab, double forces[2][3],
                       double virials[2][3][3]) {

  for (int i = 0; i < 3; i++) {
    update_all(down(i, a), down(i, b), 0.5 * a.l[i] * b.l[i], ftza, ftzb, rab,
               m1, m2, vab, pab, hab, forces, virials);
    update_all(up(i, a), down(i, b), -0.5 * ftza * b.l[i], ftza, ftzb, rab, m1,
               m2, vab, pab, hab, forces, virials);
    update_all(down(i, a), up(i, b), -0.5 * a.l[i] * ftzb, ftza, ftzb, rab, m1,
               m2, vab, pab, hab, forces, virials);
    update_all(up(i, a), up(i, b), 0.5 * ftza * ftzb, ftza, ftzb, rab, m1, m2,
               vab, pab, hab, forces, virials);
  }
}

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
    double hdab[n2][n1][3], double a_hdab[n2][n1][3][3]) {

  int la_max_local = la_max;
  int la_min_local = la_min;
  int lb_min_local = lb_min;
  int lb_max_local = lb_max;
  if (forces != NULL || hdab != NULL || virials != NULL || a_hdab != NULL) {
    la_max_local += 1; // for deriv. of gaussian, unimportant which one
    la_min_local -= 1;
    lb_min_local -= 1;
  }
  if (virials != NULL || a_hdab != NULL) {
    la_max_local += 1;
    lb_max_local += 1;
  }
  if (compute_tau) {
    la_max_local += 1;
    lb_max_local += 1;
    la_min_local -= 1;
    lb_min_local -= 1;
  }
  la_min_local = imax(la_min_local, 0);
  lb_min_local = imax(lb_min_local, 0);

  const int m1 = ncoset[la_max_local];
  const int m2 = ncoset[lb_max_local];
  double vab_mutable[m2 * m1];
  memset(vab_mutable, 0, m2 * m1 * sizeof(double));

  const double rscale = 1.0; // TODO: remove rscale from cab_to_grid
  cab_to_grid(orthorhombic, border_mask, la_max_local, la_min_local,
              lb_max_local, lb_min_local, zeta, zetb, rscale, dh, dh_inv, ra,
              rab, npts_global, npts_local, shift_local, border_width, radius,
              vab_mutable, grid);

  //  vab contains all the information needed to find the elements of hab
  //  and optionally of derivatives of these elements
  const double(*vab)[m1] = (const double(*)[m1])vab_mutable;

  const double ftza = 2.0 * zeta;
  const double ftzb = 2.0 * zetb;

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

              double *habval = &hab[o2 + index(b)][o1 + index(a)];
              const double pabval =
                  (pab == NULL) ? 0.0 : pab[o2 + index(b)][o1 + index(a)];

              // Fill hab, forces, and virials.
              if (compute_tau) {
                update_tau(a, b, ftza, ftzb, rab, m1, m2, vab, pabval, habval,
                           forces, virials);
              } else {
                update_all(a, b, 1.0, ftza, ftzb, rab, m1, m2, vab, pabval,
                           habval, forces, virials);
              }

              // Fill hdab and a_hdab.
              if (hdab != NULL) {
                assert(!compute_tau);
                update_force_a(a, b, 1.0, ftza, m1, m2, vab,
                               hdab[o2 + index(b)][o1 + index(a)]);
              }
              if (a_hdab != NULL) {
                assert(!compute_tau);
                update_virial_a(a, b, 1.0, ftza, m1, m2, vab,
                                a_hdab[o2 + index(b)][o1 + index(a)]);
              }
            }
          }
        }
      }
    }
  }
}

// EOF
