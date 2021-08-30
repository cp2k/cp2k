/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2021 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define GRID_DO_COLLOCATE 1
#include "../common/grid_common.h"
#include "grid_ref_collint.h"
#include "grid_ref_collocate.h"
#include "grid_ref_integrate.h"
#include "grid_ref_prepare_pab.h"

/*******************************************************************************
 * \brief Collocates a single product of primitiv Gaussians.
 *        See grid_ref_collocate.h for details.
 * \author Ole Schuett
 ******************************************************************************/
static void collocate_internal(
    const bool orthorhombic, const int border_mask, const enum grid_func func,
    const int la_max, const int la_min, const int lb_max, const int lb_min,
    const double zeta, const double zetb, const double rscale,
    const double dh[3][3], const double dh_inv[3][3], const double ra[3],
    const double rab[3], const int npts_global[3], const int npts_local[3],
    const int shift_local[3], const int border_width[3], const double radius,
    const int o1, const int o2, const int n1, const int n2,
    const double pab[n2][n1], double *grid) {

  int la_min_diff, la_max_diff, lb_min_diff, lb_max_diff;
  grid_ref_prepare_get_ldiffs(func, &la_min_diff, &la_max_diff, &lb_min_diff,
                              &lb_max_diff);

  const int la_min_cab = imax(la_min + la_min_diff, 0);
  const int lb_min_cab = imax(lb_min + lb_min_diff, 0);
  const int la_max_cab = la_max + la_max_diff;
  const int lb_max_cab = lb_max + lb_max_diff;
  const int n1_cab = ncoset(la_max_cab);
  const int n2_cab = ncoset(lb_max_cab);

  const size_t cab_size = n2_cab * n1_cab;
  double cab[cab_size];
  memset(cab, 0, cab_size * sizeof(double));

  grid_ref_prepare_pab(func, o1, o2, la_max, la_min, lb_max, lb_min, zeta, zetb,
                       n1, n2, pab, n1_cab, n2_cab, (double(*)[n1_cab])cab);
  cab_to_grid(orthorhombic, border_mask, la_max_cab, la_min_cab, lb_max_cab,
              lb_min_cab, zeta, zetb, rscale, dh, dh_inv, ra, rab, npts_global,
              npts_local, shift_local, border_width, radius, cab, grid);
}

/*******************************************************************************
 * \brief Writes the given arguments into a .task file.
 *        See grid_collocate_replay.h for details.
 * \author Ole Schuett
 ******************************************************************************/
static void record_collocate(
    const bool orthorhombic, const int border_mask, const enum grid_func func,
    const int la_max, const int la_min, const int lb_max, const int lb_min,
    const double zeta, const double zetb, const double rscale,
    const double dh[3][3], const double dh_inv[3][3], const double ra[3],
    const double rab[3], const int npts_global[3], const int npts_local[3],
    const int shift_local[3], const int border_width[3], const double radius,
    const int o1, const int o2, const int n1, const int n2,
    const double pab[n2][n1], const double *grid) {

  static int counter = 1;
  int my_number;

#pragma omp critical
  my_number = counter++;

  char filename[100];
  snprintf(filename, sizeof(filename), "grid_collocate_%05i.task", my_number);

  const int D = DECIMAL_DIG; // In C11 we could use DBL_DECIMAL_DIG.
  FILE *fp = fopen(filename, "w+");
  fprintf(fp, "#Grid task v10\n");
  fprintf(fp, "orthorhombic %i\n", orthorhombic);
  fprintf(fp, "border_mask %i\n", border_mask);
  fprintf(fp, "func %i\n", func);
  fprintf(fp, "la_max %i\n", la_max);
  fprintf(fp, "la_min %i\n", la_min);
  fprintf(fp, "lb_max %i\n", lb_max);
  fprintf(fp, "lb_min %i\n", lb_min);
  fprintf(fp, "zeta %.*e\n", D, zeta);
  fprintf(fp, "zetb %.*e\n", D, zetb);
  fprintf(fp, "rscale %.*e\n", D, rscale);
  for (int i = 0; i < 3; i++)
    fprintf(fp, "dh %i %.*e %.*e %.*e\n", i, D, dh[i][0], D, dh[i][1], D,
            dh[i][2]);
  for (int i = 0; i < 3; i++)
    fprintf(fp, "dh_inv %i %.*e %.*e %.*e\n", i, D, dh_inv[i][0], D,
            dh_inv[i][1], D, dh_inv[i][2]);
  fprintf(fp, "ra %.*e %.*e %.*e\n", D, ra[0], D, ra[1], D, ra[2]);
  fprintf(fp, "rab %.*e %.*e %.*e\n", D, rab[0], D, rab[1], D, rab[2]);
  fprintf(fp, "npts_global %i %i %i\n", npts_global[0], npts_global[1],
          npts_global[2]);
  fprintf(fp, "npts_local %i %i %i\n", npts_local[0], npts_local[1],
          npts_local[2]);
  fprintf(fp, "shift_local %i %i %i\n", shift_local[0], shift_local[1],
          shift_local[2]);
  fprintf(fp, "border_width %i %i %i\n", border_width[0], border_width[1],
          border_width[2]);
  fprintf(fp, "radius %.*e\n", D, radius);
  fprintf(fp, "o1 %i\n", o1);
  fprintf(fp, "o2 %i\n", o2);
  fprintf(fp, "n1 %i\n", n1);
  fprintf(fp, "n2 %i\n", n2);

  for (int i = 0; i < n2; i++) {
    for (int j = 0; j < n1; j++) {
      fprintf(fp, "pab %i %i %.*e\n", i, j, D, pab[i][j]);
    }
  }

  const int npts_local_total = npts_local[0] * npts_local[1] * npts_local[2];

  int ngrid_nonzero = 0;
  for (int i = 0; i < npts_local_total; i++) {
    if (grid[i] != 0.0) {
      ngrid_nonzero++;
    }
  }
  fprintf(fp, "ngrid_nonzero %i\n", ngrid_nonzero);

  for (int k = 0; k < npts_local[2]; k++) {
    for (int j = 0; j < npts_local[1]; j++) {
      for (int i = 0; i < npts_local[0]; i++) {
        const double val =
            grid[k * npts_local[1] * npts_local[0] + j * npts_local[0] + i];
        if (val != 0.0) {
          fprintf(fp, "grid %i %i %i %.*e\n", i, j, k, D, val);
        }
      }
    }
  }

  double hab[n2][n1];
  double forces[2][3] = {0};
  double virials[2][3][3] = {0};
  memset(hab, 0, n2 * n1 * sizeof(double));

  const bool compute_tau = (func == GRID_FUNC_DADB);

  grid_ref_integrate_pgf_product(orthorhombic, compute_tau, border_mask, la_max,
                                 la_min, lb_max, lb_min, zeta, zetb, dh, dh_inv,
                                 ra, rab, npts_global, npts_local, shift_local,
                                 border_width, radius, o1, o2, n1, n2, grid,
                                 hab, pab, forces, virials, NULL, NULL);

  for (int i = o2; i < ncoset(lb_max) + o2; i++) {
    for (int j = o1; j < ncoset(la_max) + o1; j++) {
      fprintf(fp, "hab %i %i %.*e\n", i, j, D, hab[i][j]);
    }
  }
  fprintf(fp, "force_a %.*e %.*e %.*e\n", D, forces[0][0], D, forces[0][1], D,
          forces[0][2]);
  fprintf(fp, "force_b %.*e %.*e %.*e\n", D, forces[1][0], D, forces[1][1], D,
          forces[1][2]);
  for (int i = 0; i < 3; i++)
    fprintf(fp, "virial %i %.*e %.*e %.*e\n", i, D,
            virials[0][i][0] + virials[1][i][0], D,
            virials[0][i][1] + virials[1][i][1], D,
            virials[0][i][2] + virials[1][i][2]);

  fprintf(fp, "#THE_END\n");
  fclose(fp);
  printf("Wrote %s\n", filename);
}

/*******************************************************************************
 * \brief Public entry point. A thin wrapper with the only purpose of calling
 *        record_collocate when DUMP_TASKS = true.
 * \author Ole Schuett
 ******************************************************************************/
void grid_ref_collocate_pgf_product(
    const bool orthorhombic, const int border_mask, const enum grid_func func,
    const int la_max, const int la_min, const int lb_max, const int lb_min,
    const double zeta, const double zetb, const double rscale,
    const double dh[3][3], const double dh_inv[3][3], const double ra[3],
    const double rab[3], const int npts_global[3], const int npts_local[3],
    const int shift_local[3], const int border_width[3], const double radius,
    const int o1, const int o2, const int n1, const int n2,
    const double pab[n2][n1], double *grid) {

  // Set this to true to write each task to a file.
  const bool DUMP_TASKS = false;

  double *grid_before = NULL;
  const size_t npts_local_total = npts_local[0] * npts_local[1] * npts_local[2];

  if (DUMP_TASKS) {
    const size_t sizeof_grid = sizeof(double) * npts_local_total;
    grid_before = malloc(sizeof_grid);
    memcpy(grid_before, grid, sizeof_grid);
    memset(grid, 0, sizeof_grid);
  }

  collocate_internal(orthorhombic, border_mask, func, la_max, la_min, lb_max,
                     lb_min, zeta, zetb, rscale, dh, dh_inv, ra, rab,
                     npts_global, npts_local, shift_local, border_width, radius,
                     o1, o2, n1, n2, pab, grid);

  if (DUMP_TASKS) {
    record_collocate(orthorhombic, border_mask, func, la_max, la_min, lb_max,
                     lb_min, zeta, zetb, rscale, dh, dh_inv, ra, rab,
                     npts_global, npts_local, shift_local, border_width, radius,
                     o1, o2, n1, n2, pab, grid);

    for (size_t i = 0; i < npts_local_total; i++) {
      grid[i] += grid_before[i];
    }
    free(grid_before);
  }
}

// EOF
