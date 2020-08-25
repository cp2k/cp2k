/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2020  CP2K developers group                         *
 *****************************************************************************/

#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "grid_collocate.h"
#include "grid_collocate_replay.h"
#include "ref/grid_ref_collocate.h"

//******************************************************************************
// \brief Public entry point. A thin wrapper with the only purpose of calling
//        grid_collocate_record when DUMP_TASKS = true.
// \author Ole Schuett
//******************************************************************************
void grid_collocate_pgf_product(
    const bool orthorhombic, const int border_mask, const int func,
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

  // For collocating single pgf products there's only the ref. implementation.
  grid_ref_collocate_pgf_product(
      orthorhombic, border_mask, func, la_max, la_min, lb_max, lb_min, zeta,
      zetb, rscale, dh, dh_inv, ra, rab, npts_global, npts_local, shift_local,
      border_width, radius, o1, o2, n1, n2, pab, grid);

  if (DUMP_TASKS) {
    grid_collocate_record(orthorhombic, border_mask, func, la_max, la_min,
                          lb_max, lb_min, zeta, zetb, rscale, dh, dh_inv, ra,
                          rab, npts_global, npts_local, shift_local,
                          border_width, radius, o1, o2, n1, n2, pab, grid);

    for (size_t i = 0; i < npts_local_total; i++) {
      grid[i] += grid_before[i];
    }
    free(grid_before);
  }
}

// EOF
