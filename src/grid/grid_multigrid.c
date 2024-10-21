/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2024 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>

#include "grid_multigrid.h"

/*******************************************************************************
 * \brief Allocates a multigrid which is passed to task list-based and
 *pgf_product-based routines.
 *
 * \param orthorhombic     Whether simulation box is orthorhombic.
 * \param nlevels          Number of grid levels.
 * \param npts_global     Number of global grid points in each direction.
 * \param npts_local      Number of local grid points in each direction.
 * \param shift_local     Number of points local grid is shifted wrt global grid
 * \param border_width    Width of halo region in grid points in each direction.
 * \param dh              Incremental grid matrix.
 * \param dh_inv          Inverse incremental grid matrix.
 *
 * \param multigrid        Handle to the created multigrid.
 *
 * \author Frederick Stein
 ******************************************************************************/
void grid_create_multigrid(
    const bool orthorhombic, const int nlevels,
    // const int npts_global[nlevels][3],
    // const int npts_local[nlevels][3], const int shift_local[nlevels][3],
    // const int border_width[nlevels][3], const double dh[nlevels][3][3],
    // const double dh_inv[nlevels][3][3],
    grid_multigrid **multigrid_out) {
  grid_multigrid *multigrid = NULL;

  assert(multigrid_out != NULL);

  if (*multigrid_out != NULL) {
    multigrid = *multigrid_out;
  } else {
    multigrid = calloc(1, sizeof(grid_multigrid));
  }

  multigrid->nlevels = nlevels;
  multigrid->orthorhombic = orthorhombic;

  *multigrid_out = multigrid;
}

/*******************************************************************************
 * \brief Deallocates given multigrid.
 * \author Frederick Stein
 ******************************************************************************/
void grid_free_multigrid(grid_multigrid *multigrid) {
  if (multigrid != NULL) {
    free(multigrid);
  }
}

// EOF
