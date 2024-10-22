/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2024 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include "grid_ref_multigrid.h"

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

/*******************************************************************************
 * \brief Allocates a multigrid which is passed to task list-based and
 *        pgf_product-based routines.
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
void grid_ref_create_singlegrid(
    const bool orthorhombic, const int npts_global[3], const int npts_local[3],
    const int shift_local[3], const int border_width[3], const double dh[3][3],
    const double dh_inv[3][3], const grid_mpi_comm comm,
    grid_ref_singlegrid **singlegrid_out) {
  grid_ref_singlegrid *singlegrid = NULL;

  assert(singlegrid_out != NULL);

  if (*singlegrid_out != NULL) {
    singlegrid = *singlegrid_out;
    // Always free the old communicator
    grid_mpi_comm_free(&singlegrid->comm);
  } else {
    singlegrid = calloc(1, sizeof(grid_ref_singlegrid));
  }

  singlegrid->orthorhombic = orthorhombic;
  memcpy(singlegrid->npts_global, npts_global, 3 * sizeof(int));
  memcpy(singlegrid->npts_local, npts_local, 3 * sizeof(int));
  memcpy(singlegrid->shift_local, shift_local, 3 * sizeof(int));
  memcpy(singlegrid->border_width, border_width, 3 * sizeof(int));
  memcpy(singlegrid->dh, dh, 9 * sizeof(double));
  memcpy(singlegrid->dh_inv, dh_inv, 9 * sizeof(double));
  grid_mpi_comm_dup(comm, &singlegrid->comm);

  *singlegrid_out = singlegrid;
}

/*******************************************************************************
 * \brief Deallocates given singlegrid.
 * \author Frederick Stein
 ******************************************************************************/
void grid_ref_free_singlegrid(grid_ref_singlegrid *singlegrid) {
  if (singlegrid != NULL) {
    grid_mpi_comm_free(&singlegrid->comm);
    free(singlegrid);
  }
}

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
void grid_ref_create_multigrid(
    const bool orthorhombic, const int nlevels,
    const int npts_global[nlevels][3], const int npts_local[nlevels][3],
    const int shift_local[nlevels][3], const int border_width[nlevels][3],
    const double dh[nlevels][3][3], const double dh_inv[nlevels][3][3],
    const grid_mpi_comm comm, grid_ref_multigrid **multigrid_out) {
  grid_ref_multigrid *multigrid = NULL;

  assert(multigrid_out != NULL);

  if (*multigrid_out != NULL) {
    multigrid = *multigrid_out;
    if (nlevels != multigrid->nlevels) {
      multigrid->singlegrids = realloc(multigrid->singlegrids,
                                       nlevels * sizeof(grid_ref_singlegrid *));
    }
  } else {
    multigrid = calloc(1, sizeof(grid_ref_multigrid));
    multigrid->singlegrids = calloc(nlevels, sizeof(grid_ref_singlegrid *));
  }

  multigrid->nlevels = nlevels;
  multigrid->orthorhombic = orthorhombic;

  for (int grid = 0; grid < nlevels; grid++) {
    grid_ref_create_singlegrid(orthorhombic, &npts_global[grid][0],
                               &npts_local[grid][0], &shift_local[grid][0],
                               &border_width[grid][0],
                               (const double(*)[3])(&dh[grid][0][0]),
                               (const double(*)[3])(&dh_inv[grid][0][0]), comm,
                               &(multigrid->singlegrids[grid]));
  }

  *multigrid_out = multigrid;
}

/*******************************************************************************
 * \brief Deallocates given multigrid.
 * \author Frederick Stein
 ******************************************************************************/
void grid_ref_free_multigrid(grid_ref_multigrid *multigrid) {
  if (multigrid != NULL) {
    for (int grid = 0; grid < multigrid->nlevels; grid++) {
      grid_ref_free_singlegrid(multigrid->singlegrids[grid]);
    }
    free(multigrid->singlegrids);
    free(multigrid);
  }
}

// EOF
