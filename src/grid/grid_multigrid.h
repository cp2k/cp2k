/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2024 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#ifndef GRID_MULTIGRID_H
#define GRID_MULTIGRID_H

#include <stdbool.h>

/*******************************************************************************
 * \brief Internal representation of a task list, abstracting various backends.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  bool orthorhombic;
  int nlevels;
  int (*npts_global)[3];
  int (*npts_local)[3];
  int (*shift_local)[3];
  int (*border_width)[3];
  double (*dh)[3][3];
  double (*dh_inv)[3][3];
  // more backends to be added here
} grid_multigrid;

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
void grid_create_multigrid(const bool orthorhombic, const int nlevels,
                           const int npts_global[nlevels][3],
                           const int npts_local[nlevels][3],
                           const int shift_local[nlevels][3],
                           const int border_width[nlevels][3],
                           const double dh[nlevels][3][3],
                           const double dh_inv[nlevels][3][3],
                           grid_multigrid **multigrid_out);

/*******************************************************************************
 * \brief Deallocates given multigrid.
 * \author Frederick Stein
 ******************************************************************************/
void grid_free_multigrid(grid_multigrid *multigrid);

#endif

// EOF
