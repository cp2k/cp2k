/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2024 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#ifndef GRID_MULTIGRID_H
#define GRID_MULTIGRID_H

#include "common/grid_mpi.h"
#include "cpu/grid_cpu_multigrid.h"
#include "ref/grid_ref_multigrid.h"

#include <stdbool.h>

/*******************************************************************************
 * \brief Internal representation of a multigrid, abstracting various backends.
 * \author Frederick Stein
 ******************************************************************************/
typedef struct {
  int backend;
  bool orthorhombic;
  int nlevels;
  int (*npts_global)[3];
  int (*npts_local)[3];
  int (*shift_local)[3];
  int (*border_width)[3];
  double (*dh)[3][3];
  double (*dh_inv)[3][3];
  double **grids;
  grid_mpi_comm comm;
  grid_ref_multigrid *ref;
  grid_cpu_multigrid *cpu;
  // more backends to be added here
} grid_multigrid;

/*******************************************************************************
 * \brief Wrapper araound grid_create_multigrid to deal with calls from Fortran.
 *        Compare grid_create_multigrid for information on the other arguments.
 *
 * \param fortran_comm     The integer Fortran handle
 *
 * \author Frederick Stein
 ******************************************************************************/
void grid_create_multigrid_f(
    const bool orthorhombic, const int nlevels,
    const int npts_global[nlevels][3], const int npts_local[nlevels][3],
    const int shift_local[nlevels][3], const int border_width[nlevels][3],
    const double dh[nlevels][3][3], const double dh_inv[nlevels][3][3],
    const grid_mpi_fint fortran_comm, grid_multigrid **multigrid_out);

bool grid_get_multigrid_orthorhombic(const grid_multigrid *multigrid);

int grid_get_multigrid_nlevels(const grid_multigrid *multigrid);

void grid_get_multigrid_npts_global(const grid_multigrid *multigrid,
                                    int *nlevel, int **npts_global);

void grid_get_multigrid_npts_local(const grid_multigrid *multigrid, int *nlevel,
                                   int **npts_local);

void grid_get_multigrid_shift_local(const grid_multigrid *multigrid,
                                    int *nlevel, int **shift_local);

void grid_get_multigrid_border_width(const grid_multigrid *multigrid,
                                     int *nlevel, int **border_width);

void grid_get_multigrid_dh(const grid_multigrid *multigrid, int *nlevel,
                           double **dh);

void grid_get_multigrid_dh_inv(const grid_multigrid *multigrid, int *nlevel,
                               double **dh_inv);

grid_mpi_fint grid_get_multigrid_fortran_comm(const grid_multigrid *multigrid);

grid_mpi_comm grid_get_multigrid_comm(const grid_multigrid *multigrid);

/*******************************************************************************
 * \brief Allocates a multigrid which is passed to task list-based and
 *        pgf_product-based routines.
 *
 * \param orthorhombic     Whether simulation box is orthorhombic.
 * \param nlevels          Number of grid levels.
 * \param npts_global      Number of global grid points in each direction.
 * \param npts_local       Number of local grid points in each direction.
 * \param shift_local      Number of points local grid is shifted wrt global
 *                         grid
 * \param border_width     Width of halo region in grid points in each
 *                         direction.
 * \param dh               Incremental grid matrix.
 * \param dh_inv           Inverse incremental grid matrix.
 * \param comm             The C-communicator
 *
 * \param multigrid        Handle to the created multigrid.
 *
 * \author Frederick Stein
 ******************************************************************************/
void grid_create_multigrid(
    const bool orthorhombic, const int nlevels,
    const int npts_global[nlevels][3], const int npts_local[nlevels][3],
    const int shift_local[nlevels][3], const int border_width[nlevels][3],
    const double dh[nlevels][3][3], const double dh_inv[nlevels][3][3],
    const grid_mpi_comm comm, grid_multigrid **multigrid_out);

/*******************************************************************************
 * \brief Deallocates given multigrid.
 * \author Frederick Stein
 ******************************************************************************/
void grid_free_multigrid(grid_multigrid *multigrid);

#endif

// EOF
