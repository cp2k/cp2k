/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2024 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include "grid_multigrid.h"
#include "common/grid_library.h"

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

bool grid_get_multigrid_orthorhombic(const grid_multigrid *multigrid) {
  assert(multigrid != NULL);
  return multigrid->orthorhombic;
}

int grid_get_multigrid_nlevels(const grid_multigrid *multigrid) {
  assert(multigrid != NULL);
  return multigrid->nlevels;
}

void grid_get_multigrid_npts_global(const grid_multigrid *multigrid,
                                    int *nlevels, int **npts_global) {
  assert(multigrid != NULL);
  *nlevels = multigrid->nlevels;
  *npts_global = (int *)multigrid->npts_global;
}

void grid_get_multigrid_npts_local(const grid_multigrid *multigrid,
                                   int *nlevels, int **npts_local) {
  assert(multigrid != NULL);
  *nlevels = multigrid->nlevels;
  *npts_local = (int *)multigrid->npts_local;
}

void grid_get_multigrid_shift_local(const grid_multigrid *multigrid,
                                    int *nlevels, int **shift_local) {
  assert(multigrid != NULL);
  *nlevels = multigrid->nlevels;
  *shift_local = (int *)multigrid->shift_local;
}

void grid_get_multigrid_border_width(const grid_multigrid *multigrid,
                                     int *nlevels, int **border_width) {
  assert(multigrid != NULL);
  *nlevels = multigrid->nlevels;
  *border_width = (int *)multigrid->border_width;
}

void grid_get_multigrid_dh(const grid_multigrid *multigrid, int *nlevels,
                           double **dh) {
  assert(multigrid != NULL);
  *nlevels = multigrid->nlevels;
  *dh = (double *)multigrid->dh;
}

void grid_get_multigrid_dh_inv(const grid_multigrid *multigrid, int *nlevels,
                               double **dh_inv) {
  assert(multigrid != NULL);
  *nlevels = multigrid->nlevels;
  *dh_inv = (double *)multigrid->dh_inv;
}

grid_mpi_fint grid_get_multigrid_fortran_comm(const grid_multigrid *multigrid) {
  assert(multigrid != NULL);
  return grid_mpi_comm_c2f(multigrid->comm);
}

grid_mpi_comm grid_get_multigrid_comm(const grid_multigrid *multigrid) {
  assert(multigrid != NULL);
  return multigrid->comm;
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
void grid_create_multigrid_f(
    const bool orthorhombic, const int nlevels,
    const int npts_global[nlevels][3], const int npts_local[nlevels][3],
    const int shift_local[nlevels][3], const int border_width[nlevels][3],
    const double dh[nlevels][3][3], const double dh_inv[nlevels][3][3],
    const grid_mpi_fint fortran_comm, grid_multigrid **multigrid_out) {
  grid_create_multigrid(orthorhombic, nlevels, npts_global, npts_local,
                        shift_local, border_width, dh, dh_inv,
                        grid_mpi_comm_f2c(fortran_comm), multigrid_out);
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
void grid_create_multigrid(
    const bool orthorhombic, const int nlevels,
    const int npts_global[nlevels][3], const int npts_local[nlevels][3],
    const int shift_local[nlevels][3], const int border_width[nlevels][3],
    const double dh[nlevels][3][3], const double dh_inv[nlevels][3][3],
    const grid_mpi_comm comm, grid_multigrid **multigrid_out) {

  const grid_library_config config = grid_library_get_config();

  grid_multigrid *multigrid = NULL;

  assert(multigrid_out != NULL);

  const int num_int = 3 * nlevels;
  const int num_double = 9 * nlevels;

  if (*multigrid_out != NULL) {
    multigrid = *multigrid_out;
    if (nlevels != multigrid->nlevels) {
      multigrid->npts_global =
          realloc(multigrid->npts_global, num_int * sizeof(int));
      multigrid->npts_local =
          realloc(multigrid->npts_local, num_int * sizeof(int));
      multigrid->shift_local =
          realloc(multigrid->shift_local, num_int * sizeof(int));
      multigrid->border_width =
          realloc(multigrid->border_width, num_int * sizeof(int));
      multigrid->dh = realloc(multigrid->dh, num_double * sizeof(double));
      multigrid->dh_inv =
          realloc(multigrid->dh_inv, num_double * sizeof(double));

      for (int level = 0; level < multigrid->nlevels; level++) {
        free(multigrid->grids[level]);
      }
      multigrid->grids = realloc(multigrid->grids, nlevels * sizeof(double *));
      for (int level = 0; level < nlevels; level++) {
        multigrid->grids[level] =
            malloc(npts_local[level][0] * npts_local[level][1] *
                   npts_local[level][2] * sizeof(double));
      }
    } else {
      for (int level = 0; level < nlevels; level++) {
        multigrid->grids[level] =
            realloc(multigrid->grids[level],
                    npts_local[level][0] * npts_local[level][1] *
                        npts_local[level][2] * sizeof(double));
      }
    }
    // Always free the old communicator
    grid_mpi_comm_free(&multigrid->comm);
  } else {
    multigrid = calloc(1, sizeof(grid_multigrid));
    multigrid->npts_global = calloc(num_int, sizeof(int));
    multigrid->npts_local = calloc(num_int, sizeof(int));
    multigrid->shift_local = calloc(num_int, sizeof(int));
    multigrid->border_width = calloc(num_int, sizeof(int));
    multigrid->dh = calloc(num_double, sizeof(double));
    multigrid->dh_inv = calloc(num_double, sizeof(double));
    multigrid->grids = calloc(nlevels, sizeof(double *));

    // Resolve AUTO to a concrete backend.
    if (config.backend == GRID_BACKEND_AUTO) {
#if defined(__OFFLOAD_HIP) && !defined(__NO_OFFLOAD_GRID)
      multigrid->backend = GRID_BACKEND_HIP;
#elif defined(__OFFLOAD) && !defined(__NO_OFFLOAD_GRID)
      multigrid->backend = GRID_BACKEND_GPU;
#else
      multigrid->backend = GRID_BACKEND_CPU;
#endif
    } else {
      multigrid->backend = config.backend;
    }
    for (int level = 0; level < nlevels; level++) {
      multigrid->grids[level] =
          malloc(npts_local[level][0] * npts_local[level][1] *
                 npts_local[level][2] * sizeof(double));
    }
  }

  multigrid->nlevels = nlevels;
  multigrid->orthorhombic = orthorhombic;
  memcpy(multigrid->npts_global, npts_global, num_int * sizeof(int));
  memcpy(multigrid->npts_local, npts_local, num_int * sizeof(int));
  memcpy(multigrid->shift_local, shift_local, num_int * sizeof(int));
  memcpy(multigrid->border_width, border_width, num_int * sizeof(int));
  memcpy(multigrid->dh, dh, num_double * sizeof(double));
  memcpy(multigrid->dh_inv, dh_inv, num_double * sizeof(double));
  grid_mpi_comm_dup(comm, &multigrid->comm);

  grid_ref_create_multigrid(orthorhombic, nlevels, npts_global, npts_local,
                            shift_local, border_width, dh, dh_inv, comm,
                            &(multigrid->ref));

  // We need it for collocation/integration of pgf_products
  grid_cpu_create_multigrid(orthorhombic, nlevels, npts_global, npts_local,
                            shift_local, border_width, dh, dh_inv, comm,
                            &(multigrid->cpu));

  switch (multigrid->backend) {
  case GRID_BACKEND_REF:
    break;
  case GRID_BACKEND_CPU:
    break;
  }

  for (int level = 0; level < nlevels; level++) {
    memset(multigrid->grids[level], 0,
           npts_local[level][0] * npts_local[level][1] * npts_local[level][2] *
               sizeof(double));
  }

  *multigrid_out = multigrid;

  grid_mpi_barrier(multigrid->comm);
}

/*******************************************************************************
 * \brief Deallocates given multigrid.
 * \author Frederick Stein
 ******************************************************************************/
void grid_free_multigrid(grid_multigrid *multigrid) {
  if (multigrid != NULL) {
    if (multigrid->npts_global != NULL)
      free(multigrid->npts_global);
    if (multigrid->npts_local != NULL)
      free(multigrid->npts_local);
    if (multigrid->shift_local != NULL)
      free(multigrid->shift_local);
    if (multigrid->border_width != NULL)
      free(multigrid->border_width);
    if (multigrid->dh != NULL)
      free(multigrid->dh);
    if (multigrid->dh_inv != NULL)
      free(multigrid->dh_inv);
    if (multigrid->grids != NULL) {
      for (int level = 0; level < multigrid->nlevels; level++) {
        free(multigrid->grids[level]);
      }
      free(multigrid->grids);
    }
    grid_mpi_comm_free(&multigrid->comm);
    grid_ref_free_multigrid(multigrid->ref);
    grid_cpu_free_multigrid(multigrid->cpu);
    free(multigrid);
  }
}

// EOF
