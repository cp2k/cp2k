/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2024 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include "grid_mpi.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#if defined(__parallel)
/*******************************************************************************
 * \brief Check given MPI status and upon failure abort with a nice message.
 * \author Ole Schuett
 ******************************************************************************/
#define CHECK(status)                                                          \
  if (status != MPI_SUCCESS) {                                                 \
    fprintf(stderr, "MPI error in %s:%i\n", __FILE__, __LINE__);               \
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);                                   \
  }
#endif

grid_mpi_comm grid_mpi_comm_f2c(const grid_mpi_fint fortran_comm) {
#if defined(__parallel)
  return MPI_Comm_f2c(fortran_comm);
#else
  return (grid_mpi_comm)fortran_comm;
#endif
}

grid_mpi_fint grid_mpi_comm_c2f(const grid_mpi_comm comm) {
#if defined(__parallel)
  return MPI_Comm_c2f(comm);
#else
  return (grid_mpi_fint)comm;
#endif
}

void grid_mpi_comm_dup(const grid_mpi_comm old_comm, grid_mpi_comm *new_comm) {
#if defined(__parallel)
  CHECK(MPI_Comm_dup(old_comm, new_comm));
#else
  *new_comm = old_comm;
#endif
}

void grid_mpi_comm_free(grid_mpi_comm *comm) {
#if defined(__parallel)
  CHECK(MPI_Comm_free(comm));
#else
  *comm = grid_mpi_comm_null;
#endif
}

// EOF