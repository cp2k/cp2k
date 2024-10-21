/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2024 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#ifndef GRID_MPI_H
#define GRID_MPI_H

#if defined(__parallel)
#include <mpi.h>

typedef MPI_Comm grid_mpi_comm;
typedef MPI_Fint grid_mpi_fint;

static const grid_mpi_comm grid_mpi_comm_world = MPI_COMM_WORLD;
static const grid_mpi_comm grid_mpi_comm_null = MPI_COMM_NULL;
#else
typedef int grid_mpi_comm;
typedef int grid_mpi_fint;

static const grid_mpi_comm grid_mpi_comm_world = -2;
static const grid_mpi_comm grid_mpi_comm_null = -3;
#endif

grid_mpi_comm grid_mpi_comm_f2c(const grid_mpi_fint fortran_comm);

grid_mpi_fint grid_mpi_comm_c2f(const grid_mpi_comm comm);

void grid_mpi_comm_dup(const grid_mpi_comm old_comm, grid_mpi_comm *new_comm);

void grid_mpi_comm_free(grid_mpi_comm *comm);

#endif

// EOF