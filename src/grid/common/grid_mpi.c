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

int grid_mpi_comm_size(const grid_mpi_comm comm) {
#if defined(__parallel)
  int comm_size;
  CHECK(MPI_Comm_size(comm, &comm_size));
  return comm_size;
#else
  (void)comm;
  return 1;
#endif
}

int grid_mpi_comm_rank(const grid_mpi_comm comm) {
#if defined(__parallel)
  int comm_rank;
  CHECK(MPI_Comm_rank(comm, &comm_rank));
  return comm_rank;
#else
  (void)comm;
  return 0;
#endif
}

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

void grid_mpi_barrier(const grid_mpi_comm comm) {
#if defined(__parallel)
  CHECK(MPI_Barrier(comm));
#else
  // Nothing to do in the serial case
  (void)comm;
#endif
}

bool grid_mpi_comm_is_unequal(const grid_mpi_comm comm1,
                              const grid_mpi_comm comm2) {
#if defined(__parallel)
  int result = -1;
  CHECK(MPI_Comm_compare(comm1, comm2, &result));
  return result == MPI_UNEQUAL;
#else
  return ((comm1 == grid_mpi_comm_null) && (comm2 != grid_mpi_comm_null)) ||
         ((comm1 != grid_mpi_comm_null) && (comm2 == grid_mpi_comm_null));
#endif
}

bool grid_mpi_comm_is_similar(const grid_mpi_comm comm1,
                              const grid_mpi_comm comm2) {
#if defined(__parallel)
  int result = -1;
  CHECK(MPI_Comm_compare(comm1, comm2, &result));
  return result == MPI_SIMILAR || result == MPI_CONGRUENT ||
         result == MPI_IDENT;
#else
  return ((comm1 == grid_mpi_comm_null) && (comm2 == grid_mpi_comm_null)) ||
         ((comm1 != grid_mpi_comm_null) && (comm2 != grid_mpi_comm_null));
#endif
}

bool grid_mpi_comm_is_congruent(const grid_mpi_comm comm1,
                                const grid_mpi_comm comm2) {
#if defined(__parallel)
  int result = -1;
  CHECK(MPI_Comm_compare(comm1, comm2, &result));
  return result == MPI_CONGRUENT || result == MPI_IDENT;
#else
  return ((comm1 == grid_mpi_comm_null) && (comm2 == grid_mpi_comm_null)) ||
         ((comm1 != grid_mpi_comm_null) && (comm2 != grid_mpi_comm_null));
#endif
}

bool grid_mpi_comm_is_ident(const grid_mpi_comm comm1,
                            const grid_mpi_comm comm2) {
#if defined(__parallel)
  int result = -1;
  CHECK(MPI_Comm_compare(comm1, comm2, &result));
  return result == MPI_IDENT;
#else
  return ((comm1 == grid_mpi_comm_null) && (comm2 == grid_mpi_comm_null)) ||
         ((comm1 != grid_mpi_comm_null) && (comm2 != grid_mpi_comm_null));
#endif
}

void grid_mpi_sendrecv_double(const double *sendbuffer, const int sendcount,
                              const int dest, const int sendtag,
                              double *recvbuf, const int recvcount,
                              const int source, const int recvtag,
                              const grid_mpi_comm comm) {
#if defined(__parallel)
  MPI_Status status;
  CHECK(MPI_Sendrecv(sendbuffer, sendcount, MPI_DOUBLE, dest, sendtag, recvbuf,
                     recvcount, MPI_DOUBLE, source, recvtag, comm, &status));
#else
  (void)sendbuffer;
  (void)sendcount;
  (void)dest;
  (void)sendtag;
  (void)recvbuf;
  (void)recvcount;
  (void)source;
  (void)recvtag;
  (void)comm;
  assert(false && "Communication not allowed in serial mode");
#endif
}

// EOF