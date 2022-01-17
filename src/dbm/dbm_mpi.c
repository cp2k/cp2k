/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dbm_mpi.h"

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

/*******************************************************************************
 * \brief Wrapper around MPI_Init.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mpi_init(int *argc, char ***argv) {
#if defined(__parallel)
  CHECK(MPI_Init(argc, argv));
#else
  (void)argc; // mark used
  (void)argv;
#endif
}

/*******************************************************************************
 * \brief Wrapper around MPI_Finalize.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mpi_finalize() {
#if defined(__parallel)
  CHECK(MPI_Finalize());
#endif
}

/*******************************************************************************
 * \brief Returns MPI_COMM_WORLD.
 * \author Ole Schuett
 ******************************************************************************/
dbm_mpi_comm_t dbm_mpi_get_comm_world() {
#if defined(__parallel)
  return MPI_COMM_WORLD;
#else
  return -1;
#endif
}

/*******************************************************************************
 * \brief Wrapper around MPI_Comm_f2c.
 * \author Ole Schuett
 ******************************************************************************/
dbm_mpi_comm_t dbm_mpi_comm_f2c(const int fortran_comm) {
#if defined(__parallel)
  return MPI_Comm_f2c(fortran_comm);
#else
  (void)fortran_comm; // mark used
  return -1;
#endif
}

/*******************************************************************************
 * \brief Wrapper around MPI_Comm_c2f.
 * \author Ole Schuett
 ******************************************************************************/
int dbm_mpi_comm_c2f(const dbm_mpi_comm_t comm) {
#if defined(__parallel)
  return MPI_Comm_c2f(comm);
#else
  (void)comm; // mark used
  return -1;
#endif
}

/*******************************************************************************
 * \brief Wrapper around MPI_Comm_rank.
 * \author Ole Schuett
 ******************************************************************************/
int dbm_mpi_comm_rank(const dbm_mpi_comm_t comm) {
#if defined(__parallel)
  int rank;
  CHECK(MPI_Comm_rank(comm, &rank));
  return rank;
#else
  (void)comm; // mark used
  return 0;
#endif
}

/*******************************************************************************
 * \brief Wrapper around MPI_Comm_size.
 * \author Ole Schuett
 ******************************************************************************/
int dbm_mpi_comm_size(const dbm_mpi_comm_t comm) {
#if defined(__parallel)
  int nranks;
  CHECK(MPI_Comm_size(comm, &nranks));
  return nranks;
#else
  (void)comm; // mark used
  return 1;
#endif
}

/*******************************************************************************
 * \brief Wrapper around MPI_Cart_create.
 * \author Ole Schuett
 ******************************************************************************/
dbm_mpi_comm_t dbm_mpi_cart_create(const dbm_mpi_comm_t comm_old,
                                   const int ndims, const int dims[],
                                   const int periods[], const int reorder) {
#if defined(__parallel)
  dbm_mpi_comm_t comm_cart;
  CHECK(MPI_Cart_create(comm_old, ndims, dims, periods, reorder, &comm_cart));
  return comm_cart;
#else
  (void)comm_old; // mark used
  (void)ndims;
  (void)dims;
  (void)periods;
  (void)reorder;
  return -1;
#endif
}

/*******************************************************************************
 * \brief Wrapper around MPI_Cart_get.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mpi_cart_get(const dbm_mpi_comm_t comm, int maxdims, int dims[],
                      int periods[], int coords[]) {
#if defined(__parallel)
  CHECK(MPI_Cart_get(comm, maxdims, dims, periods, coords));
#else
  (void)comm; // mark used
  for (int i = 0; i < maxdims; i++) {
    dims[i] = 1;
    periods[i] = 1;
    coords[i] = 0;
  }
#endif
}

/*******************************************************************************
 * \brief Wrapper around MPI_Cart_rank.
 * \author Ole Schuett
 ******************************************************************************/
int dbm_mpi_cart_rank(const dbm_mpi_comm_t comm, const int coords[]) {
#if defined(__parallel)
  int rank;
  CHECK(MPI_Cart_rank(comm, coords, &rank));
  return rank;
#else
  (void)comm; // mark used
  (void)coords;
  return 0;
#endif
}

/*******************************************************************************
 * \brief Wrapper around MPI_Cart_sub.
 * \author Ole Schuett
 ******************************************************************************/
dbm_mpi_comm_t dbm_mpi_cart_sub(const dbm_mpi_comm_t comm,
                                const int remain_dims[]) {
#if defined(__parallel)
  dbm_mpi_comm_t newcomm;
  CHECK(MPI_Cart_sub(comm, remain_dims, &newcomm));
  return newcomm;
#else
  (void)comm; // mark used
  (void)remain_dims;
  return -1;
#endif
}

/*******************************************************************************
 * \brief Wrapper around MPI_Comm_free.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mpi_comm_free(dbm_mpi_comm_t *comm) {
#if defined(__parallel)
  CHECK(MPI_Comm_free(comm));
#else
  (void)comm;  // mark used
#endif
}

/*******************************************************************************
 * \brief Wrapper around MPI_Comm_compare.
 * \author Ole Schuett
 ******************************************************************************/
bool dbm_mpi_comms_are_similar(const dbm_mpi_comm_t comm1,
                               const dbm_mpi_comm_t comm2) {
#if defined(__parallel)
  int res;
  CHECK(MPI_Comm_compare(comm1, comm2, &res));
  return res == MPI_IDENT || res == MPI_CONGRUENT || res == MPI_SIMILAR;
#else
  (void)comm1; // mark used
  (void)comm2;
  return true;
#endif
}

/*******************************************************************************
 * \brief Wrapper around MPI_Allreduce for op MPI_MAX and datatype MPI_INT.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mpi_max_int(int *values, const int count, const dbm_mpi_comm_t comm) {
#if defined(__parallel)
  int *recvbuf = malloc(count * sizeof(int));
  CHECK(MPI_Allreduce(values, recvbuf, count, MPI_INT, MPI_MAX, comm));
  memcpy(values, recvbuf, count * sizeof(int));
  free(recvbuf);
#else
  (void)comm; // mark used
  (void)values;
  (void)count;
#endif
}

/*******************************************************************************
 * \brief Wrapper around MPI_Allreduce for op MPI_MAX and datatype MPI_DOUBLE.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mpi_max_double(double *values, const int count,
                        const dbm_mpi_comm_t comm) {
#if defined(__parallel)
  double *recvbuf = malloc(count * sizeof(double));
  CHECK(MPI_Allreduce(values, recvbuf, count, MPI_DOUBLE, MPI_MAX, comm));
  memcpy(values, recvbuf, count * sizeof(double));
  free(recvbuf);
#else
  (void)comm; // mark used
  (void)values;
  (void)count;
#endif
}

/*******************************************************************************
 * \brief Wrapper around MPI_Allreduce for op MPI_SUM and datatype MPI_INT.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mpi_sum_int(int *values, const int count, const dbm_mpi_comm_t comm) {
#if defined(__parallel)
  int *recvbuf = malloc(count * sizeof(int));
  CHECK(MPI_Allreduce(values, recvbuf, count, MPI_INT, MPI_SUM, comm));
  memcpy(values, recvbuf, count * sizeof(int));
  free(recvbuf);
#else
  (void)comm; // mark used
  (void)values;
  (void)count;
#endif
}

/*******************************************************************************
 * \brief Wrapper around MPI_Allreduce for op MPI_SUM and datatype MPI_INT64_T.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mpi_sum_int64(int64_t *values, const int count,
                       const dbm_mpi_comm_t comm) {
#if defined(__parallel)
  int64_t *recvbuf = malloc(count * sizeof(int64_t));
  CHECK(MPI_Allreduce(values, recvbuf, count, MPI_INT64_T, MPI_SUM, comm));
  memcpy(values, recvbuf, count * sizeof(int64_t));
  free(recvbuf);
#else
  (void)comm; // mark used
  (void)values;
  (void)count;
#endif
}

/*******************************************************************************
 * \brief Wrapper around MPI_Allreduce for op MPI_SUM and datatype MPI_DOUBLE.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mpi_sum_double(double *values, const int count,
                        const dbm_mpi_comm_t comm) {
#if defined(__parallel)
  double *recvbuf = malloc(count * sizeof(double));
  CHECK(MPI_Allreduce(values, recvbuf, count, MPI_DOUBLE, MPI_SUM, comm));
  memcpy(values, recvbuf, count * sizeof(double));
  free(recvbuf);
#else
  (void)comm; // mark used
  (void)values;
  (void)count;
#endif
}

/*******************************************************************************
 * \brief Wrapper around MPI_Sendrecv for datatype MPI_BYTE.
 * \author Ole Schuett
 ******************************************************************************/
int dbm_mpi_sendrecv_byte(const void *sendbuf, const int sendcount,
                          const int dest, const int sendtag, void *recvbuf,
                          const int recvcount, const int source,
                          const int recvtag, const dbm_mpi_comm_t comm) {
#if defined(__parallel)
  MPI_Status status;
  CHECK(MPI_Sendrecv(sendbuf, sendcount, MPI_BYTE, dest, sendtag, recvbuf,
                     recvcount, MPI_BYTE, source, recvtag, comm, &status))
  int count_received;
  CHECK(MPI_Get_count(&status, MPI_BYTE, &count_received));
  return count_received;
#else
  (void)sendbuf; // mark used
  (void)sendcount;
  (void)dest;
  (void)sendtag;
  (void)recvbuf;
  (void)recvcount;
  (void)source;
  (void)recvtag;
  (void)comm;
  fprintf(stderr, "Error: dbm_mpi_sendrecv_byte not available without MPI\n");
  abort();
#endif
}

/*******************************************************************************
 * \brief Wrapper around MPI_Sendrecv for datatype MPI_DOUBLE.
 * \author Ole Schuett
 ******************************************************************************/
int dbm_mpi_sendrecv_double(const double *sendbuf, const int sendcount,
                            const int dest, const int sendtag, double *recvbuf,
                            const int recvcount, const int source,
                            const int recvtag, const dbm_mpi_comm_t comm) {
#if defined(__parallel)
  MPI_Status status;
  CHECK(MPI_Sendrecv(sendbuf, sendcount, MPI_DOUBLE, dest, sendtag, recvbuf,
                     recvcount, MPI_DOUBLE, source, recvtag, comm, &status))
  int count_received;
  CHECK(MPI_Get_count(&status, MPI_DOUBLE, &count_received));
  return count_received;
#else
  (void)sendbuf; // mark used
  (void)sendcount;
  (void)dest;
  (void)sendtag;
  (void)recvbuf;
  (void)recvcount;
  (void)source;
  (void)recvtag;
  (void)comm;
  fprintf(stderr, "Error: dbm_mpi_sendrecv_double not available without MPI\n");
  abort();
#endif
}

/*******************************************************************************
 * \brief Wrapper around MPI_Alltoall for datatype MPI_INT.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mpi_alltoall_int(const int *sendbuf, const int sendcount, int *recvbuf,
                          const int recvcount, const dbm_mpi_comm_t comm) {
#if defined(__parallel)
  CHECK(MPI_Alltoall(sendbuf, sendcount, MPI_INT, recvbuf, recvcount, MPI_INT,
                     comm));
#else
  (void)comm; // mark used
  assert(sendcount == recvcount);
  memcpy(recvbuf, sendbuf, sendcount * sizeof(int));
#endif
}

/*******************************************************************************
 * \brief Wrapper around MPI_Alltoallv for datatype MPI_BYTE.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mpi_alltoallv_byte(const void *sendbuf, const int *sendcounts,
                            const int *sdispls, void *recvbuf,
                            const int *recvcounts, const int *rdispls,
                            const dbm_mpi_comm_t comm) {
#if defined(__parallel)
  CHECK(MPI_Alltoallv(sendbuf, sendcounts, sdispls, MPI_BYTE, recvbuf,
                      recvcounts, rdispls, MPI_BYTE, comm));
#else
  (void)comm; // mark used
  assert(sendcounts[0] == recvcounts[0]);
  assert(sdispls[0] == 0 && rdispls[0] == 0);
  memcpy(recvbuf, sendbuf, sendcounts[0]);
#endif
}

/*******************************************************************************
 * \brief Wrapper around MPI_Alltoallv for datatype MPI_DOUBLE.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mpi_alltoallv_double(const double *sendbuf, const int *sendcounts,
                              const int *sdispls, double *recvbuf,
                              const int *recvcounts, const int *rdispls,
                              const dbm_mpi_comm_t comm) {
#if defined(__parallel)
  CHECK(MPI_Alltoallv(sendbuf, sendcounts, sdispls, MPI_DOUBLE, recvbuf,
                      recvcounts, rdispls, MPI_DOUBLE, comm));
#else
  (void)comm; // mark used
  assert(sendcounts[0] == recvcounts[0]);
  assert(sdispls[0] == 0 && rdispls[0] == 0);
  memcpy(recvbuf, sendbuf, sendcounts[0] * sizeof(double));
#endif
}

// EOF
