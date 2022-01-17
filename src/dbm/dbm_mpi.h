/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#ifndef DBM_MPI_H
#define DBM_MPI_H

#include <stdbool.h>
#include <stdint.h>

#if defined(__parallel)
#include <mpi.h>
typedef MPI_Comm dbm_mpi_comm_t;
#else
typedef int dbm_mpi_comm_t;
#endif

/*******************************************************************************
 * \brief Wrapper around MPI_Init.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mpi_init(int *argc, char ***argv);

/*******************************************************************************
 * \brief Wrapper around MPI_Finalize.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mpi_finalize();

/*******************************************************************************
 * \brief Returns MPI_COMM_WORLD.
 * \author Ole Schuett
 ******************************************************************************/
dbm_mpi_comm_t dbm_mpi_get_comm_world();

/*******************************************************************************
 * \brief Wrapper around MPI_Comm_f2c.
 * \author Ole Schuett
 ******************************************************************************/
dbm_mpi_comm_t dbm_mpi_comm_f2c(const int fortran_comm);

/*******************************************************************************
 * \brief Wrapper around MPI_Comm_c2f.
 * \author Ole Schuett
 ******************************************************************************/
int dbm_mpi_comm_c2f(const dbm_mpi_comm_t comm);

/*******************************************************************************
 * \brief Wrapper around MPI_Comm_rank.
 * \author Ole Schuett
 ******************************************************************************/
int dbm_mpi_comm_rank(const dbm_mpi_comm_t comm);

/*******************************************************************************
 * \brief Wrapper around MPI_Comm_size.
 * \author Ole Schuett
 ******************************************************************************/
int dbm_mpi_comm_size(const dbm_mpi_comm_t comm);

/*******************************************************************************
 * \brief Wrapper around MPI_Cart_create.
 * \author Ole Schuett
 ******************************************************************************/
dbm_mpi_comm_t dbm_mpi_cart_create(const dbm_mpi_comm_t comm_old,
                                   const int ndims, const int dims[],
                                   const int periods[], const int reorder);

/*******************************************************************************
 * \brief Wrapper around MPI_Cart_get.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mpi_cart_get(const dbm_mpi_comm_t comm, int maxdims, int dims[],
                      int periods[], int coords[]);

/*******************************************************************************
 * \brief Wrapper around MPI_Cart_rank.
 * \author Ole Schuett
 ******************************************************************************/
int dbm_mpi_cart_rank(const dbm_mpi_comm_t comm, const int coords[]);

/*******************************************************************************
 * \brief Wrapper around MPI_Cart_sub.
 * \author Ole Schuett
 ******************************************************************************/
dbm_mpi_comm_t dbm_mpi_cart_sub(const dbm_mpi_comm_t comm,
                                const int remain_dims[]);

/*******************************************************************************
 * \brief Wrapper around MPI_Comm_free.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mpi_comm_free(dbm_mpi_comm_t *comm);

/*******************************************************************************
 * \brief Wrapper around MPI_Comm_compare.
 * \author Ole Schuett
 ******************************************************************************/
bool dbm_mpi_comms_are_similar(const dbm_mpi_comm_t comm1,
                               const dbm_mpi_comm_t comm2);

/*******************************************************************************
 * \brief Wrapper around MPI_Allreduce for op MPI_MAX and datatype MPI_INT.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mpi_max_int(int *values, const int count, const dbm_mpi_comm_t comm);

/*******************************************************************************
 * \brief Wrapper around MPI_Allreduce for op MPI_MAX and datatype MPI_DOUBLE.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mpi_max_double(double *values, const int count,
                        const dbm_mpi_comm_t comm);

/*******************************************************************************
 * \brief Wrapper around MPI_Allreduce for op MPI_SUM and datatype MPI_INT.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mpi_sum_int(int *values, const int count, const dbm_mpi_comm_t comm);

/*******************************************************************************
 * \brief Wrapper around MPI_Allreduce for op MPI_SUM and datatype MPI_INT64_T.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mpi_sum_int64(int64_t *values, const int count,
                       const dbm_mpi_comm_t comm);

/*******************************************************************************
 * \brief Wrapper around MPI_Allreduce for op MPI_SUM and datatype MPI_DOUBLE.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mpi_sum_double(double *values, const int count,
                        const dbm_mpi_comm_t comm);

/*******************************************************************************
 * \brief Wrapper around MPI_Sendrecv for datatype MPI_BYTE.
 * \author Ole Schuett
 ******************************************************************************/
int dbm_mpi_sendrecv_byte(const void *sendbuf, const int sendcount,
                          const int dest, const int sendtag, void *recvbuf,
                          const int recvcount, const int source,
                          const int recvtag, const dbm_mpi_comm_t comm);

/*******************************************************************************
 * \brief Wrapper around MPI_Sendrecv for datatype MPI_DOUBLE.
 * \author Ole Schuett
 ******************************************************************************/
int dbm_mpi_sendrecv_double(const double *sendbuf, const int sendcount,
                            const int dest, const int sendtag, double *recvbuf,
                            const int recvcount, const int source,
                            const int recvtag, const dbm_mpi_comm_t comm);

/*******************************************************************************
 * \brief Wrapper around MPI_Alltoall for datatype MPI_INT.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mpi_alltoall_int(const int *sendbuf, const int sendcount, int *recvbuf,
                          const int recvcount, const dbm_mpi_comm_t comm);

/*******************************************************************************
 * \brief Wrapper around MPI_Alltoallv for datatype MPI_BYTE.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mpi_alltoallv_byte(const void *sendbuf, const int *sendcounts,
                            const int *sdispls, void *recvbuf,
                            const int *recvcounts, const int *rdispls,
                            const dbm_mpi_comm_t comm);

/*******************************************************************************
 * \brief Wrapper around MPI_Alltoallv for datatype MPI_DOUBLE.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mpi_alltoallv_double(const double *sendbuf, const int *sendcounts,
                              const int *sdispls, double *recvbuf,
                              const int *recvcounts, const int *rdispls,
                              const dbm_mpi_comm_t comm);

#endif

// EOF
