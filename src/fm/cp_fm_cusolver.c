/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2024 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#if defined(__CUSOLVERMP)

#include "../offload/offload_library.h"
#include <assert.h>
#include <cal.h>
#include <cuda_runtime.h>
#include <cusolverMp.h>
#include <mpi.h>
#include <stdlib.h>

/*******************************************************************************
 * \brief Run given CUDA command and upon failure abort with a nice message.
 * \author Ole Schuett
 ******************************************************************************/
#define CUDA_CHECK(cmd)                                                        \
  do {                                                                         \
    cudaError_t status = cmd;                                                  \
    if (status != cudaSuccess) {                                               \
      fprintf(stderr, "ERROR: %s %s %d\n", cudaGetErrorString(status),         \
              __FILE__, __LINE__);                                             \
      abort();                                                                 \
    }                                                                          \
  } while (0)

/*******************************************************************************
 * \brief Decode given cal error.
 * \author Ole Schuett
 ******************************************************************************/
static char *calGetErrorString(calError_t status) {
  switch (status) {
  case CAL_OK:
    return "CAL_OK";
  case CAL_ERROR:
    return "CAL_ERROR";
  case CAL_ERROR_INVALID_PARAMETER:
    return "CAL_ERROR_INVALID_PARAMETER";
  case CAL_ERROR_INTERNAL:
    return "CAL_ERROR_INTERNAL";
  case CAL_ERROR_CUDA:
    return "CAL_ERROR_CUDA";
  case CAL_ERROR_UCC:
    return "CAL_ERROR_UCC";
  case CAL_ERROR_NOT_SUPPORTED:
    return "CAL_ERROR_NOT_SUPPORTED";
  case CAL_ERROR_INPROGRESS:
    return "CAL_ERROR_INPROGRESS";
  default:
    return "CAL UNKNOWN ERROR";
  }
}

/*******************************************************************************
 * \brief Run given cal command and upon failure abort with a nice message.
 * \author Ole Schuett
 ******************************************************************************/
#define CAL_CHECK(cmd)                                                         \
  do {                                                                         \
    calError_t status = cmd;                                                   \
    if (status != CAL_OK) {                                                    \
      fprintf(stderr, "ERROR: %s %s %d\n", calGetErrorString(status),          \
              __FILE__, __LINE__);                                             \
      abort();                                                                 \
    }                                                                          \
  } while (0)

/*******************************************************************************
 * \brief Decode given cusolver error.
 * \author Ole Schuett
 ******************************************************************************/
static char *cusolverGetErrorString(cusolverStatus_t status) {
  switch (status) {
  case CUSOLVER_STATUS_SUCCESS:
    return "CUSOLVER_STATUS_SUCCESS";
  case CUSOLVER_STATUS_NOT_INITIALIZED:
    return "CUSOLVER_STATUS_NOT_INITIALIZED";
  case CUSOLVER_STATUS_ALLOC_FAILED:
    return "CUSOLVER_STATUS_ALLOC_FAILED";
  case CUSOLVER_STATUS_INVALID_VALUE:
    return "CUSOLVER_STATUS_INVALID_VALUE";
  case CUSOLVER_STATUS_ARCH_MISMATCH:
    return "CUSOLVER_STATUS_ARCH_MISMATCH";
  case CUSOLVER_STATUS_MAPPING_ERROR:
    return "CUSOLVER_STATUS_MAPPING_ERROR";
  case CUSOLVER_STATUS_EXECUTION_FAILED:
    return "CUSOLVER_STATUS_EXECUTION_FAILED";
  case CUSOLVER_STATUS_INTERNAL_ERROR:
    return "CUSOLVER_STATUS_INTERNAL_ERROR";
  case CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
    return "CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED";
  case CUSOLVER_STATUS_NOT_SUPPORTED:
    return "CUSOLVER_STATUS_NOT_SUPPORTED";
  case CUSOLVER_STATUS_ZERO_PIVOT:
    return "CUSOLVER_STATUS_ZERO_PIVOT";
  case CUSOLVER_STATUS_INVALID_LICENSE:
    return "CUSOLVER_STATUS_INVALID_LICENSE";
  case CUSOLVER_STATUS_IRS_PARAMS_NOT_INITIALIZED:
    return "CUSOLVER_STATUS_IRS_PARAMS_NOT_INITIALIZED";
  case CUSOLVER_STATUS_IRS_PARAMS_INVALID:
    return "CUSOLVER_STATUS_IRS_PARAMS_INVALID";
  case CUSOLVER_STATUS_IRS_PARAMS_INVALID_PREC:
    return "CUSOLVER_STATUS_IRS_PARAMS_INVALID_PREC";
  case CUSOLVER_STATUS_IRS_PARAMS_INVALID_REFINE:
    return "CUSOLVER_STATUS_IRS_PARAMS_INVALID_REFINE";
  case CUSOLVER_STATUS_IRS_PARAMS_INVALID_MAXITER:
    return "CUSOLVER_STATUS_IRS_PARAMS_INVALID_MAXITER";
  case CUSOLVER_STATUS_IRS_INTERNAL_ERROR:
    return "CUSOLVER_STATUS_IRS_INTERNAL_ERROR";
  case CUSOLVER_STATUS_IRS_NOT_SUPPORTED:
    return "CUSOLVER_STATUS_IRS_NOT_SUPPORTED";
  case CUSOLVER_STATUS_IRS_OUT_OF_RANGE:
    return "CUSOLVER_STATUS_IRS_OUT_OF_RANGE";
  case CUSOLVER_STATUS_IRS_NRHS_NOT_SUPPORTED_FOR_REFINE_GMRES:
    return "CUSOLVER_STATUS_IRS_NRHS_NOT_SUPPORTED_FOR_REFINE_GMRES";
  case CUSOLVER_STATUS_IRS_INFOS_NOT_INITIALIZED:
    return "CUSOLVER_STATUS_IRS_INFOS_NOT_INITIALIZED";
  case CUSOLVER_STATUS_IRS_INFOS_NOT_DESTROYED:
    return "CUSOLVER_STATUS_IRS_INFOS_NOT_DESTROYED";
  case CUSOLVER_STATUS_IRS_MATRIX_SINGULAR:
    return "CUSOLVER_STATUS_IRS_MATRIX_SINGULAR";
  case CUSOLVER_STATUS_INVALID_WORKSPACE:
    return "CUSOLVER_STATUS_INVALID_WORKSPACE";
  default:
    return "CUSOLVER UNKNOWN ERROR";
  }
}

/*******************************************************************************
 * \brief Run given cusolver command and upon failure abort with a nice message.
 * \author Ole Schuett
 ******************************************************************************/
#define CUSOLVER_CHECK(cmd)                                                    \
  do {                                                                         \
    cusolverStatus_t status = cmd;                                             \
    if (status != CUSOLVER_STATUS_SUCCESS) {                                   \
      fprintf(stderr, "ERROR: %s %s %d\n", cusolverGetErrorString(status),     \
              __FILE__, __LINE__);                                             \
      abort();                                                                 \
    }                                                                          \
  } while (0)

/*******************************************************************************
 * \brief Callback for cal library to initiate an allgather operation.
 * \author Ole Schuett
 ******************************************************************************/
static calError_t allgather(void *src_buf, void *recv_buf, size_t size,
                            void *data, void **req) {
  const MPI_Comm comm = *(MPI_Comm *)data;
  MPI_Request *request = malloc(sizeof(MPI_Request));
  *req = request;
  const int status = MPI_Iallgather(src_buf, size, MPI_BYTE, recv_buf, size,
                                    MPI_BYTE, comm, request);
  return (status == MPI_SUCCESS) ? CAL_OK : CAL_ERROR;
}

/*******************************************************************************
 * \brief Callback for cal library to test if a request has completed.
 * \author Ole Schuett
 ******************************************************************************/
static calError_t req_test(void *req) {
  MPI_Request *request = (MPI_Request *)(req);
  int completed;
  const int status = MPI_Test(request, &completed, MPI_STATUS_IGNORE);
  if (status != MPI_SUCCESS) {
    return CAL_ERROR;
  }
  return completed ? CAL_OK : CAL_ERROR_INPROGRESS;
}

/*******************************************************************************
 * \brief Callback for cal library to free a request.
 * \author Ole Schuett
 ******************************************************************************/
static calError_t req_free(void *req) {
  free(req);
  return CAL_OK;
}

/*******************************************************************************
 * \brief Driver routine to diagonalize a matrix with the cuSOLVERMp library.
 * \author Ole Schuett
 ******************************************************************************/
void cp_fm_diag_cusolver(const int fortran_comm, const int matrix_desc[9],
                         const int nprow, const int npcol, const int myprow,
                         const int mypcol, const int n, const double *matrix,
                         double *eigenvectors, double *eigenvalues) {

  offload_activate_chosen_device();
  const int local_device = offload_get_chosen_device();

  MPI_Comm comm = MPI_Comm_f2c(fortran_comm);
  int rank, nranks;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nranks);

  // Create CAL communicator.
  cal_comm_t cal_comm = NULL;
  cal_comm_create_params_t params;
  params.allgather = &allgather;
  params.req_test = &req_test;
  params.req_free = &req_free;
  params.data = &comm;
  params.rank = rank;
  params.nranks = nranks;
  params.local_device = local_device;
  CAL_CHECK(cal_comm_create(params, &cal_comm));

  // Create various handles.
  cudaStream_t stream = NULL;
  CUDA_CHECK(cudaStreamCreate(&stream));

  cusolverMpHandle_t cusolvermp_handle = NULL;
  CUSOLVER_CHECK(cusolverMpCreate(&cusolvermp_handle, local_device, stream));

  cusolverMpGrid_t grid = NULL;
  CUSOLVER_CHECK(cusolverMpCreateDeviceGrid(cusolvermp_handle, &grid, cal_comm,
                                            nprow, npcol,
                                            CUSOLVERMP_GRID_MAPPING_ROW_MAJOR));
  const int mb = matrix_desc[4];
  const int nb = matrix_desc[5];
  const int rsrc = matrix_desc[6];
  const int csrc = matrix_desc[7];
  const int ldA = matrix_desc[8];
  assert(rsrc == csrc);
  assert(ldA >= 1);

  const int np = cusolverMpNUMROC(n, mb, myprow, rsrc, nprow);
  const int nq = cusolverMpNUMROC(n, nb, mypcol, csrc, npcol);
  assert(np == ldA);

  const cublasFillMode_t fill_mode = CUBLAS_FILL_MODE_UPPER;
  const cudaDataType_t data_type = CUDA_R_64F; // double

  cusolverMpMatrixDescriptor_t cusolvermp_matrix_desc = NULL;
  CUSOLVER_CHECK(cusolverMpCreateMatrixDesc(
      &cusolvermp_matrix_desc, grid, data_type, n, n, mb, nb, rsrc, csrc, np));

  // Allocate workspaces.
  size_t work_dev_size, work_host_size;
  void *DUMMY = (void *)1; // Workaround to avoid crash when passing NULL.
  CUSOLVER_CHECK(cusolverMpSyevd_bufferSize(
      cusolvermp_handle, "V", fill_mode, n, DUMMY, 1, 1, cusolvermp_matrix_desc,
      NULL, NULL, 1, 1, cusolvermp_matrix_desc, data_type, &work_dev_size,
      &work_host_size));

  double *work_dev = NULL;
  CUDA_CHECK(cudaMalloc((void **)&work_dev, work_dev_size));

  double *work_host = NULL;
  CUDA_CHECK(cudaMallocHost((void **)&work_host, work_host_size));

  // Upload input matrix.
  const size_t matrix_local_size = ldA * nq * sizeof(double);
  double *matrix_dev = NULL;
  CUDA_CHECK(cudaMalloc((void **)&matrix_dev, matrix_local_size));
  CUDA_CHECK(cudaMemcpyAsync(matrix_dev, matrix, matrix_local_size,
                             cudaMemcpyHostToDevice, stream));

  // Allocate result buffers.
  int *info_dev = NULL;
  CUDA_CHECK(cudaMalloc((void **)&info_dev, sizeof(int)));

  double *eigenvectors_dev = NULL;
  CUDA_CHECK(cudaMalloc((void **)&eigenvectors_dev, matrix_local_size));

  double *eigenvalues_dev = NULL;
  CUDA_CHECK(cudaMalloc((void **)&eigenvalues_dev, n * sizeof(double)));

  // Call solver.
  CUSOLVER_CHECK(
      cusolverMpSyevd(cusolvermp_handle, "V", fill_mode, n, matrix_dev, 1, 1,
                      cusolvermp_matrix_desc, eigenvalues_dev, eigenvectors_dev,
                      1, 1, cusolvermp_matrix_desc, data_type, work_dev,
                      work_dev_size, work_host, work_host_size, info_dev));

  // Wait for solver to finish.
  CUDA_CHECK(cudaStreamSynchronize(stream));
  CAL_CHECK(cal_stream_sync(cal_comm, stream));

  // Check info.
  int info = -1;
  CUDA_CHECK(cudaMemcpy(&info, info_dev, sizeof(int), cudaMemcpyDeviceToHost));
  assert(info == 0);

  // Download results.
  CUDA_CHECK(cudaMemcpyAsync(eigenvectors, eigenvectors_dev, matrix_local_size,
                             cudaMemcpyDeviceToHost, stream));
  CUDA_CHECK(cudaMemcpyAsync(eigenvalues, eigenvalues_dev, n * sizeof(double),
                             cudaMemcpyDeviceToHost, stream));

  // Wait for download to finish.
  CUDA_CHECK(cudaStreamSynchronize(stream));

  // Free buffers.
  CUDA_CHECK(cudaFree(matrix_dev));
  CUDA_CHECK(cudaFree(info_dev));
  CUDA_CHECK(cudaFree(eigenvectors_dev));
  CUDA_CHECK(cudaFree(eigenvalues_dev));
  CUDA_CHECK(cudaFree(work_dev));
  CUDA_CHECK(cudaFreeHost(work_host));

  // Destroy handles.
  CUSOLVER_CHECK(cusolverMpDestroyMatrixDesc(cusolvermp_matrix_desc));
  CUSOLVER_CHECK(cusolverMpDestroyGrid(grid));
  CUSOLVER_CHECK(cusolverMpDestroy(cusolvermp_handle));
  CUDA_CHECK(cudaStreamDestroy(stream));
  CAL_CHECK(cal_comm_destroy(cal_comm));

  // Sync MPI ranks to include load imbalance in total timings.
  MPI_Barrier(comm);
}

#endif

// EOF
