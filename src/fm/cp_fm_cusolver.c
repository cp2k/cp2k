/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#if defined(__CUSOLVERMP)

#include "../base/cp_assert.h"
#include "../offload/offload_library.h"
#include <cuComplex.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <cusolverMp.h>
#include <math.h>
#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(__CUSOLVERMP_NCCL)
#include <nccl.h>
#else
#include <cal.h>
#endif

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

#if defined(__CUSOLVERMP_NCCL)
/*******************************************************************************
 * \brief Run given NCCL command and upon failure abort with a nice message.
 * \author Jiri Vyskocil
 ******************************************************************************/
#define NCCL_CHECK(cmd)                                                        \
  do {                                                                         \
    ncclResult_t status = cmd;                                                 \
    if (status != ncclSuccess) {                                               \
      fprintf(stderr, "ERROR: %s %s %d\n", ncclGetErrorString(status),         \
              __FILE__, __LINE__);                                             \
      abort();                                                                 \
    }                                                                          \
  } while (0)

/*******************************************************************************
 * \brief Check that the selected CUDA device is visible to this MPI rank.
 ******************************************************************************/
static void check_nccl_device_visibility(MPI_Comm comm,
                                         const int local_device) {
  int device_count = 0;
  CUDA_CHECK(cudaGetDeviceCount(&device_count));

  int rank;
  MPI_Comm_rank(comm, &rank);

  if (device_count <= 0) {
    fprintf(stderr,
            "ERROR: cuSOLVERMp with NCCL requires at least one visible GPU. "
            "No CUDA device is visible to MPI rank %d.\n",
            rank);
    MPI_Abort(comm, EXIT_FAILURE);
  }

  if (local_device < 0 || local_device >= device_count) {
    fprintf(stderr,
            "ERROR: cuSOLVERMp with NCCL selected invalid CUDA device %d on "
            "MPI rank %d; %d CUDA device(s) are visible to this rank.\n",
            local_device, rank, device_count);
    MPI_Abort(comm, EXIT_FAILURE);
  }
}

#else
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
#endif /* __CUSOLVERMP_NCCL */

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
 * \brief cuSOLVERMp execution resources owned by a BLACS environment.
 ******************************************************************************/
typedef struct {
  MPI_Comm comm;
  int device, nprow, npcol;
  cudaStream_t stream;
  cusolverMpHandle_t handle;
  cusolverDnHandle_t dense;
  cusolverMpGrid_t grid;
#if defined(__CUSOLVERMP_NCCL)
  ncclComm_t nccl_comm;
#else
  cal_comm_t cal_comm;
#endif
} cp_fm_cusolver_context;

/*******************************************************************************
 * \brief Check once per solve communicator that each rank has its own device.
 * CUDA ordinals can differ between processes even for the same physical GPU.
 ******************************************************************************/
int cp_fm_cusolver_devices_ok(int fortran_comm) {
  MPI_Comm comm = MPI_Comm_f2c(fortran_comm);
  const int device = offload_get_chosen_device();
  int ok = (device >= 0), ranks;
  MPI_Comm_size(comm, &ranks);
  if (ranks == 1)
    return ok;
  MPI_Comm node;
  MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &node);
  int nlocal, rank;
  MPI_Comm_size(node, &nlocal);
  MPI_Comm_rank(node, &rank);
  if (nlocal > 1) {
    struct cudaDeviceProp prop = {0};
    if (ok)
      CUDA_CHECK(cudaGetDeviceProperties(&prop, device));
    cudaUUID_t *ids = malloc(nlocal * sizeof(*ids));
    cp_assert(ids);
    MPI_Allgather(&prop.uuid, sizeof(*ids), MPI_BYTE, ids, sizeof(*ids),
                  MPI_BYTE, node);
    for (int i = 0; i < rank; ++i) {
      if (memcmp(&ids[i], &prop.uuid, sizeof(*ids)) == 0)
        ok = 0;
    }
    free(ids);
  }
  MPI_Comm_free(&node);
  MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, comm);
  return ok;
}

void cp_fm_cusolver_context_destroy(void **context);

/*******************************************************************************
 * \brief Lazily initialize the execution resources for a BLACS environment.
 ******************************************************************************/
static cp_fm_cusolver_context *get_context(void **context,
                                           const int fortran_comm,
                                           const int nprow, const int npcol,
                                           const int local) {
  const int device = offload_get_chosen_device();
  MPI_Comm comm = MPI_Comm_f2c(fortran_comm);
  cp_fm_cusolver_context *ctx = *context;
#if defined(__CUSOLVERMP_NCCL)
  if (!ctx) {
    check_nccl_device_visibility(comm, device);
  }
#endif
  offload_activate_chosen_device();
  if (ctx && (ctx->dense != NULL) != local) {
    cp_fm_cusolver_context_destroy(context);
    ctx = NULL;
  }
  if (ctx) {
    int comparison;
    MPI_Comm_compare(comm, ctx->comm, &comparison);
    cp_assert(comparison == MPI_IDENT || comparison == MPI_CONGRUENT);
    cp_assert(ctx->device == device && ctx->nprow == nprow &&
              ctx->npcol == npcol);
    return ctx;
  }

  ctx = calloc(1, sizeof(*ctx));
  cp_assert(ctx);
  MPI_Comm_dup(comm, &ctx->comm);
  ctx->device = device;
  ctx->nprow = nprow;
  ctx->npcol = npcol;
  CUDA_CHECK(cudaStreamCreate(&ctx->stream));
  if (local) {
    CUSOLVER_CHECK(cusolverDnCreate(&ctx->dense));
    CUSOLVER_CHECK(cusolverDnSetStream(ctx->dense, ctx->stream));
    *context = ctx;
    return ctx;
  }
  int rank, nranks;
  MPI_Comm_rank(ctx->comm, &rank);
  MPI_Comm_size(ctx->comm, &nranks);

#if defined(__CUSOLVERMP_NCCL)
  ncclUniqueId nccl_id;
  if (rank == 0) {
    NCCL_CHECK(ncclGetUniqueId(&nccl_id));
  }
  MPI_Bcast(&nccl_id, sizeof(nccl_id), MPI_BYTE, 0, ctx->comm);
  NCCL_CHECK(ncclCommInitRank(&ctx->nccl_comm, nranks, nccl_id, rank));
#else
  cal_comm_create_params_t params;
  params.allgather = &allgather;
  params.req_test = &req_test;
  params.req_free = &req_free;
  params.data = &ctx->comm;
  params.rank = rank;
  params.nranks = nranks;
  params.local_device = device;
  CAL_CHECK(cal_comm_create(params, &ctx->cal_comm));
#endif

  CUSOLVER_CHECK(cusolverMpCreate(&ctx->handle, device, ctx->stream));
#if defined(__CUSOLVERMP_NCCL)
  CUSOLVER_CHECK(cusolverMpCreateDeviceGrid(ctx->handle, &ctx->grid,
                                            ctx->nccl_comm, nprow, npcol,
                                            CUSOLVERMP_GRID_MAPPING_ROW_MAJOR));
#else
  CUSOLVER_CHECK(cusolverMpCreateDeviceGrid(ctx->handle, &ctx->grid,
                                            ctx->cal_comm, nprow, npcol,
                                            CUSOLVERMP_GRID_MAPPING_ROW_MAJOR));
#endif
  *context = ctx;
  return ctx;
}

/*******************************************************************************
 * \brief Release execution resources before their BLACS environment is freed.
 ******************************************************************************/
void cp_fm_cusolver_context_destroy(void **context) {
  cp_fm_cusolver_context *ctx = *context;
  if (!ctx) {
    return;
  }

  int device;
  CUDA_CHECK(cudaGetDevice(&device));
  CUDA_CHECK(cudaSetDevice(ctx->device));
  if (ctx->dense) {
    CUSOLVER_CHECK(cusolverDnDestroy(ctx->dense));
  } else {
    CUSOLVER_CHECK(cusolverMpDestroyGrid(ctx->grid));
    CUSOLVER_CHECK(cusolverMpDestroy(ctx->handle));
#if defined(__CUSOLVERMP_NCCL)
    NCCL_CHECK(ncclCommDestroy(ctx->nccl_comm));
#else
    CAL_CHECK(cal_comm_destroy(ctx->cal_comm));
#endif
  }
  CUDA_CHECK(cudaStreamDestroy(ctx->stream));
  MPI_Comm_free(&ctx->comm);
  free(ctx);
  *context = NULL;
  CUDA_CHECK(cudaSetDevice(device));
}

/*******************************************************************************
 * \brief Solve on a single device, overwriting A with its eigenvectors.
 ******************************************************************************/
static void solve_local(cusolverDnHandle_t h, int n, int lda, int nev, void *a,
                        void *b, double *w, int *info, cudaDataType_t dtype,
                        size_t elem, void **work, int *found) {
  const cusolverEigMode_t job = CUSOLVER_EIG_MODE_VECTOR;
  const int subset = b && nev > 0 && nev < n;
  int lwork = 0;
#define DN(r, z, ...)                                                          \
  (dtype == CUDA_R_64F ? cusolverDn##r(__VA_ARGS__)                            \
                       : cusolverDn##z(__VA_ARGS__))
  if (subset) {
    CUSOLVER_CHECK(DN(Dsygvdx_bufferSize, Zhegvdx_bufferSize, h,
                      CUSOLVER_EIG_TYPE_1, job, CUSOLVER_EIG_RANGE_I,
                      CUBLAS_FILL_MODE_LOWER, n, a, lda, b, lda, 0., 0., 1, nev,
                      found, w, &lwork));
  } else if (b) {
    CUSOLVER_CHECK(DN(Dsygvd_bufferSize, Zhegvd_bufferSize, h,
                      CUSOLVER_EIG_TYPE_1, job, CUBLAS_FILL_MODE_LOWER, n, a,
                      lda, b, lda, w, &lwork));
  } else {
    CUSOLVER_CHECK(DN(Dsyevd_bufferSize, Zheevd_bufferSize, h, job,
                      CUBLAS_FILL_MODE_UPPER, n, a, lda, w, &lwork));
  }
  cp_assert(lwork > 0);
  CUDA_CHECK(cudaMalloc(work, (size_t)lwork * elem));
  if (subset) {
    CUSOLVER_CHECK(DN(Dsygvdx, Zhegvdx, h, CUSOLVER_EIG_TYPE_1, job,
                      CUSOLVER_EIG_RANGE_I, CUBLAS_FILL_MODE_LOWER, n, a, lda,
                      b, lda, 0., 0., 1, nev, found, w, *work, lwork, info));
  } else if (b) {
    CUSOLVER_CHECK(DN(Dsygvd, Zhegvd, h, CUSOLVER_EIG_TYPE_1, job,
                      CUBLAS_FILL_MODE_LOWER, n, a, lda, b, lda, w, *work,
                      lwork, info));
  } else {
    CUSOLVER_CHECK(DN(Dsyevd, Zheevd, h, job, CUBLAS_FILL_MODE_UPPER, n, a, lda,
                      w, *work, lwork, info));
  }
#undef DN
}

/*******************************************************************************
 * \brief Solve a standard or type-1 generalized FP64 eigenproblem. Matrices
 * share a block distribution; output on the host may have a different stride.
 ******************************************************************************/
void cp_fm_cusolver_solve(void **context, int fortran_comm, const int desc[9],
                          int nprow, int npcol, int myprow, int mypcol,
                          const void *a, const void *b, int ldb, void *z,
                          int ldz, double *w, int nev, int is_complex) {
  const int n = desc[2], mb = desc[4], nb = desc[5];
  const int generalized = ldb >= 0;
  // The typed dense routines query 32-bit workspace sizes. Leave large
  // singleton problems on the distributed API with size_t workspaces.
  const int local = nprow == 1 && npcol == 1 && n <= 16384;
  const cudaDataType_t dtype = is_complex ? CUDA_C_64F : CUDA_R_64F;
  cp_fm_cusolver_context *ctx =
      get_context(context, fortran_comm, nprow, npcol, local);
  const int rsrc = desc[6], csrc = desc[7], lda = desc[8];
  cp_assert(n > 0 && desc[3] == n && nev >= 0 && nev <= n);
  cp_assert(mb > 0 && nb > 0 && rsrc == csrc);
  const int rows = cusolverMpNUMROC(n, mb, myprow, rsrc, nprow);
  const int cols = cusolverMpNUMROC(n, nb, mypcol, csrc, npcol);
  cp_assert(lda >= 1 && lda >= rows && ldz >= rows && cols >= 0);
  cp_assert(!generalized || ldb >= rows);
  const size_t elem =
      dtype == CUDA_C_64F ? sizeof(cuDoubleComplex) : sizeof(double);
  cp_assert((size_t)cols <= SIZE_MAX / elem / lda);
  const size_t bytes = (size_t)lda * cols * elem;

  cusolverMpMatrixDescriptor_t layout = NULL;
  if (!local)
    CUSOLVER_CHECK(cusolverMpCreateMatrixDesc(&layout, ctx->grid, dtype, n, n,
                                              mb, nb, rsrc, csrc, lda));

  // One allocation for A, Z, optional B, W and info. Each matrix starts at
  // a cudaMalloc-aligned address, even on a rank with no local columns.
  const size_t align = 256;
  cp_assert(bytes <= SIZE_MAX - align);
  const size_t pitch = bytes ? (bytes + align - 1) / align * align : align;
  const size_t wbytes = (size_t)n * sizeof(double);
  const int nmats = 1 + !local + generalized;
  cp_assert(pitch <= (SIZE_MAX - wbytes - sizeof(int)) / nmats);
  unsigned char *storage = NULL;
  CUDA_CHECK(
      cudaMalloc((void **)&storage, nmats * pitch + wbytes + sizeof(int)));
  void *a_dev = storage, *z_dev = local ? storage : storage + pitch;
  void *b_dev = generalized ? storage + (nmats - 1) * pitch : NULL;
  double *w_dev = (double *)(storage + nmats * pitch);
  int *info_dev = (int *)(w_dev + n);

  if (bytes && rows > 0) {
    CUDA_CHECK(
        cudaMemcpyAsync(a_dev, a, bytes, cudaMemcpyHostToDevice, ctx->stream));
    if (generalized) {
      CUDA_CHECK(cudaMemcpy2DAsync(b_dev, lda * elem, b, ldb * elem,
                                   rows * elem, cols, cudaMemcpyHostToDevice,
                                   ctx->stream));
    }
  }

  // Older cuSOLVERMp versions need valid matrix pointers for the standard
  // query.
  void *dwork = NULL, *hwork = NULL;
  int found = nev;
  if (local) {
    solve_local(ctx->dense, n, lda, nev, a_dev, b_dev, w_dev, info_dev, dtype,
                elem, &dwork, &found);
  } else {
    size_t dsize, hsize;
    if (generalized) {
      CUSOLVER_CHECK(cusolverMpSygvd_bufferSize(
          ctx->handle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR,
          CUBLAS_FILL_MODE_LOWER, n, 1, 1, layout, 1, 1, layout, 1, 1, layout,
          dtype, &dsize, &hsize));
    } else {
      CUSOLVER_CHECK(cusolverMpSyevd_bufferSize(
          ctx->handle, "V", CUBLAS_FILL_MODE_UPPER, n, a_dev, 1, 1, layout,
          w_dev, z_dev, 1, 1, layout, dtype, &dsize, &hsize));
    }
    CUDA_CHECK(cudaMalloc(&dwork, dsize ? dsize : 1));
    CUDA_CHECK(cudaMallocHost(&hwork, hsize ? hsize : 1));

    if (generalized) {
      CUSOLVER_CHECK(cusolverMpSygvd(
          ctx->handle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR,
          CUBLAS_FILL_MODE_LOWER, n, a_dev, 1, 1, layout, b_dev, 1, 1, layout,
          w_dev, z_dev, 1, 1, layout, dtype, dwork, dsize, hwork, hsize,
          info_dev));
    } else {
      CUSOLVER_CHECK(cusolverMpSyevd(ctx->handle, "V", CUBLAS_FILL_MODE_UPPER,
                                     n, a_dev, 1, 1, layout, w_dev, z_dev, 1, 1,
                                     layout, dtype, dwork, dsize, hwork, hsize,
                                     info_dev));
    }

#if !defined(__CUSOLVERMP_NCCL)
    CUDA_CHECK(cudaStreamSynchronize(ctx->stream));
    CAL_CHECK(cal_stream_sync(ctx->cal_comm, ctx->stream));
#endif
  }
  int info = -1;
  CUDA_CHECK(cudaMemcpyAsync(&info, info_dev, sizeof(info),
                             cudaMemcpyDeviceToHost, ctx->stream));
  // The generalized interface returns only the requested local columns.
  const int zcols =
      generalized ? cusolverMpNUMROC(nev, nb, mypcol, csrc, npcol) : cols;
  if (rows > 0 && zcols > 0) {
    CUDA_CHECK(cudaMemcpy2DAsync(z, ldz * elem, z_dev, lda * elem, rows * elem,
                                 zcols, cudaMemcpyDeviceToHost, ctx->stream));
  }
  if (nev > 0) {
    CUDA_CHECK(cudaMemcpyAsync(w, w_dev, nev * sizeof(double),
                               cudaMemcpyDeviceToHost, ctx->stream));
  }
  CUDA_CHECK(cudaStreamSynchronize(ctx->stream));
  if (info != 0 || found != nev) {
    fprintf(stderr,
            "ERROR: cuSOLVER %s failed: info = %d, requested %d, found %d\n",
            generalized ? "generalized solve" : "standard solve", info, nev,
            found);
    abort();
  }

  if (hwork)
    CUDA_CHECK(cudaFreeHost(hwork));
  CUDA_CHECK(cudaFree(dwork));
  CUDA_CHECK(cudaFree(storage));
  if (layout)
    CUSOLVER_CHECK(cusolverMpDestroyMatrixDesc(layout));
}

#endif

// EOF
