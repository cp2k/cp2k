/******************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations
 *  Copyright (C) 2000 - 2013  CP2K developers group
 *
 *  Authors: Benjamin G Levine, Andreas Gloess
 *
 *  2012/05/18                 Refacturing - original files:
 *                              - cuda_tools/cuda_pw_cu.cu
 *  
 *****************************************************************************/
#if defined ( __PW_CUDA )

// global dependencies
#include <cuda_runtime.h>
#include <cufft.h>
#include <cublas.h>
#include <stdio.h>

// local dependencies
#include "memory_cuda.h"
#include "error_cuda.h"
#include "fft_cuda.h"

// debug flag
#define CHECK 1

// configuration(s)
#define NTHREADS 32
#define MAXTHREADS 1024
#define MAXGRIDX 65535

// helper routine(s)
void get_grid_params(const int   ngpts,
                     const int   blocksize,
                           dim3 &threadsPerBlock,
                           dim3 &blocksPerGrid) {
  int blocks;
  if (blocksize <= MAXTHREADS) {
    threadsPerBlock.x = blocksize;
  } else {
    threadsPerBlock.x = MAXTHREADS;
    printf("WARNING: Number of threads per block (x) is too large!\n");
    printf("WARNING: Number of threads per block (x) is set to: %d.\n", MAXTHREADS);
  }
  threadsPerBlock.y = 1;
  threadsPerBlock.z = 1;
  blocks = (ngpts + threadsPerBlock.x - 1) / threadsPerBlock.x;
  blocksPerGrid.x = (int) ceil(sqrt((double) blocks));
  if (blocksPerGrid.x > MAXGRIDX) {
    printf("CUDA: Not allowed grid dimensions!\n");
    exit(1);
  }
  blocksPerGrid.y = (int) rint(sqrt((double) blocks));
  blocksPerGrid.z = 1;
}

void cudaStreamBarrier(cudaStream_t cuda_stream) {
  cudaError_t cErr;

  // stop on previous errors
  cErr = cudaGetLastError();
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  // might result in endless loop
  do {
    cErr = cudaStreamQuery(cuda_stream); 
  } while (cErr != cudaSuccess);
}

// --- CODE -------------------------------------------------------------------

/******************************************************************************
 * \brief   Performs a out-of-place copy of a double precision vector (first
 *          half filled) into a double precision complex vector on the GPU.
 *          It requires a global double precision vector 'zout' of lenght '2n'.
 *          [memory (shared):  none Byte
 *           memory (private): 4 Byte
 *           memory (global):  16*n Byte]
 *          n - size of double precision input vector
 * \author  Andreas Gloess
 * \date    2013-03-07
 * \version 0.01
 *****************************************************************************/
__global__ void pw_copy_rc_cu_z(const double *din,
                                      double *zout,
                                const int     n) {
  const int igpt = (gridDim.x * blockIdx.y + blockIdx.x) * blockDim.x + threadIdx.x;

  if (igpt < n) {
     zout[2 * igpt    ] = din[igpt];
     zout[2 * igpt + 1] = 0.0e0;
  }
}


/******************************************************************************
 * \brief   Performs a out-of-place copy of a double precision complex vector
 *          (real part) into a double precision vector on the GPU.
 *          It requires a global double precision vector 'dout' of lenght 'n'.
 *          [memory (shared):  none Byte
 *           memory (private): 4 Byte
 *           memory (global):  16*n Byte]
 *          n - size of double precision output vector
 * \author  Andreas Gloess
 * \date    2013-03-07
 * \version 0.01
 *****************************************************************************/
__global__ void pw_copy_cr_cu_z(const double *zin,
                                      double *dout,
                                const int     n) {
  const int igpt = (gridDim.x * blockIdx.y + blockIdx.x) * blockDim.x + threadIdx.x;

  if (igpt < n) {
     dout[igpt] = zin[2 * igpt];
  }
}


/******************************************************************************
 * \brief   Performs a (double precision complex) gather and scale on the GPU.
 * \author  Andreas Gloess
 * \date    2013-03-07
 * \version 0.01
 *****************************************************************************/
__global__ void pw_gather_cu_z(      double *pwcc,
                               const double *c,
                               const double  scale,
                               const int     ngpts,
                               const int    *ghatmap) {

  const int igpt = (gridDim.x * blockIdx.y + blockIdx.x) * blockDim.x + threadIdx.x;

  if (igpt < ngpts) {
    pwcc[2 * igpt    ] = scale * c[2 * ghatmap[igpt]    ];
    pwcc[2 * igpt + 1] = scale * c[2 * ghatmap[igpt] + 1];
  }
}


/******************************************************************************
 * \brief   Performs a (double precision complex) scatter and scale on the GPU.
 * \author  Andreas Gloess
 * \date    2013-03-07
 * \version 0.01
 *****************************************************************************/
__global__ void pw_scatter_cu_z(      double *c,
                                const double *pwcc,
                                const double  scale,
                                const int     ngpts,
                                const int     nmaps,
                                const int    *ghatmap) {

  const int igpt = (gridDim.x * blockIdx.y + blockIdx.x) * blockDim.x + threadIdx.x;

  if (igpt < ngpts) {
    c[2 * ghatmap[igpt]    ] = scale * pwcc[2 * igpt    ];
    c[2 * ghatmap[igpt] + 1] = scale * pwcc[2 * igpt + 1];
    if (nmaps == 2) {
      c[2 * ghatmap[igpt + ngpts]    ] =   scale * pwcc[2 * igpt    ];
      c[2 * ghatmap[igpt + ngpts] + 1] = - scale * pwcc[2 * igpt + 1];
    }
  }
}


/******************************************************************************
 * \brief   Performs a (double precision complex) FFT, followed by a (double
 *          precision complex) gather, on the GPU.
 * \author  Andreas Gloess
 * \date    2013-03-07
 * \version 0.01
 *****************************************************************************/
extern "C" void pw_cuda_cfffg_z_(const double          *din,
                                       cuDoubleComplex *zout,
                                 const int             *ghatmap,
                                 const int             *npts,
                                 const int              ngpts,
                                 const double           scale) {
  double *ptr_1, *ptr_2;
  int    *ghatmap_dev;
  int     nrpts;
  dim3    blocksPerGrid, threadsPerBlock;
  cudaStream_t *cuda_streams;
  cudaEvent_t  *cuda_events;
  cudaError_t   cErr;

  // dimensions of double and complex arrays
  nrpts = npts[0] * npts[1] * npts[2];

  // get streams
  cuda_get_streams_cu_(&cuda_streams);
  cuda_get_events_cu_(&cuda_events);

  // get device memory pointers
  cuda_device_mem_alloc_cu_(&ptr_1,       2 * nrpts);
  cuda_device_mem_alloc_cu_(&ptr_2,       2 * nrpts); //? (ngpts)
  cuda_device_mem_alloc_cu_(&ghatmap_dev,     ngpts);

  // convert the real (host) pointer 'din' into a complex (device)
  // pointer 'ptr_1'
  // copy to device (NOTE: only first half of ptr_1 is written!)
  cErr = cudaMemcpyAsync(ptr_1, din, sizeof(double) * nrpts, cudaMemcpyHostToDevice, cuda_streams[0]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  cErr = cudaEventRecord(cuda_events[0], cuda_streams[0]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // CUDA blocking for pw_copy_rc (currently only 2-D grid)
  get_grid_params(nrpts, MAXTHREADS, threadsPerBlock, blocksPerGrid);

  // real to complex blow-up
  cErr = cudaStreamWaitEvent(cuda_streams[1], cuda_events[0], 0);
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  pw_copy_rc_cu_z<<<blocksPerGrid, threadsPerBlock, 0, cuda_streams[1]>>>(ptr_1, ptr_2, nrpts);
  cErr = cudaGetLastError();
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // copy gather map array from host to the device
  cErr = cudaMemcpyAsync(ghatmap_dev, ghatmap, sizeof(int) * ngpts, cudaMemcpyHostToDevice, cuda_streams[0]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  cErr = cudaEventRecord(cuda_events[0], cuda_streams[0]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // fft on the GPU (cuda_streams[1])
  fftcu_run_3d_z_(+1, npts, 1.0e0, (cufftDoubleComplex *) ptr_2, cuda_streams[1]);

  // CUDA blocking for gather (currently only 2-D grid)
  get_grid_params(ngpts, NTHREADS, threadsPerBlock, blocksPerGrid);

  // gather on the GPU
  cErr = cudaStreamWaitEvent(cuda_streams[1], cuda_events[0], 0);
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  pw_gather_cu_z<<<blocksPerGrid, threadsPerBlock, 0, cuda_streams[1]>>>(ptr_1, ptr_2, scale, ngpts, ghatmap_dev);
  cErr = cudaGetLastError();
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  cErr = cudaEventRecord(cuda_events[1], cuda_streams[1]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // get results from device memory
  cErr = cudaStreamWaitEvent(cuda_streams[2], cuda_events[1], 0);
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  cErr = cudaMemcpyAsync(zout, ptr_1, sizeof(cuDoubleComplex) * ngpts, cudaMemcpyDeviceToHost, cuda_streams[2]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // synchronize with respect to host
  cErr = cudaStreamSynchronize(cuda_streams[2]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // release memory stack
  cuda_device_mem_free_cu_(&ptr_1);
  cuda_device_mem_free_cu_(&ptr_2);
  cuda_device_mem_free_cu_(&ghatmap_dev);
}


/******************************************************************************
 * \brief   Performs a (double precision complex) scatter, followed by a
 *          (double precision complex) FFT, on the GPU.
 * \author  Andreas Gloess
 * \date    2013-03-07
 * \version 0.01
 *****************************************************************************/
extern "C" void pw_cuda_sfffc_z_(const cuDoubleComplex *zin,
                                       double          *dout,
                                 const int             *ghatmap,
                                 const int             *npts,
                                 const int              ngpts,
                                 const int              nmaps,
                                 const double           scale) {
  double *ptr_1, *ptr_2;
  int    *ghatmap_dev;
  int    nrpts;
  dim3   blocksPerGrid, threadsPerBlock;
  cudaStream_t *cuda_streams;
  cudaEvent_t  *cuda_events;
  cudaError_t   cErr;

  // dimensions of double and complex arrays
  nrpts = npts[0] * npts[1] * npts[2];

  // get streams
  cuda_get_streams_cu_(&cuda_streams);
  cuda_get_events_cu_(&cuda_events);

  // get device memory pointers
  cuda_device_mem_alloc_cu_(&ptr_1,       2 * nrpts);
  cuda_device_mem_alloc_cu_(&ptr_2,       2 * nrpts);
  cuda_device_mem_alloc_cu_(&ghatmap_dev, nmaps * ngpts);
  
  // copy all arrays from host to the device
  cErr = cudaMemcpyAsync(ptr_1, zin, sizeof(cuDoubleComplex) * ngpts, cudaMemcpyHostToDevice, cuda_streams[0]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  cErr = cudaMemcpyAsync(ghatmap_dev, ghatmap, sizeof(int) * nmaps * ngpts, cudaMemcpyHostToDevice, cuda_streams[0]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  cErr = cudaEventRecord(cuda_events[0], cuda_streams[0]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // CUDA blocking for scatter (currently only 2-D grid)
  get_grid_params(ngpts, NTHREADS, threadsPerBlock, blocksPerGrid);

  // scatter on the GPU
  cErr = cudaStreamWaitEvent(cuda_streams[1], cuda_events[0], 0);
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  cErr = cudaMemsetAsync(ptr_2, 0, sizeof(double) * 2 * nrpts, cuda_streams[1]); // we need to do this only if spherical cut-off is used!
  if (CHECK) cuda_error_check2(cErr, __LINE__);                                  // but it turns out to be performance irrelevant
  pw_scatter_cu_z<<<blocksPerGrid, threadsPerBlock, 0, cuda_streams[1]>>>(ptr_2, ptr_1, scale, ngpts, nmaps, ghatmap_dev);
  cErr = cudaGetLastError();
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // fft on the GPU (cuda_streams[1])
  fftcu_run_3d_z_(-1, npts, 1.0e0, (cufftDoubleComplex *) ptr_2, cuda_streams[1]);

  // CUDA blocking for pw_copy_cr (currently only 2-D grid)
  get_grid_params(nrpts, MAXTHREADS, threadsPerBlock, blocksPerGrid);

  // convert the complex (device) pointer 'ptr_2' into a real (host)
  // pointer 'dout' (NOTE: Only first half of ptr_1 is written!)
  pw_copy_cr_cu_z<<<blocksPerGrid, threadsPerBlock, 0, cuda_streams[1]>>>(ptr_2, ptr_1, nrpts);
  cErr = cudaGetLastError();
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  cErr = cudaEventRecord(cuda_events[1], cuda_streams[1]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // get results from device memory
  cErr = cudaStreamWaitEvent(cuda_streams[2], cuda_events[1], 0);
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  cErr = cudaMemcpyAsync(dout, ptr_1, sizeof(double) * nrpts, cudaMemcpyDeviceToHost, cuda_streams[2]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // synchronize with respect to host
  cErr = cudaStreamSynchronize(cuda_streams[2]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // release memory stack
  cuda_device_mem_free_cu_(&ptr_1);
  cuda_device_mem_free_cu_(&ptr_2);
  cuda_device_mem_free_cu_(&ghatmap_dev);
}


/******************************************************************************
 * \brief   Performs a (double to complex double) blow-up and a (double
 *          precision complex) 2D-FFT on the GPU.
 * \author  Andreas Gloess
 * \date    2013-03-07
 * \version 0.01
 *****************************************************************************/
extern "C" void pw_cuda_cff_z_(const double          *din,
                                     cuDoubleComplex *zout,
                               const int             *npts) {
  double *ptr_1, *ptr_2;
  int     nrpts;
  dim3    blocksPerGrid, threadsPerBlock;
  cudaStream_t *cuda_streams;
  cudaEvent_t  *cuda_events;
  cudaError_t   cErr;

  // dimensions of double and complex arrays
  nrpts  = npts[0] * npts[1] * npts[2];

  // get streams
  cuda_get_streams_cu_(&cuda_streams);
  cuda_get_events_cu_(&cuda_events);

  // get device memory pointers for:
  cuda_device_mem_alloc_cu_(&ptr_1, 2 * nrpts);
  cuda_device_mem_alloc_cu_(&ptr_2, 2 * nrpts);

  // convert the real (host) pointer 'din' into a complex (device)
  // pointer 'ptr_in'
  // copy all arrays from host to the device (NOTE: Only first half of ptr_1 is written!)
  cErr = cudaMemcpyAsync(ptr_1, din, sizeof(double) * nrpts, cudaMemcpyHostToDevice, cuda_streams[0]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  cErr = cudaEventRecord(cuda_events[0], cuda_streams[0]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // CUDA blocking for pw_copy_rc (currently only 2-D grid)
  get_grid_params(nrpts, MAXTHREADS, threadsPerBlock, blocksPerGrid);

  // real to complex blow-up
  cErr = cudaStreamWaitEvent(cuda_streams[1], cuda_events[0], 0);
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  pw_copy_rc_cu_z<<<blocksPerGrid, threadsPerBlock, 0, cuda_streams[1]>>>(ptr_1, ptr_2, nrpts);
  cErr = cudaGetLastError();
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // fft on the GPU (cuda_streams[1])
  //NOTE: the following works, but CUDA does 2D-FFT in C-shaped (not optimal) order
  //fftcu_run_2dm_z_(1, npts, 1.0e0, (cufftDoubleComplex *) ptr_2, (cufftDoubleComplex *) ptr_1, cuda_streams[1]);
  fftcu_run_1dm_z_(1, npts[2], npts[0]*npts[1], 1.0e0, (cufftDoubleComplex *) ptr_2, (cufftDoubleComplex *) ptr_1, cuda_streams[1]);
  fftcu_run_1dm_z_(1, npts[1], npts[0]*npts[2], 1.0e0, (cufftDoubleComplex *) ptr_1, (cufftDoubleComplex *) ptr_2, cuda_streams[1]);
  cErr = cudaEventRecord(cuda_events[1], cuda_streams[1]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // get results from device memory
  cErr = cudaStreamWaitEvent(cuda_streams[2], cuda_events[1], 0);
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  //cErr = cudaMemcpyAsync(zout, ptr_1, sizeof(cuDoubleComplex) * nrpts, cudaMemcpyDeviceToHost, cuda_streams[2]);
  cErr = cudaMemcpyAsync(zout, ptr_2, sizeof(cuDoubleComplex) * nrpts, cudaMemcpyDeviceToHost, cuda_streams[2]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // synchronize with respect to host
  cErr = cudaStreamSynchronize(cuda_streams[2]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // release memory stack
  cuda_device_mem_free_cu_(&ptr_1);
  cuda_device_mem_free_cu_(&ptr_2);
}


/******************************************************************************
 * \brief   Performs a (double precision complex) 2D-FFT and a (double complex
 *          to double) shrink-down on the GPU.
 * \author  Andreas Gloess
 * \date    2013-03-07
 * \version 0.01
 *****************************************************************************/
extern "C" void pw_cuda_ffc_z_(const cuDoubleComplex *zin,
                                     double          *dout,
                               const int             *npts) {
  double *ptr_1, *ptr_2;
  int     nrpts;
  dim3    blocksPerGrid, threadsPerBlock;
  cudaStream_t *cuda_streams;
  cudaEvent_t  *cuda_events;
  cudaError_t   cErr;

  // dimensions of double and complex arrays
  nrpts  = npts[0] * npts[1] * npts[2];

  // get streams
  cuda_get_streams_cu_(&cuda_streams);
  cuda_get_events_cu_(&cuda_events);

  // get device memory pointers for:
  cuda_device_mem_alloc_cu_(&ptr_1, 2 * nrpts);
  cuda_device_mem_alloc_cu_(&ptr_2, 2 * nrpts);
  
  // copy input data from host to the device
  cErr = cudaMemcpyAsync(ptr_1, zin, sizeof(cuDoubleComplex) * nrpts, cudaMemcpyHostToDevice, cuda_streams[0]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  cErr = cudaEventRecord(cuda_events[0], cuda_streams[0]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // fft on the GPU (cuda_stream[1])
  cErr = cudaStreamWaitEvent(cuda_streams[1], cuda_events[0], 0);
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  //fftcu_run_2dm_z_(-1, npts, 1.0e0, (cufftDoubleComplex *) ptr_1, (cufftDoubleComplex *) ptr_2, cuda_streams[1]);
  fftcu_run_1dm_z_(-1, npts[1], npts[0]*npts[2], 1.0e0, (cufftDoubleComplex *) ptr_1, (cufftDoubleComplex *) ptr_2, cuda_streams[1]);
  fftcu_run_1dm_z_(-1, npts[2], npts[0]*npts[1], 1.0e0, (cufftDoubleComplex *) ptr_2, (cufftDoubleComplex *) ptr_1, cuda_streams[1]);

  // CUDA blocking for pw_copy_cr (currently only 2-D grid)
  get_grid_params(nrpts, MAXTHREADS, threadsPerBlock, blocksPerGrid);

  // convert the complex (device) pointer 'ptr_1' into a real (host)
  // pointer 'dout' (NOTE: Only first half of ptr_2 is written!)
  //pw_copy_cr_cu_z<<<blocksPerGrid, threadsPerBlock, 0, cuda_streams[1]>>>(ptr_2, ptr_1, nrpts);
  pw_copy_cr_cu_z<<<blocksPerGrid, threadsPerBlock, 0, cuda_streams[1]>>>(ptr_1, ptr_2, nrpts);
  cErr = cudaGetLastError();
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  cErr = cudaEventRecord(cuda_events[1], cuda_streams[1]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // get results from device memory
  cErr = cudaStreamWaitEvent(cuda_streams[2], cuda_events[1], 0);
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  //cErr = cudaMemcpyAsync(dout, ptr_1, sizeof(double) * nrpts, cudaMemcpyDeviceToHost, cuda_streams[2]);
  cErr = cudaMemcpyAsync(dout, ptr_2, sizeof(double) * nrpts, cudaMemcpyDeviceToHost, cuda_streams[2]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // synchronize with respect to host
  cErr = cudaStreamSynchronize(cuda_streams[2]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // release memory stack
  cuda_device_mem_free_cu_(&ptr_1);
  cuda_device_mem_free_cu_(&ptr_2);
}


/******************************************************************************
 * \brief   Performs a (double to complex double) blow-up and a (double
 *          precision complex) 1D-FFT on the GPU.
 * \author  Andreas Gloess
 * \date    2013-03-07
 * \version 0.01
 *****************************************************************************/
extern "C" void pw_cuda_cf_z_(const double          *din,
                                    cuDoubleComplex *zout,
                              const int             *npts) {
  double *ptr_1, *ptr_2;
  int     nrpts;
  dim3    blocksPerGrid, threadsPerBlock;
  cudaStream_t *cuda_streams;
  cudaEvent_t  *cuda_events;
  cudaError_t   cErr;

  // dimensions of double and complex arrays
  nrpts  = npts[0] * npts[1] * npts[2];

  // get streams
  cuda_get_streams_cu_(&cuda_streams);
  cuda_get_events_cu_(&cuda_events);

  // get device memory pointers for:
  cuda_device_mem_alloc_cu_(&ptr_1, 2 * nrpts);
  cuda_device_mem_alloc_cu_(&ptr_2, 2 * nrpts);

  // convert the real (host) pointer 'din' into a complex (device)
  // pointer 'ptr_2' (NOTE: Only first half of ptr_1 is written!)
  // copy all arrays from host to the device
  cErr = cudaMemcpyAsync(ptr_1, din, sizeof(double) * nrpts, cudaMemcpyHostToDevice, cuda_streams[0]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  cErr = cudaEventRecord(cuda_events[0], cuda_streams[0]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // CUDA blocking for pw_copy_rc (currently only 2-D grid)
  get_grid_params(nrpts, MAXTHREADS, threadsPerBlock, blocksPerGrid);

  // real to complex blow-up
  cErr = cudaStreamWaitEvent(cuda_streams[1], cuda_events[0], 0);
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  pw_copy_rc_cu_z<<<blocksPerGrid, threadsPerBlock, 0, cuda_streams[1]>>>(ptr_1, ptr_2, nrpts);
  cErr = cudaGetLastError();
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // fft on the GPU (cuda_streams[1])
  fftcu_run_1dm_z_(1, npts[2], npts[0]*npts[1], 1.0e0, (cufftDoubleComplex *) ptr_2, (cufftDoubleComplex *) ptr_1, cuda_streams[1]);
  cErr = cudaEventRecord(cuda_events[1], cuda_streams[1]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // get results from device memory
  cErr = cudaStreamWaitEvent(cuda_streams[2], cuda_events[1], 0);
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  cErr = cudaMemcpyAsync(zout, ptr_1, sizeof(cuDoubleComplex) * nrpts, cudaMemcpyDeviceToHost, cuda_streams[2]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // synchronize with respect to host
  cErr = cudaStreamSynchronize(cuda_streams[2]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // release memory stack
  cuda_device_mem_free_cu_(&ptr_1);
  cuda_device_mem_free_cu_(&ptr_2);
}


/******************************************************************************
 * \brief   Performs a (double precision complex) 1D-FFT and a (double complex
 *          to double) shrink-down on the GPU.
 * \author  Andreas Gloess
 * \date    2013-03-07
 * \version 0.01
 *****************************************************************************/
extern "C" void pw_cuda_fc_z_(const cuDoubleComplex *zin,
                                    double          *dout,
                              const int             *npts) {
  double *ptr_1, *ptr_2;
  int     nrpts;
  dim3    blocksPerGrid, threadsPerBlock;
  cudaStream_t *cuda_streams;
  cudaEvent_t  *cuda_events;
  cudaError_t   cErr;

  // dimensions of double and complex arrays
  nrpts  = npts[0] * npts[1] * npts[2];

  // get streams
  cuda_get_streams_cu_(&cuda_streams);
  cuda_get_events_cu_(&cuda_events);

  // get device memory pointers for:
  cuda_device_mem_alloc_cu_(&ptr_1, 2 * nrpts);
  cuda_device_mem_alloc_cu_(&ptr_2, 2 * nrpts);
  
  // copy input data from host to the device
  cErr = cudaMemcpyAsync(ptr_1, zin, sizeof(cuDoubleComplex) * nrpts, cudaMemcpyHostToDevice, cuda_streams[0]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  cErr = cudaEventRecord(cuda_events[0], cuda_streams[0]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // fft on the GPU (cuda_stream[1])
  cErr = cudaStreamWaitEvent(cuda_streams[1], cuda_events[0], 0);
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  fftcu_run_1dm_z_(-1, npts[2], npts[0]*npts[1], 1.0e0, (cufftDoubleComplex *) ptr_1, (cufftDoubleComplex *) ptr_2, cuda_streams[1]);

  // CUDA blocking for pw_copy_cr (currently only 2-D grid)
  get_grid_params(nrpts, MAXTHREADS, threadsPerBlock, blocksPerGrid);

  // convert the complex (device) pointer 'ptr_2' into a real (host)
  // pointer 'dout' (NOTE: Only first half of ptr_1 is written!)
  pw_copy_cr_cu_z<<<blocksPerGrid, threadsPerBlock, 0, cuda_streams[1]>>>(ptr_2, ptr_1, nrpts);
  cErr = cudaGetLastError();
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  cErr = cudaEventRecord(cuda_events[1], cuda_streams[1]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // get results from device memory
  cErr = cudaStreamWaitEvent(cuda_streams[2], cuda_events[1], 0);
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  cErr = cudaMemcpyAsync(dout, ptr_1, sizeof(double) * nrpts, cudaMemcpyDeviceToHost, cuda_streams[2]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // synchronize with respect to host
  cErr = cudaStreamSynchronize(cuda_streams[2]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // release memory stack
  cuda_device_mem_free_cu_(&ptr_1);
  cuda_device_mem_free_cu_(&ptr_2);
}


/******************************************************************************
 * \brief   Performs a (double precision complex) 1D-FFT on the GPU.
 * \author  Andreas Gloess
 * \date    2013-05-01
 * \version 0.01
 *****************************************************************************/
extern "C" void pw_cuda_f_z_(const cuDoubleComplex *zin,
                                   cuDoubleComplex *zout,
                             const int              dir,
                             const int              n,
                             const int              m) {
  double *ptr_1, *ptr_2;
  int     nrpts;
  cudaStream_t *cuda_streams;
  cudaEvent_t  *cuda_events;
  cudaError_t   cErr;

  // dimensions of complex arrays
  nrpts  = n * m;

  // get streams
  cuda_get_streams_cu_(&cuda_streams);
  cuda_get_events_cu_(&cuda_events);

  // get device memory pointers for:
  cuda_device_mem_alloc_cu_(&ptr_1, 2 * nrpts);
  cuda_device_mem_alloc_cu_(&ptr_2, 2 * nrpts);
  
  // copy input data from host to the device
  cErr = cudaMemcpyAsync(ptr_1, zin, sizeof(cuDoubleComplex) * nrpts, cudaMemcpyHostToDevice, cuda_streams[0]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  cErr = cudaEventRecord(cuda_events[0], cuda_streams[0]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // fft on the GPU (cuda_stream[1])
  cErr = cudaStreamWaitEvent(cuda_streams[1], cuda_events[0], 0);
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  fftcu_run_1dm_z_(dir, n, m, 1.0e0, (cufftDoubleComplex *) ptr_1, (cufftDoubleComplex *) ptr_2, cuda_streams[1]);
  cErr = cudaEventRecord(cuda_events[1], cuda_streams[1]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // get results from device memory
  cErr = cudaStreamWaitEvent(cuda_streams[2], cuda_events[1], 0);
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  cErr = cudaMemcpyAsync(zout, ptr_2, sizeof(cuDoubleComplex) * nrpts, cudaMemcpyDeviceToHost, cuda_streams[2]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // synchronize with respect to host
  cErr = cudaStreamSynchronize(cuda_streams[2]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // release memory stack
  cuda_device_mem_free_cu_(&ptr_1);
  cuda_device_mem_free_cu_(&ptr_2);
}


/******************************************************************************
 * \brief   Performs a (double precision complex) 1D-FFT, followed by a (double
 *          precision complex) gather, on the GPU.
 * \author  Andreas Gloess
 * \date    2013-03-07
 * \version 0.01
 *****************************************************************************/
extern "C" void pw_cuda_fg_z_(const cuDoubleComplex *zin,
                                    cuDoubleComplex *zout,
                              const int             *ghatmap,
                              const int             *npts,
                              const int              mmax,
                              const int              ngpts,
                              const double           scale) {
  double *ptr_1, *ptr_2;
  int    *ghatmap_dev;
  int     nrpts;
  dim3    blocksPerGrid, threadsPerBlock;
  cudaStream_t *cuda_streams;
  cudaEvent_t  *cuda_events;
  cudaError_t   cErr;

  // dimensions of double and complex arrays
  nrpts = npts[0] * mmax;

  // get streams and events
  cuda_get_streams_cu_(&cuda_streams);
  cuda_get_events_cu_(&cuda_events);

  // get device memory pointers
  cuda_device_mem_alloc_cu_(&ptr_1,      2 * nrpts);
  cuda_device_mem_alloc_cu_(&ptr_2,      2 * nrpts);
  cuda_device_mem_alloc_cu_(&ghatmap_dev,    ngpts);

  // transfer gather data from host to device
  cErr = cudaMemcpyAsync(ghatmap_dev, ghatmap, sizeof(int) * ngpts, cudaMemcpyHostToDevice, cuda_streams[0]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // transfer input data from host to device
  cErr = cudaMemcpyAsync(ptr_1, zin, sizeof(cuDoubleComplex) * nrpts, cudaMemcpyHostToDevice, cuda_streams[0]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  cErr = cudaEventRecord(cuda_events[0], cuda_streams[0]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // fft on the GPU (cuda_streams[1])
  cErr = cudaStreamWaitEvent(cuda_streams[1], cuda_events[0], 0);
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  fftcu_run_1dm_z_(1, npts[0], mmax, 1.0e0, (cufftDoubleComplex *) ptr_1, (cufftDoubleComplex *) ptr_2, cuda_streams[1]);

  // CUDA blocking for gather (currently only 2-D grid)
  get_grid_params(ngpts, NTHREADS, threadsPerBlock, blocksPerGrid);

  // gather on the GPU
  pw_gather_cu_z<<<blocksPerGrid, threadsPerBlock, 0, cuda_streams[1]>>>(ptr_1, ptr_2, scale, ngpts, ghatmap_dev);
  cErr = cudaGetLastError();
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  cErr = cudaEventRecord(cuda_events[1], cuda_streams[1]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // transfer results from device to host
  cErr = cudaStreamWaitEvent(cuda_streams[2], cuda_events[1], 0);
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  cErr = cudaMemcpyAsync(zout, ptr_1, sizeof(cuDoubleComplex) * ngpts, cudaMemcpyDeviceToHost, cuda_streams[2]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // synchronize with respect to host
  cErr = cudaStreamSynchronize(cuda_streams[2]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // release memory stack
  cuda_device_mem_free_cu_(&ptr_1);
  cuda_device_mem_free_cu_(&ptr_2);
  cuda_device_mem_free_cu_(&ghatmap_dev);
}


/******************************************************************************
 * \brief   Performs a (double precision complex) scatter, followed by a
 *          (double precision complex) 1D-FFT, on the GPU.
 * \author  Andreas Gloess
 * \date    2013-03-07
 * \version 0.01
 *****************************************************************************/
extern "C" void pw_cuda_sf_z_(const cuDoubleComplex *zin,
                                    cuDoubleComplex *zout,
                              const int             *ghatmap,
                              const int             *npts,
                              const int              mmax,
                              const int              ngpts,
                              const int              nmaps,
                              const double           scale) {
  double *ptr_1, *ptr_2;
  int    *ghatmap_dev;
  int    nrpts;
  dim3   blocksPerGrid, threadsPerBlock;
  cudaStream_t *cuda_streams;
  cudaEvent_t  *cuda_events;
  cudaError_t   cErr;

  // dimensions of double and complex arrays
  nrpts = npts[0] * mmax;

  // get streams
  cuda_get_streams_cu_(&cuda_streams);
  cuda_get_events_cu_(&cuda_events);

  // get device memory pointers
  cuda_device_mem_alloc_cu_(&ptr_1,       2 * nrpts);
  cuda_device_mem_alloc_cu_(&ptr_2,       2 * nrpts);
  cuda_device_mem_alloc_cu_(&ghatmap_dev, nmaps * ngpts);

  // transfer input data from host to the device
  cErr = cudaMemcpyAsync(ptr_1, zin, sizeof(cuDoubleComplex) * ngpts, cudaMemcpyHostToDevice, cuda_streams[0]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // transfer scatter data from host to device
  cErr = cudaMemcpyAsync(ghatmap_dev, ghatmap, sizeof(int) * nmaps * ngpts, cudaMemcpyHostToDevice, cuda_streams[0]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  cErr = cudaEventRecord(cuda_events[0], cuda_streams[0]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // CUDA blocking for scatter (currently only 2-D grid)
  get_grid_params(ngpts, NTHREADS, threadsPerBlock, blocksPerGrid);

  // scatter on the GPU
  cErr = cudaStreamWaitEvent(cuda_streams[1], cuda_events[0], 0);
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  cErr = cudaMemsetAsync(ptr_2, 0, sizeof(double) * 2 * nrpts, cuda_streams[1]); // we need to do this only if spherical cut-off is used!
  if (CHECK) cuda_error_check2(cErr, __LINE__);                                  // but it turns out to be performance irrelevant
  pw_scatter_cu_z<<<blocksPerGrid, threadsPerBlock, 0, cuda_streams[1]>>>(ptr_2, ptr_1, scale, ngpts, nmaps, ghatmap_dev);
  cErr = cudaGetLastError();
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // fft on the GPU (cuda_streams[1])
  fftcu_run_1dm_z_(-1, npts[0], mmax, 1.0e0, (cufftDoubleComplex *) ptr_2, (cufftDoubleComplex *) ptr_1, cuda_streams[1]);
  cErr = cudaGetLastError();
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  cErr = cudaEventRecord(cuda_events[1], cuda_streams[1]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // get results from device memory
  cErr = cudaStreamWaitEvent(cuda_streams[2], cuda_events[1], 0);
  if (CHECK) cuda_error_check2(cErr, __LINE__);
  cErr = cudaMemcpyAsync(zout, ptr_1, sizeof(cuDoubleComplex) * nrpts, cudaMemcpyDeviceToHost, cuda_streams[2]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // synchronize with respect to host
  cErr = cudaStreamSynchronize(cuda_streams[2]);
  if (CHECK) cuda_error_check2(cErr, __LINE__);

  // release memory stack
  cuda_device_mem_free_cu_(&ptr_1);
  cuda_device_mem_free_cu_(&ptr_2);
  cuda_device_mem_free_cu_(&ghatmap_dev);
}
#endif
