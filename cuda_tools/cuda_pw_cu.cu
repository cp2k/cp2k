#if defined ( __CUDA )

#include <cuda_runtime.h>
#include <cufft.h>
#include <cublas.h>
#include <stdio.h>

#include "cufft.h"
#include "cuda_memory_cu.h"

#define NTHREADS 64
#define NBLOCKS  64

__global__ void pw_gather_cu(float* pwcc, float* c, float scale, int* mapl, int* mapm, int* mapn, int ngpts, int* nn, int* ghat, int l1, int l2, int l3) {

  int igpt, n0, n1, limit, index;
  __shared__ int ghattmp[NTHREADS*3];
  __shared__ float tmp[NTHREADS*2];
  int nleftover;

  n0=nn[0];
  n1=nn[1];

  limit=ngpts/(NTHREADS*NBLOCKS);

  index = blockIdx.x * NTHREADS;
  for (igpt=0; igpt<limit; igpt++) {
    ghattmp[threadIdx.x] = ghat[ index * 3 + threadIdx.x ] ;
    ghattmp[threadIdx.x+NTHREADS] = ghat[ index * 3 + threadIdx.x + NTHREADS ]; 
    ghattmp[threadIdx.x+2*NTHREADS] = ghat[ index * 3 + threadIdx.x + 
                                            2 * NTHREADS ];
    __syncthreads();

    ghattmp[threadIdx.x * 3] = mapl[ ghattmp[ threadIdx.x * 3 ] - l1 ];
    ghattmp[threadIdx.x * 3 + 1] = mapm[ ghattmp[ threadIdx.x * 3 + 1 ] - l2 ];
    ghattmp[threadIdx.x * 3 + 2] = mapn[ ghattmp[ threadIdx.x * 3 + 2 ] - l3 ];

    ghattmp[threadIdx.x * 3] += n0 * ( n1 * ghattmp[threadIdx.x * 3 + 2] + 
                                ghattmp[threadIdx.x * 3 + 1]);

    tmp[ 2 * threadIdx.x ] = c[ 2 * ghattmp[ threadIdx.x * 3] ];
    tmp[ 2 * threadIdx.x + 1 ] = c[ 2 * ghattmp[ threadIdx.x * 3] + 1 ];

    __syncthreads();

    pwcc[ 2 * index + threadIdx.x ] = scale * tmp[ threadIdx.x ];
    pwcc[ 2 * index + NTHREADS + threadIdx.x ] = scale * 
           tmp[ threadIdx.x + NTHREADS ];

    index = ( igpt + 1 ) * NTHREADS * NBLOCKS + blockIdx.x * NTHREADS;
  }

  nleftover = ngpts - index;
  if (nleftover >= NTHREADS) {
    nleftover = NTHREADS;
  }
  else if (nleftover <= 0 ) {
    nleftover = 0;
  }

  if ( threadIdx.x < nleftover ) {
    ghattmp[threadIdx.x] = ghat[ index * 3 + threadIdx.x ] ;
    ghattmp[threadIdx.x+nleftover] = ghat[ index * 3 + threadIdx.x + nleftover ]; 
    ghattmp[threadIdx.x+2*nleftover] = ghat[ index * 3 + threadIdx.x + 
                                            2 * nleftover ];
  }

  __syncthreads();

  if ( threadIdx.x < nleftover ) {
    ghattmp[threadIdx.x * 3] = mapl[ ghattmp[ threadIdx.x * 3 ] - l1 ];
    ghattmp[threadIdx.x * 3 + 1] = mapm[ ghattmp[ threadIdx.x * 3 + 1 ] - l2 ];
    ghattmp[threadIdx.x * 3 + 2] = mapn[ ghattmp[ threadIdx.x * 3 + 2 ] - l3 ];

    ghattmp[threadIdx.x * 3] += n0 * ( n1 * ghattmp[threadIdx.x * 3 + 2] + 
                                ghattmp[threadIdx.x * 3 + 1]);

    tmp[ 2 * threadIdx.x ] = c[ 2 * ghattmp[ 3 * threadIdx.x ] ];
    tmp[ 2 * threadIdx.x + 1 ] = c[ 2 * ghattmp[ 3 * threadIdx.x ] + 1 ];
  }

  __syncthreads();

  if ( threadIdx.x < nleftover ) {
    pwcc[ 2 * index + threadIdx.x ] = scale * tmp[ threadIdx.x ];
    pwcc[ 2 * index + nleftover + threadIdx.x ] = scale * 
           tmp[ threadIdx.x + nleftover ];
  }
}

__global__ void pw_scatter_cu(float* pwcc, float* c, float scale, int* mapl, int* mapm, int* mapn, int ngpts, int* nn, int* ghat, int l1, int l2, int l3) {

  int igpt, n0, n1, limit, index;
  __shared__ int ghattmp[NTHREADS*3];
  __shared__ float tmp[NTHREADS*2];
  int nleftover;

  n0=nn[0];
  n1=nn[1];

  limit=ngpts/(NTHREADS*NBLOCKS);

  index = blockIdx.x * NTHREADS;
  for (igpt=0; igpt<limit; igpt++) {
    ghattmp[threadIdx.x] = ghat[ index * 3 + threadIdx.x ] ;
    ghattmp[threadIdx.x+NTHREADS] = ghat[ index * 3 + threadIdx.x + NTHREADS ]; 
    ghattmp[threadIdx.x+2*NTHREADS] = ghat[ index * 3 + threadIdx.x + 
                                            2 * NTHREADS ];
    tmp[ threadIdx.x ] = pwcc[ 2 * index + threadIdx.x ];
    tmp[ threadIdx.x + NTHREADS ] = pwcc[ 2 * index + NTHREADS + threadIdx.x ];

    __syncthreads();

    ghattmp[threadIdx.x * 3] = mapl[ ghattmp[ threadIdx.x * 3 ] - l1 ];
    ghattmp[threadIdx.x * 3 + 1] = mapm[ ghattmp[ threadIdx.x * 3 + 1 ] - l2 ];
    ghattmp[threadIdx.x * 3 + 2] = mapn[ ghattmp[ threadIdx.x * 3 + 2 ] - l3 ];

    ghattmp[threadIdx.x * 3] += n0 * ( n1 * ghattmp[threadIdx.x * 3 + 2] + 
                                ghattmp[threadIdx.x * 3 + 1]);

    __syncthreads();

    c[ 2 * ghattmp[ threadIdx.x * 3 ] ] = scale * tmp[ 2 * threadIdx.x ];
    c[ 2 * ghattmp[ threadIdx.x * 3 ] + 1 ] = scale * tmp[ 2 * threadIdx.x + 1 ];

    index = ( igpt + 1 ) * NTHREADS * NBLOCKS + blockIdx.x * NTHREADS;
  }

  nleftover = ngpts - index;
  if (nleftover >= NTHREADS) {
    nleftover = NTHREADS;
  }
  else if (nleftover <= 0 ) {
    nleftover = 0;
  }

  if ( threadIdx.x < nleftover ) {
    ghattmp[threadIdx.x] = ghat[ index * 3 + threadIdx.x ] ;
    ghattmp[threadIdx.x+nleftover] = ghat[ index * 3 + threadIdx.x + nleftover ]; 
    ghattmp[threadIdx.x+2*nleftover] = ghat[ index * 3 + threadIdx.x + 
                                            2 * nleftover ];

    tmp[ threadIdx.x ] = pwcc[ 2 * index + threadIdx.x ];
    tmp[ threadIdx.x + nleftover ] = pwcc[ 2 * index + nleftover + 
                                           threadIdx.x ];
  }

  __syncthreads();

  if ( threadIdx.x < nleftover ) {
    ghattmp[threadIdx.x * 3] = mapl[ ghattmp[ threadIdx.x * 3 ] - l1 ];
    ghattmp[threadIdx.x * 3 + 1] = mapm[ ghattmp[ threadIdx.x * 3 + 1 ] - l2 ];
    ghattmp[threadIdx.x * 3 + 2] = mapn[ ghattmp[ threadIdx.x * 3 + 2 ] - l3 ];

    ghattmp[threadIdx.x * 3] += n0 * ( n1 * ghattmp[threadIdx.x * 3 + 2] + 
                                ghattmp[threadIdx.x * 3 + 1]);
  }

  __syncthreads();

  if ( threadIdx.x < nleftover ) {
    c[ 2 * ghattmp[ 3 * threadIdx.x ] ] = tmp[ 2 * threadIdx.x ];
    c[ 2 * ghattmp[ 3 * threadIdx.x ] + 1 ] = tmp[ 2 * threadIdx.x + 1 ];
  }
}

extern "C" void pw_fft_wrap_fg_cu_(int* fsign, cufftComplex* zin, cufftComplex* zout, float* scale, int* n,  int* mapl, int* mapm, int* mapn, int* ngpts, int* ghat, int* l1, int* l2, int* l3) {
  float *ptr_in, *ptr_out;
  cufftComplex *data, *pwcc;
  int *mapl_dev, *mapm_dev, *mapn_dev, *ghat_dev, *n_dev;
  int lmem;

  lmem = 2 * n[0] * n[1] * n[2];
  cuda_device_mem_alloc_cu_(&ptr_in, lmem);
  cuda_device_mem_alloc_cu_(&ptr_out, lmem);
  cuda_device_mem_alloc_cu_(&ghat_dev, 3*(*ngpts));
  cuda_device_mem_alloc_cu_(&mapl_dev, (n[0]));
  cuda_device_mem_alloc_cu_(&mapm_dev, (n[1]));
  cuda_device_mem_alloc_cu_(&mapn_dev, (n[2]));
  cuda_device_mem_alloc_cu_(&n_dev, 3);
  
  data=(cufftComplex*)ptr_in;
  pwcc=(cufftComplex*)ptr_out;

  cudaMemcpy(data, zin, sizeof(cufftComplex)*n[0]*n[1]*n[2], cudaMemcpyHostToDevice);
  cudaMemcpy(mapl_dev, mapl, sizeof(int)*n[0], cudaMemcpyHostToDevice);
  cudaMemcpy(mapm_dev, mapm, sizeof(int)*n[1], cudaMemcpyHostToDevice);
  cudaMemcpy(mapn_dev, mapn, sizeof(int)*n[2], cudaMemcpyHostToDevice);
  cudaMemcpy(n_dev, n, sizeof(int)*3, cudaMemcpyHostToDevice);
  cudaMemcpy(ghat_dev, ghat, sizeof(int)*3*(*ngpts), cudaMemcpyHostToDevice);

  fftcu_run_3d_cu_(n, data, *fsign, 1.0f);

/*  cudaMemcpy(zin, data, sizeof(cufftComplex)*n[0]*n[1]*n[2], cudaMemcpyDeviceToHost);*/
  cudaThreadSynchronize();

  pw_gather_cu<<<NBLOCKS,NTHREADS>>>(ptr_out, ptr_in, *scale, mapl_dev, mapm_dev, mapn_dev, *ngpts, n_dev, ghat_dev, *l1, *l2, *l3);

  cudaThreadSynchronize();

  cudaMemcpy(zout, pwcc, sizeof(cufftComplex)*n[0]*n[1]*n[2], cudaMemcpyDeviceToHost);

  cuda_device_mem_free_cu_(&ptr_in);
  cuda_device_mem_free_cu_(&ptr_out);
  cuda_device_mem_free_cu_(&ghat_dev);
  cuda_device_mem_free_cu_(&mapl_dev);
  cuda_device_mem_free_cu_(&mapm_dev);
  cuda_device_mem_free_cu_(&mapn_dev);
  cuda_device_mem_free_cu_(&n_dev);
  
}

extern "C" void pw_fft_wrap_sf_cu_(int* fsign, cufftComplex* zin, cufftComplex* zout, float* scale, int* n,  int* mapl, int* mapm, int* mapn, int* ngpts, int* ghat, int* l1, int* l2, int* l3) {
  float *ptr_in, *ptr_out;
  cufftComplex *data, *pwcc;
  int *mapl_dev, *mapm_dev, *mapn_dev, *ghat_dev, *n_dev;
  int lmem;

  lmem = 2 * n[0] * n[1] * n[2];
  cuda_device_mem_alloc_cu_(&ptr_in, lmem);
  cuda_device_mem_alloc_cu_(&ptr_out, lmem);
  cuda_device_mem_alloc_cu_(&ghat_dev, 3*(*ngpts));
  cuda_device_mem_alloc_cu_(&mapl_dev, (n[0]));
  cuda_device_mem_alloc_cu_(&mapm_dev, (n[1]));
  cuda_device_mem_alloc_cu_(&mapn_dev, (n[2]));
  cuda_device_mem_alloc_cu_(&n_dev, 3);
  
  pwcc=(cufftComplex*)ptr_in;
  data=(cufftComplex*)ptr_out;

  cudaMemcpy(pwcc, zin, sizeof(cufftComplex)*n[0]*n[1]*n[2], cudaMemcpyHostToDevice);
  cudaMemcpy(mapl_dev, mapl, sizeof(int)*n[0], cudaMemcpyHostToDevice);
  cudaMemcpy(mapm_dev, mapm, sizeof(int)*n[1], cudaMemcpyHostToDevice);
  cudaMemcpy(mapn_dev, mapn, sizeof(int)*n[2], cudaMemcpyHostToDevice);
  cudaMemcpy(n_dev, n, sizeof(int)*3, cudaMemcpyHostToDevice);
  cudaMemcpy(ghat_dev, ghat, sizeof(int)*3*(*ngpts), cudaMemcpyHostToDevice);

  pw_scatter_cu<<<NBLOCKS,NTHREADS>>>(ptr_in, ptr_out, *scale, mapl_dev, mapm_dev, mapn_dev, *ngpts, n_dev, ghat_dev, *l1, *l2, *l3);

  cudaThreadSynchronize();

  fftcu_run_3d_cu_(n, data, *fsign, 1.0f);

  cudaThreadSynchronize();

  cudaMemcpy(zout, data, sizeof(cufftComplex)*n[0]*n[1]*n[2], cudaMemcpyDeviceToHost);

  cuda_device_mem_free_cu_(&ptr_in);
  cuda_device_mem_free_cu_(&ptr_out);
  cuda_device_mem_free_cu_(&ghat_dev);
  cuda_device_mem_free_cu_(&mapl_dev);
  cuda_device_mem_free_cu_(&mapm_dev);
  cuda_device_mem_free_cu_(&mapn_dev);
  cuda_device_mem_free_cu_(&n_dev);
  
}

#endif
