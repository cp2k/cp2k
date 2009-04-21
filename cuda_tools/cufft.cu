#if defined ( __FFTCU ) || defined ( __CUDA )

#include <cuda_runtime.h>
#include <cufft.h>
#include <cublas.h>
#include <stdio.h>

#include "cuda_memory_cu.h"

static const int max_plans=10;
static int n_plans=0;
static cufftHandle saved_plans[max_plans];
static int iplandims[max_plans][3];

void fftcu_plan3d(cufftHandle& plan, int* n, int& ioverflow) {
  int i;

  ioverflow=0;
  for (i=0; i<n_plans; i++) {
    if ( n[0] == iplandims[i][0] && n[1] == iplandims[i][1] && 
         n[2] == iplandims[i][2] ) {
      plan = saved_plans[i];
      return;
    }
  }
  cufftPlan3d(&plan, n[2], n[1], n[0], CUFFT_C2C);
  if ( n_plans < max_plans ) {
    fflush(stdout);
    saved_plans[n_plans] = plan;
    iplandims[n_plans][0] = n[0];
    iplandims[n_plans][1] = n[1];
    iplandims[n_plans][2] = n[2];
    n_plans++;
    return;
  }
  ioverflow=1;
}

extern "C" void fftcu_run_3d_cu_(int* n, cufftComplex* data, int fsign, float scale) {
  int ioverflow, lmem;
  cufftHandle plan;

  lmem=2*n[0]*n[1]*n[2];

  fftcu_plan3d(plan,n,ioverflow);

  //exit(0);

  if ( fsign < 0.0f  ) {
    cufftExecC2C(plan, data, data, CUFFT_INVERSE);
  }
  else {
    cufftExecC2C(plan, data, data, CUFFT_FORWARD);
  }

  if (scale /= 1.0f) {
    cublasSscal(lmem, scale, (float*)data, 1);
  }

  if (ioverflow) { cufftDestroy(plan); }

}

extern "C" void fftcu3d_cu_ (int *ifft_in_place, int *fsign, float *scale, int *n, float *zin, float *zout) {

  float *ptr;
  cufftComplex *data;
  int lmem;
  
  lmem=2*n[0]*n[1]*n[2];
  cuda_device_mem_alloc_cu_(&ptr, lmem);

  data=(cufftComplex*)ptr;

  cudaMemcpy(data, zin, sizeof(cufftComplex)*n[0]*n[1]*n[2], cudaMemcpyHostToDevice);

  fftcu_run_3d_cu_(n, data, *fsign, *scale);

  if ( *ifft_in_place == 1 ) {
    cudaMemcpy(zin, data, sizeof(cufftComplex)*n[0]*n[1]*n[2], cudaMemcpyDeviceToHost);
  }
  else {
    cudaMemcpy(zout, data, sizeof(cufftComplex)*n[0]*n[1]*n[2], cudaMemcpyDeviceToHost);
  }
  cuda_device_mem_free_cu_(&ptr);
}

extern "C" void fftcu1dm_cu_ (int *fsign, int itrans, int *n, int *m, float *zin, float *zout, float *scale) {

  float *ptr;
  cufftComplex *data;
  cufftHandle plan;
  int lmem;
  
  lmem=2*(*n)*(*m);
  cuda_device_mem_alloc_cu_(&ptr, lmem);

  data=(cufftComplex*)ptr;

  cudaMemcpy(data, zin, sizeof(cufftComplex)*(*n)*(*m), cudaMemcpyHostToDevice);

  cufftPlan1d(&plan, (*n), CUFFT_C2C, (*m));

  if ( *fsign == 1 ) {
    cufftExecC2C(plan, data, data, CUFFT_INVERSE);
  }
  else {
    cufftExecC2C(plan, data, data, CUFFT_FORWARD);
  }

  if (*scale /= 1.0f) {
    cublasSscal(lmem, *scale, ptr, 1);
  }

  cudaMemcpy(zout, data, sizeof(cufftComplex)*(*n)*(*m), cudaMemcpyDeviceToHost);

  cuda_device_mem_free_cu_(&ptr);
  cufftDestroy(plan);
}

#endif