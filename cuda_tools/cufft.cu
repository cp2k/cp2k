#if defined ( __FFTCU ) || defined ( __CUDA )

#include <cuda_runtime.h>
#include <cufft.h>
#include <stdio.h>

extern "C" void fftcu3dc_ (int *ifft_in_place, int *fsign, float *scale, int *n, float *zin, float *zout) {

  cufftComplex *data;
  cufftHandle plan;
  int i;
  
  cudaMalloc((void**)&data, sizeof(cufftComplex)*n[0]*n[1]*n[2]);

  cudaMemcpy(data, zin, sizeof(cufftComplex)*n[0]*n[1]*n[2], cudaMemcpyHostToDevice);

  cufftPlan3d(&plan, n[2], n[1], n[0], CUFFT_C2C);

  if ( *fsign == 1 ) {
    // this code does not yet use scale!
    cufftExecC2C(plan, data, data, CUFFT_INVERSE);
  }
  else {
    // this code does not yet use scale!
    cufftExecC2C(plan, data, data, CUFFT_FORWARD);
  }

  if ( *ifft_in_place == 1 ) {
    cudaMemcpy(zin, data, sizeof(cufftComplex)*n[0]*n[1]*n[2], cudaMemcpyDeviceToHost);
    if (*scale /= 1.0) {
      for (i=0; i<2*n[0]*n[1]*n[2]; i++) {
        zin[i]=zin[i]*(*scale);
      }
    }
  }
  else {
    cudaMemcpy(zout, data, sizeof(cufftComplex)*n[0]*n[1]*n[2], cudaMemcpyDeviceToHost);
    if (*scale /= 1.0) {
      for (i=0; i<2*n[0]*n[1]*n[2]; i++) {
        zin[i]=zin[i]*(*scale);
      }
    }
  }

  cufftDestroy(plan);
  cudaFree(data);
}

#endif