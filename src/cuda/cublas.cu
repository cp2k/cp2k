#if defined ( __PW_CUDA ) || defined ( __CUBLASDP )

#include <cuda_runtime.h>
#include <cublas.h>
#include <stdio.h>

#include "memory_cuda.h"

extern "C" void cublasinit_cu_() {
  cublasInit();
}

extern "C" void cublasshutdown_cu_() {
  cublasShutdown();
}

extern "C" void gpu_d_gemm_(char& transa, char& transb, int& m, int& n, int& k, double& alpha, double* a, int& lda, double* b, int& ldb, double& beta, double* c, int& ldc) {

  double *ptra;
  double *ptrb;
  double *ptrc;
  int ka, kb, la, lb, lmem;

  ka=m;
  la=k;
  if (transa=='n' || transa=='N') {
    ka=k;
    la=m;
  }

  kb=k;
  lb=n;
  if (transb=='n' || transb=='N') {
    kb=n;
    lb=k;
  }

  lmem=ka*la;
  cuda_device_mem_alloc_cu_(&ptra, lmem);
  lmem=kb*lb;
  cuda_device_mem_alloc_cu_(&ptrb, lmem);
  lmem=n*m;
  cuda_device_mem_alloc_cu_(&ptrc, lmem);

  cublasSetMatrix(la, ka, 8, a, lda, ptra, la);
  cublasSetMatrix(lb, kb, 8, b, ldb, ptrb, lb);
  cublasSetMatrix(m, n, 8, c, ldc, ptrc, m);

  cublasDgemm(transa, transb, m, n, k, alpha, ptra, la, ptrb, lb, beta, ptrc, m);

  cublasGetMatrix(m, n, 8, ptrc, m, c, ldc);

  cuda_device_mem_free_cu_(&ptra);
  cuda_device_mem_free_cu_(&ptrb);
  cuda_device_mem_free_cu_(&ptrc);
}

extern "C" void gpu_d_symm_(char& side, char& uplo, int& m, int& n, double& alpha, double* a, int& lda, double* b, int& ldb, double& beta, double* c, int& ldc) {

  double *ptra;
  double *ptrb;
  double *ptrc;
  int ka, lmem;

  ka=n;
  if (side=='L' || side=='l') {
    ka=m;
  }


  lmem=m*n;
  cuda_device_mem_alloc_cu_(&ptra, lmem);
  lmem=n*n;
  cuda_device_mem_alloc_cu_(&ptrb, lmem);
  lmem=m*n;
  cuda_device_mem_alloc_cu_(&ptrc, lmem);

  cublasSetMatrix(m, n, 8, a, lda, ptra, m);
  cublasSetMatrix(n, n, 8, b, ldb, ptrb, n);
  cublasSetMatrix(m, n, 8, c, ldc, ptrc, m);

  cublasDsymm(side, uplo, m, n, alpha, ptra, m, ptrb, n, beta, ptrc, m);

  cublasGetMatrix(m, n, 8, ptrc, m, c, ldc);

  cuda_device_mem_free_cu_(&ptra);
  cuda_device_mem_free_cu_(&ptrb);
  cuda_device_mem_free_cu_(&ptrc);
}

#endif
