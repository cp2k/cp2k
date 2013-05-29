#ifndef ERROR_CUDA_H
#define ERROR_CUDA_H
#if defined ( __PW_CUDA )
#include <cufft.h>

//extern void cuda_error_check (cudaError_t cudaError);

extern void cuda_error_check2 (cudaError_t cudaError, int line);

//extern void cufft_error_check (cufftResult_t cufftError);

extern void cufft_error_check2 (cufftResult_t cufftError, int line);

#endif
#endif
