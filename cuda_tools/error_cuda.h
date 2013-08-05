#ifndef ERROR_CUDA_H
#define ERROR_CUDA_H


extern void cuda_error_check2 (cudaError_t cudaError, int line);

#if defined ( __PW_CUDA )
#include <cufft.h>
extern void cufft_error_check2 (cufftResult_t cufftError, int line);
#endif

#endif
