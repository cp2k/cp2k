#if defined ( __FFTCU ) || defined ( __CUDA )

void fftcu_plan3d(cufftHandle& plan, int* n, int& ioverflow);

extern "C" void fftcu_run_3d_cu_(int* n, cufftComplex* data, int fsign, float scale);

extern "C" void fftcu3d_cu_ (int *ifft_in_place, int *fsign, float *scale, int *n, float *zin, float *zout);

extern "C" void fftcu1dm_cu_ (int *fsign, int itrans, int *n, int *m, float *zin, float *zout, float *scale);

#endif
