/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015 the CP2K developers group                      *
 *****************************************************************************/

typedef int (*KernelLauncher)(int *param_stack, int stack_size, cudaStream_t stream,
                              int m_max, int n_max, int k_max,
                              double *a_data, double *b_data, double *c_data);

int libcusmm_benchmark(int mat_m, int mat_n, int mat_k, int nkernel,
              KernelLauncher* launchers, char** kernel_descr, bool tune_mode);

//EOF
