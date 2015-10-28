/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015  CP2K developers group                         *
 *****************************************************************************/

typedef int (*KernelLauncher)(int *param_stack, int stack_size, cudaStream_t stream,
                              int m_max, int n_max, int k_max,
                              double *a_data, double *b_data, double *c_data);

typedef struct {
    // host-buffers
    double* mat_a;
    double* mat_b;
    double* mat_c;
    int*    stack;
    // device-buffers
    double* d_mat_a;
    double* d_mat_b;
    double* d_mat_c;
    int*    d_stack;
    int max_m, max_n, max_k, max_stack;
    cudaEvent_t t_start, t_stop;
} libcusmm_benchmark_t;

void libcusmm_benchmark_init(libcusmm_benchmark_t** handle, int max_m, int max_n, int max_k);
void libcusmm_benchmark_finalize(libcusmm_benchmark_t*);

int libcusmm_benchmark(libcusmm_benchmark_t* handle,
              int mat_m, int mat_n, int mat_k, int nkernel,
              KernelLauncher* launchers, char** kernel_descr, bool tune_mode);

//EOF
