/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2013 the CP2K developers group                      *
 *****************************************************************************/

int launch_cusmm_kernel_fallback(int *param_stack, int stack_size, cudaStream_t stream,
                      int m_max, int n_max, int k_max,
                      double *a_data, double *b_data, double *c_data);

