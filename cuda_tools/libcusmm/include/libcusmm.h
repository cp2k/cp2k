/******************************************************************************
*  CP2K: A general program to perform molecular dynamics simulations
*  Copyright (C) 2000 - 2013 the CP2K developers group
*****************************************************************************/

int libcusmm_process_d(int *param_stack, int stack_size,
    cudaStream_t stream, int m, int n, int k,
    double * a_data, double * b_data, double * c_data);

