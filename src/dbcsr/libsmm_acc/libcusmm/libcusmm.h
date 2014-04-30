/******************************************************************************
*  CP2K: A general program to perform molecular dynamics simulations
*  Copyright (C) 2000 - 2013 the CP2K developers group
*****************************************************************************/

int libcusmm_process_d(int *param_stack, int stack_size,
    cudaStream_t stream, int m, int n, int k,
    double * a_data, double * b_data, double * c_data);

int libcusmm_transpose_d(int *trs_stack, int offset, int nblks, double *buffer,
                         int m, int n, cudaStream_t * stream);

void libcusmm_list_blocksizes_d(const int **list, int *length);

//EOF
