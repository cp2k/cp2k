/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015  CP2K developers group                         *
 *****************************************************************************/

#ifndef LIBCLSMM_H
#define LIBCLSMM_H

#if defined (__ACC) && defined (__OPENCL)

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

int libclsmm_process_d (void *param_stack, int stack_size,
    void *stream, int m, int n, int k,
    void *a_data, void *b_data, void *c_data);

int libclsmm_transpose_d (void *trs_stack, int offset, int nblks, void *buffer,
                         int m, int n, void *stream);

void libclsmm_list_blocksizes_d (const int **list, int *length);

// keep compiled kernels
static cl_kernel multiply_kernel = NULL;
static cl_kernel transpose_kernel = NULL;
#endif

#endif
//EOF
