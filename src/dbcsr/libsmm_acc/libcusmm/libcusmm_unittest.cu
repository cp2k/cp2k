/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015  CP2K developers group                         *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include "libcusmm_benchmark.h"
#include "libcusmm.h"

/****************************************************************************\
 \brief Checks correctness of every libcusmm kernel and measures its performance.
\****************************************************************************/

int main(int argc, char** argv){

    KernelLauncher launcher = libcusmm_process_d;
    char buffer[1000];
    char * kernel_descr[1] = {buffer};

    const int *blocksizes;
    int n_blocksizes;
    libcusmm_list_blocksizes_d(&blocksizes, &n_blocksizes);
    printf("# Libcusmm has %d blocksizes compiled in...\n", n_blocksizes);

    int errors = 0;
    for(int i=0; i<n_blocksizes; i++){
        int m = blocksizes[3*i + 0];
        int n = blocksizes[3*i + 1];
        int k = blocksizes[3*i + 2];
        sprintf(buffer, "%d x %d x %d", m, n, k);
        errors += libcusmm_benchmark(m, n, k, 1, &launcher, kernel_descr, false);
    }

    return(errors);
}

//EOF
