/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015  CP2K developers group                         *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include "libcusmm_benchmark.h"
#include "libcusmm.h"

#if defined(__parallel)
#include <mpi.h>
#endif

/****************************************************************************\
 \brief Checks correctness of every libcusmm kernel and measures its performance.
\****************************************************************************/

int main(int argc, char** argv){
    int rank=0;

#if defined(__parallel)
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank==0) printf("# Using only MPI rank 0 for testing.\n");
#endif

    int errors = 0;
    if(rank==0){
        KernelLauncher launcher = libcusmm_process_d;
        char buffer[1000];
        char * kernel_descr[1] = {buffer};

        const int *blocksizes;
        int n_blocksizes;
        libcusmm_list_blocksizes_d(&blocksizes, &n_blocksizes);
        printf("# Libcusmm has %d blocksizes compiled in...\n", n_blocksizes);
        for(int i=0; i<n_blocksizes; i++){
            int m = blocksizes[3*i + 0];
            int n = blocksizes[3*i + 1];
            int k = blocksizes[3*i + 2];
            sprintf(buffer, "%d x %d x %d", m, n, k);
            errors += libcusmm_benchmark(m, n, k, 1, &launcher, kernel_descr, false);
        }
        printf("# Tests finished with %d errors.\n", errors);
    }

#if defined(__parallel)
    MPI_Bcast(&errors, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Finalize();
#endif

    return(errors);
}

//EOF
