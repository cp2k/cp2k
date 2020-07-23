/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2020  CP2K developers group                         *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "grid_collocate_replay.h"

int main(int argc, char *argv[]){
    if (argc != 3) {
        printf("Usage: grid_base_ref_miniapp.x <cycles> <task-file>\n");
        return 1;
    }

    int cycles;
    if (sscanf(argv[1], "%i", &cycles) != 1) {
        fprintf(stderr, "Error: Could not parse cycles.\n");
        return 1;
    }
    if (cycles <= 0) {
        fprintf(stderr, "Error: Cycles have to be greater than zero.\n");
        return 1;
    }

    const double max_diff = grid_collocate_replay(argv[2], cycles);
    if (max_diff > 1e-14 * cycles) {
        fprintf(stderr, "Error: Maximal difference is too large.\n");
        return 2;
    }

    return 0;
}

//EOF
