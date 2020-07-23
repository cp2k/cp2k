/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2020  CP2K developers group                         *
 *****************************************************************************/

#include <assert.h>
#include <stdio.h>

#include "grid_collocate_replay.h"

int main(int argc, char *argv[]){
    if (argc != 3) {
        printf("Usage: grid_base_ref_miniapp.x <cycles> <task-file>\n");
        return 1;
    }

    int cycles;
    GRID_COLLOCATE_ASSERT(sscanf(argv[1], "%i", &cycles) == 1);
    GRID_COLLOCATE_ASSERT(cycles > 0);

    const double max_diff = grid_collocate_replay(argv[2], cycles);
    GRID_COLLOCATE_ASSERT(max_diff < 1e-14 * cycles);
    return 0;
}

//EOF
