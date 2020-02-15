/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2020  CP2K developers group                         *
 *****************************************************************************/

#include <stdio.h>
#include <assert.h>
#include <string.h>
//
#include "grid_collocate_replay.h"

static int run_test(const char cp2k_root_dir[], const char task_file[]) {
    char filename[1024] = "";

    assert(strlen(cp2k_root_dir) < 512);
    assert(strcpy(filename, cp2k_root_dir) != NULL);
    if (filename[strlen(filename) - 1] != '/') {
        assert(strcat(filename, "/") != NULL);
    }

    assert(strcat(filename, "src/grid/sample_tasks/") != NULL);
    assert(strcat(filename, task_file) != NULL);

    const double max_diff = grid_collocate_replay(filename, 1);
    if (max_diff > 1e-14) {
        printf("Max diff too high, test failed.\n");
        return 1;
    } else {
        printf("Max diff looks good, test passed.\n\n");
        return 0;
    }
}

int main(int argc, char *argv[]){
    if (argc != 2) {
        printf("Usage: grid_base_ref_unittest.x <cp2k-root-dir>\n");
        return 1;
    }

    int errors = 0;
    errors += run_test(argv[1], "collocate_ortho_density_l0000.task");
    errors += run_test(argv[1], "collocate_ortho_density_l2200.task");
    errors += run_test(argv[1], "collocate_ortho_density_l3300.task");
    errors += run_test(argv[1], "collocate_ortho_density_l3333.task");
    errors += run_test(argv[1], "collocate_ortho_tau.task");
    errors += run_test(argv[1], "collocate_general_density.task");
    errors += run_test(argv[1], "collocate_general_tau.task");
    errors += run_test(argv[1], "collocate_general_subpatch0.task");
    errors += run_test(argv[1], "collocate_general_subpatch16.task");
    errors += run_test(argv[1], "collocate_ortho_non_periodic.task");

    if (errors == 0) {
        printf("All tests have passed :-)\n");
    } else {
        printf("Found %i errors :-(\n", errors);
    }
    return errors;
}

//EOF
