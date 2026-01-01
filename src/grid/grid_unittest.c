/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#include "../offload/offload_library.h"
#include "../offload/offload_mempool.h"
#include "common/grid_library.h"
#include "grid_replay.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Only used to call MPI_Init and MPI_Finalize to avoid spurious MPI error.
#if defined(__parallel)
#include <mpi.h>
#endif

/*******************************************************************************
 * \brief Wrapper for printf, passed to grid_library_print_stats.
 * \author Ole Schuett
 ******************************************************************************/
static void print_func(const char *msg, int msglen, int output_unit) {
  (void)msglen;           // mark used
  if (output_unit == 0) { // i.e. my_rank == 0
    printf("%s", msg);
  }
}

/*******************************************************************************
 * \brief Unit test for the grid code.
 * \author Ole Schuett
 ******************************************************************************/
static int run_test(const char cp2k_root_dir[], const char task_file[]) {
  if (strlen(cp2k_root_dir) > 512) {
    fprintf(stderr, "Error: cp2k_root_dir too long.\n");
    abort();
  }

  char filename[1024];
  strcpy(filename, cp2k_root_dir);
  if (filename[strlen(filename) - 1] != '/') {
    strcat(filename, "/");
  }

  strcat(filename, "src/grid/sample_tasks/");
  strcat(filename, task_file);

  const double tolerance = 1e-12;
  int errors = 0;
  for (int icol = 0; icol < 2; icol++) {
    for (int ibatch = 0; ibatch < 2; ibatch++) {
      const bool success =
          grid_replay(filename, 1, icol == 1, ibatch == 1, 1, tolerance);
      if (!success) {
        printf("Max diff too high, test failed.\n\n");
        errors++;
      }
    }
  }
  return errors;
}

int main(int argc, char *argv[]) {
#if defined(__parallel)
  MPI_Init(&argc, &argv);
#endif

  if (argc != 2) {
    printf("Usage: grid_unittest.x <cp2k-root-dir>\n");
    return 1;
  }

  offload_set_chosen_device(0);
  grid_library_init();

  int errors = 0;
  errors += run_test(argv[1], "ortho_density_l0000.task");
  errors += run_test(argv[1], "ortho_density_l0122.task");
  errors += run_test(argv[1], "ortho_density_l2200.task");
  errors += run_test(argv[1], "ortho_density_l3300.task");
  errors += run_test(argv[1], "ortho_density_l3333.task");
  errors += run_test(argv[1], "ortho_density_l0505.task");
  errors += run_test(argv[1], "ortho_non_periodic.task");
  errors += run_test(argv[1], "ortho_tau.task");
  errors += run_test(argv[1], "general_density.task");
  errors += run_test(argv[1], "general_tau.task");
  errors += run_test(argv[1], "general_subpatch0.task");
  errors += run_test(argv[1], "general_subpatch16.task");
  errors += run_test(argv[1], "general_overflow.task");

  if (errors == 0) {
    grid_library_print_stats(0 /*fortran_comm*/, &print_func, 0 /*rank*/);
    offload_mempool_stats_print(0 /*fortran_comm*/, &print_func, 0 /*rank*/);
    grid_library_finalize();
    printf("\nAll tests have passed :-)\n");
  } else {
    grid_library_finalize();
    printf("\nFound %i errors :-(\n", errors);
  }

#if defined(__parallel)
  MPI_Finalize();
#endif

  return errors;
}

// EOF
