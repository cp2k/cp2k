/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2025 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../offload/offload_library.h"
#include "common/grid_library.h"
#include "grid_replay.h"

// Only used to call MPI_Init and MPI_Finalize to avoid spurious MPI error.
#if defined(__parallel)
#include <mpi.h>
#endif

/*******************************************************************************
 * \brief Standin for mpi_sum, passed to grid_library_print_stats.
 * \author Ole Schuett
 ******************************************************************************/
static void mpi_sum_func(long *number, int mpi_comm) {
  (void)number; // mark used
  (void)mpi_comm;
}

/*******************************************************************************
 * \brief Wrapper for printf, passed to grid_library_print_stats.
 * \author Ole Schuett
 ******************************************************************************/
static void print_func(char *message, int output_unit) {
  (void)output_unit; // mark used
  printf("%s", message);
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

  grid_library_print_stats(&mpi_sum_func, 0, &print_func, 0);
  grid_library_finalize();

  if (errors == 0) {
    printf("\nAll tests have passed :-)\n");
  } else {
    printf("\nFound %i errors :-(\n", errors);
  }

#if defined(__parallel)
  MPI_Finalize();
#endif

  return errors;
}

// EOF
