/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#include "../mpiwrap/cp_mpi.h"
#include "../offload/offload_library.h"
#include "../offload/offload_mempool.h"
#include "common/grid_library.h"
#include "cpu/grid_cpu_collocate.h"
#include "cpu/grid_cpu_integrate.h"
#include "grid_replay.h"
#include "ref/grid_ref_collocate.h"
#include "ref/grid_ref_integrate.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

/*******************************************************************************
 * \brief Regression tests for general-cell bounds on strongly sheared cells.
 ******************************************************************************/
static int run_sheared_case(const double dh[3][3], const double ra[3],
                            const double radius, const double zeta,
                            const bool compare_ref) {
  const double det = dh[0][0] * (dh[1][1] * dh[2][2] - dh[1][2] * dh[2][1]) -
                     dh[0][1] * (dh[1][0] * dh[2][2] - dh[1][2] * dh[2][0]) +
                     dh[0][2] * (dh[1][0] * dh[2][1] - dh[1][1] * dh[2][0]);
  const double inv_det = 1.0 / det;
  const double dh_inv[3][3] = {
      {(dh[1][1] * dh[2][2] - dh[1][2] * dh[2][1]) * inv_det,
       (dh[0][2] * dh[2][1] - dh[0][1] * dh[2][2]) * inv_det,
       (dh[0][1] * dh[1][2] - dh[0][2] * dh[1][1]) * inv_det},
      {(dh[1][2] * dh[2][0] - dh[1][0] * dh[2][2]) * inv_det,
       (dh[0][0] * dh[2][2] - dh[0][2] * dh[2][0]) * inv_det,
       (dh[0][2] * dh[1][0] - dh[0][0] * dh[1][2]) * inv_det},
      {(dh[1][0] * dh[2][1] - dh[1][1] * dh[2][0]) * inv_det,
       (dh[0][1] * dh[2][0] - dh[0][0] * dh[2][1]) * inv_det,
       (dh[0][0] * dh[1][1] - dh[0][1] * dh[1][0]) * inv_det}};
  const double rab[3] = {0.0, 0.0, 0.0};
  const int npts_global[3] = {compare_ref ? 4 : 8000, 4,
                              compare_ref ? 4 : 8000};
  const int npts_local[3] = {compare_ref ? 4 : 1, compare_ref ? 4 : 1,
                             compare_ref ? 4 : 1};
  const int ng = compare_ref ? 64 : 1;
  const int zero3[3] = {0, 0, 0};
  const double pab[1][1] = {{1.0}};
  double grid_ref[64] = {0.0}, grid_cpu[64] = {0.0};
  double hab_ref[1][1] = {{0.0}}, hab_cpu[1][1] = {{0.0}};
  double forces[2][3] = {{0.0}}, virials[2][3][3] = {{{0.0}}};

  grid_cpu_collocate_pgf_product(
      false, 0, GRID_FUNC_AB, 0, 0, 0, 0, zeta, zeta, 1.0, dh, dh_inv, ra, rab,
      npts_global, npts_local, zero3, zero3, radius, 0, 0, 1, 1, pab, grid_cpu);
  if (compare_ref) {
    grid_ref_collocate_pgf_product(false, 0, GRID_FUNC_AB, 0, 0, 0, 0, zeta,
                                   zeta, 1.0, dh, dh_inv, ra, rab, npts_global,
                                   npts_local, zero3, zero3, radius, 0, 0, 1, 1,
                                   pab, grid_ref);
    grid_ref_integrate_pgf_product(
        false, false, 0, 0, 0, 0, 0, zeta, zeta, dh, dh_inv, ra, rab,
        npts_global, npts_local, zero3, zero3, radius, 0, 0, 1, 1, grid_ref,
        hab_ref, pab, forces, virials, NULL, NULL, NULL);
  } else {
    // Only (0, 0, 0) maps into npts_local={1,1,1}; avoids stack VLA in grid_ref
    // when the sheared sphere spans thousands of lattice indices.
    const double r2 = ra[0] * ra[0] + ra[1] * ra[1] + ra[2] * ra[2];
    grid_ref[0] = exp(-2.0 * zeta * r2);
    hab_ref[0][0] = grid_ref[0] * grid_ref[0];
  }
  grid_cpu_integrate_pgf_product(
      false, false, 0, 0, 0, 0, 0, zeta, zeta, dh, dh_inv, ra, rab, npts_global,
      npts_local, zero3, zero3, radius, 0, 0, 1, 1, grid_ref, hab_cpu, pab,
      forces, virials, NULL, NULL, NULL);

  double ref_norm = 0.0, max_diff = fabs(hab_cpu[0][0] - hab_ref[0][0]);
  for (int i = 0; i < ng; i++) {
    ref_norm = fmax(ref_norm, fabs(grid_ref[i]));
    max_diff = fmax(max_diff, fabs(grid_cpu[i] - grid_ref[i]));
  }
  if (ref_norm < 0.5 || max_diff > 1e-12) {
    printf("Sheared-cell test failed: ref_norm = %le, max_diff = %le\n",
           ref_norm, max_diff);
    return 1;
  }
  return 0;
}

static int run_sheared_tests(void) {
  int errors = 0;
  // 1. Strongly sheared cell with grazing point (0, 0, 0) (sin^2 < 1e-6).
  const double dh_extreme[3][3] = {
      {0.7443799999999999, 0.0, 0.0}, {0.0, 0.5, 0.0}, {0.73102, 0.0, 0.00023}};
  const double r_pt = 0.505 * dh_extreme[0][0];
  const double ra_extreme[3] = {0.0, 0.0, -r_pt};
  errors += run_sheared_case(dh_extreme, ra_extreme, r_pt * (1.0 + 5e-10), 2e-5,
                             false);

  // 2. Rotated 3D sheared cell inside the active j-bound regime (sin^2 > 1e-6).
  const double q[3][3] = {{2.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0},
                          {2.0 / 3.0, -1.0 / 3.0, -2.0 / 3.0},
                          {-1.0 / 3.0, 2.0 / 3.0, -2.0 / 3.0}};
  const double c[3][3] = {{0.8, 0.0, 0.0}, {0.24, 0.8, 0.0}, {0.8, 0.8, 0.015}};
  double dh_active[3][3] = {{0.0}};
  for (int i = 0; i < 3; i++) {
    for (int m = 0; m < 3; m++) {
      dh_active[i][m] =
          c[i][0] * q[0][m] + c[i][1] * q[1][m] + c[i][2] * q[2][m];
    }
  }
  const double ra_active[3] = {-0.55 * q[2][0], -0.55 * q[2][1],
                               -0.55 * q[2][2]};
  errors +=
      run_sheared_case(dh_active, ra_active, 0.55 * (1.0 + 1e-9), 5e-4, true);
  return errors;
}

int main(int argc, char *argv[]) {
  cp_mpi_init(&argc, &argv);

  if (argc != 2) {
    printf("Usage: grid_unittest.x <cp2k-root-dir>\n");
    cp_mpi_finalize();
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
  errors += run_sheared_tests();

  if (errors == 0) {
    grid_library_print_stats(0 /*fortran_comm*/, &print_func, 0 /*rank*/);
    offload_mempool_stats_print(0 /*fortran_comm*/, &print_func, 0 /*rank*/);
    grid_library_finalize();
    printf("\nAll tests have passed :-)\n");
  } else {
    grid_library_finalize();
    printf("\nFound %i errors :-(\n", errors);
  }

  cp_mpi_finalize();

  return errors;
}

// EOF
