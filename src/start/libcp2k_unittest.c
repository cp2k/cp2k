/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include "libcp2k.h"
#include "mpiwrap/cp_mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*******************************************************************************
 * \brief Unit test of the C-interface provided via libcp2k.h
 * \author Ole Schuett
 ******************************************************************************/
int main() {

  printf("Unit test starts ...\n");

  // test cp2k_get_version()
  printf("Testing cp_c_get_version(): ");
  char version_str[100];
  cp2k_get_version(version_str, 100);
  printf("%s.\n", version_str);

  cp2k_init();
  const int rank = cp_mpi_comm_rank(cp_mpi_get_comm_world());

  // create simple input file
  const char *inp_fn = "H2.inp";
  FILE *f;
  if (rank == 0) {
    f = fopen(inp_fn, "w");
    fprintf(f, "&FORCE_EVAL\n");
    fprintf(f, "  METHOD Quickstep\n");
    fprintf(f, "  STRESS_TENSOR ANALYTICAL\n");
    fprintf(f, "  &DFT\n");
    fprintf(f, "    BASIS_SET_FILE_NAME BASIS_SET\n");
    fprintf(f, "    POTENTIAL_FILE_NAME POTENTIAL\n");
    fprintf(f, "    LSD\n");
    fprintf(f, "    &MGRID\n");
    fprintf(f, "      CUTOFF 140\n");
    fprintf(f, "    &END MGRID\n");
    fprintf(f, "    &QS\n");
    fprintf(f, "      EPS_DEFAULT 1.0E-8\n");
    fprintf(f, "    &END QS\n");
    fprintf(f, "    &SCF\n");
    fprintf(f, "      EPS_DIIS 0.1\n");
    fprintf(f, "      EPS_SCF 1.0E-4\n");
    fprintf(f, "      IGNORE_CONVERGENCE_FAILURE\n");
    fprintf(f, "      MAX_DIIS 4\n");
    fprintf(f, "      MAX_SCF 3\n");
    fprintf(f, "      SCF_GUESS atomic\n");
    fprintf(f, "    &END SCF\n");
    fprintf(f, "    &XC\n");
    fprintf(f, "      &XC_FUNCTIONAL Pade\n");
    fprintf(f, "      &END XC_FUNCTIONAL\n");
    fprintf(f, "    &END XC\n");
    fprintf(f, "  &END DFT\n");
    fprintf(f, "  &SUBSYS\n");
    fprintf(f, "    &CELL\n");
    fprintf(f, "      ABC 8.0 4.0 4.0\n");
    fprintf(f, "    &END CELL\n");
    fprintf(f, "    &COORD\n");
    fprintf(f, "    H     0.000000  0.000000  0.000000\n");
    fprintf(f, "    H     1.000000  0.000000  0.000000\n");
    fprintf(f, "    &END COORD\n");
    fprintf(f, "    &KIND H\n");
    fprintf(f, "      BASIS_SET DZV-GTH-PADE\n");
    fprintf(f, "      POTENTIAL GTH-PADE-q1\n");
    fprintf(f, "    &END KIND\n");
    fprintf(f, "  &END SUBSYS\n");
    fprintf(f, "&END FORCE_EVAL\n");
    fprintf(f, "&GLOBAL\n");
    fprintf(f, "  PRINT_LEVEL SILENT\n");
    fprintf(f, "  PROJECT libcp2k_unittest_H2\n");
    fprintf(f, "&END GLOBAL\n");
    fclose(f);
  }
  cp_mpi_barrier(cp_mpi_get_comm_world());

  // use input file to create a force environment
  force_env_t force_env;
  cp2k_create_force_env(&force_env, inp_fn, "__STD_OUT__");
  int scf_status = 99;
  cp2k_get_scf_convergence(force_env, &scf_status);
  if (scf_status != -1) {
    printf("SCF status must be unavailable before calculation\n");
    return (-1);
  }
  cp2k_calc_energy_force(force_env);
  cp2k_get_scf_convergence(force_env, &scf_status);
  if (scf_status != 0) {
    printf("The deliberately truncated SCF must report non-convergence\n");
    return (-1);
  }

  // Stress is a column-major, pressure-positive potential tensor.
  double stress[9];
  int stress_available = 0;
  cp2k_get_stress_tensor(force_env, stress, &stress_available);
  if (!stress_available) {
    printf("Missing analytical stress\n");
    return (-1);
  }
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (!isfinite(stress[3 * j + i]) ||
          fabs(stress[3 * j + i] - stress[3 * i + j]) > 1e-10) {
        printf("Invalid stress tensor\n");
        return (-1);
      }
    }
  }

  // check energy
  double energy;
  cp2k_get_potential_energy(force_env, &energy);
  printf("\n ENERGY: %.12f\n", energy);
  if (fabs(-1.118912797546392 - energy) / fabs(energy) > 1e-13) {
    printf("Wrong energy\n");
    return (-1);
  }

  // Geometry and velocity updates invalidate the last SCF status.
  double positions[6], cell[9], velocities[6] = {0};
  cp2k_get_positions(force_env, positions, 6);
  cp2k_set_positions(force_env, positions, 6);
  cp2k_get_scf_convergence(force_env, &scf_status);
  if (scf_status != -1)
    return (-1);
  cp2k_calc_energy(force_env);
  cp2k_get_cell(force_env, cell);
  cp2k_set_cell(force_env, cell);
  cp2k_get_scf_convergence(force_env, &scf_status);
  if (scf_status != -1)
    return (-1);
  cp2k_calc_energy(force_env);
  cp2k_set_velocities(force_env, velocities, 6);
  cp2k_get_scf_convergence(force_env, &scf_status);
  if (scf_status != -1)
    return (-1);
  cp2k_destroy_force_env(force_env);
  // A library caller has not already opened output in the Fortran runtime.
  // run_input must create, close, and subsequently append to the named file.
  const char *run_out = "libcp2k_unittest_run.out";
  cp2k_run_input_comm(inp_fn, run_out,
                      cp_mpi_comm_c2f(cp_mpi_get_comm_world()));
  cp_mpi_barrier(cp_mpi_get_comm_world());
  long first_size = 0;
  if (rank == 0) {
    f = fopen(run_out, "r");
    if (f == NULL) {
      printf("run_input did not create its output file\n");
      return (-1);
    }
    fseek(f, 0, SEEK_END);
    first_size = ftell(f);
    fclose(f);
  }
  cp2k_run_input(inp_fn, run_out);
  cp_mpi_barrier(cp_mpi_get_comm_world());
  if (rank == 0) {
    f = fopen(run_out, "r");
    if (f == NULL) {
      printf("run_input output file disappeared\n");
      return (-1);
    }
    fseek(f, 0, SEEK_END);
    const long second_size = ftell(f);
    fclose(f);
    if (first_size <= 0 || second_size <= first_size) {
      printf("run_input did not append output\n");
      return (-1);
    }
  }

  // run_input must leave the surrounding library runtime usable.
  cp2k_create_force_env(&force_env, inp_fn, "__STD_OUT__");
  cp2k_calc_energy_force(force_env);
  cp2k_get_potential_energy(force_env, &energy);
  if (fabs(-1.118912797546392 - energy) / fabs(energy) > 1e-13) {
    printf("Wrong energy after run_input\n");
    return (-1);
  }
  cp2k_destroy_force_env(force_env);

  // clean up
  cp2k_finalize();
  if (rank == 0) {
    remove(inp_fn);
    remove(run_out);
  }

  printf("Unit test finished, found no errors\n");
  return (0);
}

// EOF
