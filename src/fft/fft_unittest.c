/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../mpiwrap/cp_mpi.h"
#include "../offload/offload_library.h"
#include "fft_api_F_test.h"
#include "fft_lib.h"
#include "fft_lib_test.h"
#include "fft_timer.h"

int run_unittests(const bool debug, const int backend, const int planning_mode,
                  const bool use_guru, const double threshold) {
  int errors = 0;
  fft_finalize_timer();
  fft_finalize_lib(NULL);

  fflush(stdout);
  cp_mpi_barrier(cp_mpi_get_comm_world());

  fft_init_timer(debug);
  fft_init_lib(backend, planning_mode, use_guru, NULL);

  errors += fft_test_local();
  fft_print_timing_report(threshold);

  fft_finalize_lib(NULL);
  fft_finalize_timer();

  return errors;
}

int main(int argc, char *argv[]) {
  cp_mpi_init(&argc, &argv);

  offload_set_chosen_device(0);

  const bool debug = false;
  const double threshold = 1e-4;
  const int backend = fft_lib_default_library();

  if (cp_mpi_comm_rank(cp_mpi_get_comm_world()) == 0) {
    printf("Number of MPI ranks: %i\n",
           cp_mpi_comm_size(cp_mpi_get_comm_world()));
    int number_of_threads = 1;
#pragma omp parallel default(none) shared(number_of_threads)
    {
#pragma omp single
      {
        number_of_threads = omp_get_num_threads();
      }
    }
    printf("Number of OpenMP threads: %i\n", number_of_threads);
    fflush(stdout);
  }
  cp_mpi_barrier(cp_mpi_get_comm_world());

  // First run tests with ESTIMATE mode, which is faster and should be
  // sufficient to catch most errors
  int errors = run_unittests(debug, backend, FFT_ESTIMATE, true, threshold);

  if (fft_lib_has_guru_interface()) {

    errors += run_unittests(debug, backend, FFT_ESTIMATE, false, threshold);
  }

  // Rerun the tests with MEASURE mode, which is slower but more accurate and
  // should catch any remaining errors
  errors += run_unittests(debug, backend, FFT_MEASURE, true, threshold);

  if (fft_lib_has_guru_interface()) {

    errors += run_unittests(debug, backend, FFT_MEASURE, false, threshold);
  }

  fft_finalize_lib(NULL);
  fft_finalize_timer();

  errors += fft_test_api_local();

  if (cp_mpi_comm_rank(cp_mpi_get_comm_world()) == 0) {
    if (errors == 0) {
      printf("\nAll tests have passed :-)\n");
    } else {
      printf("\nFound %i errors :-(\n", errors);
    }
    fflush(stdout);
  }
  cp_mpi_barrier(cp_mpi_get_comm_world());

  cp_mpi_finalize();

  return errors;
}

// EOF
