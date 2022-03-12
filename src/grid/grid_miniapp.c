/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../offload/offload_library.h"
#include "common/grid_library.h"
#include "grid_replay.h"

void mpi_sum_func(long *number, int mpi_comm) {
  *number += 0; // Nothing todo without MPI, pretend arguments are used anyways.
  mpi_comm += 0;
}

void print_func(char *message, int output_unit) {
  output_unit += 0; // Pretend argument is used.
  printf("%s", message);
}

/*******************************************************************************
 * \brief Stand-alone miniapp for running .task files.
 * \author Ole Schuett
 ******************************************************************************/
int main(int argc, char *argv[]) {
  // Parsing of optional args.
  int iarg = 1;
  int nrequired_args = 2;

  bool collocate = true;
  if (iarg < argc && strcmp(argv[iarg], "--integrate") == 0) {
    iarg++;
    collocate = false;
  }

  bool batch = false;
  if (iarg < argc && strcmp(argv[iarg], "--batch") == 0) {
    iarg++;
    batch = true;
    nrequired_args++;
  }

  // All optional args have been parsed.
  if (argc - iarg != nrequired_args) {
    fprintf(stderr, "Usage: grid_miniapp.x [--integrate] [--batch "
                    "<cycles-per-block>] <cycles> <task-file>\n");
    return 1;
  }

  int cycles_per_block = 1;
  if (batch && sscanf(argv[iarg++], "%i", &cycles_per_block) != 1) {
    fprintf(stderr, "Error: Could not parse cycles per block.\n");
    return 1;
  }

  int cycles;
  if (sscanf(argv[iarg++], "%i", &cycles) != 1) {
    fprintf(stderr, "Error: Could not parse cycles.\n");
    return 1;
  }

  offload_set_chosen_device(0);
  grid_library_init();

  const double max_diff =
      grid_replay(argv[iarg++], cycles, collocate, batch, cycles_per_block);

  grid_library_print_stats(&mpi_sum_func, 0, &print_func, 0);
  grid_library_finalize();

  if (max_diff > 1e-12 * cycles) {
    fprintf(stderr, "Error: Maximal difference is too large.\n");
    return 2;
  }

  return 0;
}

// EOF
