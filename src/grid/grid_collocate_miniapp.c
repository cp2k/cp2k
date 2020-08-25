/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2020  CP2K developers group                         *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/grid_library.h"
#include "grid_collocate_replay.h"

void mpi_sum_func(long *number) {
  *number += 0; // Nothing todo without MPI, pretend argument is used anyways.
}

void print_func(char *message) { printf("%s", message); }

//******************************************************************************
// \brief Stand-alone miniapp for running .task files.
// \author Ole Schuett
//******************************************************************************
int main(int argc, char *argv[]) {
  // Parsing of optional args.
  int iarg = 1;

  bool batch = false;
  int nrequired_args = 2;
  if (iarg < argc && strcmp(argv[iarg], "--batch") == 0) {
    iarg++;
    batch = true;
    nrequired_args++;
  }

  // All optional args have been parsed.
  if (argc - iarg != nrequired_args) {
    fprintf(stderr, "Usage: grid_base_ref_miniapp.x [--batch "
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

  grid_library_init();

  const double max_diff =
      grid_collocate_replay(argv[iarg++], cycles, batch, cycles_per_block);

  grid_library_print_stats(&mpi_sum_func, &print_func);
  grid_library_finalize();

  if (max_diff > 1e-12 * cycles) {
    fprintf(stderr, "Error: Maximal difference is too large.\n");
    return 2;
  }

  return 0;
}

// EOF
