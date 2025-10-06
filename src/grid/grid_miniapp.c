/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2025 CP2K developers group <https://cp2k.org>              */
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

  const double tolerance = 1e-12 * cycles;
  const bool success = grid_replay(argv[iarg++], cycles, collocate, batch,
                                   cycles_per_block, tolerance);

  if (success) {
    grid_library_print_stats(0 /*fortran_comm*/, &print_func, 0 /*rank*/);
    offload_mempool_stats_print(0 /*fortran_comm*/, &print_func, 0 /*rank*/);
    grid_library_finalize();
  } else {
    fprintf(stderr, "Error: Maximal difference is too large.\n");
    grid_library_finalize();
    return 2;
  }

  return 0;
}

// EOF
