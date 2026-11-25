/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#define _GNU_SOURCE
#include <mpi.h>
#include <omp.h>
#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

static double seconds(clockid_t clock) {
  struct timespec ts;
  clock_gettime(clock, &ts);
  return ts.tv_sec + 1.e-9 * ts.tv_nsec;
}

int main(int argc, char **argv) {
  int rank, size, provided, threads = 0;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (provided < MPI_THREAD_FUNNELED || omp_get_max_threads() > 256)
    MPI_Abort(MPI_COMM_WORLD, 2);
  cpu_set_t masks[256];
  double cpu_seconds[256], values[256];
  MPI_Barrier(MPI_COMM_WORLD);
  double start = seconds(CLOCK_REALTIME);
#pragma omp parallel
  {
    int t = omp_get_thread_num();
#pragma omp single
    threads = omp_get_num_threads();
    if (sched_getaffinity(0, sizeof(cpu_set_t), &masks[t]) != 0)
      abort();
    double begin = seconds(CLOCK_THREAD_CPUTIME_ID), x = 1. + t * .001;
    for (int k = 0; k < 30000000; k++)
      x = x * 1.0000000000000002 + 1.e-12;
    values[t] = x;
    cpu_seconds[t] = seconds(CLOCK_THREAD_CPUTIME_ID) - begin;
  }
  char path[128];
  snprintf(path, sizeof(path), "rank-%d.json", rank);
  FILE *out = fopen(path, "w");
  if (!out)
    MPI_Abort(MPI_COMM_WORLD, 3);
  fprintf(out,
          "{\"rank\":%d,\"ranks\":%d,\"pid\":%d,\"start\":%.9f,"
          "\"omp_threads\":%d,\"threads\":[",
          rank, size, getpid(), start, threads);
  for (int t = 0; t < threads; t++) {
    fprintf(out,
            "%s{\"thread\":%d,\"cpu_seconds\":%.9f,\"value\":%.12f,"
            "\"affinity\":[",
            t ? "," : "", t, cpu_seconds[t], values[t]);
    int comma = 0;
    for (int cpu = 0; cpu < CPU_SETSIZE; cpu++)
      if (CPU_ISSET(cpu, &masks[t])) {
        fprintf(out, "%s%d", comma ? "," : "", cpu);
        comma = 1;
      }
    fprintf(out, "]}");
  }
  fprintf(out, "]}\n");
  fclose(out);
  MPI_Barrier(MPI_COMM_WORLD);
  const char *hold = getenv("CP2K_PROBE_HOLD");
  sleep(hold ? (unsigned)atoi(hold) : 5);
  MPI_Finalize();
  return 0;
}
