/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2024 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#include "../../exts/dbcsr/src/acc/opencl/common/opencl_atomics.h"
#include "dbm_multiply_internal.h"

#define IDX(I, J, K, M, N) ((I) * (N) + (J) + (K))
#define IDT(I, J, K, M, N) IDX(J, I, K, N, M)

#define DBM_MULTIPLY_KERNEL(TASK, A, B, CV, M, N0, N1, UNROLL_N, UNROLL_K)     \
  do {                                                                         \
    UNROLL_K                                                                   \
    for (int k = 0; k < (TASK).k; ++k) {                                       \
      const int ia = IDT(M, k, (TASK).offset_a, (TASK).m, (TASK).k);           \
      const double a = (A)[ia];                                                \
      UNROLL_N                                                                 \
      for (int n = 0; n < (N1); ++n) {                                         \
        const int ib = IDX(k, n + (N0), (TASK).offset_b, (TASK).k, (TASK).n);  \
        (CV)[n] = MAD(a, (B)[ib], (CV)[n]);                                    \
      }                                                                        \
    }                                                                          \
  } while (0)

#define DBM_MULTIPLY_ACCUMULATE(ALPHA, TASK, CV, C, M, BN, N0)                 \
  do { /* flush private accumulator to global memory using atomics */          \
    UNROLL_FORCE(BN)                                                           \
    for (int n = 0; n < (BN); ++n) {                                           \
      const int ic = IDT(M, n + (N0), (TASK).offset_c, (TASK).m, (TASK).n);    \
      ACCUMULATE((C) + ic, (ALPHA) * (CV)[n]);                                 \
      cv[n] = ZERO; /* reset */                                                \
    }                                                                          \
  } while (0)

kernel void dbm_multiply(double alpha, int itask, int ntasks,
                         global const dbm_task_t *tasks,
                         global const double *restrict am,
                         global const double *restrict bm,
                         global double *restrict cm) {
  double cv[BN] = {0}; /* private accumulator */
  const int size = (int)get_global_size(0), i = (int)get_global_id(0);

  if (size != ntasks) { /* LIBDBM_TASK_SPLIT */
    const int max_m = size / ntasks, tid = i / max_m;
    const dbm_task_t task = tasks[itask + min(tid, ntasks - 1)]; /* copy */
    const int m = i - tid * max_m;
    if (m < task.m) {
      if ((BN) < task.n) {
        UNROLL_AUTO
        for (int n0 = 0; n0 < task.n; n0 += (BN)) {
          const int n1 = min(BN, task.n - n0);
          DBM_MULTIPLY_KERNEL(task, am, bm, cv, m, n0, n1, UNROLL_FORCE(BN),
                              UNROLL_AUTO);
          DBM_MULTIPLY_ACCUMULATE(alpha, task, cv, cm, m, BN, n0);
        }
      } else { /* small */
        DBM_MULTIPLY_KERNEL(task, am, bm, cv, m, 0, task.n, UNROLL_FORCE(BN),
                            UNROLL_FORCE(BN));
        DBM_MULTIPLY_ACCUMULATE(alpha, task, cv, cm, m, BN, 0);
      }
    }
  } else { /* full matrix multiplication */
    const dbm_task_t task = tasks[itask + i];
    UNROLL_OUTER(1)
    for (int m = 0; m < task.m; ++m) {
      UNROLL(1)
      for (int n0 = 0; n0 < task.n; n0 += (BN)) {
        const int n1 = min(BN, task.n - n0);
        DBM_MULTIPLY_KERNEL(task, am, bm, cv, m, n0, n1, UNROLL_AUTO,
                            UNROLL_AUTO);
        DBM_MULTIPLY_ACCUMULATE(alpha, task, cv, cm, m, BN, n0);
      }
    }
  }
}
