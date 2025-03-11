/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2025 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#if defined(DBM_MULTIPLY_OPENCL_GEN)
#include "dbm_multiply_opencl.irh"
#else
#include "../../exts/dbcsr/src/acc/opencl/common/opencl_atomics.h"
#include "dbm_internal.h"

#define SINT short

#define X(T, I) (T)->I /* task can be taken by value or by pointer */
#define XA(T) X(T, offset_a)
#define XB(T) X(T, offset_b)
#define XC(T) X(T, offset_c)
#define XM(T) (SINT) X(T, m)
#define XN(T) (SINT) X(T, n)
#define XK(T) (SINT) X(T, k)

#define DBM_MULTIPLY_STORE(ALPHA, TASK, CMAT, CVEC, M, N0, N1)                 \
  do { /* CMAT atomically accumulates CVEC */                                  \
    UNROLL_AUTO for (SINT n = 0; n < (N1); ++n) { /* flush to global */        \
      const int idx = IDT(M, n + (N0), XM(TASK), XN(TASK)) + XC(TASK);         \
      ACCUMULATE((CMAT) + idx, (ALPHA) * (CVEC)[n]);                           \
    }                                                                          \
  } while (0)

#define DBM_MULTIPLY_KERNEL(TASK, AMAT, BMAT, CVEC, M, N0, BN)                 \
  do { /* CVEC accumulates result */                                           \
    UNROLL_AUTO for (SINT k = 0; k < XK(TASK); ++k) {                          \
      const double a = (AMAT)[XA(TASK) + IDT(M, k, XM(TASK), XK(TASK))];       \
      const int idx = IDX(k, N0, XK(TASK), XN(TASK));                          \
      UNROLL_AUTO for (SINT n = 0; n < (BN); ++n) {                            \
        (CVEC)[n] = MAD(a, (BMAT)[idx + n], (CVEC)[n]);                        \
      }                                                                        \
    }                                                                          \
  } while (0)

#define DBM_MULTIPLY(ALPHA, TASK, AMAT, BMAT, CMAT, CVEC, M, BN)               \
  do { /* DBM_MULTIPLY_KERNEL specialized over N */                            \
    SINT n0 = 0, n1 = XN(TASK) - (BN);                                         \
    UNROLL_FORCE(BN) for (SINT n = 0; n < (BN); ++n) { (CVEC)[n] = ZERO; }     \
    UNROLL_OUTER(1) for (; n0 <= n1; n0 += (BN)) {                             \
      DBM_MULTIPLY_KERNEL(TASK, AMAT, BMAT, CVEC, M, n0, BN);                  \
      DBM_MULTIPLY_STORE(ALPHA, TASK, CMAT, CVEC, M, n0, BN);                  \
      UNROLL_FORCE(BN) for (SINT n = 0; n < (BN); ++n) { (CVEC)[n] = ZERO; }   \
    }                                                                          \
    DBM_MULTIPLY_KERNEL(TASK, AMAT, BMAT, CVEC, M, n0, XN(TASK) - n0);         \
    DBM_MULTIPLY_STORE(ALPHA, TASK, CMAT, CVEC, M, n0, XN(TASK) - n0);         \
  } while (0)

#if defined(WG) && (0 < WG)
__attribute__((reqd_work_group_size(WG, 1, 1)))
#if defined(SG) && (0 < SG)
__attribute__((intel_reqd_sub_group_size(SG)))
#endif
#endif
kernel void
dbm_multiply(double alpha, int itask, int ntasks, int size,
             global const dbm_task_t *tasks, global const double *restrict amat,
             global const double *restrict bmat, global double *restrict cmat) {
  const int i = (int)get_global_id(0);
#if defined(SM) && (0 < SM)
  local double tls[WG][BN + SM - 1], *const cvec = &tls[get_local_id(0)];
#else
  double cvec[BN];
#endif
#if defined(WG) && (0 < WG) && !defined(NDEBUG)
  if (i < size)
#endif
  { /* valid task */
    const int max_m = size / ntasks, tid = i / max_m;
    const SINT m = i - tid * max_m;
    global const dbm_task_t *const task = &tasks[itask + tid];
#if !defined(NDEBUG)
    if (m < XM(task))
#endif
    { /* valid slice (subtask) */
      bmat += XB(task);
      DBM_MULTIPLY(alpha, task, amat, bmat, cmat, cvec, m, BN);
    }
  }
}
#endif
