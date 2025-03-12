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
#define XC(T) X(T, offset_c)
#define XK(T) (SINT) X(T, k)

#if !defined(CLINEAR)
#define XA(T) X(T, offset_a)
#define XB(T) X(T, offset_b)
#define XM(T) (SINT) X(T, m)
#define XN(T) (SINT) X(T, n)
#define XI IDT
#else
#define XA(T) X(T, offset_b)
#define XB(T) X(T, offset_a)
#define XM(T) (SINT) X(T, n)
#define XN(T) (SINT) X(T, m)
#define XI IDX
#endif

#define DBM_MULTIPLY_STORE(ALPHA, TASK, CMAT, CVEC, M, N0, N1)                 \
  do { /* CMAT atomically accumulates CVEC */                                  \
    UNROLL_AUTO for (SINT n = 0; n < (N1); ++n) { /* flush to global */        \
      const int idx = XI(M, n + (N0), XM(TASK), XN(TASK)) + XC(TASK);          \
      ACCUMULATE((CMAT) + idx, (ALPHA) * (CVEC)[n]);                           \
    }                                                                          \
  } while (0)

#define DBM_MULTIPLY_KERNEL(TASK, AMAT, BMAT, CVEC, M, N0, BN, BK)             \
  do { /* CVEC accumulates result */                                           \
    UNROLL(BK) for (SINT k = 0; k < XK(TASK); ++k) {                           \
      const double a = (AMAT)[XA(TASK) + IDT(M, k, XM(TASK), XK(TASK))];       \
      const int idx = IDX(k, N0, XK(TASK), XN(TASK));                          \
      UNROLL_AUTO for (SINT n = 0; n < (BN); ++n) {                            \
        (CVEC)[n] = MAD(a, (BMAT)[idx + n], (CVEC)[n]);                        \
      }                                                                        \
    }                                                                          \
  } while (0)

#define DBM_MULTIPLY(ALPHA, TASK, AMAT, BMAT, CMAT, CVEC, M, BN, BK)           \
  do { /* DBM_MULTIPLY_KERNEL specialized over N */                            \
    SINT n0 = 0, n1 = XN(TASK) - (BN);                                         \
    UNROLL_FORCE(BN) for (SINT n = 0; n < (BN); ++n) { (CVEC)[n] = ZERO; }     \
    UNROLL_OUTER(1) for (; n0 <= n1; n0 += (BN)) {                             \
      DBM_MULTIPLY_KERNEL(TASK, AMAT, BMAT, CVEC, M, n0, BN, BK);              \
      DBM_MULTIPLY_STORE(ALPHA, TASK, CMAT, CVEC, M, n0, BN);                  \
      UNROLL_FORCE(BN) for (SINT n = 0; n < (BN); ++n) { (CVEC)[n] = ZERO; }   \
    }                                                                          \
    n1 = XN(TASK) - n0;                                                        \
    DBM_MULTIPLY_KERNEL(TASK, AMAT, BMAT, CVEC, M, n0, n1, BK);                \
    DBM_MULTIPLY_STORE(ALPHA, TASK, CMAT, CVEC, M, n0, n1);                    \
  } while (0)

#if defined(WG) && (0 < WG)
__attribute__((reqd_work_group_size(WG, 1, 1)))
#if defined(SG) && (0 < SG)
__attribute__((intel_reqd_sub_group_size(SG)))
#endif
#endif
kernel void
dbm_multiply(double alpha, int itask, int ntasks, int size,
             global const dbm_task_t *tasks,
#if !defined(CLINEAR)
             global const double *restrict amat,
             global const double *restrict bmat,
#else
             global const double *restrict bmat,
             global const double *restrict amat,
#endif
             global double *restrict cmat) {
  const int i = (int)get_global_id(0);
#if defined(SM) && (0 < SM)
  local double tls[WG][BN + SM - 1], *const cvec = &tls[get_local_id(0)];
#else
  double cvec[BN];
#endif
#if defined(WG) && (0 < WG)
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
      if (16 <= XK(task)) {
        DBM_MULTIPLY(alpha, task, amat, bmat, cmat, cvec, m, BN, 16);
      } else if (8 <= XK(task)) {
        DBM_MULTIPLY(alpha, task, amat, bmat, cmat, cvec, m, BN, 8);
      } else if (4 <= XK(task)) {
        DBM_MULTIPLY(alpha, task, amat, bmat, cmat, cvec, m, BN, 4);
      } else {
        DBM_MULTIPLY(alpha, task, amat, bmat, cmat, cvec, m, BN, 1);
      }
    }
  }
}
#endif
