/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2024 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#if defined(DBM_MULTIPLY_OPENCL_GEN)
#include "dbm_multiply_opencl.ir.h"
#else
#include "../../exts/dbcsr/src/acc/opencl/common/opencl_atomics.h"
#include "dbm_multiply_internal.h"

#if defined(GPU) && defined(WG) && (0 < WG) && (200 <= ACC_OPENCL_VERSION)
#if defined(SG) && (0 < SG)
#define BCST_WG(V) sub_group_broadcast(V, 0)
#else
#define BCST_WG(V) work_group_broadcast(V, 0)
#endif
#endif
#define BCST_NO(V) (V)

#define SINT short

#define X(T, I) (T)->I /* task can be taken by value or by pointer */
#define XA(T) X(T, offset_a)
#define XB(T) X(T, offset_b)
#define XC(T) X(T, offset_c)
#define XM(T) (SINT) X(T, m)
#define XN(T) (SINT) X(T, n)
#define XK(T) (SINT) X(T, k)

#define DBM_MULTIPLY_SHM(ALPHA, TASK, AMAT, BMAT, CMAT, SHM, WG, BM, BN)       \
  do { /* matrix multiplication per work-group using shared memory */          \
    local double *restrict const ashm = (SHM);                                 \
    local double *restrict const bshm = (SHM) + (WG);                          \
    const int mk = XM(TASK) * XK(TASK), kn = XK(TASK) * XN(TASK);              \
    const SINT tid = (SINT)get_local_id(0);                                    \
    /* y/s can exceed BN/BM (up to BK), and x/t is fast index (up to BM/BN) */ \
    const SINT y = tid / (BM), x = tid - y * (BM), bk = (WG) / MAX(BM, BN);    \
    const SINT s = tid / (BN), t = tid - s * (BN);                             \
    for (SINT m0 = 0; m0 < XM(TASK); m0 += (BM)) {                             \
      for (SINT n0 = 0; n0 < XN(TASK); n0 += (BN)) {                           \
        double r = ZERO;                                                       \
        UNROLL_AUTO for (SINT k0 = 0; k0 < XK(TASK); k0 += bk) {               \
          if (x < (BM) && y < bk) { /* load A-tile */                          \
            const int idx = IDT(m0 + x, k0 + y, XM(TASK), XK(TASK));           \
            ashm[y * (BM) + x] = (idx < mk ? (AMAT)[XA(TASK) + idx] : ZERO);   \
          }                                                                    \
          if (s < bk && t < (BN)) { /* load B-tile */                          \
            const int idx = IDX(k0 + s, n0 + t, XK(TASK), XN(TASK));           \
            bshm[s * (BN) + t] = (idx < kn ? (BMAT)[XB(TASK) + idx] : ZERO);   \
          }                                                                    \
          BARRIER(CLK_LOCAL_MEM_FENCE);                                        \
          if (x < (BM) && y < (BN)) { /* multiply tiles */                     \
            UNROLL_AUTO for (SINT z = 0; z < bk; ++z) {                        \
              r = MAD(ashm[z * (BM) + x], bshm[z * (BN) + y], r);              \
            }                                                                  \
          }                                                                    \
          BARRIER(CLK_LOCAL_MEM_FENCE);                                        \
        }                                                                      \
        if (x < (BM) && y < (BN)) { /* flush to global */                      \
          const SINT m = m0 + x, n = n0 + y;                                   \
          if (m < XM(TASK) && n < XN(TASK)) {                                  \
            const int idx = IDT(m, n, XM(TASK), XN(TASK));                     \
            ACCUMULATE((CMAT) + XC(TASK) + idx, (ALPHA) * r);                  \
          }                                                                    \
        }                                                                      \
      }                                                                        \
    }                                                                          \
  } while (0)

#define DBM_MULTIPLY_KERNEL(ALPHA, TASK, AMAT, BMAT, CMAT, CVEC, M, N0, N1, K, \
                            BCST)                                              \
  UNROLL_AUTO for (SINT k = 0; k < (K); ++k) {                                 \
    const double a = (AMAT)[XA(TASK) + IDT(M, k, XM(TASK), K)];                \
    UNROLL_AUTO for (SINT n = 0; n < (N1); ++n) {                              \
      const double b = (BMAT)[XB(TASK) + IDX(k, n + (N0), K, XN(TASK))];       \
      (CVEC)[n] = MAD(a, BCST(b), (CVEC)[n]);                                  \
    }                                                                          \
  }                                                                            \
  UNROLL_AUTO for (SINT n = 0; n < (N1); ++n) { /* flush to global */          \
    const int idx = IDT(M, n + (N0), XM(TASK), XN(TASK));                      \
    ACCUMULATE((CMAT) + XC(TASK) + idx, (ALPHA) * (CVEC)[n]);                  \
    (CVEC)[n] = ZERO; /* reset */                                              \
  }

#define DBM_MULTIPLY(ALPHA, TASK, AMAT, BMAT, CMAT, M, BN, BCST)               \
  do { /* DBM_MULTIPLY_KERNEL unrolled/specialized over N and K */             \
    double cvec[BN];                                                           \
    SINT n0 = 0;                                                               \
    UNROLL_AUTO for (SINT n = 0; n < (BN); ++n) { cvec[n] = ZERO; }            \
    if ((BN) <= XN(TASK)) {                                                    \
      if (1 < XK(TASK)) {                                                      \
        UNROLL_OUTER(1) for (; (n0 + (BN)) <= XN(TASK); n0 += (BN)) {          \
          DBM_MULTIPLY_KERNEL(ALPHA, TASK, AMAT, BMAT, CMAT, cvec, M, n0, BN,  \
                              XK(TASK), BCST);                                 \
        }                                                                      \
      } else { /* K = 1 */                                                     \
        UNROLL_OUTER(1) for (; (n0 + (BN)) <= XN(TASK); n0 += (BN)) {          \
          DBM_MULTIPLY_KERNEL(ALPHA, TASK, AMAT, BMAT, CMAT, cvec, M, n0, BN,  \
                              1, BCST);                                        \
        }                                                                      \
      }                                                                        \
    } else if (1 != XK(TASK)) { /* N < BN */                                   \
      DBM_MULTIPLY_KERNEL(ALPHA, TASK, AMAT, BMAT, CMAT, cvec, M, 0, 1,        \
                          XK(TASK), BCST);                                     \
      n0 = 1;                                                                  \
    } else { /* N < BN, K = 1 */                                               \
      DBM_MULTIPLY_KERNEL(ALPHA, TASK, AMAT, BMAT, CMAT, cvec, M, 0, 1, 1,     \
                          BCST);                                               \
      n0 = 1;                                                                  \
    }                                                                          \
    /*if (n0 < XN(TASK))*/ { /* handle remainder */                            \
      DBM_MULTIPLY_KERNEL(ALPHA, TASK, AMAT, BMAT, CMAT, cvec, M, n0,          \
                          XN(TASK) - n0, XK(TASK), BCST);                      \
    }                                                                          \
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
#if defined(SPLIT) && (1 < SPLIT) && defined(WG) && (0 < WG)
  local double shm[WG * 2];
  global const dbm_task_t *const task = &tasks[itask + get_group_id(0)];
  const SINT rmin = MIN(XM(task), XN(task)), rmax = MAX(XM(task), XN(task));
  if ((rmax - rmin) <= BN) {
    if ((rmin * 4) < BN) {
      DBM_MULTIPLY_SHM(alpha, task, amat, bmat, cmat, shm, WG, BN / 4, BN / 4);
    } else if ((rmin * 2) < BN) {
      DBM_MULTIPLY_SHM(alpha, task, amat, bmat, cmat, shm, WG, BN / 2, BN / 2);
    } else {
      DBM_MULTIPLY_SHM(alpha, task, amat, bmat, cmat, shm, WG, BN, BN);
    }
  } else if (XM(task) <= XN(task)) {
    const SINT r1 = BLR(XM(task), BN);
    const SINT r2 = BLR(XM(task), BN / 2) * 2;
    const SINT r3 = BLR(XM(task), BN / 4) * 4;
    if (r1 <= r2) {
      if (r1 <= r3) {
        DBM_MULTIPLY_SHM(alpha, task, amat, bmat, cmat, shm, WG, BN, BN);
      } else {
        DBM_MULTIPLY_SHM(alpha, task, amat, bmat, cmat, shm, WG, BN / 4,
                         BN * 4);
      }
    } else if (r2 <= r3) {
      DBM_MULTIPLY_SHM(alpha, task, amat, bmat, cmat, shm, WG, BN / 2, BN * 2);
    } else {
      DBM_MULTIPLY_SHM(alpha, task, amat, bmat, cmat, shm, WG, BN / 4, BN * 4);
    }
  } else {
    const SINT r1 = BLR(XN(task), BN);
    const SINT r2 = BLR(XN(task), BN / 2) * 2;
    const SINT r3 = BLR(XN(task), BN / 4) * 4;
    if (r1 <= r2) {
      if (r1 <= r3) {
        DBM_MULTIPLY_SHM(alpha, task, amat, bmat, cmat, shm, WG, BN, BN);
      } else {
        DBM_MULTIPLY_SHM(alpha, task, amat, bmat, cmat, shm, WG, BN * 4,
                         BN / 4);
      }
    } else if (r2 <= r3) {
      DBM_MULTIPLY_SHM(alpha, task, amat, bmat, cmat, shm, WG, BN * 2, BN / 2);
    } else {
      DBM_MULTIPLY_SHM(alpha, task, amat, bmat, cmat, shm, WG, BN * 4, BN / 4);
    }
  }
#elif defined(SPLIT) && (0 != SPLIT)
  const int i = (int)get_global_id(0);
#if defined(BCST_WG)
  if (i < size)
#endif
  { /* DBM_MULTIPLY_SPLIT */
    const int max_m = size / ntasks, tid = i / max_m;
    const SINT m = i - tid * max_m;
    global const dbm_task_t *const task = &tasks[itask + tid];
    if (m < XM(task)) { /* valid task */
#if defined(BCST_WG)
      if (XM(task) <= XN(task)) { /* BCST_WG to broadcast B-values */
        DBM_MULTIPLY(alpha, task, amat, bmat, cmat, m, BN, BCST_WG);
      } else
#endif
      {
        DBM_MULTIPLY(alpha, task, amat, bmat, cmat, m, BN, BCST_NO);
      }
    }
  }
#else
#if defined(BCST_WG)
  if (get_global_id(0) < size)
#endif
  { /* full matrix multiplication per work-item (thread) */
    global const dbm_task_t *const task = &tasks[itask + get_global_id(0)];
    UNROLL_OUTER(1) for (SINT m = 0; m < XM(task); ++m) {
      DBM_MULTIPLY(alpha, task, amat, bmat, cmat, m, BN, BCST_NO);
    }
  }
#endif
}
#endif
