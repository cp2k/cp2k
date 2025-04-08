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

#if !defined(CLINEAR)
#define XM(T) T[0]
#define XN(T) T[1]
#define XI IDT
#else
#define XM(T) T[1]
#define XN(T) T[0]
#define XI IDX
#endif

#define XK(T) T[2]
#define XA(T, IBASE) (XM(T) - IBASE)
#define XB(T, IBASE) (XN(T) - IBASE)
#define XC(T, IBASE) (XK(T) - IBASE)

#define DBM_MULTIPLY_STORE(ALPHA, IBASE, SHIFT, SHAPE, C, CVEC, M, N0, N1)     \
  do { /* C atomically accumulates CVEC */                                     \
    UNROLL_AUTO for (SINT n = 0; n < (N1); ++n) { /* flush to global */        \
      const int im = XI(M, n + (N0), XM(SHAPE), XN(SHAPE));                    \
      ACCUMULATE((C) + XC(SHIFT, IBASE) + im, (ALPHA) * (CVEC)[n]);            \
    }                                                                          \
  } while (0)

#define DBM_MULTIPLY_KERNEL(IBASE, SHIFT, SHAPE, A, B, CVEC, M, N0, BN, BK)    \
  do { /* CVEC accumulates result */                                           \
    UNROLL(BK) for (SINT k = 0; k < XK(SHAPE); ++k) {                          \
      const int ik = IDX(k, N0, XK(SHAPE), XN(SHAPE));                         \
      const int ia = IDT(M, k, XM(SHAPE), XK(SHAPE));                          \
      const double ak = (A)[XA(SHIFT, IBASE) + ia];                            \
      UNROLL_AUTO for (SINT n = 0; n < (BN); ++n) {                            \
        (CVEC)[n] = MAD(ak, (B)[ik + n], (CVEC)[n]);                           \
      }                                                                        \
    }                                                                          \
  } while (0)

#define DBM_MULTIPLY(ALPHA, IBASE, SHIFT, SHAPE, A, B, C, CVEC, M, BN, BK)     \
  do { /* DBM_MULTIPLY_KERNEL specialized over N */                            \
    SINT n0 = 0, n1 = XN(SHAPE) - (BN);                                        \
    UNROLL_FORCE(BN) for (SINT n = 0; n < (BN); ++n) { (CVEC)[n] = ZERO; }     \
    UNROLL_OUTER(1) for (; n0 <= n1; n0 += (BN)) {                             \
      DBM_MULTIPLY_KERNEL(IBASE, SHIFT, SHAPE, A, B, CVEC, M, n0, BN, BK);     \
      DBM_MULTIPLY_STORE(ALPHA, IBASE, SHIFT, SHAPE, C, CVEC, M, n0, BN);      \
      UNROLL_FORCE(BN) for (SINT n = 0; n < (BN); ++n) { (CVEC)[n] = ZERO; }   \
    }                                                                          \
    n1 = XN(SHAPE) - n0;                                                       \
    DBM_MULTIPLY_KERNEL(IBASE, SHIFT, SHAPE, A, B, CVEC, M, n0, n1, BK);       \
    DBM_MULTIPLY_STORE(ALPHA, IBASE, SHIFT, SHAPE, C, CVEC, M, n0, n1);        \
  } while (0)

#if defined(WG) && (0 < WG)
__attribute__((reqd_work_group_size(WG, 1, 1)))
#if defined(SG) && (0 < SG)
__attribute__((intel_reqd_sub_group_size(SG)))
#endif
#endif
kernel void
dbm_multiply(double alpha, int itask, int ntasks, int size, int param_format,
             global const int *params,
#if !defined(CLINEAR)
             global const double *restrict a, global const double *restrict b,
#else
             global const double *restrict b, global const double *restrict a,
#endif
             global double *restrict c) {
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
    SINT shape[3], ibase = 0, m;
    int tid = i;
    shape[0] = size / ntasks;
    shape[1] = 0xFF & (param_format >> 8);
    shape[2] = 0xFF & (param_format >> 16);
    tid /= shape[0];
    m = i - tid * shape[0];
    if (0 == param_format) {
      const int task = (itask + tid) * 6;
      shape[0] = params[task + 0];
      shape[1] = params[task + 1];
      shape[2] = params[task + 2];
      params += task + 3;
    } else {
      params += (itask + tid) * 3;
      ibase = 1;
    }
#if !defined(NDEBUG)
    if (m < XM(shape))
#endif
    { /* valid slice (subtask) */
      b += XB(params, ibase);
      if (16 <= XK(shape)) {
        DBM_MULTIPLY(alpha, ibase, params, shape, a, b, c, cvec, m, BN, 16);
      } else if (8 <= XK(shape)) {
        DBM_MULTIPLY(alpha, ibase, params, shape, a, b, c, cvec, m, BN, 8);
      } else if (4 <= XK(shape)) {
        DBM_MULTIPLY(alpha, ibase, params, shape, a, b, c, cvec, m, BN, 4);
      } else {
        DBM_MULTIPLY(alpha, ibase, params, shape, a, b, c, cvec, m, BN, 1);
      }
    }
  }
}
#endif
