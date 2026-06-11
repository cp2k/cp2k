/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#include <opencl/libxstream_atomics.h>

#define SINT short

#if defined(PRECISION) && (1 == PRECISION)
#define CVT(A) convert_float(A)
#define TYPE float
#else
#define TYPE double
#define CVT(A) A
#endif

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

/* When exact shape is known (homogeneous batches), use K for unrolling
   but cap at 8 to avoid instruction-cache pressure for larger K values.
   Otherwise fall back to the BK threshold from the host. */
#if defined(DBM_M) && defined(DBM_N) && defined(DBM_K)
#undef BK
#if (DBM_K <= 8)
#define BK DBM_K
#else
#define BK 8
#endif
#elif !defined(BK) || (0 >= BK)
#define BK 1
#endif

/* K-block for BLKRD_P: cap at 8 to limit register pressure
   (a_reg[BKP] + c_acc[SG] must fit in available GRF). */
#if (BK <= 8)
#define BKP BK
#else
#define BKP 8
#endif

/* Broadcast tile size: use exact M when known and fits in one tile,
   avoiding wasted broadcast iterations when M < BN. */
#if defined(DBM_M) && (0 < DBM_M) && (DBM_M <= BN)
#define BM DBM_M
#else
#define BM BN
#endif

/* Skip zero contributions (NZ): avoids atomic overhead for zero values. */
#if defined(NZ) && (0 != NZ)
#define DBM_ACCUMULATE(PTR, VAL)                                               \
  do {                                                                         \
    const TYPE dbm_nz_v = (VAL);                                               \
    if ((TYPE)0 != dbm_nz_v) {                                                 \
      ACCUMULATE(PTR, dbm_nz_v);                                               \
    }                                                                          \
  } while (0)
#else
#define DBM_ACCUMULATE(PTR, VAL) ACCUMULATE(PTR, VAL)
#endif

/* When block reads are enabled, override SG for optimal tile geometry.
   For BLKRD_A (homogeneous): SG = DBM_M so each sub-group = one task.
   For BLKRD_P (heterogeneous): SG = 16 for efficient sub-group ops. */
#if defined(BLKRD_A) && defined(DBM_M) && defined(SG) && (DBM_M != SG)
#undef SG
#define SG DBM_M
#elif defined(BLKRD_P) && defined(SG) && (16 < SG)
#undef SG
#define SG 16
#endif

/* Sub-group block read for A: each lane loads one row-element from a
   contiguous column of A. Requires DBM_M == SG (flat dispatch). */
#if (defined(BLKRD_A) || defined(BLKRD_P)) && defined(SG) && (0 < SG)
#if defined(PRECISION) && (1 == PRECISION)
#define A_BLOCK_READ(PTR)                                                      \
  as_float(intel_sub_group_block_read((const global uint *)(PTR)))
#else
#define A_BLOCK_READ(PTR)                                                      \
  as_double(intel_sub_group_block_read_ul((const global ulong *)(PTR)))
#endif
#endif

#if defined(A_BLOCK_READ)
#define LOAD_A(A, IBASE, SHIFT, SHAPE, M, K)                                   \
  CVT(A_BLOCK_READ((A) + XA(SHIFT, IBASE) + (K) * XM(SHAPE)))
#else
#define LOAD_A(A, IBASE, SHIFT, SHAPE, M, K)                                   \
  CVT((A)[XA(SHIFT, IBASE) + IDT(M, K, XM(SHAPE), XK(SHAPE))])
#endif

/* Row-oriented store: flush CVEC[0..N1-1] to C[M, N0..N0+N1-1] */
#define DBM_MUL_STORE(ALPHA, IBASE, SHIFT, SHAPE, C, CVEC, M, N0, N1)          \
  do {                                                                         \
    UNROLL_AUTO for (SINT n = 0; n < (N1); ++n) {                              \
      DBM_ACCUMULATE((C) + XC(SHIFT, IBASE) +                                  \
                         XI(M, n + (N0), XM(SHAPE), XN(SHAPE)),                \
                     (ALPHA) * (CVEC)[n]);                                     \
    }                                                                          \
  } while (0)

/* Broadcast B values when all lanes in a sub-group handle the same task.
   Requires WG > 0 (guarantees intel_reqd_sub_group_size) and stride >= SG
   so sub-group boundaries align with task boundaries.  Only valid for
   homogeneous batches (DBM_M) where all lanes are unconditionally active. */
#if defined(WG) && (0 < WG) && defined(BCST_SG) && defined(DBM_M) &&           \
    (DBM_M >= SG)
#define LOAD_B(V) BCST_SG(V, 0)
#else
#define LOAD_B(V) (V)
#endif

/* Row-oriented kernel: one row M, columns N0..N0+BN-1 */
#define DBM_MUL_KERNEL(IBASE, SHIFT, SHAPE, A, B, CVEC, M, N0, BN)             \
  do {                                                                         \
    UNROLL_AUTO for (SINT k = 0; k < XK(SHAPE); ++k) {                         \
      const int ik = IDX(k, N0, XK(SHAPE), XN(SHAPE));                         \
      const TYPE ak = LOAD_A(A, IBASE, SHIFT, SHAPE, M, k);                    \
      UNROLL_AUTO for (SINT n = 0; n < (BN); ++n) {                            \
        (CVEC)[n] = MAD(ak, LOAD_B(CVT((B)[ik + n])), (CVEC)[n]);              \
      }                                                                        \
    }                                                                          \
  } while (0)

/* Row-oriented multiply: tiles N by BN, compile-time BK */
#define DBM_MUL(ALPHA, IBASE, SHIFT, SHAPE, A, B, C, CVEC, M, BN)              \
  do {                                                                         \
    SINT n0 = 0, n1 = XN(SHAPE) - (BN);                                        \
    UNROLL_FORCE(BN) for (SINT n = 0; n < (BN); ++n) { (CVEC)[n] = ZERO; }     \
    UNROLL_AUTO for (; n0 <= n1; n0 += (BN)) {                                 \
      DBM_MUL_KERNEL(IBASE, SHIFT, SHAPE, A, B, CVEC, M, n0, BN);              \
      DBM_MUL_STORE(ALPHA, IBASE, SHIFT, SHAPE, C, CVEC, M, n0, BN);           \
      UNROLL_FORCE(BN) for (SINT n = 0; n < (BN); ++n) { (CVEC)[n] = ZERO; }   \
    }                                                                          \
    n1 = XN(SHAPE) - n0;                                                       \
    DBM_MUL_KERNEL(IBASE, SHIFT, SHAPE, A, B, CVEC, M, n0, n1);                \
    DBM_MUL_STORE(ALPHA, IBASE, SHIFT, SHAPE, C, CVEC, M, n0, n1);             \
  } while (0)

#if defined(SGBCST) && defined(BCST_SG)
/* Column-oriented store: flush CVEC[0..ME-1] to C[MB..MB+ME-1, N] */
#define DBM_BCST_STORE(ALPHA, IBASE, SHIFT, SHAPE, C, CVEC, N, MB, ME)         \
  do {                                                                         \
    UNROLL_AUTO for (SINT m = 0; m < (ME); ++m) {                              \
      DBM_ACCUMULATE((C) + XC(SHIFT, IBASE) +                                  \
                         XI(m + (MB), N, XM(SHAPE), XN(SHAPE)),                \
                     (ALPHA) * (CVEC)[m]);                                     \
    }                                                                          \
  } while (0)

/* Column-oriented kernel: one column N, rows MB..MB+BN-1 via broadcast */
#define DBM_BCST_KERNEL(IBASE, SHIFT, SHAPE, A, B, CVEC, SID, N, MB, BN)       \
  do {                                                                         \
    UNROLL_AUTO for (SINT k = 0; k < XK(SHAPE); ++k) {                         \
      const TYPE ak =                                                          \
          CVT((A)[XA(SHIFT, IBASE) +                                           \
                  IDT(MIN((SINT)(SID) + (MB), (SINT)(XM(SHAPE) - 1)), k,       \
                      XM(SHAPE), XK(SHAPE))]);                                 \
      const TYPE bv = CVT((B)[IDX(k, N, XK(SHAPE), XN(SHAPE))]);               \
      UNROLL_FORCE(BN) for (SINT m = 0; m < (BN); ++m) {                       \
        (CVEC)[m] = MAD(BCST_SG(ak, (uint)m), bv, (CVEC)[m]);                  \
      }                                                                        \
    }                                                                          \
  } while (0)

/* Column-oriented multiply: tiles M by BN, compile-time BK.
   ACTIVE=0 for inactive lanes (participate in broadcasts, skip C writes). */
#define DBM_BCST(ALPHA, IBASE, SHIFT, SHAPE, A, B, C, CVEC, SID, N, ACTIVE,    \
                 BN)                                                           \
  do {                                                                         \
    const SINT xm = XM(SHAPE);                                                 \
    UNROLL_AUTO for (SINT mb = 0; mb < xm; mb += (BN)) {                       \
      UNROLL_FORCE(BN) for (SINT m = 0; m < (BN); ++m) { (CVEC)[m] = ZERO; }   \
      DBM_BCST_KERNEL(IBASE, SHIFT, SHAPE, A, B, CVEC, SID, N, mb, BN);        \
      if (ACTIVE) {                                                            \
        DBM_BCST_STORE(ALPHA, IBASE, SHIFT, SHAPE, C, CVEC, N, mb,             \
                       MIN(xm - mb, (SINT)(BN)));                              \
      }                                                                        \
    }                                                                          \
  } while (0)
#endif

/* Decode task parameters: shape, ibase, params offset */
#define DBM_TASK_DECODE(ITASK, TID, PFMT, PARAMS, SHAPE, IBASE)                \
  do {                                                                         \
    if (0 == (PFMT)) {                                                         \
      const int task = ((ITASK) + (TID)) * 6;                                  \
      (SHAPE)[0] = (PARAMS)[task + 0];                                         \
      (SHAPE)[1] = (PARAMS)[task + 1];                                         \
      (SHAPE)[2] = (PARAMS)[task + 2];                                         \
      (PARAMS) += task + 3;                                                    \
    } else {                                                                   \
      (SHAPE)[0] = (SINT)(0xFF & (PFMT));                                      \
      (SHAPE)[1] = 0xFF & ((PFMT) >> 8);                                       \
      (SHAPE)[2] = 0xFF & ((PFMT) >> 16);                                      \
      (PARAMS) += ((ITASK) + (TID)) * 3;                                       \
      (IBASE) = 1;                                                             \
    }                                                                          \
  } while (0)

/* Override shape with compile-time constants (homogeneous batches).
   Applied after DBM_TASK_DECODE to enable constant folding. */
#if defined(DBM_M) && defined(DBM_N) && defined(DBM_K)
#if !defined(CLINEAR)
#define DBM_SHAPE_OVERRIDE(SHAPE)                                              \
  do {                                                                         \
    (SHAPE)[0] = DBM_M;                                                        \
    (SHAPE)[1] = DBM_N;                                                        \
    (SHAPE)[2] = DBM_K;                                                        \
  } while (0)
#else
#define DBM_SHAPE_OVERRIDE(SHAPE)                                              \
  do {                                                                         \
    (SHAPE)[0] = DBM_N;                                                        \
    (SHAPE)[1] = DBM_M;                                                        \
    (SHAPE)[2] = DBM_K;                                                        \
  } while (0)
#endif
#else
#define DBM_SHAPE_OVERRIDE(SHAPE)                                              \
  do {                                                                         \
    (void)(SHAPE);                                                             \
  } while (0)
#endif

#if defined(WG) && (0 < WG)
__attribute__((reqd_work_group_size(WG, 1, 1)))
#if defined(SG) && (0 < SG) && defined(INTEL) && (0 != INTEL)
__attribute__((intel_reqd_sub_group_size(SG)))
#endif
#endif
kernel void
dbm_multiply(double alpha, int itask, int ntasks, int size, int param_format,
             CONSTANT const int *restrict params,
/* CLINEAR swaps a/b in the signature so that the host's
   fixed arg order (adata=arg6, bdata=arg7) transposes
   the access pattern for coalesced memory reads. */
#if !defined(CLINEAR)
             CONSTANT const double *restrict a,
             CONSTANT const double *restrict b,
#else
             CONSTANT const double *restrict b,
             CONSTANT const double *restrict a,
#endif
             global double *restrict c) {
#if defined(SM) && (0 < SM)
  local TYPE tls[WG][BN + SM - 1];
  local TYPE *restrict const cvec = &tls[get_local_id(0)][0];
#else
  TYPE cvec[BN];
#endif
#if defined(BLKRD_P) && defined(A_BLOCK_READ) && defined(BCST_SG)
  /* per-task dispatch: block-read A, B distributed across lanes.
     Each lane holds one M-row (via block read).  N is tiled by BN
     with sub_group_broadcast fanning B out to BN columns per K-step. */
  const int tid = (int)get_group_id(0);
  const SINT sid = (SINT)get_sub_group_local_id();
  SINT shape[3], ibase = 0;
  DBM_TASK_DECODE(itask, tid, param_format, params, shape, ibase);
  DBM_SHAPE_OVERRIDE(shape);
  {
    CONSTANT const double *restrict al = a + XA(params, ibase);
    CONSTANT const double *restrict bl = b + XB(params, ibase);
    const int c0 = XC(params, ibase);
    const SINT xm = XM(shape), xn = XN(shape), xk = XK(shape);
    const SINT nsg = (SINT)get_num_sub_groups();
    TYPE c_acc[SG];
    /* M-tiling: each sub-group handles SG consecutive rows */
    UNROLL_OUTER(1)
    for (SINT mb = (SINT)get_sub_group_id() * SG; mb < xm; mb += nsg * SG) {
      const SINT m = mb + sid;
      /* N-tiling by SG columns: lane sid holds B-column sid+n0 */
      UNROLL_AUTO for (SINT n0 = 0; n0 < xn; n0 += SG) {
        SINT k = 0;
        UNROLL_FORCE(SG) for (SINT i = 0; i < SG; ++i) c_acc[i] = ZERO;
        UNROLL_AUTO for (; k + BKP <= xk; k += BKP) {
          TYPE a_reg[BKP];
          if (mb + SG <= xm) {
            UNROLL_FORCE(BKP) for (SINT kb = 0; kb < BKP; ++kb) {
              a_reg[kb] = CVT(A_BLOCK_READ(al + (k + kb) * xm + mb));
            }
          } else {
            UNROLL_FORCE(BKP) for (SINT kb = 0; kb < BKP; ++kb) {
              a_reg[kb] = (m < xm) ? CVT(al[IDT(m, k + kb, xm, xk)]) : ZERO;
            }
          }
          UNROLL_FORCE(BKP) for (SINT kb = 0; kb < BKP; ++kb) {
            const TYPE bv =
                (sid + n0 < xn) ? CVT(bl[IDX(k + kb, sid + n0, xk, xn)]) : ZERO;
            UNROLL_FORCE(SG) for (SINT n = 0; n < SG; ++n) {
              c_acc[n] = MAD(a_reg[kb], BCST_SG(bv, (uint)n), c_acc[n]);
            }
          }
        }
        /* K remainder */
        UNROLL_AUTO for (; k < xk; ++k) {
          const TYPE ak = (mb + SG <= xm)
                              ? CVT(A_BLOCK_READ(al + k * xm + mb))
                              : ((m < xm) ? CVT(al[IDT(m, k, xm, xk)]) : ZERO);
          const TYPE bv =
              (sid + n0 < xn) ? CVT(bl[IDX(k, sid + n0, xk, xn)]) : ZERO;
          UNROLL_FORCE(SG) for (SINT n = 0; n < SG; ++n) {
            c_acc[n] = MAD(ak, BCST_SG(bv, (uint)n), c_acc[n]);
          }
        }
        /* store: only active M-rows and valid N-columns */
        if (m < xm) {
          const SINT ncols = MIN((SINT)SG, xn - n0);
          UNROLL_AUTO for (SINT n = 0; n < ncols; ++n) {
            DBM_ACCUMULATE(c + c0 + XI(m, n0 + n, xm, xn), alpha * c_acc[n]);
          }
        }
      }
    }
  }
#elif defined(SGBCST) && defined(BCST_SG) && !defined(BLKRD_A)
  /* per-task dispatch: broadcast shares A */
  const int tid = (int)get_group_id(0);
  const SINT sid = (SINT)get_local_id(0);
  SINT shape[3], ibase = 0;
  DBM_TASK_DECODE(itask, tid, param_format, params, shape, ibase);
  DBM_SHAPE_OVERRIDE(shape);
  b += XB(params, ibase);
  UNROLL_AUTO for (SINT nb0 = 0; nb0 < XN(shape); nb0 += WG) { /* all lanes */
    const SINT nb = nb0 + sid;
    const int active = (nb < XN(shape));
    DBM_BCST(alpha, ibase, params, shape, a, b, c, cvec, sid, active ? nb : 0,
             active, BM);
  }
#else
  /* flat dispatch: global work-item maps to (task, row) */
  const int i = (int)get_global_id(0);
#if defined(WG) && (0 < WG)
  if (i < size)
#endif
  {
    SINT shape[3], ibase = 0, m;
    int tid = i;
#if defined(DBM_M) && defined(DBM_N) && defined(DBM_K)
    shape[0] = DBM_M;
    shape[1] = DBM_N;
    shape[2] = DBM_K;
#elif defined(MAX_M)
    shape[0] = MAX_M;
#else
    shape[0] = (0 != param_format ? (SINT)(0xFF & param_format)
                                  : (SINT)(size / ntasks));
    shape[1] = 0xFF & (param_format >> 8);
    shape[2] = 0xFF & (param_format >> 16);
#endif
    tid /= shape[0];
    m = i - tid * shape[0];
    DBM_TASK_DECODE(itask, tid, param_format, params, shape, ibase);
    DBM_SHAPE_OVERRIDE(shape);
    if (m < XM(shape)) {
      b += XB(params, ibase);
      DBM_MUL(alpha, ibase, params, shape, a, b, c, cvec, m, BN);
    }
  }
#endif
}
