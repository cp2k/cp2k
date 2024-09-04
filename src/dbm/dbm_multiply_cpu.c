/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2024 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <stddef.h>
#include <string.h>

#if defined(__LIBXSMM)
#include <libxsmm.h>
#if !defined(DBM_LIBXSMM_PREFETCH)
// #define DBM_LIBXSMM_PREFETCH LIBXSMM_GEMM_PREFETCH_AL2_AHEAD
#define DBM_LIBXSMM_PREFETCH LIBXSMM_GEMM_PREFETCH_NONE
#endif
#if LIBXSMM_VERSION4(1, 17, 0, 3710) > LIBXSMM_VERSION_NUMBER
#define libxsmm_dispatch_gemm libxsmm_dispatch_gemm_v2
#endif
#endif

#include "dbm_hyperparams.h"
#include "dbm_multiply_cpu.h"

/*******************************************************************************
 * \brief Prototype for BLAS dgemm.
 * \author Ole Schuett
 ******************************************************************************/
void dgemm_(const char *transa, const char *transb, const int *m, const int *n,
            const int *k, const double *alpha, const double *a, const int *lda,
            const double *b, const int *ldb, const double *beta, double *c,
            const int *ldc);

/*******************************************************************************
 * \brief Private convenient wrapper to hide Fortran nature of dgemm_.
 * \author Ole Schuett
 ******************************************************************************/
static inline void dbm_dgemm(const char transa, const char transb, const int m,
                             const int n, const int k, const double alpha,
                             const double *a, const int lda, const double *b,
                             const int ldb, const double beta, double *c,
                             const int ldc) {

  dgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c,
         &ldc);
}

/*******************************************************************************
 * \brief Private hash function based on Szudzik's elegant pairing.
 *        Using unsigned int to return a positive number even after overflow.
 *        https://en.wikipedia.org/wiki/Pairing_function#Other_pairing_functions
 *        https://stackoverflow.com/a/13871379
 *        http://szudzik.com/ElegantPairing.pdf
 * \author Ole Schuett
 ******************************************************************************/
static inline unsigned int hash(const dbm_task_t task) {
  const unsigned int m = task.m, n = task.n, k = task.k;
  const unsigned int mn = (m >= n) ? m * m + m + n : m + n * n;
  const unsigned int mnk = (mn >= k) ? mn * mn + mn + k : mn + k * k;
  return mnk;
}

/*******************************************************************************
 * \brief Internal routine for executing the tasks in given batch on the CPU.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply_cpu_process_batch(const int ntasks, dbm_task_t batch[ntasks],
                                    const double alpha,
                                    const dbm_pack_t *pack_a,
                                    const dbm_pack_t *pack_b,
                                    dbm_shard_t *shard_c) {

  if (0 >= ntasks) { // nothing to do
    return;
  }
  dbm_shard_allocate_promised_blocks(shard_c);

#if defined(__LIBXSMM)

  // Sort tasks approximately by m,n,k via bucket sort.
  int buckets[BATCH_NUM_BUCKETS];
  memset(buckets, 0, BATCH_NUM_BUCKETS * sizeof(int));
  for (int itask = 0; itask < ntasks; ++itask) {
    const int i = hash(batch[itask]) % BATCH_NUM_BUCKETS;
    ++buckets[i];
  }
  for (int i = 1; i < BATCH_NUM_BUCKETS; ++i) {
    buckets[i] += buckets[i - 1];
  }
  assert(buckets[BATCH_NUM_BUCKETS - 1] == ntasks);
  int batch_order[ntasks];
  for (int itask = 0; itask < ntasks; ++itask) {
    const int i = hash(batch[itask]) % BATCH_NUM_BUCKETS;
    --buckets[i];
    batch_order[buckets[i]] = itask;
  }

  // Prepare arguments for libxsmm's kernel-dispatch.
  const int flags = LIBXSMM_GEMM_FLAG_TRANS_B; // transa = "N", transb = "T"
  const int prefetch = DBM_LIBXSMM_PREFETCH;
  int kernel_m = 0, kernel_n = 0, kernel_k = 0;
  dbm_task_t task_next = batch[batch_order[0]];

#if (LIBXSMM_GEMM_PREFETCH_NONE != DBM_LIBXSMM_PREFETCH)
  double *data_a_next = NULL, *data_b_next = NULL, *data_c_next = NULL;
#endif
#if LIBXSMM_VERSION2(1, 17) < LIBXSMM_VERSION_NUMBER
  libxsmm_gemmfunction kernel_func = NULL;
#else
  libxsmm_dmmfunction kernel_func = NULL;
  const double beta = 1.0;
#endif

  // Loop over tasks.
  for (int itask = 0; itask < ntasks; ++itask) {
    const dbm_task_t task = task_next;
    task_next = batch[batch_order[(itask + 1) < ntasks ? (itask + 1) : itask]];

    if (task.m != kernel_m || task.n != kernel_n || task.k != kernel_k) {
#if LIBXSMM_VERSION2(1, 17) < LIBXSMM_VERSION_NUMBER
      const libxsmm_gemm_shape shape = libxsmm_create_gemm_shape(
          task.m, task.n, task.k, task.m /*lda*/, task.n /*ldb*/,
          task.m /*ldc*/, LIBXSMM_DATATYPE_F64 /*aprec*/,
          LIBXSMM_DATATYPE_F64 /*bprec*/, LIBXSMM_DATATYPE_F64 /*cprec*/,
          LIBXSMM_DATATYPE_F64 /*calcp*/);
      kernel_func = (LIBXSMM_FEQ(1.0, alpha)
                         ? libxsmm_dispatch_gemm(shape, (libxsmm_bitfield)flags,
                                                 (libxsmm_bitfield)prefetch)
                         : NULL);
#else
      kernel_func = libxsmm_dmmdispatch(task.m, task.n, task.k, NULL /*lda*/,
                                        NULL /*ldb*/, NULL /*ldc*/, &alpha,
                                        &beta, &flags, &prefetch);
#endif
      kernel_m = task.m;
      kernel_n = task.n;
      kernel_k = task.k;
    }

    // gemm_param wants non-const data even for A and B
    double *const data_a = pack_a->data + task.offset_a;
    double *const data_b = pack_b->data + task.offset_b;
    double *const data_c = shard_c->data + task.offset_c;

    if (kernel_func != NULL) {
#if LIBXSMM_VERSION2(1, 17) < LIBXSMM_VERSION_NUMBER
      libxsmm_gemm_param gemm_param;
      gemm_param.a.primary = data_a;
      gemm_param.b.primary = data_b;
      gemm_param.c.primary = data_c;
#if (LIBXSMM_GEMM_PREFETCH_NONE != DBM_LIBXSMM_PREFETCH)
      gemm_param.a.quaternary = pack_a->data + task_next.offset_a;
      gemm_param.b.quaternary = pack_b->data + task_next.offset_b;
      gemm_param.c.quaternary = shard_c->data + task_next.offset_c;
#endif
      kernel_func(&gemm_param);
#elif (LIBXSMM_GEMM_PREFETCH_NONE != DBM_LIBXSMM_PREFETCH)
      kernel_func(data_a, data_b, data_c, pack_a->data + task_next.offset_a,
                  pack_b->data + task_next.offset_b,
                  shard_c->data + task_next.offset_c);
#else
      kernel_func(data_a, data_b, data_c);
#endif
    } else {
      dbm_dgemm('N', 'T', task.m, task.n, task.k, alpha, data_a, task.m, data_b,
                task.n, 1.0, data_c, task.m);
    }
  }
#else
  // Fallback to BLAS when libxsmm is not available.
  for (int itask = 0; itask < ntasks; ++itask) {
    const dbm_task_t task = batch[itask];
    const double *data_a = &pack_a->data[task.offset_a];
    const double *data_b = &pack_b->data[task.offset_b];
    double *data_c = &shard_c->data[task.offset_c];
    dbm_dgemm('N', 'T', task.m, task.n, task.k, alpha, data_a, task.m, data_b,
              task.n, 1.0, data_c, task.m);
  }
#endif
}

// EOF
