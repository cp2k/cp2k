/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#include "dbm_multiply_cpu.h"
#include "dbm_hyperparams.h"

#include <assert.h>
#include <stddef.h>

#if defined(__LIBXSMM)
#include <libxsmm.h>
#endif
#if defined(__LIBXS)
#include <libxs/libxs_gemm.h>
#endif

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
void dbm_multiply_cpu_process_batch(int ntasks, const dbm_task_t batch[ntasks],
                                    double alpha, const dbm_pack_t *pack_a,
                                    const dbm_pack_t *pack_b,
                                    dbm_shard_t *shard_c, int options) {
  if (0 >= ntasks) { // nothing to do
    return;
  }
  dbm_shard_allocate_promised_blocks(shard_c);

  int batch_order[ntasks];
  if (DBM_MULTIPLY_TASK_REORDER & options) {
    // Sort tasks approximately by m,n,k via bucket sort.
    int buckets[DBM_BATCH_NUM_BUCKETS] = {0};
    for (int itask = 0; itask < ntasks; ++itask) {
      const int i = hash(batch[itask]) % DBM_BATCH_NUM_BUCKETS;
      ++buckets[i];
    }
    for (int i = 1; i < DBM_BATCH_NUM_BUCKETS; ++i) {
      buckets[i] += buckets[i - 1];
    }
    assert(buckets[DBM_BATCH_NUM_BUCKETS - 1] == ntasks);
    for (int itask = 0; itask < ntasks; ++itask) {
      const int i = hash(batch[itask]) % DBM_BATCH_NUM_BUCKETS;
      --buckets[i];
      batch_order[buckets[i]] = itask;
    }
  } else {
    for (int itask = 0; itask < ntasks; ++itask) {
      batch_order[itask] = itask;
    }
  }

#if defined(__LIBXS)
  const libxs_gemm_config_t *gemm_config = NULL;
  int kernel_m = 0, kernel_n = 0, kernel_k = 0;
#endif

  // Loop over tasks.
  dbm_task_t task_next = batch[batch_order[0]];
  for (int itask = 0; itask < ntasks; ++itask) {
    const dbm_task_t task = task_next;
    task_next = batch[batch_order[(itask + 1) < ntasks ? (itask + 1) : itask]];

#if defined(__LIBXS)
    if (0 == (DBM_MULTIPLY_BLAS_LIBRARY & options) &&
        (task.m != kernel_m || task.n != kernel_n || task.k != kernel_k)) {
      const double beta = 1.0;
      gemm_config = libxs_gemm_dispatch(LIBXS_DATATYPE_F64, 'N', 'T', task.m,
                                        task.n, task.k, task.m, task.n, task.m,
                                        &alpha, &beta, NULL);
      kernel_m = task.m;
      kernel_n = task.n;
      kernel_k = task.k;
    }
#endif

    double *const data_a = pack_a->data + task.offset_a;
    double *const data_b = pack_b->data + task.offset_b;
    double *const data_c = shard_c->data + task.offset_c;

#if defined(__LIBXS)
    if (NULL != gemm_config) {
      libxs_gemm_call(gemm_config, data_a, data_b, data_c);
    } else
#endif
    {
      dbm_dgemm('N', 'T', task.m, task.n, task.k, alpha, data_a, task.m, data_b,
                task.n, 1.0, data_c, task.m);
    }
  }
}

// EOF
