/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <stddef.h>
#include <string.h>

#if defined(__LIBXSMM)
#include <libxsmm.h>
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
static inline void dbm_dgemm(const bool transa, const bool transb, const int m,
                             const int n, const int k, const double alpha,
                             const double *a, const int lda, const double *b,
                             const int ldb, const double beta, double *c,
                             const int ldc) {
  const char transa_char = (transa) ? 'T' : 'N';
  const char transb_char = (transb) ? 'T' : 'N';

#if defined(__LIBXSMM)
  libxsmm_dgemm(&transa_char, &transb_char, &m, &n, &k, &alpha, a, &lda, b,
                &ldb, &beta, c, &ldc);
#else
  dgemm_(&transa_char, &transb_char, &m, &n, &k, &alpha, a, &lda, b, &ldb,
         &beta, c, &ldc);
#endif
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
                                    const bool transa, const bool transb,
                                    const double alpha,
                                    const dbm_pack_t *pack_a,
                                    const dbm_pack_t *pack_b,
                                    dbm_shard_t *shard_c) {

  dbm_shard_allocate_promised_blocks(shard_c);

#if defined(__LIBXSMM)

  // Sort tasks approximately by m,n,k via bucket sort.
  int buckets[BATCH_NUM_BUCKETS];
  memset(buckets, 0, BATCH_NUM_BUCKETS * sizeof(int));
  for (int itask = 0; itask < ntasks; itask++) {
    const int i = hash(batch[itask]) % BATCH_NUM_BUCKETS;
    buckets[i]++;
  }
  for (int i = 1; i < BATCH_NUM_BUCKETS; i++) {
    buckets[i] += buckets[i - 1];
  }
  assert(buckets[BATCH_NUM_BUCKETS - 1] == ntasks);
  int batch_order[ntasks];
  for (int itask = 0; itask < ntasks; itask++) {
    const int i = hash(batch[itask]) % BATCH_NUM_BUCKETS;
    buckets[i]--;
    batch_order[buckets[i]] = itask;
  }

  // Prepare arguments for libxsmm_dmmdispatch().
  const double beta = 1.0;
  int flags = 0;
  flags |= (transa) ? LIBXSMM_GEMM_FLAG_TRANS_A : 0;
  flags |= (transb) ? LIBXSMM_GEMM_FLAG_TRANS_B : 0;
  const int prefetch = LIBXSMM_PREFETCH_NONE; // somehow prefetching is slower

  // Loop over tasks.
  libxsmm_dmmfunction kernel_func = NULL;
  int kernel_m = 0, kernel_n = 0, kernel_k = 0;
  for (int itask = 0; itask < ntasks; itask++) {
    const dbm_task_t task = batch[batch_order[itask]];

    if (task.m != kernel_m || task.n != kernel_n || task.k != kernel_k) {
      kernel_func = libxsmm_dmmdispatch(task.m, task.n, task.k, NULL /*lda*/,
                                        NULL /*ldb*/, NULL /*ldc*/, &alpha,
                                        &beta, &flags, &prefetch);
      kernel_m = task.m;
      kernel_n = task.n;
      kernel_k = task.k;
    }

    const double *data_a = &pack_a->data[task.offset_a];
    const double *data_b = &pack_b->data[task.offset_b];
    double *data_c = &shard_c->data[task.offset_c];

    if (kernel_func != NULL) {
      kernel_func(data_a, data_b, data_c);
    } else {
      const int lda = (transa) ? task.k : task.m;
      const int ldb = (transb) ? task.n : task.k;
      const int ldc = task.m;
      dbm_dgemm(transa, transb, task.m, task.n, task.k, alpha, data_a, lda,
                data_b, ldb, 1.0, data_c, ldc);
    }
  }

#else

  // Fallback to BLAS when libxsmm is not available.
  for (int itask = 0; itask < ntasks; itask++) {
    const dbm_task_t task = batch[itask];
    const int lda = (transa) ? task.k : task.m;
    const int ldb = (transb) ? task.n : task.k;
    const int ldc = task.m;
    const double *data_a = &pack_a->data[task.offset_a];
    const double *data_b = &pack_b->data[task.offset_b];
    double *data_c = &shard_c->data[task.offset_c];
    dbm_dgemm(transa, transb, task.m, task.n, task.k, alpha, data_a, lda,
              data_b, ldb, 1.0, data_c, ldc);
  }

#endif
}

// EOF
