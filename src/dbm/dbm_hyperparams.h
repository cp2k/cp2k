/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2025 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#ifndef DBM_HYPERPARAMS_H
#define DBM_HYPERPARAMS_H

#define DBM_OMP_SCHEDULE schedule(dynamic, 1)

#define DBM_MEMPOOL_DEVICE 1
#define DBM_MEMPOOL_HOST 1

#define DBM_ALLOC_OFFLOAD 1
#define DBM_ALLOC_OPENMP 1
#define DBM_ALLOC_MPI 0

#define DBM_OVERCOMMIT_DEVICE 1.3
#define DBM_OVERCOMMIT_HOST 1.4
#define DBM_HASHTABLE_FACTOR 3.0

#define DBM_SHARDS_PER_THREAD 1.0
#define DBM_MAX_BATCH_SIZE 30000
#define DBM_BATCH_NUM_BUCKETS 1000

#endif

// EOF
