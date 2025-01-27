/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2025 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#ifndef DBM_HYPERPARAMS_H
#define DBM_HYPERPARAMS_H

// TODO: Tune parameters, check if dynamic OpenMP scheduling is really faster?

#define HASHTABLE_FACTOR 3.0
#define ALLOCATION_FACTOR 1.5
#define SHARDS_PER_THREAD 1.0
#define MAX_BATCH_SIZE 30000
#define BATCH_NUM_BUCKETS 1000

#endif

// EOF
