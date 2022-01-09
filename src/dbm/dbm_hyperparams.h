/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#ifndef DBM_HYPERPARAMS_H

// TODO: Tune parameters, check if dynamic OpenMP scheduling is really faster?

static const float HASHTABLE_FACTOR = 3.0;
static const float ALLOCATION_FACTOR = 1.5;
static const float SHARDS_PER_THREAD = 1.0;
static const int MAX_BATCH_SIZE = 10000;
static const int BATCH_NUM_BUCKETS = 1000;
static const int INITIAL_NBLOCKS_ALLOCATED = 100;
static const int INITIAL_DATA_ALLOCATED = 1024;

#endif

// EOF
