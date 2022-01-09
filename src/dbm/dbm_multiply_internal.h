/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#ifndef DBM_MULTIPLY_INTERNAL_H
#define DBM_MULTIPLY_INTERNAL_H

#include "dbm_shard.h"

/*******************************************************************************
 * \brief Internal struct for storing a pack - essentially a shard for MPI.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  int nblocks;
  int data_size;
  dbm_block_t *blocks;
  double *data;
} dbm_pack_t;

/*******************************************************************************
 * \brief Internal struct for storing a task, ie. a single block multiplication.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  int m;
  int n;
  int k;
  int offset_a;
  int offset_b;
  int offset_c;
} dbm_task_t;

#endif

// EOF
