/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2024 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#ifndef DBM_MULTIPLY_INTERNAL_H
#define DBM_MULTIPLY_INTERNAL_H

/*******************************************************************************
 * \brief Internal struct for storing a dbm_block_t plus its norm.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  int free_index; // Free index in Einstein notation of matrix multiplication.
  int sum_index;  // Summation index - also called dummy index.
  int offset;
  float norm;
} dbm_pack_block_t;

/*******************************************************************************
 * \brief Internal struct for storing a pack - essentially a shard for MPI.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  int nblocks;
  int data_size;
  dbm_pack_block_t *blocks;
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
