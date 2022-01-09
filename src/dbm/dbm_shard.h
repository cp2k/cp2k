/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#ifndef DBM_SHARD_H
#define DBM_SHARD_H

#include <omp.h>

/*******************************************************************************
 * \brief Internal struct for storing a block's metadata.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  int row;
  int col;
  int offset;
  float norm;
} dbm_block_t;

/*******************************************************************************
 * \brief Internal struct for storing a matrix shard.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  int nblocks;
  int nblocks_allocated;
  dbm_block_t *blocks;

  int hashtable_size;
  int *hashtable; // maps row/col to block numbers

  int data_promised;  // referenced by a dbm_block_t.offset, but not yet
                      // allocated
  int data_allocated; // we over allocate to amortized the resizing cost
  int data_size;      // actually allocated and initialized
  double *data;

  omp_lock_t lock; // used by dbm_put_block
} dbm_shard_t;

/*******************************************************************************
 * \brief Internal routine for initializing a shard.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_shard_init(dbm_shard_t *shard);

/*******************************************************************************
 * \brief Internal routine for releasing a shard.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_shard_release(dbm_shard_t *shard);

/*******************************************************************************
 * \brief Internal routine for looking up a block from a shard.
 * \author Ole Schuett
 ******************************************************************************/
dbm_block_t *dbm_shard_lookup(const dbm_shard_t *shard, const int row,
                              const int col);

/*******************************************************************************
 * \brief Internal routine for allocating the metadata of a new block.
 * \author Ole Schuett
 ******************************************************************************/
dbm_block_t *dbm_shard_promise_new_block(dbm_shard_t *shard, const int row,
                                         const int col, const int block_size);

/*******************************************************************************
 * \brief Internal routine for allocating and zeroing any promised block's data.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_shard_allocate_promised_blocks(dbm_shard_t *shard);

/*******************************************************************************
 * \brief Internal routine for getting block or promising a new one.
 * \author Ole Schuett
 ******************************************************************************/
dbm_block_t *dbm_shard_get_or_promise_block(dbm_shard_t *shard, const int row,
                                            const int col,
                                            const int block_size);

/*******************************************************************************
 * \brief Internal routine for getting block or allocating a new one.
 * \author Ole Schuett
 ******************************************************************************/
dbm_block_t *dbm_shard_get_or_allocate_block(dbm_shard_t *shard, const int row,
                                             const int col,
                                             const int block_size);

#endif

// EOF
