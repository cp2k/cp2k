/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#ifndef DBM_SHARD_H
#define DBM_SHARD_H

#include <omp.h>

#ifdef __cplusplus
extern "C" {
#endif

/*******************************************************************************
 * \brief Internal struct for storing a block's metadata.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  int row;
  int col;
  int offset;
} dbm_block_t;

/*******************************************************************************
 * \brief Internal struct for storing a matrix shard.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  int nblocks;
  int nblocks_allocated;
  dbm_block_t *blocks;

  int hashtable_size;  // should be a power of two
  int hashtable_prime; // should be a coprime of hashtable_size
  int *hashtable;      // maps row/col to block numbers

  int data_promised;  // ref'd by dbm_block_t.offset (not yet allocated)
  int data_allocated; // actually allocated (capacity of data buffer)
  int data_size;      // actually allocated and initialized
  double *data;       // potentially over-allocated (to amortize resizing)

  omp_lock_t lock; // used by dbm_put_block
} dbm_shard_t;

/*******************************************************************************
 * \brief Internal routine for initializing a shard.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_shard_init(dbm_shard_t *shard);

/*******************************************************************************
 * \brief Internal routine for copying content of shard_b into shard_a.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_shard_copy(dbm_shard_t *shard_a, const dbm_shard_t *shard_b);

/*******************************************************************************
 * \brief Internal routine for releasing a shard.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_shard_release(dbm_shard_t *shard);

/*******************************************************************************
 * \brief Internal hash function based on Cantor pairing function.
 *        https://en.wikipedia.org/wiki/Pairing_function#Cantor_pairing_function
 *        Szudzik's elegant pairing proved to be too asymmetric wrt. row / col.
 *        Using unsigned int to return a positive number even after overflow.
 * \author Ole Schuett
 ******************************************************************************/
static inline unsigned int dbm_shard_hash(const unsigned int row,
                                          const unsigned int col) {
  return (row + col) * (row + col + 1) / 2 + row; // Division by 2 is cheap.
}

/*******************************************************************************
 * \brief Internal routine for the slot of a block in the shard's hashtable.
 * \author Hans Pabst
 ******************************************************************************/
static inline int dbm_shard_slot(const dbm_shard_t *shard, const int row,
                                 const int col) {
  return (shard->hashtable_prime * dbm_shard_hash(row, col)) &
         (shard->hashtable_size - 1);
}

/*******************************************************************************
 * \brief Internal routine for the slot following the given one (probing).
 * \author Hans Pabst
 ******************************************************************************/
static inline int dbm_shard_nextslot(const dbm_shard_t *shard, const int slot) {
  return (slot + 1) & (shard->hashtable_size - 1);
}

/*******************************************************************************
 * \brief Internal routines for prefetching an upcoming lookup in two stages:
 *        the hashtable slot first, then (once cached) the block it refers to.
 *        A lookup otherwise waits for two dependent cache misses.
 * \author Hans Pabst
 ******************************************************************************/
static inline void dbm_shard_prefetch_slot(const dbm_shard_t *shard,
                                           const int row, const int col) {
  __builtin_prefetch(&shard->hashtable[dbm_shard_slot(shard, row, col)]);
}
static inline void dbm_shard_prefetch_block(const dbm_shard_t *shard,
                                            const int row, const int col) {
  const int block_idx = shard->hashtable[dbm_shard_slot(shard, row, col)];
  if (0 < block_idx) { // first probe only, 1-based, 0 means empty
    __builtin_prefetch(&shard->blocks[block_idx - 1]);
  }
}

/*******************************************************************************
 * \brief Internal routine for looking up a block from a shard.
 * \author Ole Schuett
 ******************************************************************************/
dbm_block_t *dbm_shard_lookup(const dbm_shard_t *shard, const int row,
                              const int col);

/*******************************************************************************
 * \brief Internal routine for growing the block metadata of a shard such that
 *        it holds nblocks blocks in total without growing again.
 * \author Ole Schuett and Hans Pabst
 ******************************************************************************/
void dbm_shard_reserve(dbm_shard_t *shard, int nblocks);

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

#ifdef __cplusplus
}
#endif

#endif

// EOF
