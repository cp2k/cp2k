/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <omp.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "dbm_hyperparams.h"
#include "dbm_shard.h"

/*******************************************************************************
 * \brief Internal routine for initializing a shard.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_shard_init(dbm_shard_t *shard) {
  shard->nblocks = 0;
  shard->nblocks_allocated = INITIAL_NBLOCKS_ALLOCATED;
  shard->blocks = malloc(shard->nblocks_allocated * sizeof(dbm_block_t));

  shard->hashtable_size = HASHTABLE_FACTOR * shard->nblocks_allocated;
  shard->hashtable = calloc(shard->hashtable_size, sizeof(int));

  shard->data_size = 0;
  shard->data_promised = 0;
  shard->data_allocated = INITIAL_DATA_ALLOCATED;
  shard->data = malloc(shard->data_allocated * sizeof(double));

  omp_init_lock(&shard->lock);
}

/*******************************************************************************
 * \brief Internal routine for releasing a shard.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_shard_release(dbm_shard_t *shard) {
  free(shard->blocks);
  free(shard->hashtable);
  free(shard->data);
  omp_destroy_lock(&shard->lock);
}

/*******************************************************************************
 * \brief Private hash function based on Szudzik's elegant pairing.
 *        Using unsigned int to return a positive number even after overflow.
 *        https://en.wikipedia.org/wiki/Pairing_function#Other_pairing_functions
 *        https://stackoverflow.com/a/13871379
 *        http://szudzik.com/ElegantPairing.pdf
 * \author Ole Schuett
 ******************************************************************************/
static inline unsigned int hash(const unsigned int row,
                                const unsigned int col) {
  return (row >= col) ? row * row + row + col : row + col * col;
}

/*******************************************************************************
 * \brief Private routine for inserting a block into a shard's hashtable.
 * \author Ole Schuett
 ******************************************************************************/
static void hashtable_insert(dbm_shard_t *shard, const int block_idx) {
  assert(0 <= block_idx && block_idx < shard->nblocks);
  const dbm_block_t *blk = &shard->blocks[block_idx];
  int slot = hash(blk->row, blk->col) % shard->hashtable_size;
  while (true) {
    if (shard->hashtable[slot] == 0) {
      shard->hashtable[slot] = block_idx + 1; // 1-based because 0 means empty
      return;
    }
    // linear probing
    slot = (slot + 1) % shard->hashtable_size;
  }
}

/*******************************************************************************
 * \brief Internal routine for looking up a block from a shard.
 * \author Ole Schuett
 ******************************************************************************/
dbm_block_t *dbm_shard_lookup(const dbm_shard_t *shard, const int row,
                              const int col) {
  int slot = hash(row, col) % shard->hashtable_size;
  while (true) {
    const int block_idx = shard->hashtable[slot] - 1; // 1-based, 0 means empty.
    if (block_idx < 0) {
      return NULL; // block not found
    }
    assert(0 <= block_idx && block_idx < shard->nblocks);
    dbm_block_t *blk = &shard->blocks[block_idx];
    if (blk->row == row && blk->col == col) {
      return blk;
    }
    // linear probing
    slot = (slot + 1) % shard->hashtable_size;
  }
}

/*******************************************************************************
 * \brief Internal routine for allocating the metadata of a new block.
 * \author Ole Schuett
 ******************************************************************************/
dbm_block_t *dbm_shard_promise_new_block(dbm_shard_t *shard, const int row,
                                         const int col, const int block_size) {
  // Grow blocks array if nessecary.
  if (shard->nblocks_allocated < shard->nblocks + 1) {
    dbm_block_t *old_blocks = shard->blocks;
    shard->nblocks_allocated = ALLOCATION_FACTOR * (shard->nblocks + 1);
    shard->blocks = malloc(shard->nblocks_allocated * sizeof(dbm_block_t));
    memcpy(shard->blocks, old_blocks, shard->nblocks * sizeof(dbm_block_t));
    free(old_blocks);

    // rebuild hashtable
    free(shard->hashtable);
    shard->hashtable_size = HASHTABLE_FACTOR * shard->nblocks_allocated;
    shard->hashtable = calloc(shard->hashtable_size, sizeof(int));
    for (int i = 0; i < shard->nblocks; i++) {
      hashtable_insert(shard, i);
    }
  }

  const int new_block_idx = shard->nblocks;
  shard->nblocks++;
  dbm_block_t *new_block = &shard->blocks[new_block_idx];
  new_block->row = row;
  new_block->col = col;
  new_block->offset = shard->data_promised;
  new_block->norm = 0.0; // initially data will be zeroed
  shard->data_promised += block_size;
  // The data_size will be increase after the memory is allocated and zeroed.
  hashtable_insert(shard, new_block_idx);
  return new_block;
}

/*******************************************************************************
 * \brief Internal routine for allocating and zeroing any promised block's data.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_shard_allocate_promised_blocks(dbm_shard_t *shard) {

  // Reallocate data array if nessecary.
  if (shard->data_promised > shard->data_allocated) {
    double *old_data = shard->data;
    shard->data_allocated = ALLOCATION_FACTOR * shard->data_promised;
    shard->data = malloc(shard->data_allocated * sizeof(double));
    memcpy(shard->data, old_data, shard->data_size * sizeof(double));
    free(old_data);
  }

  // Zero new blocks.
  // The following memset is usually the first touch of the memory, which leads
  // to frequent page faults. The executing thread determines the NUMA location
  if (shard->data_promised > shard->data_size) {
    const int tail = shard->data_promised - shard->data_size;
    memset(&shard->data[shard->data_size], 0, tail * sizeof(double));
    shard->data_size = shard->data_promised;
  }
}

/*******************************************************************************
 * \brief Internal routine for getting block or promising a new one.
 * \author Ole Schuett
 ******************************************************************************/
dbm_block_t *dbm_shard_get_or_promise_block(dbm_shard_t *shard, const int row,
                                            const int col,
                                            const int block_size) {
  dbm_block_t *existing_blk = dbm_shard_lookup(shard, row, col);
  if (existing_blk != NULL) {
    return existing_blk;
  } else {
    return dbm_shard_promise_new_block(shard, row, col, block_size);
  }
}

/*******************************************************************************
 * \brief Internal routine for getting block or allocating a new one.
 * \author Ole Schuett
 ******************************************************************************/
dbm_block_t *dbm_shard_get_or_allocate_block(dbm_shard_t *shard, const int row,
                                             const int col,
                                             const int block_size) {
  dbm_block_t *existing_blk = dbm_shard_lookup(shard, row, col);
  if (existing_blk != NULL) {
    return existing_blk;
  }

  // Create a new block.
  dbm_block_t *new_blk =
      dbm_shard_promise_new_block(shard, row, col, block_size);
  dbm_shard_allocate_promised_blocks(shard);

  return new_blk;
}

// EOF
