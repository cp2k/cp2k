/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2025 CP2K developers group <https://cp2k.org>              */
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
 * \brief Internal routine for finding a power of two greater than given number.
 * \author Ole Schuett
 ******************************************************************************/
static int next_power2(const int start) {
  int candidate = 2;
  while (candidate < start) {
    candidate *= 2;
  }
  return candidate;
}

/*******************************************************************************
 * \brief Internal routine for finding a prime greater equal than given number.
 * \author Ole Schuett
 ******************************************************************************/
static int next_prime(const int start) {
  int candidate = start, divisor = 0;
  while (divisor < candidate) {
    for (divisor = 2; divisor < candidate; divisor++) {
      if (candidate % divisor == 0) {
        candidate++;
        break;
      }
    }
  }
  return candidate;
}

/*******************************************************************************
 * \brief Internal routine for initializing a shard's hashtable.
 * \author Ole Schuett
 ******************************************************************************/
static void hashtable_init(dbm_shard_t *shard) {
  // Choosing size as power of two allows to replace modulo with bitwise AND.
  shard->hashtable_size =
      next_power2(HASHTABLE_FACTOR * shard->nblocks_allocated);
  shard->hashtable_mask = shard->hashtable_size - 1;
  shard->hashtable_prime = next_prime(shard->hashtable_size);
  shard->hashtable = calloc(shard->hashtable_size, sizeof(int));
  assert(shard->hashtable != NULL);
}

/*******************************************************************************
 * \brief Internal routine for initializing a shard.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_shard_init(dbm_shard_t *shard) {
  shard->nblocks = 0;
  shard->nblocks_allocated = 0;
  shard->blocks = NULL;
  hashtable_init(shard);
  shard->data_size = 0;
  shard->data_promised = 0;
  shard->data_allocated = 0;
  shard->data = NULL;
  omp_init_lock(&shard->lock);
}

/*******************************************************************************
 * \brief Internal routine for copying content of shard_b into shard_a.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_shard_copy(dbm_shard_t *shard_a, const dbm_shard_t *shard_b) {
  assert(shard_a != NULL && shard_b != NULL);

  if (shard_a->nblocks_allocated < shard_b->nblocks) {
    free(shard_a->blocks);
    shard_a->blocks = malloc(shard_b->nblocks * sizeof(dbm_block_t));
    shard_a->nblocks_allocated = shard_b->nblocks;
    assert(shard_a->blocks != NULL);
  }
  shard_a->nblocks = shard_b->nblocks;

  if (shard_a->hashtable_size < shard_b->hashtable_size) {
    free(shard_a->hashtable);
    shard_a->hashtable = malloc(shard_b->hashtable_size * sizeof(int));
    assert(shard_a->hashtable != NULL);
  }
  shard_a->hashtable_size = shard_b->hashtable_size;
  shard_a->hashtable_mask = shard_b->hashtable_mask;
  shard_a->hashtable_prime = shard_b->hashtable_prime;

  if (shard_a->data_allocated < shard_b->data_size) {
    free(shard_a->data);
    shard_a->data = malloc(shard_b->data_size * sizeof(double));
    shard_a->data_allocated = shard_b->data_size;
    assert(shard_a->data != NULL);
  }
  shard_a->data_size = shard_b->data_size;

  if (shard_a->blocks != NULL) {
    memcpy(shard_a->blocks, shard_b->blocks,
           shard_b->nblocks * sizeof(dbm_block_t));
  }
  if (shard_a->hashtable != NULL) {
    memcpy(shard_a->hashtable, shard_b->hashtable,
           shard_b->hashtable_size * sizeof(int));
  }
  if (shard_a->data != NULL) {
    memcpy(shard_a->data, shard_b->data, shard_b->data_size * sizeof(double));
  }
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
 * \brief Private hash function based on Cantor pairing function.
 *        https://en.wikipedia.org/wiki/Pairing_function#Cantor_pairing_function
 *        Szudzik's elegant pairing proved to be too asymmetric wrt. row / col.
 *        Using unsigned int to return a positive number even after overflow.
 * \author Ole Schuett
 ******************************************************************************/
static inline unsigned int hash(const unsigned int row,
                                const unsigned int col) {
  return (row + col) * (row + col + 1) / 2 + row; // Division by 2 is cheap.
}

/*******************************************************************************
 * \brief Private routine for inserting a block into a shard's hashtable.
 * \author Ole Schuett
 ******************************************************************************/
static void hashtable_insert(dbm_shard_t *shard, const int block_idx) {
  assert(0 <= block_idx && block_idx < shard->nblocks);
  const dbm_block_t *blk = &shard->blocks[block_idx];
  const int row = blk->row, col = blk->col;
  int slot = (shard->hashtable_prime * hash(row, col)) & shard->hashtable_mask;
  while (true) {
    if (shard->hashtable[slot] == 0) {
      shard->hashtable[slot] = block_idx + 1; // 1-based because 0 means empty
      return;
    }
    // linear probing
    slot = (slot + 1) & shard->hashtable_mask;
  }
}

/*******************************************************************************
 * \brief Internal routine for looking up a block from a shard.
 * \author Ole Schuett
 ******************************************************************************/
dbm_block_t *dbm_shard_lookup(const dbm_shard_t *shard, const int row,
                              const int col) {
  int slot = (shard->hashtable_prime * hash(row, col)) & shard->hashtable_mask;
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
    slot = (slot + 1) & shard->hashtable_mask;
  }
}

/*******************************************************************************
 * \brief Internal routine for allocating the metadata of a new block.
 * \author Ole Schuett
 ******************************************************************************/
dbm_block_t *dbm_shard_promise_new_block(dbm_shard_t *shard, const int row,
                                         const int col, const int block_size) {
  // Grow blocks array if necessary.
  if (shard->nblocks_allocated < shard->nblocks + 1) {
    shard->nblocks_allocated = ALLOCATION_FACTOR * (shard->nblocks + 1);
    shard->blocks =
        realloc(shard->blocks, shard->nblocks_allocated * sizeof(dbm_block_t));
    assert(shard->blocks != NULL);

    // rebuild hashtable
    free(shard->hashtable);
    hashtable_init(shard);
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

  // Reallocate data array if necessary.
  if (shard->data_promised > shard->data_allocated) {
    shard->data_allocated = ALLOCATION_FACTOR * shard->data_promised;
    shard->data = realloc(shard->data, shard->data_allocated * sizeof(double));
    assert(shard->data != NULL);
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
