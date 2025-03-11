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
#include <stdio.h>
#include <stdlib.h>

#include "../offload/offload_library.h"
#include "../offload/offload_runtime.h"
#include "dbm_hyperparams.h"
#include "dbm_mempool.h"
#include "dbm_mpi.h"

/*******************************************************************************
 * \brief Private struct for storing a chunk of memory.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct dbm_memchunk {
  void *mem; // first: allows to cast memchunk into memptr...
  struct dbm_memchunk *next;
  size_t size;
} dbm_memchunk_t;

/*******************************************************************************
 * \brief Private linked list of memory chunks that are available.
 * \author Ole Schuett
 ******************************************************************************/
static dbm_memchunk_t *mempool_available_head = NULL;

/*******************************************************************************
 * \brief Private linked list of memory chunks that are in use.
 * \author Ole Schuett
 ******************************************************************************/
static dbm_memchunk_t *mempool_allocated_head = NULL;

/*******************************************************************************
 * \brief Private statistics (survives dbm_mempool_clear).
 * \author Hans Pabst
 ******************************************************************************/
static dbm_memstats_t mempool_stats = {0};

/*******************************************************************************
 * \brief Internal routine for allocating device memory from the pool
 * \author Ole Schuett
 ******************************************************************************/
void *dbm_mempool_device_malloc(const size_t size) {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_DBM)
  if (size == 0) {
    return NULL;
  }

  dbm_memchunk_t *chunk = NULL;
#pragma omp critical(dbm_mempool_modify)
  {
    // Find a suitable chunk in mempool_available.
    dbm_memchunk_t **indirect = &mempool_available_head;
    dbm_memchunk_t **hit = NULL, **fallback = NULL;
    for (; NULL != *indirect; indirect = &(*indirect)->next) {
      if ((*indirect)->size < size) {
        if (NULL == fallback || (*fallback)->size < (*indirect)->size) {
          fallback = indirect;
        }
      } else if (NULL == hit || (*indirect)->size < (*hit)->size) {
        hit = indirect;
        if (size == (*hit)->size) {
          break;
        }
      }
    }
    if (NULL == hit) {
      hit = fallback;
    }

    // If a chunck was found, remove it from mempool_available.
    if (NULL != hit) {
      chunk = *hit;
      *hit = chunk->next;
    } else { // Allocate a new chunk.
      assert(chunk == NULL);
      chunk = malloc(sizeof(dbm_memchunk_t));
      assert(chunk != NULL);
      chunk->size = 0;
      chunk->mem = NULL;
    }

    // Insert chunk into mempool_allocated.
    chunk->next = mempool_allocated_head;
    mempool_allocated_head = chunk;

    // Update statistics
    if (chunk->size < size) {
      mempool_stats.size += size - chunk->size;
      ++mempool_stats.nmallocs;
    }
  }

  // Resize chunk if needed (outside of critical section).
  if (chunk->size < size) {
    offload_activate_chosen_device();
    offloadFree(chunk->mem);
    offloadMalloc(&chunk->mem, size);
    assert(chunk->mem != NULL);
    chunk->size = size; // update
  }

  return chunk->mem;
#else
  (void)size; // mark used
  return NULL;
#endif
}

/*******************************************************************************
 * \brief Internal routine for releasing memory back to the pool.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mempool_free(void *mem) {
  if (mem == NULL) {
    return;
  }

#pragma omp critical(dbm_mempool_modify)
  {
    // Find chunk in mempool_allocated.
    dbm_memchunk_t **indirect = &mempool_allocated_head;
    while (*indirect != NULL && (*indirect)->mem != mem) {
      indirect = &(*indirect)->next;
    }
    dbm_memchunk_t *chunk = *indirect;
    assert(chunk != NULL && chunk->mem == mem);

    // Remove chunk from mempool_allocated.
    *indirect = chunk->next;

    // Add chunk to mempool_available.
    chunk->next = mempool_available_head;
    mempool_available_head = chunk;
  }
}

/*******************************************************************************
 * \brief Internal routine for freeing all memory in the pool (not thread-safe).
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mempool_clear(void) {
  assert(omp_get_num_threads() == 1);
  assert(mempool_allocated_head == NULL); // check for memory leak
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_DBM)
  offload_activate_chosen_device();
#if 0
  while (mempool_allocated_head != NULL) {
    dbm_memchunk_t *chunk = mempool_allocated_head;
    mempool_allocated_head = chunk->next;
    printf("Found alloacted memory chunk of size: %lu\n", chunk->size);
    offloadFree(chunk->mem);
    free(chunk);
  }
#endif
  // Free chunks in mempool_available.
  while (mempool_available_head != NULL) {
    dbm_memchunk_t *chunk = mempool_available_head;
    mempool_available_head = chunk->next;
    offloadFree(chunk->mem);
    free(chunk);
  }
#endif
}

/*******************************************************************************
 * \brief Internal routine to query statistics (not thread-safe).
 * \author Hans Pabst
 ******************************************************************************/
void dbm_mempool_statistics(dbm_memstats_t *memstats) {
  assert(NULL != memstats);
  *memstats = mempool_stats;
}

// EOF
