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
#include <string.h>

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
  void *mem; // first: allows to cast memchunk into mem-ptr...
  struct dbm_memchunk *next;
  size_t size, used;
} dbm_memchunk_t;

/*******************************************************************************
 * \brief Private linked list of memory chunks that are available.
 * \author Ole Schuett
 ******************************************************************************/
static dbm_memchunk_t *mempool_device_available_head = NULL;
static dbm_memchunk_t *mempool_host_available_head = NULL;

/*******************************************************************************
 * \brief Private linked list of memory chunks that are in use.
 * \author Ole Schuett
 ******************************************************************************/
static dbm_memchunk_t *mempool_device_allocated_head = NULL;
static dbm_memchunk_t *mempool_host_allocated_head = NULL;

/*******************************************************************************
 * \brief Private statistics.
 * \author Hans Pabst
 ******************************************************************************/
static dbm_memstats_t mempool_stats = {0};

/*******************************************************************************
 * \brief Private routine for actually allocating system memory.
 * \author Ole Schuett
 ******************************************************************************/
static void *actual_malloc(size_t size, bool on_device) {
  void *memory = NULL;

  if (0 != size) {
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_DBM)
    if (on_device) {
      offload_activate_chosen_device();
      offloadMalloc(&memory, size);
    } else {
#if (0 != DBM_OFFLOAD_ALLOC)
      offloadMallocHost(&memory, size);
#else
      memory = dbm_mpi_alloc_mem(size);
#endif
    }
#else
    (void)on_device; // mark used
    memory = dbm_mpi_alloc_mem(size);
#endif
    assert(memory != NULL);

    // Update statistics.
    if (on_device) {
#pragma omp atomic
      ++mempool_stats.device_mallocs;
    } else {
#pragma omp atomic
      ++mempool_stats.host_mallocs;
    }
  }

  return memory;
}

/*******************************************************************************
 * \brief Private routine for actually freeing system memory.
 * \author Ole Schuett
 ******************************************************************************/
static void actual_free(const void *memory, bool on_device) {
  if (NULL != memory) {
    void *mem = (void *)(uintptr_t)memory;
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_DBM)
    if (on_device) {
      offload_activate_chosen_device();
      offloadFree(mem);
    } else {
#if (0 != DBM_OFFLOAD_ALLOC)
      offloadFreeHost(mem);
#else
      dbm_mpi_free_mem(mem);
#endif
    }
#else
    (void)on_device; // mark used
    dbm_mpi_free_mem(mem);
#endif
  }
}

/*******************************************************************************
 * \brief Private routine for allocating host or device memory from the pool.
 * \author Ole Schuett
 ******************************************************************************/
#if (0 != DBM_MEMPOOL_DEVICE) || (0 != DBM_MEMPOOL_HOST)
static void *internal_mempool_malloc(dbm_memchunk_t **available_head,
                                     dbm_memchunk_t **allocated_head,
                                     size_t size) {
  if (size == 0) {
    return NULL;
  }

  dbm_memchunk_t *chunk = NULL;
  const bool on_device = (&mempool_device_available_head == available_head);
  assert(on_device || &mempool_host_available_head == available_head);
  assert(on_device || &mempool_host_allocated_head == allocated_head);

#pragma omp critical(dbm_mempool_modify)
  {
    // Find a suitable chunk in mempool_available.
    dbm_memchunk_t **reuse = NULL, **reclaim = NULL;
    for (; NULL != *available_head; available_head = &(*available_head)->next) {
      const size_t s = (*available_head)->size;
      if (size <= s && (NULL == reuse || s < (*reuse)->size)) {
        reuse = available_head;
        if (size == (*reuse)->size) {
          break; // exact match
        }
      } else if (NULL != reclaim) {
        if (s > (*reclaim)->size) {
          reclaim = available_head;
        }
      } else {
        reclaim = available_head;
      }
    }
    if (NULL == reuse) {
      reuse = reclaim;
    }

    // Remove chunk from mempool_available.
    if (NULL != reuse) {
      chunk = *reuse;
      *reuse = chunk->next;
    } else { // Allocate a new chunk.
      chunk = calloc(1, sizeof(dbm_memchunk_t));
      assert(chunk != NULL);
    }

    // Insert chunk into mempool_allocated.
    chunk->next = *allocated_head;
    *allocated_head = chunk;
  }

  // Resize chunk (not in critical section).
  if (chunk->size < size) {
    void *memory = chunk->mem;
    chunk->mem = NULL; // race ok (free and stats)
    actual_free(memory, on_device);
    chunk->mem = actual_malloc(size, on_device);
    chunk->size = size;
  }
  chunk->used = size; // stats

  return chunk->mem;
}
#endif

/*******************************************************************************
 * \brief Internal routine for allocating host memory from the pool.
 * \author Ole Schuett
 ******************************************************************************/
void *dbm_mempool_host_malloc(size_t size) {
#if (0 != DBM_MEMPOOL_HOST)
  return internal_mempool_malloc(&mempool_host_available_head,
                                 &mempool_host_allocated_head, size);
#else
  return actual_malloc(size, false);
#endif
}

/*******************************************************************************
 * \brief Internal routine for allocating device memory from the pool
 * \author Ole Schuett
 ******************************************************************************/
void *dbm_mempool_device_malloc(size_t size) {
#if (0 != DBM_MEMPOOL_DEVICE)
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_DBM)
  return internal_mempool_malloc(&mempool_device_available_head,
                                 &mempool_device_allocated_head, size);
#else
  return dbm_mempool_host_malloc(size);
#endif
#else
  return actual_malloc(size, true);
#endif
}

/*******************************************************************************
 * \brief Private routine for releasing memory back to the pool.
 * \author Ole Schuett
 ******************************************************************************/
#if (0 != DBM_MEMPOOL_DEVICE) || (0 != DBM_MEMPOOL_HOST)
static void internal_mempool_free(dbm_memchunk_t **available_head,
                                  dbm_memchunk_t **allocated_head,
                                  const void *mem) {
  if (NULL != mem) {
#pragma omp critical(dbm_mempool_modify)
    {
      // Find chunk in allocated chunks.
      while (NULL != *allocated_head && (*allocated_head)->mem != mem) {
        allocated_head = &(*allocated_head)->next;
      }
      dbm_memchunk_t *chunk = *allocated_head;
      assert(NULL != chunk && chunk->mem == mem);

      // Remove chunk from mempool_allocated.
      *allocated_head = chunk->next;

      // Add chunk to mempool_available.
      chunk->next = *available_head;
      *available_head = chunk;
    }
  }
}
#endif

/*******************************************************************************
 * \brief Internal routine for releasing memory back to the pool.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mempool_host_free(const void *memory) {
#if (0 != DBM_MEMPOOL_HOST)
  internal_mempool_free(&mempool_host_available_head,
                        &mempool_host_allocated_head, memory);
#else
  actual_free(memory, false);
#endif
}

/*******************************************************************************
 * \brief Internal routine for releasing memory back to the pool.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mempool_device_free(const void *memory) {
#if (0 != DBM_MEMPOOL_DEVICE)
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_DBM)
  internal_mempool_free(&mempool_device_available_head,
                        &mempool_device_allocated_head, memory);
#else
  dbm_mempool_host_free(memory);
#endif
#else
  actual_free(memory, true);
#endif
}

/*******************************************************************************
 * \brief Private routine for freeing all memory in the pool.
 * \author Ole Schuett
 ******************************************************************************/
static void internal_mempool_clear(dbm_memchunk_t **available_head) {
  const bool on_device = (&mempool_device_available_head == available_head);
  assert(on_device || &mempool_host_available_head == available_head);

  // Free chunks in mempool_available.
  while (NULL != *available_head) {
    dbm_memchunk_t *chunk = *available_head;
    *available_head = chunk->next; // remove chunk
    actual_free(chunk->mem, on_device);
    free(chunk);
  }
}

/*******************************************************************************
 * \brief Internal routine for freeing all memory in the pool.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mempool_clear(void) {
  // check for memory leak
  assert(mempool_device_allocated_head == NULL);
  assert(mempool_host_allocated_head == NULL);

#pragma omp critical(dbm_mempool_modify)
  {
    internal_mempool_clear(&mempool_device_available_head);
    internal_mempool_clear(&mempool_host_available_head);
  }
}

/*******************************************************************************
 * \brief Internal routine to query statistics.
 * \author Hans Pabst
 ******************************************************************************/
void dbm_mempool_statistics(dbm_memstats_t *memstats) {
  assert(NULL != memstats);
  memset(memstats, 0, sizeof(dbm_memstats_t));
#pragma omp critical(dbm_mempool_modify)
  {
    dbm_memchunk_t *chunk;
    for (chunk = mempool_device_available_head; NULL != chunk;
         chunk = chunk->next) {
      memstats->device_used += chunk->used;
      memstats->device_size += chunk->size;
    }
    for (chunk = mempool_device_allocated_head; NULL != chunk;
         chunk = chunk->next) {
      memstats->device_used += chunk->used;
      memstats->device_size += chunk->size;
    }
    for (chunk = mempool_host_available_head; NULL != chunk;
         chunk = chunk->next) {
      memstats->host_used += chunk->used;
      memstats->host_size += chunk->size;
    }
    for (chunk = mempool_host_allocated_head; NULL != chunk;
         chunk = chunk->next) {
      memstats->host_used += chunk->used;
      memstats->host_size += chunk->size;
    }
    if (mempool_stats.device_used < memstats->device_used) {
      mempool_stats.device_used = memstats->device_used;
    }
    if (mempool_stats.device_size < memstats->device_size) {
      mempool_stats.device_size = memstats->device_size;
    }
    if (mempool_stats.host_used < memstats->host_used) {
      mempool_stats.host_used = memstats->host_used;
    }
    if (mempool_stats.host_size < memstats->host_size) {
      mempool_stats.host_size = memstats->host_size;
    }
    memcpy(memstats, &mempool_stats, sizeof(dbm_memstats_t));
  }
}

// EOF
