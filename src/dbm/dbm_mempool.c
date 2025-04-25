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

#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_DBM)
#define DBM_MEMPOOL_OFFLOAD_ENABLED 1
#else
#define DBM_MEMPOOL_OFFLOAD_ENABLED 0
#endif
#define DBM_MEMPOOL_DEVICE_ENABLED                                             \
  (DBM_MEMPOOL_DEVICE && DBM_MEMPOOL_OFFLOAD_ENABLED)
#define DBM_MEMPOOL_HOST_ENABLED                                               \
  ((DBM_MEMPOOL_HOST && DBM_ALLOC_OFFLOAD && DBM_MEMPOOL_OFFLOAD_ENABLED) ||   \
   (1 < DBM_MEMPOOL_HOST))

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
 * \brief Private single-linked lists of memory chunks available and allocated.
 * \author Ole Schuett
 ******************************************************************************/
#if DBM_MEMPOOL_DEVICE_ENABLED
static dbm_memchunk_t *mempool_device_available_head = NULL;
static dbm_memchunk_t *mempool_device_allocated_head = NULL;
#endif

/*******************************************************************************
 * \brief Private single-linked lists of memory chunks available and allocated.
 * \author Ole Schuett
 ******************************************************************************/
#if DBM_MEMPOOL_HOST_ENABLED
static dbm_memchunk_t *mempool_host_available_head = NULL;
static dbm_memchunk_t *mempool_host_allocated_head = NULL;
#endif

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
#if DBM_MEMPOOL_OFFLOAD_ENABLED
    if (on_device) {
      offload_activate_chosen_device();
      offloadMalloc(&memory, size);
    } else {
#if DBM_ALLOC_OFFLOAD
      offload_activate_chosen_device();
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
#if DBM_MEMPOOL_OFFLOAD_ENABLED
    if (on_device) {
      offload_activate_chosen_device();
      offloadFree(mem);
    } else {
#if DBM_ALLOC_OFFLOAD
      offload_activate_chosen_device();
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
#if DBM_MEMPOOL_DEVICE_ENABLED || DBM_MEMPOOL_HOST_ENABLED
static void *internal_mempool_malloc(dbm_memchunk_t **available_head,
                                     dbm_memchunk_t **allocated_head,
                                     size_t size) {
  if (size == 0) {
    return NULL;
  }

  dbm_memchunk_t *chunk = NULL;
#if DBM_MEMPOOL_DEVICE_ENABLED
  const bool on_device = (&mempool_device_available_head == available_head);
#else
  const bool on_device = false;
#endif
#if DBM_MEMPOOL_HOST_ENABLED
  assert(on_device || &mempool_host_available_head == available_head);
  assert(on_device || &mempool_host_allocated_head == allocated_head);
#endif

#pragma omp critical(dbm_mempool_modify)
  {
    // Find a suitable chunk in mempool_available.
    dbm_memchunk_t **reuse = NULL, **reclaim = NULL;
    for (; NULL != *available_head; available_head = &(*available_head)->next) {
      const size_t s = (*available_head)->size;
      if (size <= s && (NULL == reuse || s < (*reuse)->size)) {
        reuse = available_head;
        if (size == (*reuse)->size) {
          break; // early exit
        }
      } else if (NULL == reclaim || s > (*reclaim)->size) {
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
#if DBM_MEMPOOL_HOST_ENABLED
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
#if DBM_MEMPOOL_DEVICE_ENABLED
  return internal_mempool_malloc(&mempool_device_available_head,
                                 &mempool_device_allocated_head, size);
#elif DBM_MEMPOOL_DEVICE
  return dbm_mempool_host_malloc(size);
#else
  return actual_malloc(size, true);
#endif
}

/*******************************************************************************
 * \brief Private routine for releasing memory back to the pool.
 * \author Ole Schuett
 ******************************************************************************/
#if DBM_MEMPOOL_DEVICE_ENABLED || DBM_MEMPOOL_HOST_ENABLED
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
#if DBM_MEMPOOL_HOST_ENABLED
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
#if DBM_MEMPOOL_DEVICE_ENABLED
  internal_mempool_free(&mempool_device_available_head,
                        &mempool_device_allocated_head, memory);
#elif DBM_MEMPOOL_DEVICE
  dbm_mempool_host_free(memory);
#else
  actual_free(memory, true);
#endif
}

/*******************************************************************************
 * \brief Private routine for freeing all memory in the pool.
 * \author Ole Schuett
 ******************************************************************************/
#if DBM_MEMPOOL_DEVICE_ENABLED || DBM_MEMPOOL_HOST_ENABLED
static void internal_mempool_clear(dbm_memchunk_t **available_head) {
#if DBM_MEMPOOL_DEVICE_ENABLED
  const bool on_device = (&mempool_device_available_head == available_head);
#else
  const bool on_device = false;
#endif
#if DBM_MEMPOOL_HOST_ENABLED
  assert(on_device || &mempool_host_available_head == available_head);
#endif

  // Free chunks in mempool_available.
  while (NULL != *available_head) {
    dbm_memchunk_t *chunk = *available_head;
    *available_head = chunk->next; // remove chunk
    actual_free(chunk->mem, on_device);
    free(chunk);
  }
}
#endif

/*******************************************************************************
 * \brief Internal routine for freeing all memory in the pool.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mempool_clear(void) {
#pragma omp critical(dbm_mempool_modify)
  {
#if DBM_MEMPOOL_DEVICE_ENABLED
    assert(mempool_device_allocated_head == NULL); // check for leak
    internal_mempool_clear(&mempool_device_available_head);
#endif
#if DBM_MEMPOOL_HOST_ENABLED
    assert(mempool_host_allocated_head == NULL); // check for leak
    internal_mempool_clear(&mempool_host_available_head);
#endif
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
#if DBM_MEMPOOL_DEVICE_ENABLED
    for (dbm_memchunk_t *chunk = mempool_device_available_head; NULL != chunk;
         chunk = chunk->next) {
      memstats->device_used += chunk->used;
      memstats->device_size += chunk->size;
    }
    for (dbm_memchunk_t *chunk = mempool_device_allocated_head; NULL != chunk;
         chunk = chunk->next) {
      memstats->device_used += chunk->used;
      memstats->device_size += chunk->size;
    }
    if (mempool_stats.device_used < memstats->device_used) {
      mempool_stats.device_used = memstats->device_used;
    }
    if (mempool_stats.device_size < memstats->device_size) {
      mempool_stats.device_size = memstats->device_size;
    }
#endif
#if DBM_MEMPOOL_HOST_ENABLED
    for (dbm_memchunk_t *chunk = mempool_host_available_head; NULL != chunk;
         chunk = chunk->next) {
      memstats->host_used += chunk->used;
      memstats->host_size += chunk->size;
    }
    for (dbm_memchunk_t *chunk = mempool_host_allocated_head; NULL != chunk;
         chunk = chunk->next) {
      memstats->host_used += chunk->used;
      memstats->host_size += chunk->size;
    }
    if (mempool_stats.host_used < memstats->host_used) {
      mempool_stats.host_used = memstats->host_used;
    }
    if (mempool_stats.host_size < memstats->host_size) {
      mempool_stats.host_size = memstats->host_size;
    }
#endif
    memcpy(memstats, &mempool_stats, sizeof(dbm_memstats_t));
  }
}

// EOF
