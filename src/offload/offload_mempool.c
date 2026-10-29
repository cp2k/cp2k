/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#include "offload_mempool.h"
#include "../mpiwrap/cp_mpi.h"
#include "offload_library.h"
#include "offload_runtime.h"

#include <assert.h>
#include <inttypes.h>
#include <limits.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(__parallel)
#include <mpi.h>
#endif

#define OFFLOAD_MEMPOOL_PRINT(FN, MSG, OUTPUT_UNIT)                            \
  ((FN)(MSG, (int)strlen(MSG), OUTPUT_UNIT))
#define OFFLOAD_MEMPOOL_OMPALLOC 1

enum { MEMPOOL_NBINS = 4 * sizeof(size_t) * CHAR_BIT };
#define MEMPOOL_MIN_CACHE (UINT64_C(64) * 1024 * 1024)

/*******************************************************************************
 * \brief Private struct for storing a chunk of memory.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct offload_memchunk {
  void *mem; // first: allows to cast memchunk into mem-ptr...
  struct offload_memchunk *next, **prev;  // hash chain or size-bin links
  struct offload_memchunk *older, *newer; // cached chunks, oldest first
  size_t size, used;
} offload_memchunk_t;

/*******************************************************************************
 * \brief Private struct for storing a memory pool.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct offload_mempool {
  offload_memchunk_t *available[MEMPOOL_NBINS], *oldest, *newest;
  offload_memchunk_t **allocated; // pointer hash table
  size_t nbuckets, nallocated;
  uint64_t used, live_size, cached_size, peak_live_size, peak_size;
} offload_mempool_t;

/*******************************************************************************
 * \brief Private pools for host and device memory.
 * \author Ole Schuett
 ******************************************************************************/
static offload_mempool_t mempool_host = {0}, mempool_device = {0};

/*******************************************************************************
 * \brief Private counters for statistics.
 * \author Hans Pabst
 ******************************************************************************/
static uint64_t host_malloc_counter = 0, device_malloc_counter = 0;

/*******************************************************************************
 * \brief Returns the larger of two given integer (missing from the C standard)
 * \author Ole Schuett
 ******************************************************************************/
static inline uint64_t imax(uint64_t x, uint64_t y) { return (x > y ? x : y); }

/*******************************************************************************
 * \brief Private routine for actually allocating system memory.
 * \author Ole Schuett
 ******************************************************************************/
static void *actual_malloc(const size_t size, const bool on_device) {
  if (size == 0) {
    return NULL;
  }

  void *memory = NULL;

// This pool must not sit on top of LIBXSTREAM's pool (double caching).
#if defined(__OFFLOAD) && defined(__LIBXSTREAM)
  offload_activate_chosen_device();
  if (on_device) {
    OFFLOAD_CHECK(libxstream_mem_dev_allocate_hint(&memory, size,
                                                   LIBXSTREAM_MALLOC_CP2K));
  } else {
    OFFLOAD_CHECK(libxstream_mem_host_allocate_hint(&memory, size,
                                                    LIBXSTREAM_MALLOC_CP2K));
  }
#elif defined(__OFFLOAD)
  if (on_device) {
    offload_activate_chosen_device();
    offloadMalloc(&memory, size);
  } else {
    offload_activate_chosen_device();
    offloadMallocHost(&memory, size);
  }
#elif OFFLOAD_MEMPOOL_OMPALLOC && (201811 /*v5.0*/ <= _OPENMP)
  memory = omp_alloc(size, omp_null_allocator);
#elif defined(__parallel) && !OFFLOAD_MEMPOOL_OMPALLOC
  if (MPI_SUCCESS != MPI_Alloc_mem((MPI_Aint)size, MPI_INFO_NULL, &memory)) {
    fprintf(stderr, "ERROR: MPI_Alloc_mem failed at %s:%i\n", name, __FILE__,
            __LINE__);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
#else
  memory = malloc(size);
#endif

  // Update statistics.
  if (on_device) {
#pragma omp atomic
    ++device_malloc_counter;
  } else {
#pragma omp atomic
    ++host_malloc_counter;
  }

  assert(memory != NULL);
  return memory;
}

/*******************************************************************************
 * \brief Private routine for actually freeing system memory.
 * \author Ole Schuett
 ******************************************************************************/
static void actual_free(void *memory, const bool on_device) {
  if (NULL == memory) {
    return;
  }

// This pool must not sit on top of LIBXSTREAM's pool (double caching).
#if defined(__OFFLOAD) && defined(__LIBXSTREAM)
  offload_activate_chosen_device();
  if (on_device) {
    OFFLOAD_CHECK(
        libxstream_mem_dev_deallocate_hint(memory, LIBXSTREAM_MALLOC_CP2K));
  } else {
    OFFLOAD_CHECK(
        libxstream_mem_host_deallocate_hint(memory, LIBXSTREAM_MALLOC_CP2K));
  }
#elif defined(__OFFLOAD)
  if (on_device) {
    offload_activate_chosen_device();
    offloadFree(memory);
  } else {
    offload_activate_chosen_device();
    offloadFreeHost(memory);
  }
#elif OFFLOAD_MEMPOOL_OMPALLOC && (201811 /*v5.0*/ <= _OPENMP)
  (void)on_device; // mark used
  omp_free(memory, omp_null_allocator);
#elif defined(__parallel) && !OFFLOAD_MEMPOOL_OMPALLOC
  (void)on_device; // mark used
  if (MPI_SUCCESS != MPI_Free_mem(memory)) {
    fprintf(stderr, "ERROR: MPI_Free_mem failed at %s:%i\n", name, __FILE__,
            __LINE__);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
#else
  (void)on_device; // mark used
  free(memory);
#endif
}

/*******************************************************************************
 * \brief Four size classes per power of two, with less than 25% rounding slack.
 ******************************************************************************/
static int mempool_bin(const size_t size, size_t *capacity) {
  size_t n = size - 1;
  int shift = 0;
  while (n > 7) {
    n >>= 1;
    ++shift;
  }
  // Sizes 1..8 are exact. Larger ranges have four evenly spaced capacities.
  // Requests that cannot be rounded without overflow bypass the cache.
  const int bin =
      n + 1 > (SIZE_MAX >> shift) ? MEMPOOL_NBINS : 4 * shift + (int)n;
  if (capacity != NULL) {
    *capacity = bin < MEMPOOL_NBINS ? (n + 1) << shift : size;
  }
  return bin;
}

/*******************************************************************************
 * \brief Mix aligned host/device addresses before masking the hash bucket.
 ******************************************************************************/
static size_t mempool_hash(const void *mem, const size_t nbuckets) {
  uint64_t key = (uintptr_t)mem;
  key ^= key >> 32;
  key *= UINT64_C(11400714819323198485);
  key ^= key >> 32;
  return key & (nbuckets - 1);
}

/*******************************************************************************
 * \brief Register a borrowed chunk. All metadata helpers require the pool lock.
 ******************************************************************************/
static void mempool_register(offload_mempool_t *pool, offload_memchunk_t *chunk,
                             const size_t used) {
  if (pool->nallocated == pool->nbuckets) {
    const size_t n = pool->nbuckets ? 2 * pool->nbuckets : 256;
    offload_memchunk_t **buckets = calloc(n, sizeof(*buckets));
    assert(buckets != NULL);
    for (size_t i = 0; i < pool->nbuckets; i++) {
      offload_memchunk_t *c = pool->allocated[i];
      while (c != NULL) {
        offload_memchunk_t *next = c->next;
        const size_t h = mempool_hash(c->mem, n);
        c->next = buckets[h];
        buckets[h] = c;
        c = next;
      }
    }
    free(pool->allocated);
    pool->allocated = buckets;
    pool->nbuckets = n;
  }
  const size_t h = mempool_hash(chunk->mem, pool->nbuckets);
  chunk->next = pool->allocated[h];
  pool->allocated[h] = chunk;
  pool->nallocated++;
  chunk->used = used;
  pool->used += used;
  pool->live_size += chunk->size;
  pool->peak_live_size = imax(pool->peak_live_size, pool->live_size);
  pool->peak_size = imax(pool->peak_size, pool->live_size + pool->cached_size);
}

/*******************************************************************************
 * \brief Remove a cached chunk from both its size bin and the eviction queue.
 ******************************************************************************/
static void mempool_uncache(offload_mempool_t *pool,
                            offload_memchunk_t *chunk) {
  *chunk->prev = chunk->next;
  if (chunk->next != NULL) {
    chunk->next->prev = chunk->prev;
  }
  if (chunk->older != NULL) {
    chunk->older->newer = chunk->newer;
  } else {
    pool->oldest = chunk->newer;
  }
  if (chunk->newer != NULL) {
    chunk->newer->older = chunk->older;
  } else {
    pool->newest = chunk->older;
  }
  pool->cached_size -= chunk->size;
}

/*******************************************************************************
 * \brief Release detached chunks without holding the pool lock.
 ******************************************************************************/
static void mempool_dispose(offload_memchunk_t *chunk, const bool on_device) {
  while (chunk != NULL) {
    offload_memchunk_t *next = chunk->newer;
    actual_free(chunk->mem, on_device);
    free(chunk);
    chunk = next;
  }
}

/*******************************************************************************
 * \brief Private routine for allocating host or device memory from the pool.
 * \author Ole Schuett and Hans Pabst
 ******************************************************************************/
static void *internal_mempool_malloc(offload_mempool_t *pool, const size_t size,
                                     const bool on_device) {
  if (size == 0) {
    return NULL;
  }

  size_t capacity;
  const int bin = mempool_bin(size, &capacity);
  offload_memchunk_t *chunk = NULL;
#pragma omp critical(offload_mempool_modify)
  {
    // All chunks in a size class have the same capacity. Reuse only that class.
    if (bin < MEMPOOL_NBINS) {
      chunk = pool->available[bin];
    }
    if (chunk != NULL) {
      mempool_uncache(pool, chunk);
      mempool_register(pool, chunk, size);
    }
  }

  if (chunk == NULL) {
    // Backend allocation and memory pinning can be slow; keep them off the
    // lock.
    chunk = malloc(sizeof(*chunk));
    assert(chunk != NULL);
    chunk->size = capacity;
    chunk->mem = actual_malloc(chunk->size, on_device);
#pragma omp critical(offload_mempool_modify)
    {
      mempool_register(pool, chunk, size);
    }
  }
  return chunk->mem;
}

/*******************************************************************************
 * \brief Internal routine for allocating host memory from the pool.
 * \author Ole Schuett
 ******************************************************************************/
void *offload_mempool_host_malloc(const size_t size) {
  return internal_mempool_malloc(&mempool_host, size, false);
}

/*******************************************************************************
 * \brief Internal routine for allocating device memory from the pool
 * \author Ole Schuett
 ******************************************************************************/
void *offload_mempool_device_malloc(const size_t size) {
  return internal_mempool_malloc(&mempool_device, size, true);
}

/*******************************************************************************
 * \brief Private routine for releasing memory back to the pool.
 * \author Ole Schuett
 ******************************************************************************/
static void internal_mempool_free(offload_mempool_t *pool, const void *mem,
                                  const bool on_device) {
  if (mem == NULL) {
    return;
  }

  offload_memchunk_t *discard = NULL;
#pragma omp critical(offload_mempool_modify)
  {
    assert(pool->nallocated != 0);
    const size_t h = mempool_hash(mem, pool->nbuckets);
    offload_memchunk_t **indirect = &pool->allocated[h];
    while (*indirect != NULL && (*indirect)->mem != mem) {
      indirect = &(*indirect)->next;
    }
    offload_memchunk_t *chunk = *indirect;
    assert(chunk != NULL && chunk->mem == mem);
    *indirect = chunk->next;
    pool->nallocated--;
    pool->used -= chunk->used;
    pool->live_size -= chunk->size;
    chunk->used = 0;

    const int bin = mempool_bin(chunk->size, NULL);
    if (bin < MEMPOOL_NBINS) {
      chunk->next = pool->available[bin];
      chunk->prev = &pool->available[bin];
      if (chunk->next != NULL) {
        chunk->next->prev = &chunk->next;
      }
      pool->available[bin] = chunk;
      chunk->older = pool->newest;
      chunk->newer = NULL;
      if (pool->newest != NULL) {
        pool->newest->newer = chunk;
      } else {
        pool->oldest = chunk;
      }
      pool->newest = chunk;
      pool->cached_size += chunk->size;
    } else {
      chunk->newer = discard;
      discard = chunk;
    }

    // Keep at most one peak live working set, with a small floor for repeated
    // short-lived work. Evict the oldest cached chunks, retaining recent sizes.
    const uint64_t limit = imax(MEMPOOL_MIN_CACHE, pool->peak_live_size);
    while (pool->cached_size > limit) {
      offload_memchunk_t *victim = pool->oldest;
      mempool_uncache(pool, victim);
      victim->newer = discard;
      discard = victim;
    }
  }
  mempool_dispose(discard, on_device);
}

/*******************************************************************************
 * \brief Internal routine for releasing memory back to the pool.
 * \author Ole Schuett
 ******************************************************************************/
void offload_mempool_host_free(const void *memory) {
  internal_mempool_free(&mempool_host, memory, false);
}

/*******************************************************************************
 * \brief Internal routine for releasing memory back to the pool.
 * \author Ole Schuett
 ******************************************************************************/
void offload_mempool_device_free(const void *memory) {
  internal_mempool_free(&mempool_device, memory, true);
}

/*******************************************************************************
 * \brief Private routine for freeing all memory in the pool.
 * \author Ole Schuett
 ******************************************************************************/
static void internal_mempool_clear(offload_mempool_t *pool,
                                   const bool on_device) {
  offload_memchunk_t *discard;
#pragma omp critical(offload_mempool_modify)
  {
    assert(pool->nallocated == 0);
    discard = pool->oldest;
    memset(pool->available, 0, sizeof(pool->available));
    pool->oldest = pool->newest = NULL;
    pool->cached_size = pool->peak_live_size = 0;
    free(pool->allocated);
    pool->allocated = NULL;
    pool->nbuckets = 0;
  }
  mempool_dispose(discard, on_device);
}

/*******************************************************************************
 * \brief Internal routine for freeing all memory in the pool.
 * \author Ole Schuett and Hans Pabst
 ******************************************************************************/
void offload_mempool_clear(void) {
  internal_mempool_clear(&mempool_host, false);
  internal_mempool_clear(&mempool_device, true);
}

/*******************************************************************************
 * \brief Internal routine to query statistics.
 * \author Hans Pabst
 ******************************************************************************/
void offload_mempool_stats_get(offload_mempool_stats_t *memstats) {
  assert(NULL != memstats);
#pragma omp atomic read
  memstats->host_mallocs = host_malloc_counter;
#pragma omp atomic read
  memstats->device_mallocs = device_malloc_counter;
#pragma omp critical(offload_mempool_modify)
  {
    memstats->host_used = mempool_host.used;
    memstats->host_size = mempool_host.live_size + mempool_host.cached_size;
    memstats->host_peak = mempool_host.peak_size;
    memstats->device_used = mempool_device.used;
    memstats->device_size =
        mempool_device.live_size + mempool_device.cached_size;
    memstats->device_peak = mempool_device.peak_size;
  }
}

/*******************************************************************************
 * \brief Print allocation statistics..
 * \author Hans Pabst
 ******************************************************************************/
void offload_mempool_stats_print(int fortran_comm,
                                 void (*print_func)(const char *, int, int),
                                 int output_unit) {
  assert(omp_get_num_threads() == 1);

  char buffer[100];
  const cp_mpi_comm_t comm = cp_mpi_comm_f2c(fortran_comm);
  offload_mempool_stats_t memstats;
  offload_mempool_stats_get(&memstats);
  cp_mpi_max_uint64(&memstats.device_mallocs, 1, comm);
  cp_mpi_max_uint64(&memstats.host_mallocs, 1, comm);

  if (0 != memstats.device_mallocs || 0 != memstats.host_mallocs) {
    OFFLOAD_MEMPOOL_PRINT(print_func, "\n", output_unit);
    OFFLOAD_MEMPOOL_PRINT(
        print_func,
        " ----------------------------------------------------------------"
        "---------------\n",
        output_unit);
    OFFLOAD_MEMPOOL_PRINT(
        print_func,
        " -                                                               "
        "              -\n",
        output_unit);

    OFFLOAD_MEMPOOL_PRINT(
        print_func,
        " -                          OFFLOAD MEMPOOL STATISTICS           "
        "              -\n",
        output_unit);
    OFFLOAD_MEMPOOL_PRINT(
        print_func,
        " -                                                               "
        "              -\n",
        output_unit);
    OFFLOAD_MEMPOOL_PRINT(
        print_func,
        " ----------------------------------------------------------------"
        "---------------\n",
        output_unit);
    OFFLOAD_MEMPOOL_PRINT(print_func,
                          " Memory consumption               "
                          " Number of allocations  Used [MiB]  Size [MiB]\n",
                          output_unit);
  }
#if defined(__OFFLOAD)
  if (0 < memstats.device_mallocs) {
    cp_mpi_max_uint64(&memstats.device_peak, 1, comm);
    snprintf(buffer, sizeof(buffer),
             " Device                            "
             " %20" PRIuPTR "  %10" PRIuPTR "  %10" PRIuPTR "\n",
             (uintptr_t)memstats.device_mallocs,
             (uintptr_t)((memstats.device_used + (512U << 10)) >> 20),
             (uintptr_t)((memstats.device_peak + (512U << 10)) >> 20));
    OFFLOAD_MEMPOOL_PRINT(print_func, buffer, output_unit);
  }
#endif
  if (0 < memstats.host_mallocs) {
    cp_mpi_max_uint64(&memstats.host_peak, 1, comm);
    snprintf(buffer, sizeof(buffer),
             " Host                              "
             " %20" PRIuPTR "  %10" PRIuPTR "  %10" PRIuPTR "\n",
             (uintptr_t)memstats.host_mallocs,
             (uintptr_t)((memstats.host_used + (512U << 10)) >> 20),
             (uintptr_t)((memstats.host_peak + (512U << 10)) >> 20));
    OFFLOAD_MEMPOOL_PRINT(print_func, buffer, output_unit);
  }
  if (0 < memstats.device_mallocs || 0 < memstats.host_mallocs) {
    OFFLOAD_MEMPOOL_PRINT(
        print_func,
        " ----------------------------------------------------------------"
        "---------------\n",
        output_unit);
  }
}

// EOF
