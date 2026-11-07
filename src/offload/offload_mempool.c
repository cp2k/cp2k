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
// A reused chunk exceeds the requested size by at most this fraction, and the
// pool holds at most the most ever requested at once plus this fraction of it:
// what a workload needs again stays cached, whereas sizes that do not fit
// cannot let the pool grow without bound (the bound covers the reuse waste).
#define OFFLOAD_MEMPOOL_SLACK 4 /* 1/4 */
enum { MEMPOOL_INITIAL_BUCKETS = 64 };

/*******************************************************************************
 * \brief Private struct for storing a chunk of memory.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct offload_memchunk {
  void *mem; // first: allows to cast memchunk into mem-ptr...
  struct offload_memchunk *next; // hash chain or idle stack of one size
  struct offload_memchunk *left;
  union {
    struct offload_memchunk *right; // idle size tree
    size_t used;                    // requested bytes while on loan
  };
  size_t size;
} offload_memchunk_t;

/*******************************************************************************
 * \brief Private struct for storing a memory pool.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct offload_mempool {
  offload_memchunk_t *available;  // one tree node per distinct idle capacity
  offload_memchunk_t **allocated; // allocated chunks indexed by address
  offload_memchunk_t *initial_buckets[MEMPOOL_INITIAL_BUCKETS];
  size_t nbuckets, nallocated;
  uint64_t mallocs; // backend allocations, retained across clear
  uint64_t available_size, allocated_size;
  uint64_t allocated_used, used_peak; // requested sizes in use, most ever
  uint64_t peak_size;                 // peak total capacity before trimming
} offload_mempool_t;

/*******************************************************************************
 * \brief Private pools for host and device memory.
 * \author Ole Schuett
 ******************************************************************************/
static offload_mempool_t mempool_host = {.allocated =
                                             mempool_host.initial_buckets,
                                         .nbuckets = MEMPOOL_INITIAL_BUCKETS};
static offload_mempool_t mempool_device = {.allocated =
                                               mempool_device.initial_buckets,
                                           .nbuckets = MEMPOOL_INITIAL_BUCKETS};

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
  (void)on_device; // used only by offload backends

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
 * \brief Bucket for an aligned host or device address (power-of-two table).
 ******************************************************************************/
static size_t mempool_bucket(const void *mem, size_t nbuckets) {
  uint64_t key = (uintptr_t)mem;
  key ^= key >> 32;
  key *= UINT64_C(11400714819323198485);
  return (key ^ (key >> 32)) & (nbuckets - 1);
}

/*******************************************************************************
 * \brief Grow the address index. Chaining still works if resizing fails.
 ******************************************************************************/
static void mempool_grow_table(offload_mempool_t *pool) {
  const size_t n = 2 * pool->nbuckets;
  offload_memchunk_t **buckets = calloc(n, sizeof(*buckets));
  if (buckets == NULL) {
    return;
  }
  for (size_t i = 0; i < pool->nbuckets; ++i) {
    offload_memchunk_t *entry = pool->allocated[i];
    while (entry != NULL) {
      offload_memchunk_t *next = entry->next;
      const size_t bucket = mempool_bucket(entry->mem, n);
      entry->next = buckets[bucket];
      buckets[bucket] = entry;
      entry = next;
    }
  }
  if (pool->allocated != pool->initial_buckets) {
    free(pool->allocated);
  }
  pool->allocated = buckets;
  pool->nbuckets = n;
}

/*******************************************************************************
 * \brief Register a loan and account for the requested and actual chunk sizes.
 *        All index operations on the pool use offload_mempool_modify.
 ******************************************************************************/
static void mempool_register(offload_mempool_t *pool, offload_memchunk_t *chunk,
                             size_t requested) {
  if (pool->nallocated >= pool->nbuckets) {
    mempool_grow_table(pool);
  }
  const size_t bucket = mempool_bucket(chunk->mem, pool->nbuckets);
  chunk->next = pool->allocated[bucket];
  pool->allocated[bucket] = chunk;
  ++pool->nallocated;
  chunk->used = requested;
  pool->allocated_size += chunk->size;
  pool->allocated_used += requested;
  pool->used_peak = imax(pool->used_peak, pool->allocated_used);
}

/*******************************************************************************
 * \brief Bring an idle capacity (or the last node on its search path) to the
 *        root. Top-down splaying needs no parent pointers or recursion.
 ******************************************************************************/
static offload_memchunk_t *mempool_splay(offload_memchunk_t *root,
                                         size_t size) {
  offload_memchunk_t head, *less = &head, *greater = &head;
  head.left = head.right = NULL;
  for (;;) {
    if (size < root->size && root->left != NULL) {
      offload_memchunk_t *child = root->left;
      if (size < child->size) {
        root->left = child->right;
        child->right = root;
        root = child;
        if (root->left == NULL) {
          break;
        }
      }
      greater->left = root;
      greater = root;
      root = root->left;
    } else if (size > root->size && root->right != NULL) {
      offload_memchunk_t *child = root->right;
      if (size > child->size) {
        root->right = child->left;
        child->left = root;
        root = child;
        if (root->right == NULL) {
          break;
        }
      }
      less->right = root;
      less = root;
      root = root->right;
    } else {
      break;
    }
  }
  less->right = root->left;
  greater->left = root->right;
  root->left = head.right;
  root->right = head.left;
  return root;
}

/*******************************************************************************
 * \brief Insert a returned chunk by its actual capacity. Within a size, the
 *        most recently returned chunk is the first one reused or trimmed.
 ******************************************************************************/
static void mempool_cache(offload_mempool_t *pool, offload_memchunk_t *chunk) {
  offload_memchunk_t *root = pool->available;
  if (root != NULL && root->size != chunk->size) {
    root = mempool_splay(root, chunk->size);
  }
  if (root != NULL && root->size == chunk->size) {
    chunk->next = root->next;
    root->next = chunk;
    pool->available = root;
  } else {
    chunk->next = chunk->left = chunk->right = NULL;
    if (root != NULL) {
      if (chunk->size < root->size) {
        chunk->left = root->left;
        chunk->right = root;
        root->left = NULL;
      } else {
        chunk->right = root->right;
        chunk->left = root;
        root->right = NULL;
      }
    }
    pool->available = chunk;
  }
  pool->available_size += chunk->size;
}

/*******************************************************************************
 * \brief Take a chunk from the size currently at the root of the idle tree.
 ******************************************************************************/
static offload_memchunk_t *mempool_take_root(offload_mempool_t *pool) {
  offload_memchunk_t *root = pool->available;
  offload_memchunk_t *chunk = root->next;
  if (chunk != NULL) {
    root->next = chunk->next;
  } else {
    chunk = root;
    if (root->left == NULL) {
      pool->available = root->right;
    } else if (root->right == NULL) {
      pool->available = root->left;
    } else {
      pool->available = mempool_splay(root->left, root->size);
      pool->available->right = root->right;
    }
  }
  pool->available_size -= chunk->size;
  return chunk;
}

/*******************************************************************************
 * \brief Take the smallest idle chunk within the existing reuse tolerance.
 ******************************************************************************/
static offload_memchunk_t *mempool_take_fit(offload_mempool_t *pool,
                                            size_t size) {
  offload_memchunk_t *root = pool->available;
  if (root == NULL) {
    return NULL;
  }
  if (root->size != size) {
    root = mempool_splay(root, size);
    if (root->size < size && root->right != NULL) {
      // The search ended at the predecessor. Promote its successor.
      offload_memchunk_t *successor = mempool_splay(root->right, size);
      assert(successor->left == NULL);
      root->right = NULL;
      successor->left = root;
      root = successor;
    }
    pool->available = root;
  }
  if (root->size < size || root->size - size > size / OFFLOAD_MEMPOOL_SLACK) {
    return NULL;
  }
  return mempool_take_root(pool);
}

/*******************************************************************************
 * \brief Trim the largest idle chunks after a backend allocation. A cache hit
 *        or return only transfers capacity between the idle and live sets.
 ******************************************************************************/
static offload_memchunk_t *internal_mempool_trim(offload_mempool_t *pool) {
  offload_memchunk_t *release = NULL;
  uint64_t total = pool->available_size + pool->allocated_size;
  while (total > pool->used_peak &&
         total - pool->used_peak > pool->used_peak / OFFLOAD_MEMPOOL_SLACK) {
    assert(pool->available != NULL);
    if (pool->available->right != NULL) {
      pool->available = mempool_splay(pool->available, SIZE_MAX);
    }
    offload_memchunk_t *chunk = mempool_take_root(pool);
    total -= chunk->size;
    chunk->next = release;
    release = chunk;
  }
  return release;
}

/*******************************************************************************
 * \brief Private routine for freeing the chunks taken out of the pool.
 ******************************************************************************/
static void internal_mempool_release(offload_memchunk_t *release,
                                     const bool on_device) {
  while (release != NULL) {
    offload_memchunk_t *const chunk = release;
    release = chunk->next;
    actual_free(chunk->mem, on_device);
    free(chunk);
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

  offload_memchunk_t *chunk;
#pragma omp critical(offload_mempool_modify)
  {
    chunk = mempool_take_fit(pool, size);
    if (chunk != NULL) {
      mempool_register(pool, chunk, size);
    }
  }

  if (chunk == NULL) {
    chunk = malloc(sizeof(*chunk));
    if (chunk == NULL) {
      return NULL;
    }
    chunk->mem = actual_malloc(size, on_device);
    if (chunk->mem == NULL) {
      free(chunk);
      return NULL;
    }
    chunk->size = size;
    offload_memchunk_t *release;
#pragma omp critical(offload_mempool_modify)
    {
      mempool_register(pool, chunk, size);
      ++pool->mallocs;
      pool->peak_size =
          imax(pool->peak_size, pool->available_size + pool->allocated_size);
      release = internal_mempool_trim(pool);
    }
    internal_mempool_release(release, on_device);
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
static void internal_mempool_free(offload_mempool_t *pool, const void *mem) {
  if (mem == NULL) {
    return;
  }
#pragma omp critical(offload_mempool_modify)
  {
    assert(pool->nallocated > 0);
    const size_t bucket = mempool_bucket(mem, pool->nbuckets);
    offload_memchunk_t **link = &pool->allocated[bucket];
    while (*link != NULL && (*link)->mem != mem) {
      link = &(*link)->next;
    }
    offload_memchunk_t *chunk = *link;
    assert(chunk != NULL);
    *link = chunk->next;
    --pool->nallocated;
    pool->allocated_size -= chunk->size;
    pool->allocated_used -= chunk->used;
    mempool_cache(pool, chunk);
  }
}

/*******************************************************************************
 * \brief Internal routine for releasing memory back to the pool.
 * \author Ole Schuett
 ******************************************************************************/
void offload_mempool_host_free(const void *memory) {
  internal_mempool_free(&mempool_host, memory);
}

/*******************************************************************************
 * \brief Internal routine for releasing memory back to the pool.
 * \author Ole Schuett
 ******************************************************************************/
void offload_mempool_device_free(const void *memory) {
  internal_mempool_free(&mempool_device, memory);
}

/*******************************************************************************
 * \brief Private routine for freeing all memory in the pool.
 * \author Ole Schuett
 ******************************************************************************/
static void internal_mempool_clear(offload_mempool_t *pool,
                                   const bool on_device) {
  offload_memchunk_t *available, **buckets = NULL;
#pragma omp critical(offload_mempool_modify)
  {
    assert(pool->nallocated == 0);
    if (pool->nallocated == 0 && pool->allocated != pool->initial_buckets) {
      buckets = pool->allocated;
      memset(pool->initial_buckets, 0, sizeof(pool->initial_buckets));
      pool->allocated = pool->initial_buckets;
      pool->nbuckets = MEMPOOL_INITIAL_BUCKETS;
    }
    available = pool->available;
    pool->available = NULL;
    pool->available_size = 0;
    pool->used_peak = pool->allocated_used;
  }
  free(buckets);
  while (available != NULL) {
    if (available->right != NULL) {
      available = mempool_splay(available, SIZE_MAX);
    }
    offload_memchunk_t *next = available->left;
    internal_mempool_release(available, on_device);
    available = next;
  }
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
#pragma omp critical(offload_mempool_modify)
  {
    memstats->host_mallocs = mempool_host.mallocs;
    memstats->host_used = mempool_host.allocated_used;
    memstats->host_size =
        mempool_host.available_size + mempool_host.allocated_size;
    memstats->host_peak = mempool_host.peak_size;

    memstats->device_mallocs = mempool_device.mallocs;
    memstats->device_used = mempool_device.allocated_used;
    memstats->device_size =
        mempool_device.available_size + mempool_device.allocated_size;
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
