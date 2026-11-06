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

/*******************************************************************************
 * \brief Private struct for storing a chunk of memory.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct offload_memchunk {
  void *mem; // first: allows to cast memchunk into mem-ptr...
  struct offload_memchunk *next;
  size_t size, used;
} offload_memchunk_t;

/*******************************************************************************
 * \brief Private struct for storing a memory pool.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct offload_mempool {
  offload_memchunk_t *available_head, *allocated_head; // single-linked lists
  uint64_t available_size, allocated_size; // sums of chunk sizes in the lists
  uint64_t allocated_used, used_peak;      // requested sizes in use, most ever
  uint64_t peak_size;                      // for statistics
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
 * \brief Private routine for taking the largest available chunks out of the
 *        pool while it holds more than its bound (caller in critical region).
 *        Returns the chunks to be freed outside of the critical region.
 ******************************************************************************/
static offload_memchunk_t *internal_mempool_trim(offload_mempool_t *pool) {
  offload_memchunk_t *release = NULL;
  const uint64_t bound =
      pool->used_peak + pool->used_peak / OFFLOAD_MEMPOOL_SLACK;
  while (pool->available_head != NULL &&
         bound < pool->available_size + pool->allocated_size) {
    offload_memchunk_t **largest = &pool->available_head;
    for (offload_memchunk_t **indirect = &(*largest)->next; *indirect != NULL;
         indirect = &(*indirect)->next) {
      if ((*largest)->size < (*indirect)->size) {
        largest = indirect;
      }
    }
    offload_memchunk_t *const chunk = *largest;
    *largest = chunk->next; // remove chunk from available list
    pool->available_size -= chunk->size;
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

  offload_memchunk_t *chunk = NULL, *release = NULL;

#pragma omp critical(offload_mempool_modify)
  {
    // Find the smallest available chunk that fits (not overly large).
    offload_memchunk_t **reuse = NULL; // ** for easy list removal
    offload_memchunk_t **indirect = &pool->available_head;
    while (*indirect != NULL) {
      const size_t s = (*indirect)->size;
      if (size <= s && s - size <= size / OFFLOAD_MEMPOOL_SLACK &&
          (reuse == NULL || s < (*reuse)->size)) {
        reuse = indirect;
        if (s == size) {
          break; // perfect match, exit early
        }
      }
      indirect = &(*indirect)->next;
    }
    if (reuse != NULL) {
      chunk = *reuse;
      *reuse = chunk->next; // remove chunk from available list
      pool->available_size -= chunk->size;
    }
  }

  // Allocate a new chunk outside of critical region (none fits).
  if (chunk == NULL) {
    chunk = calloc(1, sizeof(offload_memchunk_t));
    assert(chunk != NULL);
    chunk->mem = actual_malloc(size, on_device);
    chunk->size = size;
  }

  chunk->used = size; // for statistics

  // Insert chunk into allocated list, and keep the pool within its bound.
#pragma omp critical(offload_mempool_modify)
  {
    chunk->next = pool->allocated_head;
    pool->allocated_head = chunk;
    pool->allocated_size += chunk->size;
    pool->allocated_used += size;
    pool->used_peak = imax(pool->used_peak, pool->allocated_used);
    pool->peak_size =
        imax(pool->peak_size, pool->available_size + pool->allocated_size);
    release = internal_mempool_trim(pool);
  }

  internal_mempool_release(release, on_device);
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

  offload_memchunk_t *release = NULL;

#pragma omp critical(offload_mempool_modify)
  {
    // Find chunk in allocated list.
    offload_memchunk_t **indirect = &pool->allocated_head;
    while (*indirect != NULL && (*indirect)->mem != mem) {
      indirect = &(*indirect)->next;
    }
    offload_memchunk_t *chunk = *indirect;
    assert(chunk != NULL && chunk->mem == mem);

    // Remove chunk from allocated list.
    *indirect = chunk->next;
    pool->allocated_size -= chunk->size;
    pool->allocated_used -= chunk->used;

    // Add chunk to available list, and keep the pool within its bound.
    chunk->next = pool->available_head;
    pool->available_head = chunk;
    pool->available_size += chunk->size;
    release = internal_mempool_trim(pool);
  }

  internal_mempool_release(release, on_device);
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

#pragma omp critical(offload_mempool_modify)
  {
    uint64_t pool_size = 0;

    // Check for leaks, i.e. that the allocated list is empty.
    assert(pool->allocated_head == NULL);

    // Free all chunks in available list.
    while (pool->available_head != NULL) {
      offload_memchunk_t *chunk = pool->available_head;
      pool->available_head = chunk->next; // remove chunk
      actual_free(chunk->mem, on_device);
      pool_size += chunk->size;
      free(chunk);
    }
    pool->available_size = 0;
    pool->used_peak = 0;

    // Update stats.
    pool->peak_size = imax(pool->peak_size, pool_size);
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
 * \brief Private routine for summing used sizes of all chunks in given list.
 * \author Ole Schuett
 ******************************************************************************/
static uint64_t sum_chunks_used(const offload_memchunk_t *head) {
  uint64_t used_sum = 0;
  for (const offload_memchunk_t *chunk = head; chunk != NULL;
       chunk = chunk->next) {
    used_sum += chunk->used;
  }
  return used_sum;
}

/*******************************************************************************
 * \brief Internal routine to query statistics.
 * \author Hans Pabst
 ******************************************************************************/
void offload_mempool_stats_get(offload_mempool_stats_t *memstats) {
  assert(NULL != memstats);
#pragma omp critical(offload_mempool_modify)
  {
    memstats->host_mallocs = host_malloc_counter;
    memstats->host_used = sum_chunks_used(mempool_host.available_head) +
                          sum_chunks_used(mempool_host.allocated_head);
    memstats->host_size =
        mempool_host.available_size + mempool_host.allocated_size;
    memstats->host_peak = imax(mempool_host.peak_size, memstats->host_size);

    memstats->device_mallocs = device_malloc_counter;
    memstats->device_used = sum_chunks_used(mempool_device.available_head) +
                            sum_chunks_used(mempool_device.allocated_head);
    memstats->device_size =
        mempool_device.available_size + mempool_device.allocated_size;
    memstats->device_peak =
        imax(mempool_device.peak_size, memstats->device_size);
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
