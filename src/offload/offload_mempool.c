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

#if defined(__LIBXSTREAM)
#include <libxstream/libxstream.h>
#include <libxstream/libxstream_opencl.h>
#elif defined(__LIBXS)
#include <libxs/libxs_malloc.h>
#endif

#define OFFLOAD_MEMPOOL_PRINT(FN, MSG, OUTPUT_UNIT)                            \
  ((FN)(MSG, (int)strlen(MSG), OUTPUT_UNIT))
#define OFFLOAD_MEMPOOL_OMPALLOC 1

#if !defined(__LIBXSTREAM)
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
} offload_mempool_t;

/*******************************************************************************
 * \brief Private pools for host and device memory.
 * \author Ole Schuett
 ******************************************************************************/
static offload_mempool_t mempool_host = {0};
static offload_mempool_t mempool_device = {0};

/*******************************************************************************
 * \brief Private counters for statistics.
 * \author Hans Pabst
 ******************************************************************************/
static struct {
  uint64_t mallocs, mempeak;
} host_stats = {0, 0};
static struct {
  uint64_t mallocs, mempeak;
} device_stats = {0, 0};

/*******************************************************************************
 * \brief Private routine for actually allocating system memory.
 * \author Ole Schuett
 ******************************************************************************/
static void *actual_malloc(const size_t size, const bool on_device) {
  if (size == 0) {
    return NULL;
  }

  void *memory = NULL;
#if defined(__OFFLOAD)
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
    ++device_stats.mallocs;
  } else {
#pragma omp atomic
    ++host_stats.mallocs;
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

#if defined(__OFFLOAD)
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
    // Find a possible chunk to reuse or reclaim in available list.
    offload_memchunk_t **reuse = NULL,
                       **reclaim = NULL; // ** for easy list removal
    offload_memchunk_t **indirect = &pool->available_head;
    while (*indirect != NULL) {
      const size_t s = (*indirect)->size;
      if (size <= s && (reuse == NULL || s < (*reuse)->size)) {
        reuse = indirect; // reuse smallest suitable chunk
        if (s == size) {
          break; // perfect match, exit early
        }
      } else if (reclaim == NULL || (*reclaim)->size < s) {
        reclaim = indirect; // reclaim largest unsuitable chunk
      }
      indirect = &(*indirect)->next;
    }

    // Select an existing chunk or allocate a new one.
    if (reuse != NULL) {
      // Reusing an exising chunk that's already large enough.
      chunk = *reuse;
      *reuse = chunk->next; // remove chunk from available list.
    } else if (reclaim != NULL) {
      // Reclaiming an existing chunk (resize will happen outside crit. region).
      chunk = *reclaim;
      *reclaim = chunk->next; // remove chunk from available list.
    } else {
      // Found no available chunk, allocate a new one.
      chunk = calloc(1, sizeof(offload_memchunk_t));
      assert(chunk != NULL);
    }
  }

  // Resize chunk outside of critical region before adding it to allocated list.
  if (chunk->size < size) {
    actual_free(chunk->mem, on_device);
    chunk->mem = actual_malloc(size, on_device);
    chunk->size = size;
  }

  chunk->used = size; // for statistics

  // Insert chunk into allocated list.
#pragma omp critical(offload_mempool_modify)
  {
    chunk->next = pool->allocated_head;
    pool->allocated_head = chunk;
  }

  return chunk->mem;
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
    offload_memchunk_t **indirect = &pool->allocated_head;
    while (*indirect != NULL && (*indirect)->mem != mem) {
      indirect = &(*indirect)->next;
    }
    offload_memchunk_t *chunk = *indirect;
    assert(chunk != NULL && chunk->mem == mem);
    *indirect = chunk->next;
    chunk->next = pool->available_head;
    pool->available_head = chunk;
  }
}

/*******************************************************************************
 * \brief Private routine for freeing all memory in the pool.
 * \author Ole Schuett and Hans Pabst
 ******************************************************************************/
static void internal_mempool_clear(offload_mempool_t *pool,
                                   const bool on_device) {
#pragma omp critical(offload_mempool_modify)
  {
    assert(pool->allocated_head == NULL);
    while (pool->available_head != NULL) {
      offload_memchunk_t *chunk = pool->available_head;
      pool->available_head = chunk->next;
      actual_free(chunk->mem, on_device);
      free(chunk);
    }
  }
}

/*******************************************************************************
 * \brief Private routine for summing alloc sizes of all chunks in given list.
 * \author Ole Schuett and Hans Pabst
 ******************************************************************************/
static uint64_t sum_chunks_size(const offload_memchunk_t *head, size_t offset) {
  uint64_t result = 0;
  for (const offload_memchunk_t *chunk = head; chunk != NULL;
       chunk = chunk->next) {
    result += *(const size_t *)((const char *)chunk + offset);
  }
  return result;
}
#endif /* !defined(__LIBXSTREAM) */

/*******************************************************************************
 * \brief Internal routine for allocating host memory from the pool.
 * \author Ole Schuett
 ******************************************************************************/
void *offload_mempool_host_malloc(const size_t size) {
#if defined(__LIBXSTREAM)
  return libxs_malloc(libxstream_opencl_config.pool_hst, size,
                      LIBXS_MALLOC_AUTO);
#else
  return internal_mempool_malloc(&mempool_host, size, false);
#endif
}

/*******************************************************************************
 * \brief Internal routine for allocating device memory from the pool
 * \author Ole Schuett
 ******************************************************************************/
void *offload_mempool_device_malloc(const size_t size) {
#if defined(__LIBXSTREAM)
  void *memory = NULL;
  const int result = libxstream_mem_allocate(&memory, size);
  assert(EXIT_SUCCESS == result);
  return memory;
#else
  return internal_mempool_malloc(&mempool_device, size, true);
#endif
}

/*******************************************************************************
 * \brief Internal routine for releasing memory back to the pool.
 * \author Ole Schuett
 ******************************************************************************/
void offload_mempool_host_free(const void *memory) {
#if defined(__LIBXSTREAM)
  libxs_free((void *)memory);
#else
  internal_mempool_free(&mempool_host, memory);
#endif
}

/*******************************************************************************
 * \brief Internal routine for releasing memory back to the pool.
 * \author Ole Schuett
 ******************************************************************************/
void offload_mempool_device_free(const void *memory) {
#if defined(__LIBXSTREAM)
  const int result = libxstream_mem_deallocate((void *)memory);
  assert(EXIT_SUCCESS == result);
#else
  internal_mempool_free(&mempool_device, memory);
#endif
}

/*******************************************************************************
 * \brief Internal routine for freeing all memory in the pool.
 * \author Ole Schuett
 ******************************************************************************/
void offload_mempool_clear(void) {
#if defined(__LIBXSTREAM)
  (void)0;
#else
  {
    const uint64_t hsize = sum_chunks_size(mempool_host.available_head,
                                           offsetof(offload_memchunk_t, size)) +
                           sum_chunks_size(mempool_host.allocated_head,
                                           offsetof(offload_memchunk_t, size));
    const uint64_t dsize = sum_chunks_size(mempool_device.available_head,
                                           offsetof(offload_memchunk_t, size)) +
                           sum_chunks_size(mempool_device.allocated_head,
                                           offsetof(offload_memchunk_t, size));
    if (host_stats.mempeak < hsize)
      host_stats.mempeak = hsize;
    if (device_stats.mempeak < dsize)
      device_stats.mempeak = dsize;
  }
  internal_mempool_clear(&mempool_host, false);
  internal_mempool_clear(&mempool_device, true);
#endif
}

/*******************************************************************************
 * \brief Internal routine to query statistics.
 * \author Hans Pabst
 ******************************************************************************/
void offload_mempool_stats_get(offload_mempool_stats_t *memstats) {
  assert(NULL != memstats);
#pragma omp critical(offload_mempool_modify)
  {
#if defined(__LIBXSTREAM)
    if (NULL != libxstream_opencl_config.pool_hst) {
      libxs_malloc_pool_info_t info;
      libxs_malloc_pool_info(libxstream_opencl_config.pool_hst, &info);
      memstats->host_mallocs = info.nmallocs;
      memstats->host_used = info.used;
      memstats->host_size = info.size;
      memstats->host_peak = info.peak;
    } else {
      memstats->host_mallocs = 0;
      memstats->host_used = 0;
      memstats->host_size = 0;
      memstats->host_peak = 0;
    }
    if (NULL != libxstream_opencl_config.pool_dev) {
      libxs_malloc_pool_info_t info;
      libxs_malloc_pool_info(libxstream_opencl_config.pool_dev, &info);
      memstats->device_mallocs = info.nmallocs;
      memstats->device_used = info.used;
      memstats->device_size = info.size;
      memstats->device_peak = info.peak;
    } else {
      memstats->device_mallocs = 0;
      memstats->device_used = 0;
      memstats->device_size = 0;
      memstats->device_peak = 0;
    }
#else
    memstats->host_mallocs = host_stats.mallocs;
    memstats->host_used = sum_chunks_size(mempool_host.available_head,
                                          offsetof(offload_memchunk_t, used)) +
                          sum_chunks_size(mempool_host.allocated_head,
                                          offsetof(offload_memchunk_t, used));
    memstats->host_size = sum_chunks_size(mempool_host.available_head,
                                          offsetof(offload_memchunk_t, size)) +
                          sum_chunks_size(mempool_host.allocated_head,
                                          offsetof(offload_memchunk_t, size));
    memstats->host_peak = memstats->host_size < host_stats.mempeak
                              ? host_stats.mempeak
                              : memstats->host_size;
    memstats->device_mallocs = device_stats.mallocs;
    memstats->device_used =
        sum_chunks_size(mempool_device.available_head,
                        offsetof(offload_memchunk_t, used)) +
        sum_chunks_size(mempool_device.allocated_head,
                        offsetof(offload_memchunk_t, used));
    memstats->device_size =
        sum_chunks_size(mempool_device.available_head,
                        offsetof(offload_memchunk_t, size)) +
        sum_chunks_size(mempool_device.allocated_head,
                        offsetof(offload_memchunk_t, size));
    memstats->device_peak = memstats->device_size < device_stats.mempeak
                                ? device_stats.mempeak
                                : memstats->device_size;
#endif
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
