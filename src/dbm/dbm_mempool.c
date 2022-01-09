/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <omp.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#if defined(__DBM_CUDA)
#include <cuda_runtime.h>
#endif

#include "../offload/offload_library.h"
#include "dbm_mempool.h"

/*******************************************************************************
 * \brief Check given Cuda status and upon failure abort with a nice message.
 * \author Ole Schuett
 ******************************************************************************/
#define CHECK(status)                                                          \
  if (status != cudaSuccess) {                                                 \
    fprintf(stderr, "ERROR: %s %s %d\n", cudaGetErrorString(status), __FILE__, \
            __LINE__);                                                         \
    abort();                                                                   \
  }

/*******************************************************************************
 * \brief Private routine for actually allocating system memory.
 * \author Ole Schuett
 ******************************************************************************/
static void *actual_malloc(const size_t size, const bool on_device) {
  (void)on_device; // mark used

#if defined(__DBM_CUDA)
  if (on_device) {
    void *memory;
    offload_set_device();
    CHECK(cudaMalloc(&memory, size));
    assert(memory != NULL);
    return memory;
  }
#else
  (void)on_device; // mark used
#endif

  void *memory = malloc(size);
  assert(memory != NULL);
  return memory;
}

/*******************************************************************************
 * \brief Private routine for actually freeing system memory.
 * \author Ole Schuett
 ******************************************************************************/
static void actual_free(void *memory, const bool on_device) {
  if (memory == NULL) {
    return;
  }

#if defined(__DBM_CUDA)
  if (on_device) {
    offload_set_device();
    CHECK(cudaFree(memory));
    return;
  }
#else
  (void)on_device; // mark used
#endif

  free(memory);
}

/*******************************************************************************
 * \brief Private struct for storing a chunk of memory.
 * \author Ole Schuett
 ******************************************************************************/
struct dbm_memchunk {
  bool on_device;
  size_t size;
  void *mem;
  struct dbm_memchunk *next;
};
typedef struct dbm_memchunk dbm_memchunk_t;

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
 * \brief Private routine for allocating host or device memory from the pool.
 * \author Ole Schuett
 ******************************************************************************/
static void *internal_mempool_malloc(const size_t size, const bool on_device) {
  if (size == 0) {
    return NULL;
  }

  dbm_memchunk_t *chunk;

#pragma omp critical(dbm_mempool_modify)
  {
    // Find a suitable chuck in mempool_available.
    dbm_memchunk_t **indirect = &mempool_available_head;
    while (*indirect != NULL && (*indirect)->on_device != on_device) {
      indirect = &(*indirect)->next;
    }
    chunk = *indirect;

    // If a chunck was found, remove it from mempool_available.
    if (chunk != NULL) {
      assert(chunk->on_device == on_device);
      *indirect = chunk->next;
    }

    // If no chunk was found, allocate a new one.
    if (chunk == NULL) {
      chunk = malloc(sizeof(dbm_memchunk_t));
      chunk->on_device = on_device;
      chunk->size = 0;
      chunk->mem = NULL;
    }

    // Resize chunk if needed.
    if (chunk->size < size) {
      actual_free(chunk->mem, chunk->on_device);
      chunk->mem = actual_malloc(size, chunk->on_device);
      chunk->size = size;
    }

    // Insert chunk into mempool_allocated.
    chunk->next = mempool_allocated_head;
    mempool_allocated_head = chunk;
  }

  return chunk->mem;
}

/*******************************************************************************
 * \brief Internal routine for allocating host memory from the pool.
 * \author Ole Schuett
 ******************************************************************************/
void *dbm_mempool_host_malloc(const size_t size) {
  return internal_mempool_malloc(size, false);
}

/*******************************************************************************
 * \brief Internal routine for allocating device memory from the pool
 * \author Ole Schuett
 ******************************************************************************/
void *dbm_mempool_device_malloc(const size_t size) {
  return internal_mempool_malloc(size, true);
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
    // Find chuck in mempool_allocated.
    dbm_memchunk_t **indirect = &mempool_allocated_head;
    while (*indirect != NULL && (*indirect)->mem != mem) {
      indirect = &(*indirect)->next;
    }
    dbm_memchunk_t *chunk = *indirect;
    assert(chunk != NULL && chunk->mem == mem);

    // Remove chuck from mempool_allocated.
    *indirect = chunk->next;

    // Add chuck to mempool_available.
    chunk->next = mempool_available_head;
    mempool_available_head = chunk;
  }
}

/*******************************************************************************
 * \brief Internal routine for freeing all memory in the pool.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mempool_clear(void) {
  assert(omp_get_num_threads() == 1);
  assert(mempool_allocated_head == NULL); // check for memory leak

  // while (mempool_allocated_head != NULL) {
  //  dbm_memchunk_t *chunk = mempool_allocated_head;
  //  mempool_allocated_head = chunk->next;
  //  printf("Found alloacted memory chunk of size: %lu\n", chunk->size);
  //  actual_free(chunk->mem, chunk->on_device);
  //  free(chunk);
  //}

  // Free chunks in mempool_avavailable.
  while (mempool_available_head != NULL) {
    dbm_memchunk_t *chunk = mempool_available_head;
    mempool_available_head = chunk->next;
    actual_free(chunk->mem, chunk->on_device);
    free(chunk);
  }
}

// EOF
