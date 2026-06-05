/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#ifndef OFFLOAD_MEMPOOL_H
#define OFFLOAD_MEMPOOL_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/*******************************************************************************
 * \brief Internal routine for allocating host memory from the pool.
 * \author Ole Schuett
 ******************************************************************************/
void *offload_mempool_host_malloc(size_t size);

/*******************************************************************************
 * \brief Internal routine for allocating device memory from the pool.
 * \author Ole Schuett
 ******************************************************************************/
void *offload_mempool_device_malloc(size_t size);

/*******************************************************************************
 * \brief Internal routine for releasing memory back to the pool.
 * \author Ole Schuett
 ******************************************************************************/
void offload_mempool_host_free(const void *memory);

/*******************************************************************************
 * \brief Internal routine for releasing memory back to the pool.
 * \author Ole Schuett
 ******************************************************************************/
void offload_mempool_device_free(const void *memory);

/*******************************************************************************
 * \brief Internal routine for freeing all memory in the pool.
 * \author Ole Schuett
 ******************************************************************************/
void offload_mempool_clear(void);

/*******************************************************************************
 * \brief Internal struct for pool statistics.
 * \author Hans Pabst
 ******************************************************************************/
typedef struct offload_mempool_stats_t {
  // Memory consumption (entire pool).
  uint64_t host_size, device_size;
  // Memory currently used.
  uint64_t host_used, device_used;
  // Memory peak consumption.
  uint64_t host_peak, device_peak;
  // Number of allocations.
  uint64_t host_mallocs, device_mallocs;
} offload_mempool_stats_t;

/*******************************************************************************
 * \brief Internal routine to query statistics.
 * \author Hans Pabst
 ******************************************************************************/
void offload_mempool_stats_get(offload_mempool_stats_t *memstats);

/*******************************************************************************
 * \brief Print allocation statistics.
 * \author Hans Pabst
 ******************************************************************************/
void offload_mempool_stats_print(int fortran_comm,
                                 void (*print_func)(const char *, int, int),
                                 int output_unit);

#ifdef __cplusplus
}
#endif

#endif

// EOF
