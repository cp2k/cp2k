/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2025 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#ifndef DBM_MEMPOOL_H
#define DBM_MEMPOOL_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

/*******************************************************************************
 * \brief Internal routine for allocating host memory from the pool.
 * \author Ole Schuett
 ******************************************************************************/
void *dbm_mempool_host_malloc(const size_t size);

/*******************************************************************************
 * \brief Internal routine for allocating device memory from the pool.
 * \author Ole Schuett
 ******************************************************************************/
void *dbm_mempool_device_malloc(const size_t size);

/*******************************************************************************
 * \brief Internal routine for releasing memory back to the pool.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mempool_free(void *memory);

/*******************************************************************************
 * \brief Internal routine for freeing all memory in the pool (not thread-safe).
 * \author Ole Schuett
 ******************************************************************************/
void dbm_mempool_clear(void);

/*******************************************************************************
 * \brief Internal struct for pool statistics.
 * \author Hans Pabst
 ******************************************************************************/
typedef struct dbm_memstats_t {
  // Memory consumption (maximum).
  uint64_t host_size, device_size;
  // Number of allocations.
  uint64_t host_mallocs, device_mallocs;
} dbm_memstats_t;

/*******************************************************************************
 * \brief Internal routine to query statistics (not thread-safe).
 * \author Hans Pabst
 ******************************************************************************/
void dbm_mempool_statistics(dbm_memstats_t *memstats);

#ifdef __cplusplus
}
#endif

#endif

// EOF
