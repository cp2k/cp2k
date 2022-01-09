/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#ifndef DBM_MULTIPLY_COMM_H
#define DBM_MULTIPLY_COMM_H

#include "dbm_distribution.h"
#include "dbm_matrix.h"
#include "dbm_multiply_internal.h"

#include <stdbool.h>

/*******************************************************************************
 * \brief Internal struct for storing a packed matrix.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  const dbm_dist_1d_t *dist_indices;
  const dbm_dist_1d_t *dist_ticks;
  int nsend_packs;
  dbm_pack_t *send_packs;
  dbm_pack_t recv_pack;
  int max_nblocks; // Max across all ranks in dist_ticks.
  int max_data_size;
} dbm_packed_matrix_t;

/*******************************************************************************
 * \brief Internal struct for storing a communication iterator.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  int nticks;
  int itick;
  dbm_distribution_t *dist;
  dbm_packed_matrix_t packed_a;
  dbm_packed_matrix_t packed_b;
} dbm_comm_iterator_t;

/*******************************************************************************
 * \brief Internal routine for creating a communication iterator.
 * \author Ole Schuett
 ******************************************************************************/
dbm_comm_iterator_t *dbm_comm_iterator_start(const bool transa,
                                             const bool transb,
                                             const dbm_matrix_t *matrix_a,
                                             const dbm_matrix_t *matrix_b,
                                             const dbm_matrix_t *matrix_c);

/*******************************************************************************
 * \brief Internal routine for retriving next pair of packs from given iterator.
 * \author Ole Schuett
 ******************************************************************************/
bool dbm_comm_iterator_next(dbm_comm_iterator_t *iter, dbm_pack_t **pack_a,
                            dbm_pack_t **pack_b);

/*******************************************************************************
 * \brief Internal routine for releasing the given communication iterator.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_comm_iterator_stop(dbm_comm_iterator_t *iter);

#endif

// EOF
