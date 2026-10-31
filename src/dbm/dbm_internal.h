/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#ifndef DBM_INTERNAL_H
#define DBM_INTERNAL_H

#include <stdbool.h>
#include <stdint.h>

/*******************************************************************************
 * \brief Returns the larger of two given integers (missing from the C standard)
 * \author Ole Schuett
 ******************************************************************************/
static inline int imax(int x, int y) { return (x > y ? x : y); }

/*******************************************************************************
 * \brief Internal struct for storing a dbm_block_t plus its norm.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  int free_index; // Free index in Einstein notation of matrix multiplication.
  int sum_index;  // Summation index - also called dummy index.
  int offset;
  float norm;
} dbm_pack_block_t;

/*******************************************************************************
 * \brief Internal struct for storing a pack - essentially a shard for MPI.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  int nblocks;
  int data_size;
  dbm_pack_block_t *blocks;
  double *data;
} dbm_pack_t;

/*******************************************************************************
 * \brief Internal struct for storing a task, ie. a single block multiplication.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  int m;
  int n;
  int k;
  int offset_a;
  int offset_b;
  int offset_c;
} dbm_task_t;

/*******************************************************************************
 * \brief Internal struct for the shape of a batch, which is tracked while the
 *        batch is generated such that a backend needs no pass over the tasks.
 *        The batch is homogeneous if flops == 2*ntasks*max_m*max_n*max_k,
 *        otherwise their ratio is the fill of a kernel padded to the maxima.
 * \author Hans Pabst
 ******************************************************************************/
typedef struct {
  int max_m;
  int max_n;
  int max_k;
  int64_t flops;
} dbm_batch_shape_t;

#endif

// EOF
