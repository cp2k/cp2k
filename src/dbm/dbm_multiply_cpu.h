/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#ifndef DBM_MULTIPLY_CPU_H
#define DBM_MULTIPLY_CPU_H

#include "dbm_internal.h"
#include "dbm_shard.h"

enum dbm_multiply_cpu_options {
  DBM_MULTIPLY_BLAS_LIBRARY = 1,
  DBM_MULTIPLY_TASK_REORDER = 2
};

// Approximate order of the tasks within a batch (bucket sort).
enum dbm_multiply_task_order {
  DBM_TASK_ORDER_NONE = 0, // as generated
  DBM_TASK_ORDER_SHAPE = 1 // by m,n,k
};

/*******************************************************************************
 * \brief Internal routine for ordering the tasks of a batch approximately by
 *        shape. The order is returned as permutation of the task indexes.
 * \author Hans Pabst
 ******************************************************************************/
void dbm_multiply_cpu_task_order(int ntasks, const dbm_task_t batch[ntasks],
                                 int order[ntasks]);

/*******************************************************************************
 * \brief Internal routine for executing the tasks in given batch on the CPU.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply_cpu_process_batch(int ntasks, const dbm_task_t batch[ntasks],
                                    double alpha, const dbm_pack_t *pack_a,
                                    const dbm_pack_t *pack_b,
                                    dbm_shard_t *shard_c, int options);

#endif

// EOF
