/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#ifndef DBM_MULTIPLY_CPU_H
#define DBM_MULTIPLY_CPU_H

#include <stdbool.h>

#include "dbm_multiply_internal.h"
#include "dbm_shard.h"

/*******************************************************************************
 * \brief Internal routine for executing the tasks in given batch on the CPU.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply_cpu_process_batch(const int ntasks, dbm_task_t batch[ntasks],
                                    const bool transa, const bool transb,
                                    const double alpha,
                                    const dbm_pack_t *pack_a,
                                    const dbm_pack_t *pack_b,
                                    dbm_shard_t *shard_c);

#endif

// EOF
