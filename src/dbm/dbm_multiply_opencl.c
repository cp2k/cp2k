/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#include "../offload/offload_runtime.h"
#if defined(__OFFLOAD_OPENCL) && !defined(__NO_OFFLOAD_DBM)

#include "dbm_multiply_gpu_kernel.h"

/* OpenCL backend is maintained in LIBXSTREAM (samples/dbm) */
#include <dbm_opencl.h>

_Static_assert(sizeof(dbm_task_t) == DBM_OPENCL_TASK_SIZE * sizeof(int),
               "dbm_task_t must match the native task format of LIBXSTREAM");

void dbm_multiply_gpu_launch_kernel(offloadStream_t stream, double alpha,
                                    int ntasks, const dbm_batch_shape_t *shape,
                                    const dbm_task_t *tasks_host,
                                    const dbm_task_t *tasks,
                                    const double *pack_a_data,
                                    const double *pack_b_data,
                                    double *shard_c_data) {
  // Fill of a kernel padded to the maxima in percent, 100 if homogeneous.
  const int64_t padded =
      2LL * ntasks * shape->max_m * shape->max_n * shape->max_k;
  const int shape_info[] = {shape->max_m, shape->max_n, shape->max_k,
                            (int)(100 * shape->flops / padded)};
  const int result = dbm_multiply_opencl_launch_kernel(
      stream, alpha, ntasks, 0 /*param_format*/, shape_info, &tasks_host->m,
      &tasks->m, pack_a_data, pack_b_data, shard_c_data);
  OFFLOAD_CHECK(result);
}

int dbm_multiply_gpu_task_order(void) {
  return dbm_multiply_opencl_task_order();
}

#endif /* defined(__OFFLOAD_OPENCL) && !defined(__NO_OFFLOAD_DBM) */

/* EOF */
