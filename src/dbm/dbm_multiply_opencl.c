/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2024 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include "../offload/offload_runtime.h"
#if defined(__OFFLOAD_OPENCL) && !defined(__NO_OFFLOAD_DBM)

#include "dbm_multiply_gpu_kernel.h"
#include "dbm_multiply_opencl.cl.h"

void dbm_multiply_gpu_launch_kernel(const offloadStream_t stream,
                                    const int mnk_range[3][2], double alpha,
                                    int ntasks, const dbm_task_t *tasks,
                                    const double *pack_a_data,
                                    const double *pack_b_data,
                                    double *shard_c_data) {
  static cl_kernel kernel = NULL;
  static int task_split = 0;
  int result = EXIT_SUCCESS, verbosity = c_dbcsr_acc_opencl_config.verbosity;
  cl_event event, *const perf_event =
                      ((0 <= verbosity && 2 >= verbosity) ? NULL : &event);
  const c_dbcsr_acc_opencl_stream_t *const str = ACC_OPENCL_STREAM(stream);
  const int max_m = mnk_range[0][1], max_n = mnk_range[1][1];
  const size_t work_tasks = ntasks, wgsize = 0;
  size_t ibatch = 0, iadata = 0, ibdata = 0, icdata = 0, work_size = 0;
  c_dbcsr_acc_opencl_info_memptr_t adata, bdata, cdata, batch;
  assert(NULL != pack_a_data && NULL != pack_b_data && NULL != shard_c_data);
  assert(0 < mnk_range[0][0] && 0 < max_m && mnk_range[0][0] <= max_m);
  assert(0 < mnk_range[1][0] && 0 < max_n && mnk_range[1][0] <= max_n);
  assert(0 < mnk_range[2][0] && 0 < mnk_range[2][1] &&
         mnk_range[2][0] <= mnk_range[2][1]);
  assert(NULL != str && NULL != str->queue);
  assert(0 < ntasks && NULL != tasks);
  /* creating/calling kernel must be consistent across threads */
  ACC_OPENCL_ACQUIRE(c_dbcsr_acc_opencl_config.lock_main);
#if defined(OPENCL_DBM_SOURCE_MULTIPLY_OPENCL)
  if (NULL == kernel) { /* first-time check if kernel is present */
    char params[ACC_OPENCL_BUFFERSIZE] = "-DLU=0 -DBN=4";
    const char *const flags = "-cl-fast-relaxed-math -cl-denorms-are-zero";
    const char *extensions[] = {NULL, NULL};
    const size_t nextensions = sizeof(extensions) / sizeof(*extensions);
    const size_t offset = strlen(params);
    const libxsmm_timer_tickint start = libxsmm_timer_tick();
    const int nchar = c_dbcsr_acc_opencl_flags_atomics(
        &c_dbcsr_acc_opencl_config.device, c_dbcsr_acc_opencl_atomic_fp_64,
        extensions, nextensions, params + offset, sizeof(params) - offset);
    if (0 < nchar && sizeof(params) > (offset + nchar)) {
      const char *const task_split_env = getenv("LIBDBM_TASK_SPLIT");
      task_split = (NULL == task_split_env ? 1 : atoi(task_split_env));
      result |= c_dbcsr_acc_opencl_kernel(
          0 /*source_is_file*/, OPENCL_DBM_SOURCE_MULTIPLY_OPENCL,
          "dbm_multiply", params,
          0 == c_dbcsr_acc_opencl_config.debug ? flags : NULL, NULL /*try*/,
          NULL /*try_ok*/, extensions, nextensions, &kernel);
      if (2 <= verbosity || 0 > verbosity) {
        if (EXIT_SUCCESS == result) {
          fprintf(stderr, "INFO ACC/LIBDBM: DBM-kernel ms=%.1f\n",
                  1E3 * libxsmm_timer_duration(start, libxsmm_timer_tick()));
        } else {
          fprintf(stderr, "INFO ACC/LIBDBM: DBM-kernel failed to generate\n");
        }
      }
    }
  }
#else
#error "OpenCL kernel code not found!"
#endif
  result |= c_dbcsr_acc_opencl_info_devptr_lock(&adata, NULL /*lock*/,
                                                pack_a_data, 1 /*esize*/,
                                                NULL /*amount*/, &iadata);
  result |= c_dbcsr_acc_opencl_info_devptr_lock(&bdata, NULL /*lock*/,
                                                pack_b_data, 1 /*esize*/,
                                                NULL /*amount*/, &ibdata);
  result |= c_dbcsr_acc_opencl_info_devptr_lock(&cdata, NULL /*lock*/,
                                                shard_c_data, 1 /*esize*/,
                                                NULL /*amount*/, &icdata);
  result |= c_dbcsr_acc_opencl_info_devptr_lock(
      &batch, NULL /*lock*/, tasks /*batch*/, sizeof(dbm_task_t), &work_tasks,
      &ibatch);
  assert(0 == iadata && 0 == ibdata && 0 == icdata);
  work_size = (0 != task_split ? (work_tasks * max_m) : work_tasks);
  result |= clSetKernelArg(kernel, 0, sizeof(cl_double), &alpha);
  result |= clSetKernelArg(kernel, 1, sizeof(cl_int), &ibatch);
  result |= clSetKernelArg(kernel, 2, sizeof(cl_int), &ntasks);
  result |= c_dbcsr_acc_opencl_set_kernel_ptr(kernel, 3, batch.memory);
  result |= c_dbcsr_acc_opencl_set_kernel_ptr(kernel, 4, adata.memory);
  result |= c_dbcsr_acc_opencl_set_kernel_ptr(kernel, 5, bdata.memory);
  result |= c_dbcsr_acc_opencl_set_kernel_ptr(kernel, 6, cdata.memory);
  result |= clEnqueueNDRangeKernel(str->queue, kernel, 1 /*work_dim*/,
                                   NULL /*offset*/, &work_size,
                                   0 != wgsize ? &wgsize : NULL, 0 /*num_wait*/,
                                   NULL /*wait_list*/, perf_event);
  if (NULL != perf_event && EXIT_SUCCESS == result) {
    cl_ulong begin = 0, end = 0;
    clWaitForEvents(1, perf_event);
    result |= clGetEventProfilingInfo(*perf_event, CL_PROFILING_COMMAND_START,
                                      sizeof(cl_ulong), &begin, NULL);
    result |= clGetEventProfilingInfo(*perf_event, CL_PROFILING_COMMAND_END,
                                      sizeof(cl_ulong), &end, NULL);
    if (EXIT_SUCCESS == result) {
      const double duration_ms = 1E-6 * LIBXSMM_DELTA(begin, end);
      const double gflops =
          (1E-6 * max_m * max_n * mnk_range[2][1] * ntasks) / duration_ms;
      fprintf(stderr,
              "INFO ACC/LIBDBM: DBM-kernel mnk=%ix%ix%i "
              "ntasks=%i gflops=%.1f ms=%.2g\n",
              max_m, max_n, mnk_range[2][1], ntasks, gflops, duration_ms);
    }
  }
  ACC_OPENCL_RELEASE(c_dbcsr_acc_opencl_config.lock_main);
  OFFLOAD_CHECK(result);
}

#endif // defined(__OFFLOAD_OPENCL) && !defined(__NO_OFFLOAD_DBM)

// EOF
