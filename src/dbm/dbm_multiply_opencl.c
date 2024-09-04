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
  static int ndims = 1, split = 0;
  static size_t wgsize[] = {0, 0, 0};
  int result = EXIT_SUCCESS, verbosity = c_dbcsr_acc_opencl_config.verbosity;
  cl_event event, *const perf_event =
                      ((0 <= verbosity && 2 >= verbosity) ? NULL : &event);
  const c_dbcsr_acc_opencl_stream_t *const str = ACC_OPENCL_STREAM(stream);
  const size_t max_m = mnk_range[0][1], work_tasks = ntasks;
  size_t work_size[] = {1, 1, 1}, ibatch = 0, iadata = 0, ibdata = 0,
         icdata = 0;
  c_dbcsr_acc_opencl_info_memptr_t adata, bdata, cdata, batch;
  assert(NULL != pack_a_data && NULL != pack_b_data && NULL != shard_c_data);
  assert(0 < mnk_range[0][0] && 0 < mnk_range[0][1] &&
         mnk_range[0][0] <= mnk_range[0][1]);
  assert(0 < mnk_range[1][0] && 0 < mnk_range[1][1] &&
         mnk_range[1][0] <= mnk_range[1][1]);
  assert(0 < mnk_range[2][0] && 0 < mnk_range[2][1] &&
         mnk_range[2][0] <= mnk_range[2][1]);
  assert(NULL != str && NULL != str->queue);
  assert(0 < ntasks && NULL != tasks);
  /* creating/calling kernel must be consistent across threads */
  ACC_OPENCL_ACQUIRE(c_dbcsr_acc_opencl_config.lock_main);
#if defined(OPENCL_DBM_SOURCE_MULTIPLY)
  if (NULL == kernel) { /* first-time check if kernel is present */
    const libxsmm_timer_tickint start = libxsmm_timer_tick();
    char params[ACC_OPENCL_BUFFERSIZE] =
        "-cl-fast-relaxed-math -cl-denorms-are-zero";
    const char *const gen_env = getenv("DBM_MULTIPLY_GEN");
    const char *const xf_env = getenv("DBM_MULTIPLY_XF");
    const char *const lu_env = getenv("DBM_MULTIPLY_LU");
    const char *const bn_env = getenv("DBM_MULTIPLY_BN");
    const int gpu =
        (CL_DEVICE_TYPE_GPU == c_dbcsr_acc_opencl_config.device.type);
    const int gen = (NULL == gen_env ? 0 /*default*/ : atoi(gen_env));
    const int xf = (NULL == xf_env ? -1 /*default*/ : atoi(xf_env));
    const int lu = LIBXSMM_CLMP(NULL == lu_env ? 0 : atoi(lu_env), -2, 1);
    int bn = (NULL == bn_env ? 8 : atoi(bn_env));
    const char *extensions[] = {NULL, NULL}, *flags = NULL;
    size_t nextensions = sizeof(extensions) / sizeof(*extensions);
    const size_t wgsize0 = c_dbcsr_acc_opencl_config.device.wgsize[0];
    const size_t wgsize1 = c_dbcsr_acc_opencl_config.device.wgsize[1];
    size_t wgsize2 = c_dbcsr_acc_opencl_config.device.wgsize[2];
    size_t offset = (0 == c_dbcsr_acc_opencl_config.debug ? strlen(params) : 0);
    offset += (size_t)c_dbcsr_acc_opencl_flags_atomics(
        &c_dbcsr_acc_opencl_config.device, c_dbcsr_acc_opencl_atomic_fp_64,
        extensions, &nextensions, params + offset, sizeof(params) - offset);
    if (2 <= gen || (0 != gen && 0 != wgsize2 /*subgroups*/ &&
                     2 <= *c_dbcsr_acc_opencl_config.device.std_level &&
                     NULL != extensions[1] &&
                     NULL != strstr(extensions[1], "cl_ext_float_atomics"))) {
      offset +=
          (size_t)LIBXSMM_SNPRINTF(params + offset, sizeof(params) - offset,
                                   " -DDBM_MULTIPLY_OPENCL_GEN");
      wgsize[1] = wgsize[2] = 1;
      wgsize[0] = 16;
      ndims = 3;
    } else {
      const char *const split_env = getenv("DBM_MULTIPLY_SPLIT");
      const char *const wg_env = getenv("DBM_MULTIPLY_WG");
      split = (NULL == split_env ? 1 /*default*/ : atoi(split_env));
      wgsize[0] =
          (NULL == wg_env ? (1 != split ? (wgsize1 * LIBXSMM_ABS(split)) : 0)
                          : strtoul(wg_env, NULL, 10));
      if (0 != split && 1 != split && (bn * bn) > (int)wgsize[0]) {
        wgsize[0] = bn * bn;
      }
      if (0 != split && 0 != wgsize2 && 0 < wgsize[0]) { /* subgroups */
        if (LIBXSMM_DELTA(wgsize[0], wgsize1) <=
            LIBXSMM_DELTA(wgsize[0], wgsize2)) { /* select SG-size */
          wgsize2 = wgsize1;
        }
        wgsize[0] = LIBXSMM_UP(wgsize[0], wgsize2);
      } else {
        wgsize[0] = LIBXSMM_UP(wgsize[0], wgsize1);
        wgsize2 = 0;
      }
      wgsize[0] = LIBXSMM_CLMP(wgsize[0], 0, wgsize0);
      if (NULL == bn_env && 0 != split && 1 != split &&
          (bn * bn) < (int)wgsize[0]) {
        bn = libxsmm_isqrt2_u32(wgsize[0]);
      }
      bn = LIBXSMM_CLMP(bn, 4, 32);
      offset += (size_t)LIBXSMM_SNPRINTF(
          params + offset, sizeof(params) - offset,
          " %s -DSPLIT=%i -DBN=%i -DWG=%i -DSG=%i -DLU=%i",
          0 != gpu ? "-DGPU" : "", split, bn, (int)wgsize[0], (int)wgsize2, lu);
    }
    if (0 != c_dbcsr_acc_opencl_config.device.intel && 0 < xf) {
      flags = "-cl-intel-256-GRF-per-thread";
    }
    result |= (sizeof(params) > offset ? EXIT_SUCCESS : EXIT_FAILURE);
    result |= c_dbcsr_acc_opencl_kernel(
        0 /*source_is_file*/, OPENCL_DBM_SOURCE_MULTIPLY, "dbm_multiply",
        params, flags, NULL /*try*/, NULL /*try_ok*/, extensions, nextensions,
        &kernel);
    if (2 <= verbosity || 0 > verbosity) {
      if (EXIT_SUCCESS == result) {
        const double d = libxsmm_timer_duration(start, libxsmm_timer_tick());
        fprintf(stderr, "INFO ACC/LIBDBM: DBM-kernel gpu=%i", gpu);
        if (0 == gen) {
          fprintf(stderr, " split=%i lu=%i bn=%i", split, lu, bn);
        } else { /* generated kernel */
          fprintf(stderr, " gen=%i", gen);
        }
        fprintf(stderr, " wg=%i sg=%i ms=%.1f\n", (int)wgsize[0], (int)wgsize2,
                1E3 * d);
      } else {
        fprintf(stderr, "INFO ACC/LIBDBM: DBM-kernel failed to generate\n");
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
  result |= clSetKernelArg(kernel, 0, sizeof(cl_double), &alpha);
  result |= clSetKernelArg(kernel, 1, sizeof(cl_int), &ibatch);
  if (1 < ndims) { /* generated kernel */
    const cl_uint zero = 0;
    assert(0 != wgsize[1] && 0 != wgsize[1] && 0 != wgsize[2]);
    work_size[0] = 16;
    assert(1 == work_size[1]);
    work_size[2] = work_tasks;
    result |= c_dbcsr_acc_opencl_set_kernel_ptr(kernel, 2, batch.memory);
    result |= clSetKernelArg(kernel, 3, sizeof(cl_uint), &zero /*shape*/);
    result |= c_dbcsr_acc_opencl_set_kernel_ptr(kernel, 4, adata.memory);
    result |= clSetKernelArg(kernel, 5, sizeof(cl_uint), &zero /*A_shape0*/);
    result |= c_dbcsr_acc_opencl_set_kernel_ptr(kernel, 6, bdata.memory);
    result |= clSetKernelArg(kernel, 7, sizeof(cl_uint), &zero /*B_shape0*/);
    result |= c_dbcsr_acc_opencl_set_kernel_ptr(kernel, 8, cdata.memory);
    result |= clSetKernelArg(kernel, 9, sizeof(cl_uint), &zero /*C_shape0*/);
  } else {
    result |= clSetKernelArg(kernel, 2, sizeof(cl_int), &ntasks);
    if (0 != split) {
      if (1 == split || 0 == wgsize[0]) {
        work_size[0] = work_tasks * max_m;
        result |= clSetKernelArg(kernel, 3, sizeof(cl_int), work_size);
        if (0 < wgsize[0]) { /* fixup to be a multiple of the WG-size */
          work_size[0] = LIBXSMM_UP(work_size[0], wgsize[0]);
        }
      } else {
        work_size[0] = work_tasks * wgsize[0];
        result |= clSetKernelArg(kernel, 3, sizeof(cl_int), work_size);
      }
    } else {
      work_size[0] = work_tasks;
      result |= clSetKernelArg(kernel, 3, sizeof(cl_int), work_size);
    }
    result |= c_dbcsr_acc_opencl_set_kernel_ptr(kernel, 4, batch.memory);
    result |= c_dbcsr_acc_opencl_set_kernel_ptr(kernel, 5, adata.memory);
    result |= c_dbcsr_acc_opencl_set_kernel_ptr(kernel, 6, bdata.memory);
    result |= c_dbcsr_acc_opencl_set_kernel_ptr(kernel, 7, cdata.memory);
  }
  result |= clEnqueueNDRangeKernel(
      str->queue, kernel, ndims, NULL, work_size, 0 < wgsize[0] ? wgsize : NULL,
      0 /*num_wait*/, NULL /*wait_list*/, perf_event);
  if (NULL != perf_event && EXIT_SUCCESS == result) {
    cl_ulong begin = 0, end = 0;
    clWaitForEvents(1, perf_event);
    result |= clGetEventProfilingInfo(*perf_event, CL_PROFILING_COMMAND_START,
                                      sizeof(cl_ulong), &begin, NULL);
    result |= clGetEventProfilingInfo(*perf_event, CL_PROFILING_COMMAND_END,
                                      sizeof(cl_ulong), &end, NULL);
    if (EXIT_SUCCESS == result) {
      const double duration_ns = LIBXSMM_DELTA(begin, end);
      const double gflops =
          (max_m * mnk_range[1][1] * mnk_range[2][1] * ntasks) / duration_ns;
      fprintf(stderr,
              "INFO ACC/LIBDBM: DBM-kernel mnk=%ix%ix%i "
              "ntasks=%i gflops=%.1f ms=%.2g\n",
              mnk_range[0][1], mnk_range[1][1], mnk_range[2][1], ntasks, gflops,
              1E-6 * duration_ns);
    }
  }
  ACC_OPENCL_RELEASE(c_dbcsr_acc_opencl_config.lock_main);
  OFFLOAD_CHECK(result);
}

#endif // defined(__OFFLOAD_OPENCL) && !defined(__NO_OFFLOAD_DBM)

// EOF
