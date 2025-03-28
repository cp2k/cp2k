/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2025 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include "../offload/offload_runtime.h"
#if defined(__OFFLOAD_OPENCL) && !defined(__NO_OFFLOAD_DBM)

#include "dbm_multiply_gpu_kernel.h"
#include "dbm_multiply_opencl.cl.h"

#if defined(__DBCSR_ACC)
#include <smm/opencl_libsmm.h>
#endif

#define DBM_TIMER_DIFF(A, B) libxsmm_timer_duration(A, B)
#define DBM_TIMER_TICK() libxsmm_timer_tick()
#define DBM_TIMER_TICKINT libxsmm_timer_tickint

typedef struct {
  int max_m, max_n, max_k, mnk_count;
} dbm_multiply_gpu_launch_info_t;

static void dbm_multiply_gpu_launch_info(dbm_multiply_gpu_launch_info_t *info,
                                         const dbm_task_t *tasks, int ntasks) {
  int avg_m = tasks[0].m, avg_n = tasks[0].n, avg_k = tasks[0].k, i = 1;
  info->max_m = avg_m;
  info->max_n = avg_n;
  info->max_k = avg_k;
  for (info->mnk_count = 0; i < ntasks; ++i) {
    const int m = tasks[i].m, n = tasks[i].n, k = tasks[i].k;
    info->max_m = imax(info->max_m, m);
    info->max_n = imax(info->max_n, n);
    info->max_k = imax(info->max_k, k);
    if (m != avg_m || n != avg_n || k != avg_k) {
      avg_m = (avg_m + m) / 2;
      avg_n = (avg_n + n) / 2;
      avg_k = (avg_k + k) / 2;
      ++info->mnk_count;
    }
  }
}

void dbm_multiply_gpu_launch_kernel(const offloadStream_t stream, double alpha,
                                    int ntasks, const dbm_task_t *tasks_host,
                                    const dbm_task_t *tasks,
                                    const double *pack_a_data,
                                    const double *pack_b_data,
                                    double *shard_c_data) {
  const DBM_TIMER_TICKINT start = DBM_TIMER_TICK();
  const c_dbcsr_acc_opencl_config_t *const config = &c_dbcsr_acc_opencl_config;
#if defined(OPENCL_LIBSMM_PFORMAT)
  const char *const smm_env = getenv("DBM_MULTIPLY_SMM");
  int max_kernel_dim = (NULL == smm_env ? 0 /*default*/ : atoi(smm_env));
#endif
  const int verbosity = config->verbosity;
  int result = EXIT_SUCCESS;
  cl_event e = NULL, *const event =
                         ((0 <= verbosity && 2 >= verbosity) ? NULL : &e);
  dbm_multiply_gpu_launch_info_t info = {0};
  dbm_multiply_gpu_launch_info(&info, tasks_host, ntasks);
  assert(NULL != pack_a_data && NULL != pack_b_data && NULL != shard_c_data);
  assert(0 < ntasks && NULL != tasks);
#if defined(OPENCL_LIBSMM_PFORMAT)
  if (0 != info.mnk_count || 1 != alpha ||
      (max_kernel_dim * max_kernel_dim) < (info.max_m * info.max_n))
#endif
  {
#if defined(OPENCL_DBM_SOURCE_MULTIPLY)
    /* creating/calling kernel must be consistent across threads */
    static cl_kernel kernel_global = NULL;
    static LIBXSMM_TLS cl_kernel kernel = NULL;
    static int ndims = 1, clinear = 0;
    static size_t wgsize[] = {0, 0, 0};
    const c_dbcsr_acc_opencl_stream_t *const str = ACC_OPENCL_STREAM(stream);
    c_dbcsr_acc_opencl_info_memptr_t adata, bdata, cdata, batch;
    size_t work_size[] = {1, 1, 1}, ibatch = 0;
    size_t iadata = 0, ibdata = 0, icdata = 0;
    const size_t work_tasks = ntasks;
    assert(NULL != str && NULL != str->queue);
    if (NULL == kernel_global) { /* initial check if kernel is present */
      ACC_OPENCL_ACQUIRE(config->lock_main);
      if (NULL == kernel_global) {
        char params[ACC_OPENCL_BUFFERSIZE] =
            "-cl-fast-relaxed-math -cl-denorms-are-zero";
        const char *const gen_env = getenv("DBM_MULTIPLY_GEN");
        const char *const lin_env = getenv("DBM_MULTIPLY_LIN");
        const char *const bn_env = getenv("DBM_MULTIPLY_BN");
        const char *const sm_env = getenv("DBM_MULTIPLY_SM");
        const char *const wg_env = getenv("DBM_MULTIPLY_WG");
        const char *const lu_env = getenv("DBM_MULTIPLY_LU");
        const char *const xf_env = getenv("DBM_MULTIPLY_XF");
        const c_dbcsr_acc_opencl_device_t *const devinfo = &config->device;
        int sm = (NULL == sm_env ? 0 /*default*/ : atoi(sm_env));
        const int bn0 = (0 == devinfo->nv ? (0 == devinfo->amd ? 4 : 8) : 2);
        const int bn1 = ((0 == sm && 0 == clinear) ? bn0 : (bn0 * 2));
        int bn = LIBXSMM_CLMP(NULL == bn_env ? bn1 : atoi(bn_env), 1, 32);
        int lu = LIBXSMM_CLMP(NULL == lu_env ? 0 : atoi(lu_env), -2, 1);
        int gen = ((NULL == bn_env && NULL == sm_env && NULL == wg_env &&
                    NULL == lu_env && NULL == lin_env)
                       ? (NULL == gen_env ? 1 /*default*/ : atoi(gen_env))
                       : 0);
        const int gpu = (CL_DEVICE_TYPE_GPU == devinfo->type);
        const int xf = (NULL == xf_env ? -1 /*default*/ : atoi(xf_env));
        const char *extensions[] = {NULL, NULL}, *flags = NULL;
        size_t nextensions = sizeof(extensions) / sizeof(*extensions);
        const size_t wgsize0 = devinfo->wgsize[0], wgsize1 = devinfo->wgsize[1];
        size_t wgsize2 = devinfo->wgsize[2];
        size_t offset =
            ((0 == config->debug && 0 == config->dump) ? strlen(params) : 0);
        offset += (size_t)c_dbcsr_acc_opencl_flags_atomics(
            devinfo, c_dbcsr_acc_opencl_atomic_fp_64, extensions, &nextensions,
            params + offset, sizeof(params) - offset);
        if (2 <= gen ||
            (0 != gen && 0 != wgsize2 /*subgroups*/ &&
             2 <= *devinfo->std_level && NULL != extensions[1] &&
             NULL != strstr(extensions[1], "cl_ext_float_atomics"))) {
          offset +=
              (size_t)LIBXSMM_SNPRINTF(params + offset, sizeof(params) - offset,
                                       " -DDBM_MULTIPLY_OPENCL_GEN");
          wgsize[1] = wgsize[2] = 1;
          wgsize[0] = 16;
          lu = bn = 0;
          ndims = 3;
        } else {
          wgsize[0] = (NULL == wg_env ? (unsigned long int)LIBXSMM_ABS(sm)
                                      : strtoul(wg_env, NULL, 10));
          if (0 != wgsize2 && 0 < wgsize[0]) { /* subgroups */
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
          sm = ((0 != sm && 0 != wgsize[0])
                    ? (LIBXSMM_ISPOT(bn * sizeof(double)) + 1)
                    : 0);
          clinear = (NULL == lin_env ? 0 /*default*/ : atoi(lin_env));
          offset += (size_t)LIBXSMM_SNPRINTF(
              params + offset, sizeof(params) - offset,
              " %s %s -DBN=%i -DSM=%i -DLU=%i -DWG=%i -DSG=%i",
              0 != gpu ? "-DGPU" : "", 0 == clinear ? "" : "-DCLINEAR", bn, sm,
              lu, (int)wgsize[0], (int)wgsize2);
          gen = 0;
        }
        if (0 != devinfo->intel && 0 < xf) {
          flags = "-cl-intel-256-GRF-per-thread";
        }
        result |= (sizeof(params) > offset ? EXIT_SUCCESS : EXIT_FAILURE);
        result |= c_dbcsr_acc_opencl_kernel(
            0 /*source_is_file*/, OPENCL_DBM_SOURCE_MULTIPLY, "dbm_multiply",
            params, flags, NULL /*try*/, NULL /*try_ok*/, extensions,
            nextensions, &kernel_global);
        if (2 <= verbosity || 0 > verbosity) {
          if (EXIT_SUCCESS == result) {
            const double ds = DBM_TIMER_DIFF(start, DBM_TIMER_TICK());
            fprintf(stderr, "INFO ACC/LIBDBM: DBM-kernel gpu=%i", gpu);
            if (0 != gen) { /* generated kernel */
              fprintf(stderr, " gen=%i", gen);
            }
            if (0 != clinear) {
              fprintf(stderr, " lin=%i", clinear);
            }
            if (0 != bn) {
              fprintf(stderr, " bn=%i", bn);
            }
            if (0 != sm) {
              fprintf(stderr, " sm=%i", sm);
            }
            if (0 != wgsize[0]) {
              fprintf(stderr, " wg=%i", (int)wgsize[0]);
            }
            if (0 != wgsize2) {
              fprintf(stderr, " sg=%i", (int)wgsize2);
            }
            if (0 != lu) {
              fprintf(stderr, " lu=%i", lu);
            }
            fprintf(stderr, " ms=%.1f\n", 1E3 * ds);
          } else {
            fprintf(stderr, "INFO ACC/LIBDBM: DBM-kernel failed to generate\n");
          }
        }
      }
      kernel = clCloneKernel(kernel_global, &result); /* always clone */
      ACC_OPENCL_RELEASE(config->lock_main);
    } else if (NULL == kernel) {
      kernel = clCloneKernel(kernel_global, &result);
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
    if (1 < ndims) { /* DBM_MULTIPLY_GEN */
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
      size_t size = work_tasks;
      size *= (0 == clinear ? info.max_m : info.max_n);
      /* fixup to be a multiple of the WG-size */
      work_size[0] = (0 < wgsize[0] ? LIBXSMM_UP(size, wgsize[0]) : size);
      result |= clSetKernelArg(kernel, 2, sizeof(cl_int), &ntasks);
      result |= clSetKernelArg(kernel, 3, sizeof(cl_int), &size);
      result |= c_dbcsr_acc_opencl_set_kernel_ptr(kernel, 4, batch.memory);
      result |= c_dbcsr_acc_opencl_set_kernel_ptr(kernel, 5, adata.memory);
      result |= c_dbcsr_acc_opencl_set_kernel_ptr(kernel, 6, bdata.memory);
      result |= c_dbcsr_acc_opencl_set_kernel_ptr(kernel, 7, cdata.memory);
    }
    result |= clEnqueueNDRangeKernel(str->queue, kernel, ndims, NULL, work_size,
                                     0 < wgsize[0] ? wgsize : NULL,
                                     0 /*num_wait*/, NULL /*wait_list*/, event);
  }
#if defined(OPENCL_LIBSMM_PFORMAT)
  else { /* homogeneous */
    const int pzero = 0, pbase = 3, pnext = 6;
    const int param_format = pzero | (pbase << 8) | (pnext << 16);
    result |= opencl_libsmm_acc_process(
        NULL /*tasks_host*/, &tasks->m, ntasks, dbcsr_type_real_8, pack_a_data,
        pack_b_data, shard_c_data, info.max_m, info.max_n, info.max_k,
        max_kernel_dim, 1 /*homogeneous*/, stream, NULL /*c_stream*/,
        param_format, event);
  }
#endif
  if (NULL != event && NULL != *event && EXIT_SUCCESS == result &&
      EXIT_SUCCESS == clWaitForEvents(1, event)) {
    static LIBXSMM_TLS DBM_TIMER_TICKINT start2 = 0;
    const DBM_TIMER_TICKINT stop = DBM_TIMER_TICK();
    const double dhost = DBM_TIMER_DIFF(start, stop);
    const double diter = (0 < start2 ? DBM_TIMER_DIFF(start, start2) : dhost);
#if defined(OPENCL_LIBSMM_PFORMAT)
    const char *const kind = (0 == max_kernel_dim ? "DBM" : "SMM");
#else
    const char *const kind = "DBM";
#endif
    const int pure = (100 * (ntasks - info.mnk_count) + ntasks - 1) / ntasks;
    double dkrnl = dhost, dtotl;
    if (c_dbcsr_acc_opencl_timer_host == config->timer) {
      cl_ulong begin = 0, end = 0;
      const int r0 = clGetEventProfilingInfo(*event, CL_PROFILING_COMMAND_START,
                                             sizeof(cl_ulong), &begin, NULL);
      const int r1 = clGetEventProfilingInfo(*event, CL_PROFILING_COMMAND_END,
                                             sizeof(cl_ulong), &end, NULL);
      if (EXIT_SUCCESS == r0 && EXIT_SUCCESS == r1) {
        dkrnl = 1E-9 * LIBXSMM_DELTA(begin, end);
      }
    }
    start2 = stop;
    dtotl = LIBXSMM_MIN(LIBXSMM_MIN(diter, dhost), dkrnl);
    fprintf(stderr,
            "INFO ACC/LIBDBM: %s-kernel mnk=%ix%ix%i pure=%i%% ntasks=%i "
            "ims=%.1f hms=%.1f kms=%.1f gflops=%.1f\n",
            kind, info.max_m, info.max_n, info.max_k, pure, ntasks,
            1E+3 * diter, 1E+3 * dhost, 1E+3 * dkrnl,
            1E-9 * info.max_m * info.max_n * info.max_k * ntasks / dtotl);
  }
  OFFLOAD_CHECK(result);
}

#endif // defined(__OFFLOAD_OPENCL) && !defined(__NO_OFFLOAD_DBM)

// EOF
