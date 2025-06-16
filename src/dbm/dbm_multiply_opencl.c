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

int dbm_multiply_opencl_launch_kernel(void *stream, double alpha, int ntasks,
                                      int param_format, const int *params_host,
                                      const int *params,
                                      const double *pack_a_data,
                                      const double *pack_b_data,
                                      double *shard_c_data);

#if defined(OPENCL_LIBSMM_PFORMAT) && (0 < OPENCL_LIBSMM_PFORMAT)
int dbm_multiply_opencl_initialized /*= 0*/;
int dbm_multiply_opencl_smm /*= 0*/;

LIBXSMM_ATTRIBUTE_CTOR static void dbm_multiply_opencl_initialize(void) {
  const char *const smm_env = getenv("DBM_MULTIPLY_SMM");
  const int smm = (NULL == smm_env ? 0 /*default*/ : atoi(smm_env));
  dbm_multiply_opencl_smm =
      LIBXSMM_MIN(1 != smm ? smm : 64, (1 << (OPENCL_LIBSMM_PFORMAT - 1)) - 1);
  if (0 > dbm_multiply_opencl_smm) {
    opencl_libsmm_acc_set_dbm_launch_fn(dbm_multiply_opencl_launch_kernel);
  }
  ++dbm_multiply_opencl_initialized;
}
#endif

typedef struct {
  int max_m, max_n, max_k, mnk_changes;
} dbm_multiply_gpu_launch_info_t;

static void dbm_multiply_gpu_launch_info(dbm_multiply_gpu_launch_info_t *info,
                                         const int *params, int ntasks,
                                         int param_format) {
  if (0 == param_format) { /* native */
    const int stride = sizeof(dbm_task_t) / sizeof(int);
    int avg_m = params[0], avg_n = params[1], avg_k = params[2], i = stride;
    info->max_m = avg_m;
    info->max_n = avg_n;
    info->max_k = avg_k;
    for (info->mnk_changes = 0; i < (ntasks * stride); i += stride) {
      const int m = params[i + 0], n = params[i + 1], k = params[i + 2];
      info->max_m = imax(info->max_m, m);
      info->max_n = imax(info->max_n, n);
      info->max_k = imax(info->max_k, k);
      if (m != avg_m || n != avg_n || k != avg_k) { /* approximation */
        avg_m = (avg_m + m) / 2;
        avg_n = (avg_n + n) / 2;
        avg_k = (avg_k + k) / 2;
        ++info->mnk_changes;
      }
    }
  } else {
#if defined(OPENCL_LIBSMM_PFORMAT) && (0 < OPENCL_LIBSMM_PFORMAT)
    const int mask = (1 << OPENCL_LIBSMM_PFORMAT) - 1;
    info->max_m = mask & (param_format);
    info->max_n = mask & (param_format >> (OPENCL_LIBSMM_PFORMAT));
    info->max_k = mask & (param_format >> (OPENCL_LIBSMM_PFORMAT * 2));
    info->mnk_changes = 0; /* homogeneous */
#else
    assert(0);
#endif
  }
}

static void dbm_multiply_opencl_print(FILE *stream, const char *name, int val) {
  if (0 != val) {
    fprintf(stream, " %s=%i", name, val);
  }
}

int dbm_multiply_opencl_launch_kernel(void *stream, double alpha, int ntasks,
                                      int param_format, const int *params_host,
                                      const int *params,
                                      const double *pack_a_data,
                                      const double *pack_b_data,
                                      double *shard_c_data) {
  const DBM_TIMER_TICKINT start = DBM_TIMER_TICK();
  const c_dbcsr_acc_opencl_config_t *const config = &c_dbcsr_acc_opencl_config;
  const int verbosity = config->verbosity,
            info = (0 > verbosity || 2 < verbosity);
  int result = EXIT_SUCCESS;
  dbm_multiply_gpu_launch_info_t task = {0};
  assert(NULL != pack_a_data && NULL != pack_b_data && NULL != shard_c_data);
  assert(NULL != params_host || 0 == ntasks);
  assert(NULL != params || 0 == ntasks);
  if (0 == ntasks) {
    return result;
  }
#if defined(OPENCL_LIBSMM_PFORMAT) && (0 < OPENCL_LIBSMM_PFORMAT)
  if (0 == dbm_multiply_opencl_initialized) {
    dbm_multiply_opencl_initialize();
  }
  if (0 != dbm_multiply_opencl_smm || 0 != info) {
    dbm_multiply_gpu_launch_info(&task, params_host, ntasks, param_format);
  }
  if (0 > dbm_multiply_opencl_smm || dbm_multiply_opencl_smm < task.max_m ||
      dbm_multiply_opencl_smm < task.max_n ||
      dbm_multiply_opencl_smm < task.max_k || 0 == task.max_k || 1 != alpha)
#endif
  {
#if defined(OPENCL_DBM_SOURCE_MULTIPLY)
    /* creating/calling kernel must be consistent across threads */
    static cl_kernel kernel_global = NULL;
    static LIBXSMM_TLS cl_kernel kernel = NULL;
    static int ndims = 1, clinear = 0;
    static size_t wgsize[] = {0, 0, 0};
    const c_dbcsr_acc_opencl_stream_t *const str = ACC_OPENCL_STREAM(stream);
    const c_dbcsr_acc_opencl_device_t *const devinfo = &config->device;
    ACC_OPENCL_LOCKTYPE *const lock_memory =
        (NULL != devinfo->clSetKernelArgMemPointerINTEL ? NULL
                                                        : config->lock_memory);
    c_dbcsr_acc_opencl_info_memptr_t adata, bdata, cdata, batch;
    const int stride = (0 == param_format ? 6 : 3);
    size_t work_size[] = {1, 1, 1}, ibatch = 0;
    size_t iadata = 0, ibdata = 0, icdata = 0;
    const size_t work_tasks = ntasks;
    assert(NULL != str && NULL != str->queue);
    if (NULL == kernel_global) { /* initial check if kernel is present */
      ACC_OPENCL_ACQUIRE(config->lock_main);
      if (NULL == kernel_global) {
        char flags[ACC_OPENCL_BUFFERSIZE] =
            "-cl-fast-relaxed-math -cl-denorms-are-zero";
        const char *const gen_env = getenv("DBM_MULTIPLY_GEN");
        const char *const lin_env = getenv("DBM_MULTIPLY_LIN");
        const char *const bn_env = getenv("DBM_MULTIPLY_BN");
        const char *const sm_env = getenv("DBM_MULTIPLY_SM");
        const char *const wg_env = getenv("DBM_MULTIPLY_WG");
        const char *const lu_env = getenv("DBM_MULTIPLY_LU");
        const char *const xf_env = getenv("DBM_MULTIPLY_XF");
        int sm = (NULL == sm_env ? 0 /*default*/ : atoi(sm_env));
        const int bn0 = (0 == devinfo->nv ? (0 == devinfo->amd ? 4 : 8) : 2);
        const int bn1 = ((0 == sm && 0 == clinear) ? bn0 : (bn0 * 2));
        int bn = LIBXSMM_CLMP(NULL == bn_env ? bn1 : atoi(bn_env), 1, 32);
        int lu = LIBXSMM_CLMP(NULL == lu_env ? 0 : atoi(lu_env), -2, 1);
        int gen = ((NULL == bn_env && NULL == sm_env && NULL == wg_env &&
                    NULL == lu_env && NULL == lin_env && 0 == param_format)
                       ? (NULL == gen_env ? 1 /*default*/ : atoi(gen_env))
                       : 0);
        const int gpu = (CL_DEVICE_TYPE_GPU == devinfo->type);
        const int xf = (NULL == xf_env ? -1 /*default*/ : atoi(xf_env));
        const char *extensions[] = {NULL, NULL}, *options = NULL;
        size_t nextensions = sizeof(extensions) / sizeof(*extensions);
        const size_t wgsize0 = devinfo->wgsize[0], wgsize1 = devinfo->wgsize[1];
        size_t wgsize2 = devinfo->wgsize[2];
        size_t offset =
            ((0 == config->debug && 0 == config->dump) ? strlen(flags) : 0);
        offset += (size_t)c_dbcsr_acc_opencl_flags_atomics(
            devinfo, c_dbcsr_acc_opencl_atomic_fp_64, extensions, &nextensions,
            flags + offset, sizeof(flags) - offset);
        if (2 <= gen ||
            (0 != gen && 0 != wgsize2 /*subgroups*/ &&
             2 <= *devinfo->std_level && NULL != extensions[1] &&
             NULL != strstr(extensions[1], "cl_ext_float_atomics"))) {
          offset +=
              (size_t)LIBXSMM_SNPRINTF(flags + offset, sizeof(flags) - offset,
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
              flags + offset, sizeof(flags) - offset,
              " %s %s -DBN=%i -DSM=%i -DLU=%i -DWG=%i -DSG=%i",
              0 != gpu ? "-DGPU" : "", 0 == clinear ? "" : "-DCLINEAR", bn, sm,
              lu, (int)wgsize[0], (int)wgsize2);
          gen = 0;
        }
        if (0 != devinfo->intel && 0 < xf) {
          options = "-cl-intel-256-GRF-per-thread";
        }
        result |= (sizeof(flags) > offset ? EXIT_SUCCESS : EXIT_FAILURE);
        result |= c_dbcsr_acc_opencl_kernel(
            0 /*source_is_file*/, OPENCL_DBM_SOURCE_MULTIPLY, "dbm_multiply",
            flags, options, NULL /*try*/, NULL /*try_ok*/, extensions,
            nextensions, &kernel_global);
        if (2 <= verbosity || 0 > verbosity) {
          if (EXIT_SUCCESS == result) {
            const double ds = DBM_TIMER_DIFF(start, DBM_TIMER_TICK());
            fprintf(stderr, "INFO ACC/LIBDBM: DBM-kernel gpu=%i", gpu);
            dbm_multiply_opencl_print(stderr, "gen", gen); /* generated */
            dbm_multiply_opencl_print(stderr, "lin", clinear);
            dbm_multiply_opencl_print(stderr, "bn", bn);
            dbm_multiply_opencl_print(stderr, "sm", sm);
            dbm_multiply_opencl_print(stderr, "wg", (int)wgsize[0]);
            dbm_multiply_opencl_print(stderr, "sg", (int)wgsize2);
            dbm_multiply_opencl_print(stderr, "lu", lu);
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
    if (NULL != lock_memory) {
      ACC_OPENCL_ACQUIRE(lock_memory);
    }
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
        &batch, NULL /*lock*/, params /*batch*/, sizeof(int) * stride,
        &work_tasks, &ibatch);
    if (NULL != lock_memory) {
      ACC_OPENCL_RELEASE(lock_memory);
    }
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
#if !(defined(OPENCL_LIBSMM_PFORMAT) && (0 < OPENCL_LIBSMM_PFORMAT))
      if (0 != info) {
        dbm_multiply_gpu_launch_info(&task, params_host, ntasks, param_format);
      }
#endif
    } else {
      size_t size = work_tasks;
#if defined(OPENCL_LIBSMM_PFORMAT) && (0 < OPENCL_LIBSMM_PFORMAT)
      if (0 == dbm_multiply_opencl_smm && 0 == info)
#endif
      {
        dbm_multiply_gpu_launch_info(&task, params_host, ntasks, param_format);
      }
      size *= (0 == clinear ? task.max_m : task.max_n);
      /* fixup to be a multiple of the WG-size */
      work_size[0] = (0 < wgsize[0] ? LIBXSMM_UP(size, wgsize[0]) : size);
      result |= clSetKernelArg(kernel, 2, sizeof(cl_int), &ntasks);
      result |= clSetKernelArg(kernel, 3, sizeof(cl_int), &size);
      result |= clSetKernelArg(kernel, 4, sizeof(cl_int), &param_format);
      result |= c_dbcsr_acc_opencl_set_kernel_ptr(kernel, 5, batch.memory);
      result |= c_dbcsr_acc_opencl_set_kernel_ptr(kernel, 6, adata.memory);
      result |= c_dbcsr_acc_opencl_set_kernel_ptr(kernel, 7, bdata.memory);
      result |= c_dbcsr_acc_opencl_set_kernel_ptr(kernel, 8, cdata.memory);
    }
    result |= clEnqueueNDRangeKernel(str->queue, kernel, ndims, NULL, work_size,
                                     0 < wgsize[0] ? wgsize : NULL,
                                     0 /*num_wait*/, NULL /*wait_list*/, NULL);
  }
#if defined(OPENCL_LIBSMM_PFORMAT) && (0 < OPENCL_LIBSMM_PFORMAT)
  else { /* homogeneous */
    result |= opencl_libsmm_acc_process(
        params_host, params, ntasks, dbcsr_type_real_8, pack_a_data,
        pack_b_data, shard_c_data, task.max_m, task.max_n, task.max_k,
        dbm_multiply_opencl_smm, 1 /*homogeneous*/, stream, NULL /*c_stream*/,
        task.max_m | task.max_n << OPENCL_LIBSMM_PFORMAT |
            (task.max_k << (OPENCL_LIBSMM_PFORMAT * 2)),
        NULL);
  }
#endif
  if (0 != info && EXIT_SUCCESS == result) {
    static LIBXSMM_TLS DBM_TIMER_TICKINT start2 = 0;
    const DBM_TIMER_TICKINT stop = DBM_TIMER_TICK();
    const double dhost = DBM_TIMER_DIFF(start, stop);
    const double diter = (0 < start2 ? DBM_TIMER_DIFF(start, start2) : dhost);
#if defined(OPENCL_LIBSMM_PFORMAT) && (0 < OPENCL_LIBSMM_PFORMAT)
    const char *const kind = (0 >= dbm_multiply_opencl_smm ? "DBM" : "SMM");
#else
    const char *const kind = "DBM";
#endif
    const int pure = (100 * (ntasks - task.mnk_changes) + ntasks - 1) / ntasks;
    const double dtotl = LIBXSMM_MAX(diter, dhost);
    start2 = stop;
    fprintf(stderr,
            "INFO ACC/LIBDBM: %s-kernel mnk=%ix%ix%i "
            "pure=%i%% ntasks=%i ms=%.1f\n",
            kind, task.max_m, task.max_n, task.max_k, pure, ntasks,
            1E+3 * dtotl);
  }
  return result;
}

void dbm_multiply_gpu_launch_kernel(offloadStream_t stream, double alpha,
                                    int ntasks, const dbm_task_t *tasks_host,
                                    const dbm_task_t *tasks,
                                    const double *pack_a_data,
                                    const double *pack_b_data,
                                    double *shard_c_data) {
  const int result = dbm_multiply_opencl_launch_kernel(
      stream, alpha, ntasks, 0 /*param_format*/, &tasks_host->m, &tasks->m,
      pack_a_data, pack_b_data, shard_c_data);
  OFFLOAD_CHECK(result);
}

#endif // defined(__OFFLOAD_OPENCL) && !defined(__NO_OFFLOAD_DBM)

// EOF
