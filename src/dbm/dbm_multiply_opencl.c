/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#include "../offload/offload_runtime.h"
#if defined(__OFFLOAD_OPENCL) && !defined(__NO_OFFLOAD_DBM)

#include "dbm_multiply_gpu_kernel.h"
#include "dbm_multiply_opencl.cl.h"
#include <libxs/libxs_reg.h>
#include <libxs/libxs_timer.h>
#include <libxstream/libxstream_opencl.h>

#if !defined(OPENCL_DBM_SOURCE_MULTIPLY)
#error "OpenCL kernel source code not found!"
#endif

#if defined(OPENCL_LIBSMM_PFORMAT)
#define DBM_OPENCL_LIBSMM_PFORMAT OPENCL_LIBSMM_PFORMAT
#else
#define DBM_OPENCL_LIBSMM_PFORMAT 8
#endif

#define DBM_OPENCL_CMEM LIBXSTREAM_CMEM
#define DBM_TIMER_DIFF(A, B) libxs_timer_duration(A, B)
#define DBM_TIMER_TICK() libxs_timer_tick()
#define DBM_TIMER_TICKINT libxs_timer_tick_t

#if 0 < DBM_OPENCL_LIBSMM_PFORMAT
#define DBM_OPENCL_DBCSR_TYPE_REAL_8 3
typedef int (*opencl_libsmm_acc_dbm_launch_fn_t)(
  void *stream, double alpha, int ntasks, int param_format,
  const int *params_host, const int *params, const double *pack_a_data,
  const double *pack_b_data, double *shard_c_data);
void opencl_libsmm_acc_set_dbm_launch_fn(
  opencl_libsmm_acc_dbm_launch_fn_t launch_fn);
int opencl_libsmm_acc_process(
  const int *host_param_stack, const int *dev_param_stack, int stack_size,
  int datatype, const void *dev_a_data, const void *dev_b_data,
  void *dev_c_data, int m_max, int n_max, int k_max, int max_kernel_dim,
  int def_mnk, void *stream, void *c_stream, int param_format, void *event);
LIBXS_PRAGMA_WEAK(opencl_libsmm_acc_set_dbm_launch_fn)
LIBXS_PRAGMA_WEAK(opencl_libsmm_acc_process)
#endif

/* Dispatch key for specialized kernels (must be memset-initialized).
   For homogeneous batches: m,n,k hold the exact shape (bk=k).
   For heterogeneous batches: m=n=k=0, bk holds the unroll threshold,
   max_m enables compile-time division in the flat dispatch mapping. */
typedef struct {
  int bk;
  int m, n, k;
  int max_m;
} dbm_multiply_opencl_key_t;

typedef struct {
  int max_m, max_n, max_k, mnk_changes;
} dbm_multiply_gpu_launch_info_t;

#if 0 < DBM_OPENCL_LIBSMM_PFORMAT
int dbm_multiply_opencl_initialized /*= 0*/;
int dbm_multiply_opencl_smm /*= 0*/;
#endif

int dbm_multiply_opencl_launch_kernel(void *stream, double alpha, int ntasks,
                                      int param_format, const int *params_host,
                                      const int *params,
                                      const double *pack_a_data,
                                      const double *pack_b_data,
                                      double *shard_c_data);

#if 0 < DBM_OPENCL_LIBSMM_PFORMAT
LIBXS_ATTRIBUTE_CTOR static void dbm_multiply_opencl_initialize(void) {
  const char *const smm_env = getenv("DBM_MULTIPLY_SMM");
  const int smm = (NULL == smm_env ? 0 /*default*/ : atoi(smm_env));
  dbm_multiply_opencl_smm =
      LIBXS_MIN(1 != smm ? smm : 64,
                (1 << (DBM_OPENCL_LIBSMM_PFORMAT - 1)) - 1);
  if (0 > dbm_multiply_opencl_smm &&
      NULL != opencl_libsmm_acc_set_dbm_launch_fn) {
    opencl_libsmm_acc_set_dbm_launch_fn(dbm_multiply_opencl_launch_kernel);
  }
  LIBXS_ATOMIC_STORE(&dbm_multiply_opencl_initialized, 1, LIBXS_ATOMIC_SEQ_CST);
}
#endif

static int dbm_multiply_gpu_launch_info(dbm_multiply_gpu_launch_info_t *info,
                                        const int *params, int ntasks,
                                        int param_format, int stop_at_impure) {
  assert(0 < ntasks);
  if (0 == param_format) { /* native */
    const int stride = sizeof(dbm_task_t) / sizeof(int);
    const int first_m = params[0], first_n = params[1], first_k = params[2];
    int i = stride;
    info->max_m = first_m;
    info->max_n = first_n;
    info->max_k = first_k;
    for (info->mnk_changes = 0; i < (ntasks * stride); i += stride) {
      const int m = params[i + 0], n = params[i + 1], k = params[i + 2];
      info->max_m = imax(info->max_m, m);
      info->max_n = imax(info->max_n, n);
      info->max_k = imax(info->max_k, k);
      if (m != first_m || n != first_n || k != first_k) {
        info->mnk_changes = 1;
        if (0 != stop_at_impure) {
          return 0;
        }
      }
    }
  } else {
#if 0 < DBM_OPENCL_LIBSMM_PFORMAT
    const int mask = (1 << DBM_OPENCL_LIBSMM_PFORMAT) - 1;
    info->max_m = mask & (param_format);
    info->max_n = mask & (param_format >> (DBM_OPENCL_LIBSMM_PFORMAT));
    info->max_k = mask & (param_format >> (DBM_OPENCL_LIBSMM_PFORMAT * 2));
    info->mnk_changes = 0; /* homogeneous */
#else
    assert(0);
#endif
  }
  return 1;
}

static void dbm_multiply_opencl_print(FILE *stream, const char *name, int val) {
  if (0 != val) {
    fprintf(stream, " %s=%i", name, val);
  }
}

static int dbm_multiply_opencl_bk(int max_k) {
  int result = 1;
  if (16 <= max_k) {
    result = 16;
  } else if (8 <= max_k) {
    result = 8;
  } else if (4 <= max_k) {
    result = 4;
  }
  return result;
}

int dbm_multiply_opencl_launch_kernel(void *stream, double alpha, int ntasks,
                                      int param_format, const int *params_host,
                                      const int *params,
                                      const double *pack_a_data,
                                      const double *pack_b_data,
                                      double *shard_c_data) {
  const DBM_TIMER_TICKINT start = DBM_TIMER_TICK();
  const libxstream_opencl_config_t *const config = &libxstream_opencl_config;
  const int verbosity = config->verbosity,
            trace = (0 > verbosity || 2 < verbosity);
  int result = EXIT_SUCCESS;
#if 0 < DBM_OPENCL_LIBSMM_PFORMAT
  int dbcsr = 0;
#endif
  dbm_multiply_gpu_launch_info_t task = {0};
  int task_complete = 0;
  assert(NULL != pack_a_data && NULL != pack_b_data && NULL != shard_c_data);
  assert(NULL != params_host || 0 == ntasks);
  assert(NULL != params || 0 == ntasks);
  if (0 < ntasks) {
#if 0 < DBM_OPENCL_LIBSMM_PFORMAT
    if (0 == LIBXS_ATOMIC_LOAD(&dbm_multiply_opencl_initialized,
                               LIBXS_ATOMIC_SEQ_CST)) {
      dbm_multiply_opencl_initialize();
    }
    if (0 != dbm_multiply_opencl_smm || 0 != trace) {
      task_complete = dbm_multiply_gpu_launch_info(
          &task, params_host, ntasks, param_format, 0 == trace);
    }
    if (NULL == opencl_libsmm_acc_process || 0 > dbm_multiply_opencl_smm ||
        0 != task.mnk_changes ||
        dbm_multiply_opencl_smm < task.max_m ||
        dbm_multiply_opencl_smm < task.max_n ||
        dbm_multiply_opencl_smm < task.max_k || 0 == task.max_k || 1 != alpha)
#endif
    { /* base init state: computed once, shared across all specializations */
      static int clinear = 0, sgbcst = 0, bk_max = 0;
      static int nz = 0, blkrd = 0, base_ready = 0;
      static size_t wgsize[] = {1, 1, 1}, sgsize_s = 0;
      static char base_flags[LIBXSTREAM_BUFFERSIZE];
      static const char *base_options /*= NULL*/;
      static const char *base_source /*= NULL*/;
      static const char *base_exts[2] /*= {NULL, NULL}*/;
      static size_t base_nexts, base_source_kind;
      /* registry of BK-specialized kernels */
      static libxs_registry_t *kernel_registry /*= NULL*/;
      /* serializes kernel compilation (get-compile-set) */
      static libxs_lock_t compile_lock /*= LIBXS_LOCK_INITIALIZER*/;
      /* serializes kernel arg setting + enqueue (no TLS clones needed) */
      static libxs_lock_t kernel_lock /*= LIBXS_LOCK_INITIALIZER*/;
      const libxstream_opencl_stream_t *const str =
          (const libxstream_opencl_stream_t *)(stream);
      const libxstream_opencl_device_t *const devinfo = &config->device;
      libxs_lock_t *const lock_memory =
          (NULL != devinfo->clSetKernelArgMemPointerINTEL
               ? NULL
               : config->lock_memory);
      libxstream_opencl_info_memptr_t adata, bdata, cdata, batch;
      const int stride = (0 == param_format ? 6 : 3);
      size_t work_size[] = {1, 1, 1}, ibatch = 0;
      const size_t work_tasks = ntasks;
      cl_kernel kernel = NULL;
      int bk, bk0, use_blkrd = 0;
      cl_int size;
      assert(NULL != str && NULL != str->queue);
      /* base init (once): env vars, flags, source, extensions */
      if (0 == LIBXS_ATOMIC_LOAD(&base_ready, LIBXS_ATOMIC_SEQ_CST)) {
        LIBXS_LOCK_ACQUIRE(LIBXS_LOCK, config->lock_main);
        if (0 == base_ready) {
          const char *const krn_env = getenv("DBM_MULTIPLY_KERNEL");
          const char *const sgb_env = getenv("DBM_MULTIPLY_SGB");
          const char *const blk_env = getenv("DBM_MULTIPLY_BLK");
          const char *const lin_env = getenv("DBM_MULTIPLY_LIN");
          const char *const fp_env = getenv("DBM_MULTIPLY_FP");
          const char *const bn_env = getenv("DBM_MULTIPLY_BN");
          const char *const bk_env = getenv("DBM_MULTIPLY_BK");
          const char *const sm_env = getenv("DBM_MULTIPLY_SM");
          const char *const wg_env = getenv("DBM_MULTIPLY_WG");
          const char *const lu_env = getenv("DBM_MULTIPLY_LU");
          const char *const ro_env = getenv("DBM_MULTIPLY_RO");
          const char *const xf_env = getenv("DBM_MULTIPLY_XF");
          const char *const nz_env = getenv("DBM_MULTIPLY_NZ");
          const char *options = NULL;
          const int dd = (0 != config->debug && 0 != config->dump);
          const int ro = (NULL == ro_env ? -1 /*default*/ : atoi(ro_env));
          const int xf = (NULL == xf_env ? -1 /*default*/ : atoi(xf_env));
          const int sm0 = (NULL == sm_env ? 0 : atoi(sm_env));
          int source_kind = 0, sm = LIBXS_ABS(sm0);
          const int bn0 = (0 == devinfo->nv ? 8 : 2);
          const int bn1 = ((0 == sm && 0 == clinear) ? bn0 : (bn0 * sm * 2));
          const int gpu = (CL_DEVICE_TYPE_GPU == devinfo->type);
          const int precision = (NULL == fp_env ? 0 /*default*/ : atoi(fp_env));
          int bn = LIBXS_CLMP(NULL == bn_env ? bn1 : atoi(bn_env), 1, 32);
          int lu = LIBXS_CLMP(NULL == lu_env ? 0 : atoi(lu_env), -2, 1);
          size_t sgsize = devinfo->wgsize[2];
          size_t offset;
          const char *source = OPENCL_DBM_SOURCE_MULTIPLY, *cmem = NULL;
          LIBXS_MEMZERO(base_flags);
          LIBXS_SNPRINTF(base_flags, sizeof(base_flags),
                         "-cl-fast-relaxed-math -cl-denorms-are-zero");
          offset = (0 == dd ? strlen(base_flags) : 0);
          base_exts[0] = base_exts[1] = NULL;
          base_nexts = sizeof(base_exts) / sizeof(*base_exts);
          offset += (size_t)libxstream_opencl_flags_atomics(
              devinfo, libxstream_opencl_atomic_fp_64, base_exts, &base_nexts,
              base_flags + offset, sizeof(base_flags) - offset);
          if (NULL != krn_env) {
            FILE *const krn_file = fopen(krn_env, "rb");
            if (NULL != krn_file) {
              fclose(krn_file);
              source = krn_env;
              source_kind = 1;
            }
          }
          wgsize[0] = (NULL == wg_env ? LIBXS_MAX((unsigned long int)sm,
                                                  devinfo->wgsize[1])
                                      : strtoul(wg_env, NULL, 10));
          if (1 < sgsize && 0 < wgsize[0]) { /* subgroups */
            if (LIBXS_DELTA(wgsize[0], devinfo->wgsize[1]) <=
                LIBXS_DELTA(wgsize[0], sgsize)) {
              sgsize = devinfo->wgsize[1];
            }
            wgsize[0] = LIBXS_UP(wgsize[0], sgsize);
          } else {
            wgsize[0] = LIBXS_UP(wgsize[0], devinfo->wgsize[1]);
            sgsize = 0;
          }
          { /* 256-GRF */
            const int biggrf = (0 <= xf ? xf : devinfo->biggrf);
            const size_t max_wgs =
                (0 != biggrf) ? devinfo->wgsize[0] / 2 : devinfo->wgsize[0];
            wgsize[0] = LIBXS_CLMP(wgsize[0], 0, max_wgs);
            if (0 != biggrf && 0 != devinfo->intel && 0 == devinfo->biggrf) {
              options = "-cl-intel-256-GRF-per-thread";
            }
          }
          sm = ((0 != sm && 0 != wgsize[0])
                    ? (LIBXS_ISPOT(bn * sizeof(double)) + 1)
                    : 0);
          clinear = (NULL == lin_env ? 0 /*default*/ : atoi(lin_env));
          sgbcst = (0 != gpu && 0 < sgsize && 0 < wgsize[0] &&
                    2 <= devinfo->std_level[0] &&
                    (NULL == sgb_env ? 1 /*default*/ : (0 != atoi(sgb_env))));
          cmem =
#if defined(DBM_OPENCL_CMEM)
              (0 > ro && EXIT_SUCCESS == libxstream_opencl_use_cmem(devinfo))
                  ? "constant"
                  :
#endif
                  (0 >= ro ? "global" : "constant");
          blkrd = (0 == clinear && 0 != devinfo->intel && 0 < (int)sgsize &&
                   (NULL == blk_env ? 1 /*default*/ : (0 != atoi(blk_env))));
          if (0 != blkrd && 'g' != cmem[0]) {
            cmem = "global"; /* block reads require global address space */
          }
          offset += (size_t)LIBXS_SNPRINTF(
              base_flags + offset, sizeof(base_flags) - offset,
              " %s %s -DCONSTANT=%s"
              " -DBN=%i -DSM=%i -DLU=%i -DSG=%i -DINTEL=%i",
              0 != gpu ? "-DGPU" : "", 0 == clinear ? "" : "-DCLINEAR", cmem,
              bn, sm, lu, (int)sgsize, (int)(0 != devinfo->intel));
          if (0 != precision) {
            offset += (size_t)LIBXS_SNPRINTF(base_flags + offset,
                                             sizeof(base_flags) - offset,
                                             " -DPRECISION=%i", precision);
          }
          bk_max = (NULL == bk_env ? 0 /*default*/ : atoi(bk_env));
          nz = (NULL == nz_env ? 0 /*default*/ : atoi(nz_env));
          if (0 != nz) {
            offset += (size_t)LIBXS_SNPRINTF(base_flags + offset,
                                             sizeof(base_flags) - offset,
                                             " -DNZ=%i", nz);
          }
          sgsize_s = sgsize;
          base_source = source;
          base_source_kind = source_kind;
          base_options = options;
          kernel_registry = libxs_registry_create();
          if (2 <= verbosity || 0 > verbosity) {
            fprintf(stderr, "INFO ACC/LIBDBM: DBM-kernel gpu=%i", gpu);
            dbm_multiply_opencl_print(stderr, "sgb", sgbcst);
            dbm_multiply_opencl_print(stderr, "lin", clinear);
            dbm_multiply_opencl_print(stderr, "fp", precision);
            dbm_multiply_opencl_print(stderr, "bn", bn);
            dbm_multiply_opencl_print(stderr, "sm", sm);
            dbm_multiply_opencl_print(stderr, "wg", (int)wgsize[0]);
            dbm_multiply_opencl_print(stderr, "sg", (int)sgsize);
            dbm_multiply_opencl_print(stderr, "lu", lu);
            dbm_multiply_opencl_print(stderr, "nz", nz);
            dbm_multiply_opencl_print(stderr, "blk", blkrd);
            fprintf(stderr, " -> %.1f ms\n",
                    1E3 * DBM_TIMER_DIFF(start, DBM_TIMER_TICK()));
          }
          LIBXS_ATOMIC_STORE(&base_ready, 1, LIBXS_ATOMIC_SEQ_CST);
        }
        LIBXS_LOCK_RELEASE(LIBXS_LOCK, config->lock_main);
      }
      /* per-launch: compute task info and dispatch key */
#if 0 < DBM_OPENCL_LIBSMM_PFORMAT
      if (0 == task_complete)
#endif
      {
        task_complete = dbm_multiply_gpu_launch_info(&task, params_host, ntasks,
                                                     param_format, 0);
      }
      bk0 = dbm_multiply_opencl_bk(task.max_k);
      bk = (0 < bk_max ? LIBXS_MIN(bk0, bk_max) : bk0);
      { /* per-shape kernel lookup/compile */
        dbm_multiply_opencl_key_t key;
        cl_kernel *kptr;
        LIBXS_MEMZERO(&key);
        key.bk = bk;
        if (0 == task.mnk_changes) { /* homogeneous: exact shape */
          key.m = task.max_m;
          key.n = task.max_n;
          key.k = task.max_k;
          key.bk = task.max_k; /* exact K for full unrolling */
        } else { /* heterogeneous: max_m for compile-time division.
                    With BLKRD_P (per-task dispatch), max_m is unused
                    so omit it to consolidate into fewer kernels. */
          if (0 == blkrd) {
            key.max_m = (0 == clinear ? task.max_m : task.max_n);
          }
        }
        use_blkrd = (0 != blkrd && 0 != key.m && 16 <= key.m &&
                     key.m <= (int)sgsize_s && 0 == (key.m & (key.m - 1)));
        kptr = (cl_kernel *)libxs_registry_get(
            kernel_registry, &key, sizeof(key),
            libxs_registry_lock(kernel_registry));
        if (NULL == kptr || NULL == *kptr) { /* compile specialization */
          LIBXS_LOCK_ACQUIRE(LIBXS_LOCK, &compile_lock);
          kptr = (cl_kernel *)libxs_registry_get(
              kernel_registry, &key, sizeof(key),
              libxs_registry_lock(kernel_registry));
          if (NULL == kptr || NULL == *kptr) {
            char flags[LIBXSTREAM_BUFFERSIZE];
            const cl_device_id device_id = config->devices[config->device_id];
            cl_kernel kernel_new = NULL;
            size_t wgs[3];
            if (0 != key.m) { /* homogeneous: add shape defines */
              const int use_wg = (0 != use_blkrd || 0 != sgbcst);
              const int n = LIBXS_SNPRINTF(
                  flags, sizeof(flags),
                  "%s -DWG=%i -DBK=%i -DDBM_M=%i -DDBM_N=%i -DDBM_K=%i%s%s",
                  base_flags, use_wg ? (int)wgsize[0] : 0, key.bk, key.m, key.n,
                  key.k, 0 != use_blkrd ? " -DBLKRD_A" : "",
                  (0 != sgbcst && 0 == use_blkrd) ? " -DSGBCST" : "");
              assert(0 < n && (size_t)n < sizeof(flags));
              LIBXS_UNUSED(n);
            } else if (0 < key.max_m) { /* heterogeneous with known max_m */
              const int n = LIBXS_SNPRINTF(
                  flags, sizeof(flags), "%s -DWG=%i -DBK=%i -DMAX_M=%i%s%s",
                  base_flags, (int)wgsize[0], bk, key.max_m,
                  0 != blkrd ? " -DBLKRD_P" : "",
                  (0 != sgbcst && 0 == blkrd) ? " -DSGBCST" : "");
              assert(0 < n && (size_t)n < sizeof(flags));
              LIBXS_UNUSED(n);
            } else { /* heterogeneous: BK only */
              const int n = LIBXS_SNPRINTF(
                  flags, sizeof(flags), "%s -DWG=%i -DBK=%i%s%s", base_flags,
                  (int)wgsize[0], bk, 0 != blkrd ? " -DBLKRD_P" : "",
                  (0 != sgbcst && 0 == blkrd) ? " -DSGBCST" : "");
              assert(0 < n && (size_t)n < sizeof(flags));
              LIBXS_UNUSED(n);
            }
            result |= libxstream_opencl_kernel(
                base_source_kind, base_source, "dbm_multiply", flags,
                base_options, NULL /*try*/, NULL /*try_ok*/, base_exts,
                base_nexts, &kernel_new);
            if (EXIT_SUCCESS == result &&
                EXIT_SUCCESS ==
                    clGetKernelWorkGroupInfo(kernel_new, device_id,
                                             CL_KERNEL_COMPILE_WORK_GROUP_SIZE,
                                             sizeof(wgs), wgs, NULL) &&
                0 != wgs[0] && 0 != wgs[1]) {
              wgsize[0] = wgs[0];
              wgsize[1] = wgs[1];
            }
            kptr = (cl_kernel *)libxs_registry_set(
                kernel_registry, &key, sizeof(key), &kernel_new,
                sizeof(kernel_new), libxs_registry_lock(kernel_registry));
            if (2 <= verbosity || 0 > verbosity || EXIT_SUCCESS != result) {
              const char *const kind =
                  (EXIT_SUCCESS == result ? "INFO" : "ERROR");
              fprintf(stderr, "%s ACC/LIBDBM: DBM-kernel bk=%i", kind, key.bk);
              if (0 != key.m) {
                fprintf(stderr, " mnk=%ix%ix%i", key.m, key.n, key.k);
              }
              fprintf(stderr, " -> ");
              if (EXIT_SUCCESS == result) {
                fprintf(stderr, "%.1f ms\n",
                        1E3 * DBM_TIMER_DIFF(start, DBM_TIMER_TICK()));
              } else {
                fprintf(stderr, "FAILED!\n");
              }
            }
          }
          LIBXS_LOCK_RELEASE(LIBXS_LOCK, &compile_lock);
        }
        kernel = (NULL != kptr) ? *kptr : NULL;
      }
      LIBXS_LOCK_ACQUIRE(LIBXS_LOCK, &kernel_lock);
      if (NULL != lock_memory) {
        LIBXS_LOCK_ACQUIRE(LIBXS_LOCK, lock_memory);
      }
      { /* assume A, B, and C do not carry an offset */
        size_t iadata = 0, ibdata = 0, icdata = 0;
        result |= libxstream_opencl_info_devptr_lock(&adata, NULL /*lock*/,
                                                     pack_a_data, 1 /*esize*/,
                                                     NULL /*amount*/, &iadata);
        result |= libxstream_opencl_info_devptr_lock(&bdata, NULL /*lock*/,
                                                     pack_b_data, 1 /*esize*/,
                                                     NULL /*amount*/, &ibdata);
        result |= libxstream_opencl_info_devptr_lock(&cdata, NULL /*lock*/,
                                                     shard_c_data, 1 /*esize*/,
                                                     NULL /*amount*/, &icdata);
        assert(0 == iadata && 0 == ibdata && 0 == icdata);
      }
      result |= libxstream_opencl_info_devptr_lock(
          &batch, NULL /*lock*/, params /*batch*/, sizeof(int) * stride,
          &work_tasks, &ibatch);
      if (NULL != lock_memory) {
        LIBXS_LOCK_RELEASE(LIBXS_LOCK, lock_memory);
      }
      if (EXIT_SUCCESS ==
          result) { /* determine dispatch mode: per-task vs flat */
        const int per_task = (0 != sgbcst && 0 == use_blkrd) ||
                             (0 != blkrd && 0 != task.mnk_changes);
        const int use_wg =
            (0 != task.mnk_changes || 0 != use_blkrd || 0 != sgbcst);
        size = (cl_int)(work_tasks * (0 == clinear ? task.max_m : task.max_n));
        if (per_task) {
          work_size[0] = work_tasks * wgsize[0];
        } else if (use_wg) {
          work_size[0] = LIBXS_UP((size_t)size, wgsize[0]);
        } else {
          work_size[0] = (size_t)size;
        }
        result |= clSetKernelArg(kernel, 0, sizeof(cl_double), &alpha);
        result |= clSetKernelArg(kernel, 1, sizeof(cl_int), &ibatch);
        result |= clSetKernelArg(kernel, 2, sizeof(cl_int), &ntasks);
        result |= clSetKernelArg(kernel, 3, sizeof(cl_int), &size);
        result |= clSetKernelArg(kernel, 4, sizeof(cl_int), &param_format);
        result |= libxstream_opencl_set_kernel_ptr(kernel, 5, batch.memory);
        result |= libxstream_opencl_set_kernel_ptr(kernel, 6, adata.memory);
        result |= libxstream_opencl_set_kernel_ptr(kernel, 7, bdata.memory);
        result |= libxstream_opencl_set_kernel_ptr(kernel, 8, cdata.memory);
        result |=
            clEnqueueNDRangeKernel(str->queue, kernel, 1, NULL, work_size,
                                   0 < wgsize[0] ? wgsize : NULL,
                                   0 /*num_wait*/, NULL /*wait_list*/, NULL);
      }
      LIBXS_LOCK_RELEASE(LIBXS_LOCK, &kernel_lock);
    }
#if 0 < DBM_OPENCL_LIBSMM_PFORMAT
    else { /* homogeneous */
      result |= opencl_libsmm_acc_process(
          params_host, params, ntasks, DBM_OPENCL_DBCSR_TYPE_REAL_8,
          pack_a_data, pack_b_data, shard_c_data, task.max_m, task.max_n,
          task.max_k, dbm_multiply_opencl_smm, 1 /*homogeneous*/, stream,
          NULL /*c_stream*/,
          task.max_m | task.max_n << DBM_OPENCL_LIBSMM_PFORMAT |
              (task.max_k << (DBM_OPENCL_LIBSMM_PFORMAT * 2)),
          NULL);
      dbcsr = 1;
    }
#endif
    if (0 != trace && EXIT_SUCCESS == result) {
      static LIBXS_TLS DBM_TIMER_TICKINT start2 = 0;
      const DBM_TIMER_TICKINT stop = DBM_TIMER_TICK();
      const double dhost = DBM_TIMER_DIFF(start, stop);
      const double diter = (0 < start2 ? DBM_TIMER_DIFF(start, start2) : dhost);
#if 0 < DBM_OPENCL_LIBSMM_PFORMAT
      const char *const kind = (0 == dbcsr ? "DBM" : "SMM");
#else
      const char *const kind = "DBM";
#endif
      const int pure = (0 == task.mnk_changes ? 100 : 0);
      const double dtotl = LIBXS_MAX(diter, dhost);
      start2 = stop;
      fprintf(stderr,
              "INFO ACC/LIBDBM: %s-kernel mnk=%ix%ix%i "
              "pure=%i%% ntasks=%i ms=%.1f\n",
              kind, task.max_m, task.max_n, task.max_k, pure, ntasks,
              1E+3 * dtotl);
    }
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

#endif /* defined(__OFFLOAD_OPENCL) && !defined(__NO_OFFLOAD_DBM) */

/* EOF */
