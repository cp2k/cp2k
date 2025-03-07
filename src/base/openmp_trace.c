/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2025 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

static int openmp_trace_issues_n;
static int openmp_trace_level;

int openmp_trace_issues(void);
int openmp_trace_issues(void) { /* routine is exposed in Fortran interface */
  return 0 != openmp_trace_level ? openmp_trace_issues_n : -1 /*disabled*/;
}

#if defined(_OPENMP)
/* #include <omp.h>: avoid functionality being traced */
#include <assert.h>
#include <ctype.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Simple compile-time check if OMPT is available (omp/iomp, not gomp).
 * __clang__: omp and iomp/icx, __INTEL_COMPILER: iomp/icc
 * __INTEL_LLVM_COMPILER: already covered by __clang__
 */
#if (defined(__clang__) || defined(__INTEL_COMPILER))
#include <omp-tools.h>
#else
typedef struct ompt_frame_t ompt_frame_t;
typedef void *ompt_initialize_t;
typedef void *ompt_finalize_t;
typedef void *ompt_callback_t;
typedef union ompt_data_t {
  uint64_t value;
  void *ptr;
} ompt_data_t;
typedef struct ompt_start_tool_result_t {
  ompt_initialize_t initialize;
  ompt_finalize_t finalize;
  ompt_data_t tool_data;
} ompt_start_tool_result_t;
typedef enum ompt_scope_endpoint_t {
  ompt_scope_begin = 1,
  ompt_scope_end,
  ompt_scope_beginend
} ompt_scope_endpoint_t;
typedef enum ompt_set_result_t {
  ompt_set_never = 1,
} ompt_set_result_t;
typedef enum ompt_callbacks_t {
  ompt_callback_parallel_begin = 3,
  ompt_callback_parallel_end = 4,
  ompt_callback_work = 20,
  ompt_callback_master = 21,
  ompt_callback_sync_region = 23,
} ompt_callbacks_t;
typedef enum ompt_parallel_flag_t {
  ompt_parallel_team = 0x80000000
} ompt_parallel_flag_t;
typedef enum ompt_sync_region_t {
  ompt_sync_region_barrier = 1,
  ompt_sync_region_barrier_implicit,
  ompt_sync_region_barrier_explicit,
  ompt_sync_region_barrier_implementation
} ompt_sync_region_t;
typedef enum ompt_work_t {
  ompt_work_loop = 1,
  ompt_work_sections,
  ompt_work_single_executor,
  ompt_work_single_other,
  ompt_work_workshare,
} ompt_work_t;

typedef void (*ompt_interface_fn_t)(void);
typedef int (*ompt_get_parallel_info_t)(int, ompt_data_t **, int *);
typedef ompt_interface_fn_t (*ompt_function_lookup_t)(const char *);
typedef ompt_set_result_t (*ompt_set_callback_t)(ompt_callbacks_t,
                                                 ompt_callback_t);
#endif

#if !defined(_WIN32) && !defined(__CYGWIN__) && !defined(OPENMP_TRACE_SYMBOL)
#define OPENMP_TRACE_SYMBOL
#include <execinfo.h>
#include <unistd.h>
#endif

#define OPENMP_TRACE_PTR_KIND(PTR) (int)(((uintptr_t)(PTR)) >> 56)
#define OPENMP_TRACE_PTR_SYMBOL(PTR)                                           \
  (const void *)(0x0FFFFFFFFFFFFFFF & ((uintptr_t)(PTR)))
#define OPENMP_TRACE_PTR(PTR, KIND)                                            \
  (const void *)((((uintptr_t)(0xF & (KIND))) << 56) |                         \
                 (uintptr_t)OPENMP_TRACE_PTR_SYMBOL(PTR))
#define OPENMP_TRACE_SET_CALLBACK(PREFIX, NAME)                                \
  if (ompt_set_never ==                                                        \
      set_callback(ompt_callback_##NAME, (ompt_callback_t)PREFIX##_##NAME)) {  \
    ++openmp_trace_issues_n;                                                   \
  }
#define OPENMP_TRACE_PRINT(KIND, FORMAT, ...)                                  \
  fprintf(stderr, "OMP/TRACE %s: " FORMAT, KIND, __VA_ARGS__)
#define OPENMP_TRACE_UNUSED(VAR) (void)VAR
#if 0
#define OPENMP_TRACE_ENABLE(FEATURE) (FEATURE)
#else
#define OPENMP_TRACE_ENABLE(FEATURE) 0
#endif

enum {
  openmp_trace_level_deflt = 2,
  openmp_trace_level_high = 4,
  openmp_trace_level_warn,
  openmp_trace_level_info
};

static int openmp_trace_parallel_n;
static int openmp_trace_sync_n;

static const void *openmp_trace_parallel;
static ompt_data_t *openmp_trace_sync;

static ompt_get_parallel_info_t openmp_trace_get_parallel_info;

/* translate debug symbol (address) to character string */
static void openmp_trace_symbol(const void *symbol, char *str, size_t size,
                                int cleanup) {
  if (NULL != str && 0 < size) {
#if !defined(OPENMP_TRACE_SYMBOL)
    OPENMP_TRACE_UNUSED(symbol);
#else
    int pipefd[2];
    if (NULL != symbol && 0 == pipe(pipefd)) {
      void *const backtrace[] = {(void *)(uintptr_t)symbol};
      backtrace_symbols_fd(backtrace, 1, pipefd[1]);
      close(pipefd[1]);
      if (0 < read(pipefd[0], str, size)) {
        char *s = (char *)(0 != cleanup ? memchr(str, '(', size) : NULL);
        char *t =
            (char *)(NULL != s ? memchr(s + 1, '+', size - (s - str)) : NULL);
        if (NULL != t) {
          *t = '\0';
          memmove(str, s + 1, t - s);
        }
        s = (char *)memchr(str, '\n', size);
        if (NULL != s) {
          *s = '\0';
        }
        for (s = str; s < (str + size) && '\0' != *s; ++s) {
          if (0 == isprint(*s)) {
            *str = '\0';
            break;
          }
        }
      } else {
        *str = '\0';
      }
      close(pipefd[0]);
    } else
#endif
    { *str = '\0'; }
  }
}

/* give a name to a kind of synchronization construct */
static const char *openmp_trace_sync_name(int kind) {
  static const char *kinds[] = {
      "master",   "barrier", "implicit barrier", "explicit barrier",
      "sections", "single",  "single",           "workshare"};
  return (kind * sizeof(*kinds)) < sizeof(kinds) ? kinds[kind]
                                                 : "synchronization";
}

/* https://www.openmp.org/spec-html/5.0/openmpsu187.html */
static void openmp_trace_parallel_begin(
    ompt_data_t *encountering_task_data,
    const ompt_frame_t *encountering_task_frame, ompt_data_t *parallel_data,
    unsigned int requested_parallelism, int flags, const void *codeptr_ra) {
  OPENMP_TRACE_UNUSED(encountering_task_data);
  OPENMP_TRACE_UNUSED(encountering_task_frame);
  OPENMP_TRACE_UNUSED(parallel_data);
  OPENMP_TRACE_UNUSED(requested_parallelism);
  if (ompt_parallel_team & flags) {
    const ompt_data_t *sync;
#pragma omp atomic read
    sync = openmp_trace_sync;
    if (NULL != sync) {
      const int kind = OPENMP_TRACE_PTR_KIND(sync->ptr);
      if (ompt_sync_region_barrier_implementation > kind) {
        ++openmp_trace_issues_n;
      }
      if (1 /*assert*/ < openmp_trace_level || 0 > openmp_trace_level) {
        const char *const type =
            (ompt_sync_region_barrier_implementation > kind ? "ERROR" : "WARN");
        if ('E' == *type || openmp_trace_level_warn <= openmp_trace_level) {
          const char *const name = openmp_trace_sync_name(kind);
          char symbol[1024], symbol2[1024];
          openmp_trace_symbol(codeptr_ra, symbol, sizeof(symbol),
                              1 /*cleanup*/);
          openmp_trace_symbol(OPENMP_TRACE_PTR_SYMBOL(sync->ptr), symbol2,
                              sizeof(symbol2), 1 /*cleanup*/);
          if ('\0' != *symbol) {
            if ('\0' != *symbol2) {
              OPENMP_TRACE_PRINT(type,
                                 "parallel region \"%s\" opened in %s \"%s\"\n",
                                 symbol, name, symbol2);
            } else {
              OPENMP_TRACE_PRINT(type, "parallel region \"%s\" opened in %s\n",
                                 symbol, name);
            }
          } else {
            if ('\0' != *symbol2) {
              OPENMP_TRACE_PRINT(type, "parallel region opened in %s \"%s\"\n",
                                 name, symbol2);
            } else {
              OPENMP_TRACE_PRINT(type, "parallel region opened in %s\n", name);
            }
          }
        }
      } else {
        assert(0);
      }
    }
  }
}

/* https://www.openmp.org/spec-html/5.0/openmpsu187.html */
static void openmp_trace_parallel_end(ompt_data_t *parallel_data,
                                      ompt_data_t *encountering_task_data,
                                      int flags, const void *codeptr_ra) {
  OPENMP_TRACE_UNUSED(parallel_data);
  OPENMP_TRACE_UNUSED(encountering_task_data);
  if (0 != (ompt_parallel_team & flags) &&
      0 != openmp_trace_get_parallel_info(openmp_trace_parallel_n + 1, NULL,
                                          NULL)) {
    openmp_trace_parallel = codeptr_ra;
    ++openmp_trace_parallel_n;
  }
}

/* https://www.openmp.org/spec-html/5.0/openmpsu187.html */
static void openmp_trace_master(ompt_scope_endpoint_t endpoint,
                                ompt_data_t *parallel_data,
                                ompt_data_t *task_data,
                                const void *codeptr_ra) {
  OPENMP_TRACE_UNUSED(task_data);
  if (NULL != parallel_data) {
    int sync_n;
    switch ((int)endpoint) {
    case OPENMP_TRACE_ENABLE(ompt_scope_beginend):
    case ompt_scope_begin: {
      if (OPENMP_TRACE_ENABLE(ompt_scope_beginend) != endpoint) {
#pragma omp atomic capture
        sync_n = openmp_trace_sync_n++;
      } else {
#pragma omp atomic read
        sync_n = openmp_trace_sync_n;
      }
      if (0 == sync_n) {
        assert(OPENMP_TRACE_PTR(codeptr_ra, 0) == codeptr_ra);
        parallel_data->ptr = (void *)(uintptr_t)codeptr_ra;
        openmp_trace_sync = parallel_data;
      }
    } break;
    case ompt_scope_end: {
#pragma omp atomic capture
      sync_n = --openmp_trace_sync_n;
      if (0 == sync_n) {
        openmp_trace_sync = NULL;
      }
    } break;
    }
  }
}

/* https://www.openmp.org/spec-html/5.0/openmpsu187.html */
void openmp_trace_sync_region(ompt_sync_region_t kind,
                              ompt_scope_endpoint_t endpoint,
                              ompt_data_t *parallel_data,
                              ompt_data_t *task_data, const void *codeptr_ra) {
  OPENMP_TRACE_UNUSED(task_data);
  assert(0 < kind);
  if (NULL != parallel_data && ompt_sync_region_barrier_implementation > kind) {
    int sync_n;
    switch ((int)endpoint) {
    case OPENMP_TRACE_ENABLE(ompt_scope_beginend):
    case ompt_scope_begin: {
      if (OPENMP_TRACE_ENABLE(ompt_scope_beginend) != endpoint) {
#pragma omp atomic capture
        sync_n = openmp_trace_sync_n++;
      } else {
#pragma omp atomic read
        sync_n = openmp_trace_sync_n;
      }
      if (0 == sync_n) {
        assert(OPENMP_TRACE_PTR(codeptr_ra, 0) == codeptr_ra);
        parallel_data->ptr =
            (void *)(uintptr_t)OPENMP_TRACE_PTR(codeptr_ra, kind);
        openmp_trace_sync = parallel_data;
      } else if (openmp_trace_level_warn <= openmp_trace_level ||
                 0 > openmp_trace_level) {
        const ompt_data_t *sync;
#pragma omp atomic read
        sync = openmp_trace_sync;
        if (NULL != sync && parallel_data != sync) {
          const char *const name = openmp_trace_sync_name(kind);
          char symbol[1024], symbol2[1024];
          openmp_trace_symbol(codeptr_ra, symbol, sizeof(symbol),
                              1 /*cleanup*/);
          openmp_trace_symbol(OPENMP_TRACE_PTR_SYMBOL(sync->ptr), symbol2,
                              sizeof(symbol2), 1 /*cleanup*/);
          if ('\0' != *symbol) {
            if ('\0' != *symbol2) {
              OPENMP_TRACE_PRINT("WARN",
                                 "potential deadlock at \"%s\" in %s \"%s\"\n",
                                 symbol2, name, symbol);
            } else {
              OPENMP_TRACE_PRINT("WARN", "potential deadlock in %s \"%s\"\n",
                                 name, symbol);
            }
          } else {
            if ('\0' != *symbol2) {
              OPENMP_TRACE_PRINT("WARN", "potential deadlock at \"%s\" in %s\n",
                                 symbol2, name);
            } else {
              OPENMP_TRACE_PRINT("WARN", "potential deadlock in %s\n", name);
            }
          }
        }
      }
    } break;
    case ompt_scope_end: {
#pragma omp atomic capture
      sync_n = --openmp_trace_sync_n;
      if (0 == sync_n) {
        openmp_trace_sync = NULL;
      }
    } break;
    }
  }
}

/* https://www.openmp.org/spec-html/5.0/openmpsu187.html */
static void openmp_trace_work(ompt_work_t wstype,
                              ompt_scope_endpoint_t endpoint,
                              ompt_data_t *parallel_data,
                              ompt_data_t *task_data, uint64_t count,
                              const void *codeptr_ra) {
  OPENMP_TRACE_UNUSED(task_data);
  OPENMP_TRACE_UNUSED(count);
  assert(0 < wstype);
  if (NULL != parallel_data && ompt_work_sections <= wstype &&
      wstype <= ompt_work_workshare) {
    int sync_n;
    switch ((int)endpoint) {
    case OPENMP_TRACE_ENABLE(ompt_scope_beginend):
    case ompt_scope_begin: {
      if (OPENMP_TRACE_ENABLE(ompt_scope_beginend) != endpoint) {
#pragma omp atomic capture
        sync_n = openmp_trace_sync_n++;
      } else {
#pragma omp atomic read
        sync_n = openmp_trace_sync_n;
      }
      if (0 == sync_n) {
        const int kind = wstype - ompt_work_sections +
                         ompt_sync_region_barrier_implementation;
        assert(OPENMP_TRACE_PTR(codeptr_ra, 0) == codeptr_ra);
        parallel_data->ptr =
            (void *)(uintptr_t)OPENMP_TRACE_PTR(codeptr_ra, kind);
        openmp_trace_sync = parallel_data;
      }
    } break;
    case ompt_scope_end: {
#pragma omp atomic capture
      sync_n = --openmp_trace_sync_n;
      if (0 == sync_n) {
        openmp_trace_sync = NULL;
      }
    } break;
    }
  }
}

/* initially, events of interest are registered */
static int openmp_trace_initialize(ompt_function_lookup_t lookup,
                                   int initial_device_num,
                                   ompt_data_t *tool_data) {
  const ompt_set_callback_t set_callback =
      (ompt_set_callback_t)lookup("ompt_set_callback");
  openmp_trace_get_parallel_info =
      (ompt_get_parallel_info_t)lookup("ompt_get_parallel_info");
  OPENMP_TRACE_UNUSED(initial_device_num);
  OPENMP_TRACE_UNUSED(tool_data);
  OPENMP_TRACE_SET_CALLBACK(openmp_trace, parallel_begin);
  OPENMP_TRACE_SET_CALLBACK(openmp_trace, parallel_end);
  OPENMP_TRACE_SET_CALLBACK(openmp_trace, master);
  if (openmp_trace_level_deflt < openmp_trace_level || 0 > openmp_trace_level) {
    OPENMP_TRACE_SET_CALLBACK(openmp_trace, sync_region);
  }
  if (openmp_trace_level_high <= openmp_trace_level || 0 > openmp_trace_level) {
    OPENMP_TRACE_SET_CALLBACK(openmp_trace, work);
  }
  assert(NULL != openmp_trace_get_parallel_info);
  return 0 == openmp_trace_issues();
}

/* here tool_data might be freed and analysis concludes */
static void openmp_trace_finalize(ompt_data_t *tool_data) {
  OPENMP_TRACE_UNUSED(tool_data);
  if (openmp_trace_level_info <= openmp_trace_level || 0 > openmp_trace_level) {
    if (1 < openmp_trace_parallel_n) { /* nested */
      char symbol[1024];
      openmp_trace_symbol(openmp_trace_parallel, symbol, sizeof(symbol),
                          1 /*cleanup*/);
      if ('\0' != *symbol) {
        OPENMP_TRACE_PRINT("INFO", "parallelism in \"%s\" is nested (%i)\n",
                           symbol, openmp_trace_parallel_n);
      } else {
        OPENMP_TRACE_PRINT("INFO", "parallelism is nested (%i)\n",
                           openmp_trace_parallel_n);
      }
    }
  }
}

/* entry point which is automatically called by the OpenMP runtime */
ompt_start_tool_result_t *ompt_start_tool(unsigned int omp_version,
                                          const char *runtime_version) {
  static ompt_start_tool_result_t openmp_start_tool;
  const char *const enabled_env = getenv("CP2K_OMP_TRACE");
  ompt_start_tool_result_t *result = NULL;
  openmp_trace_level = (NULL == enabled_env ? 0 : atoi(enabled_env));
  OPENMP_TRACE_UNUSED(omp_version);
  OPENMP_TRACE_UNUSED(runtime_version);
  if (0 != openmp_trace_level) { /* trace OpenMP constructs */
    openmp_start_tool.initialize = (ompt_initialize_t)openmp_trace_initialize;
    openmp_start_tool.finalize = (ompt_finalize_t)openmp_trace_finalize;
    openmp_start_tool.tool_data.ptr = NULL;
    result = &openmp_start_tool;
#if defined(NDEBUG)
    if (1 == openmp_trace_level) {
      openmp_trace_level = 2; /* adjust trace level */
    }
#endif
  }
  return result;
}

#endif
