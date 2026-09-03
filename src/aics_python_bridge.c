/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include "aics_python_loader.h"
#include <Python.h>
#include <dlfcn.h>
#include <mpi.h>

#include <errno.h>
#include <limits.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static Cp2kPythonApi python_api;
static int python_api_ready;
static int python_owned_by_cp2k;
static int python_normal_output_enabled = 1;
static char *python_loaded_library;

#define CP2K_AICS_ERROR_SIZE 1024
static char cp2k_aics_last_error[CP2K_AICS_ERROR_SIZE];

/* Keep Python.h's types and inline reference-counting rules, but route every
 * dynamic CPython entry point and exported data object through the exact-handle
 * API table. */
#define PyCallable_Check python_api.callable_check
#define PyConfig_Clear python_api.config_clear
#define PyConfig_InitPythonConfig python_api.config_init_python
#define PyConfig_Read python_api.config_read
#define PyConfig_SetBytesString python_api.config_set_bytes_string
#undef Py_DECREF
#undef Py_XDECREF
#define Py_DECREF(object) python_api.dec_ref((PyObject *)(object))
#define Py_XDECREF(object)                                                     \
  do {                                                                         \
    if ((object) != NULL)                                                      \
      python_api.dec_ref((PyObject *)(object));                                \
  } while (0)
#define PyErr_Format python_api.err_format
#define PyErr_Clear python_api.err_clear
#define PyErr_Fetch python_api.err_fetch
#define PyErr_NoMemory python_api.err_no_memory
#define PyErr_NormalizeException python_api.err_normalize_exception
#define PyErr_Occurred python_api.err_occurred
#define PyErr_Print python_api.err_print
#define PyErr_Restore python_api.err_restore
#define PyErr_SetString python_api.err_set_string
#define PyExc_RuntimeError (*python_api.exc_runtime_error)
#define PyExc_TypeError (*python_api.exc_type_error)
#define PyFloat_FromDouble python_api.float_from_double
#define PyImport_ImportModule python_api.import_module
#define PyList_Insert python_api.list_insert
#define PyLong_FromLong python_api.long_from_long
#define PyObject_CallFunctionObjArgs python_api.object_call_function_obj_args
#define PyObject_CallObject python_api.object_call_object
#define PyObject_GetAttrString python_api.object_get_attr_string
#define PyObject_Str python_api.object_str
#undef PyRun_SimpleString
#define PyRun_SimpleString python_api.run_simple_string
#define PySequence_Contains python_api.sequence_contains
#define PyStatus_Exception python_api.status_exception
#define PyStatus_IsExit python_api.status_is_exit
#define PySys_GetObject python_api.sys_get_object
#define PyUnicode_DecodeFSDefault python_api.unicode_decode_fs_default
#define PyUnicode_AsUTF8 python_api.unicode_as_utf8
#define PyUnicode_FromString python_api.unicode_from_string
#define Py_FinalizeEx python_api.finalize_ex
#define Py_GetPrefix python_api.get_prefix
#define Py_GetProgramFullPath python_api.get_program_full_path
#define Py_GetVersion python_api.get_version
#define Py_InitializeFromConfig python_api.initialize_from_config
#define Py_IsInitialized python_api.is_initialized

static const char *cp2k_aics_python_relative_dir =
    CP2K_AICS_PYTHON_RELATIVE_DIR;

_Static_assert(sizeof(MPI_Fint) == sizeof(int),
               "MPI_Fint must match Fortran integer(c_int) communicator ABI");

typedef struct {
  PyObject *module;
  PyObject *prepare;
  PyObject *initialize;
  PyObject *solve;
  PyObject *finalize;
  PyObject *pending_request;
  int initialized;
  int current_scope_depth;
  int owner_scope_depth;
  int python_owned;
} Cp2kAicsBridge;

static Cp2kAicsBridge cp2k_aics_bridge = {NULL, NULL, NULL, NULL, NULL,
                                          NULL, 0,    0,    0,    0};

static int runtime_diagnostics_enabled(void) {
  return python_normal_output_enabled &&
         getenv("CP2K_AICS_RUNTIME_DIAGNOSTICS") != NULL;
}

static void clear_last_error(void) { cp2k_aics_last_error[0] = '\0'; }

static void set_last_error(const char *message) {
  if (message == NULL)
    message = "CP2K-AICS failed without an available diagnostic.";
  snprintf(cp2k_aics_last_error, sizeof(cp2k_aics_last_error), "%s", message);
}

static void set_last_errorf(const char *format, ...) {
  va_list arguments;

  va_start(arguments, format);
  vsnprintf(cp2k_aics_last_error, sizeof(cp2k_aics_last_error), format,
            arguments);
  va_end(arguments);
}

static void set_last_error_if_empty(const char *message) {
  if (cp2k_aics_last_error[0] == '\0')
    set_last_error(message);
}

/* Collapse newlines and redact absolute path tokens from normal output. */
static void set_concise_text(const char *raw) {
  char concise[CP2K_AICS_ERROR_SIZE];
  size_t input = 0;
  size_t output = 0;
  int previous_was_space = 0;

  while (raw[input] != '\0' && output + 1 < sizeof(concise)) {
    unsigned char current = (unsigned char)raw[input];
    int path_boundary = input == 0 || raw[input - 1] == ' ' ||
                        raw[input - 1] == '\t' || raw[input - 1] == '\n' ||
                        raw[input - 1] == '\'' || raw[input - 1] == '"' ||
                        raw[input - 1] == '(' || raw[input - 1] == '[' ||
                        raw[input - 1] == '=';

    if (current == '/' && path_boundary) {
      static const char redacted[] = "<path>";
      size_t redacted_length = sizeof(redacted) - 1;
      size_t available = sizeof(concise) - output - 1;
      size_t copy_length =
          redacted_length < available ? redacted_length : available;

      memcpy(concise + output, redacted, copy_length);
      output += copy_length;
      input++;
      while (raw[input] != '\0' && raw[input] != ' ' && raw[input] != '\t' &&
             raw[input] != '\r' && raw[input] != '\n' && raw[input] != '\'' &&
             raw[input] != '"' && raw[input] != ')' && raw[input] != ']')
        input++;
      previous_was_space = 0;
      continue;
    }
    if (current == ' ' || current == '\t' || current == '\r' ||
        current == '\n') {
      if (!previous_was_space && output > 0)
        concise[output++] = ' ';
      previous_was_space = 1;
      input++;
      continue;
    }
    concise[output++] = (char)current;
    previous_was_space = 0;
    input++;
  }
  while (output > 0 && concise[output - 1] == ' ')
    output--;
  concise[output] = '\0';
  set_last_error(concise);
}

/* Copy one exception type/value into the process-local fixed buffer. */
static void set_concise_python_error(const char *type_name,
                                     const char *value_text,
                                     const char *operation) {
  char raw[CP2K_AICS_ERROR_SIZE * 2];

  if (type_name != NULL && value_text != NULL && value_text[0] != '\0')
    snprintf(raw, sizeof(raw), "%s: %s", type_name, value_text);
  else if (value_text != NULL && value_text[0] != '\0')
    snprintf(raw, sizeof(raw), "%s", value_text);
  else if (type_name != NULL)
    snprintf(raw, sizeof(raw), "%s while %s", type_name, operation);
  else
    snprintf(raw, sizeof(raw), "CP2K-AICS Python failure while %s", operation);
  set_concise_text(raw);
}

void cp2k_aics_get_last_error(char *buffer, size_t buffer_size) {
  size_t length;

  if (buffer == NULL || buffer_size == 0)
    return;
  length = strlen(cp2k_aics_last_error);
  if (length >= buffer_size)
    length = buffer_size - 1;
  memcpy(buffer, cp2k_aics_last_error, length);
  buffer[length] = '\0';
}

static void report_python_initialization_status(const char *stage,
                                                PyStatus status) {
  char message[CP2K_AICS_ERROR_SIZE * 2];

  snprintf(message, sizeof(message),
           "CP2K-AICS could not initialize Python during %s%s%s.", stage,
           status.err_msg != NULL ? ": " : "",
           status.err_msg != NULL ? status.err_msg : "");
  set_concise_text(message);
  if (runtime_diagnostics_enabled()) {
    fprintf(stderr,
            "CP2K_AICS_DIAGNOSTIC Python initialization failed during %s",
            stage);
    if (status.func != NULL)
      fprintf(stderr, " in %s", status.func);
    if (status.err_msg != NULL)
      fprintf(stderr, ": %s", status.err_msg);
    if (PyStatus_IsExit(status))
      fprintf(stderr, " (exit code %d)", status.exitcode);
    fprintf(stderr, ".\n");
  }
}

int cp2k_python_initialize(void) {
  const char *executable_override = getenv("CP2K_PYTHON_EXECUTABLE");
  const char *library_override = getenv("CP2K_PYTHON_LIBRARY");
  const int has_executable_override =
      executable_override != NULL && executable_override[0] != '\0';
  const int has_library_override =
      library_override != NULL && library_override[0] != '\0';
  const char *selected_executable = CP2K_CONFIGURED_PYTHON_EXECUTABLE;
  const char *selected_library = CP2K_CONFIGURED_PYTHON_LIBRARY;
  char runtime_version[64];
  char error[1024];
  char *external_library = NULL;
  int preinitialized;
  PyConfig config;
  PyStatus status;

  if (python_api_ready && python_owned_by_cp2k && Py_IsInitialized())
    return 0;

  if (has_executable_override != has_library_override) {
    set_last_error("CP2K_AICS requires CP2K_PYTHON_EXECUTABLE and "
                   "CP2K_PYTHON_LIBRARY to be set together.");
    return -1;
  }
  if (has_executable_override) {
    selected_executable = executable_override;
    selected_library = library_override;
  }

  if (cp2k_python_validate_executable(selected_executable, error,
                                      sizeof(error)) != 0) {
    set_last_error("CP2K-AICS configured Python executable is unavailable or "
                   "not executable.");
    if (runtime_diagnostics_enabled())
      fprintf(stderr, "CP2K_AICS_DIAGNOSTIC %s.\n", error);
    return -1;
  }

  /* AICS deliberately supports only an interpreter initialized and owned by
   * CP2K.  This inspection uses only Py_IsInitialized through a temporary
   * loader handle; no CPython object API or cached API state is entered. */
  preinitialized =
      cp2k_python_is_preinitialized(&external_library, error, sizeof(error));
  if (preinitialized < 0) {
    set_last_error("CP2K-AICS could not inspect the active Python runtime.");
    if (runtime_diagnostics_enabled())
      fprintf(stderr, "CP2K_AICS_DIAGNOSTIC %s.\n", error);
    return -1;
  }
  if (preinitialized > 0) {
    set_last_error(
        "CP2K-AICS does not support an externally initialized Python "
        "interpreter; run AICS in a CP2K-owned process.");
    if (runtime_diagnostics_enabled())
      fprintf(stderr, "CP2K_AICS_DIAGNOSTIC active external runtime '%s'.\n",
              external_library);
    free(external_library);
    return -1;
  }
  if (cp2k_python_api_load_exact(selected_library, &python_api, error,
                                 sizeof(error)) != 0) {
    set_last_error("CP2K-AICS could not load the configured Python runtime.");
    if (runtime_diagnostics_enabled())
      fprintf(stderr, "CP2K_AICS_DIAGNOSTIC %s.\n", error);
    return -1;
  } else {
    python_api_ready = 1;
    python_loaded_library = strdup(selected_library);
  }

  snprintf(runtime_version, sizeof(runtime_version), "%s", Py_GetVersion());
  runtime_version[strcspn(runtime_version, " \r\n")] = '\0';
  if (strcmp(runtime_version, CP2K_CONFIGURED_PYTHON_VERSION) != 0) {
    set_last_errorf("CP2K-AICS Python build/runtime version mismatch "
                    "(configured %s, loaded %s, library %s).",
                    CP2K_CONFIGURED_PYTHON_VERSION, runtime_version,
                    selected_library);
    return -1;
  }

  PyConfig_InitPythonConfig(&config);
  config.parse_argv = 0;
  config.site_import = 1;

  status = PyConfig_SetBytesString(&config, &config.program_name,
                                   selected_executable);
  if (PyStatus_Exception(status)) {
    report_python_initialization_status("setting program_name", status);
    PyConfig_Clear(&config);
    return -1;
  }
  status = PyConfig_Read(&config);
  if (PyStatus_Exception(status)) {
    report_python_initialization_status("reading path configuration", status);
    PyConfig_Clear(&config);
    return -1;
  }
  status = Py_InitializeFromConfig(&config);
  if (PyStatus_Exception(status)) {
    report_python_initialization_status("initializing the interpreter", status);
    PyConfig_Clear(&config);
    return -1;
  }
  PyConfig_Clear(&config);

  if (!Py_IsInitialized()) {
    set_last_error("CP2K-AICS Python initialization returned without an active "
                   "interpreter.");
    return -1;
  }
  python_owned_by_cp2k = 1;
  if (Py_GetPrefix() == NULL || Py_GetPrefix()[0] == L'\0' ||
      Py_GetProgramFullPath() == NULL || Py_GetProgramFullPath()[0] == L'\0') {
    set_last_error("CP2K-AICS initialized Python with inconsistent path "
                   "configuration.");
    if (Py_FinalizeEx() != 0)
      set_last_error("CP2K-AICS could not clean up a failed Python "
                     "initialization.");
    python_owned_by_cp2k = 0;
    return -1;
  }
  if (runtime_diagnostics_enabled())
    fprintf(stderr,
            "CP2K-AICS loaded exact Python runtime '%s' (version %s) "
            "for executable '%s'.\n",
            selected_library, runtime_version, selected_executable);
  return 0;
}

static int cp2k_aics_api_is_cached(void) {
  return cp2k_aics_bridge.module != NULL || cp2k_aics_bridge.prepare != NULL ||
         cp2k_aics_bridge.initialize != NULL ||
         cp2k_aics_bridge.solve != NULL || cp2k_aics_bridge.finalize != NULL ||
         cp2k_aics_bridge.pending_request != NULL;
}

// This helper may only be called while the CPython interpreter is alive.
static void clear_cp2k_aics_api_while_python_alive(void) {
  Py_XDECREF(cp2k_aics_bridge.pending_request);
  Py_XDECREF(cp2k_aics_bridge.finalize);
  Py_XDECREF(cp2k_aics_bridge.solve);
  Py_XDECREF(cp2k_aics_bridge.initialize);
  Py_XDECREF(cp2k_aics_bridge.prepare);
  Py_XDECREF(cp2k_aics_bridge.module);
  cp2k_aics_bridge.pending_request = NULL;
  cp2k_aics_bridge.finalize = NULL;
  cp2k_aics_bridge.solve = NULL;
  cp2k_aics_bridge.initialize = NULL;
  cp2k_aics_bridge.prepare = NULL;
  cp2k_aics_bridge.module = NULL;
  cp2k_aics_bridge.initialized = 0;
  cp2k_aics_bridge.owner_scope_depth = 0;
}

static int cp2k_aics_scope_state_is_valid(void) {
  if (cp2k_aics_bridge.current_scope_depth < 0 ||
      cp2k_aics_bridge.owner_scope_depth < 0 ||
      cp2k_aics_bridge.owner_scope_depth >
          cp2k_aics_bridge.current_scope_depth ||
      (!cp2k_aics_bridge.initialized &&
       cp2k_aics_bridge.owner_scope_depth != 0)) {
    set_last_error("CP2K-AICS Python runtime state is inconsistent.");
    return 0;
  }
  return 1;
}

static PyObject *get_callable(PyObject *module, const char *name) {
  PyObject *callable = PyObject_GetAttrString(module, name);

  if (callable != NULL && !PyCallable_Check(callable)) {
    PyErr_Format(PyExc_TypeError, "cp2k_aics.%s is not callable", name);
    Py_DECREF(callable);
    callable = NULL;
  }
  return callable;
}

static int add_cp2k_aics_package_directory(void) {
  Dl_info library_info;
  const char *separator;
  size_t library_dir_length;
  size_t package_dir_length;
  char *package_dir = NULL;
  PyObject *python_package_dir = NULL;
  PyObject *sys_path = NULL;
  int contains;
  int status = -1;

  if (dladdr((const void *)&add_cp2k_aics_package_directory, &library_info) ==
          0 ||
      library_info.dli_fname == NULL) {
    PyErr_SetString(PyExc_RuntimeError,
                    "CP2K-AICS could not locate its loaded native library");
    return -1;
  }

  separator = strrchr(library_info.dli_fname, '/');
  library_dir_length =
      separator == NULL ? 1 : (size_t)(separator - library_info.dli_fname);
  package_dir_length =
      library_dir_length + 1 + strlen(cp2k_aics_python_relative_dir) + 1;
  package_dir = malloc(package_dir_length);
  if (package_dir == NULL) {
    PyErr_NoMemory();
    goto fail;
  }
  if (separator == NULL) {
    package_dir[0] = '.';
  } else {
    memcpy(package_dir, library_info.dli_fname, library_dir_length);
  }
  package_dir[library_dir_length] = '/';
  memcpy(package_dir + library_dir_length + 1, cp2k_aics_python_relative_dir,
         strlen(cp2k_aics_python_relative_dir) + 1);

  python_package_dir = PyUnicode_DecodeFSDefault(package_dir);
  if (python_package_dir == NULL)
    goto fail;
  sys_path = PySys_GetObject("path");
  if (sys_path == NULL || !PyList_Check(sys_path)) {
    PyErr_SetString(PyExc_RuntimeError, "sys.path is unavailable");
    goto fail;
  }
  contains = PySequence_Contains(sys_path, python_package_dir);
  if (contains < 0)
    goto fail;
  if (contains == 0 && PyList_Insert(sys_path, 0, python_package_dir) < 0)
    goto fail;
  status = 0;

fail:
  Py_XDECREF(python_package_dir);
  free(package_dir);
  return status;
}

static PyObject *get_loaded_cp2k_library_path(void) {
  Dl_info library_info;

  if (dladdr((const void *)&get_loaded_cp2k_library_path, &library_info) == 0 ||
      library_info.dli_fname == NULL || library_info.dli_fname[0] == '\0') {
    PyErr_SetString(PyExc_RuntimeError,
                    "CP2K-AICS could not identify its loaded native library");
    return NULL;
  }
  return PyUnicode_DecodeFSDefault(library_info.dli_fname);
}

static int load_cp2k_aics_api(void) {
  PyObject *module = NULL;
  PyObject *prepare = NULL;
  PyObject *initialize = NULL;
  PyObject *solve = NULL;
  PyObject *finalize = NULL;

  if (cp2k_aics_bridge.module != NULL)
    return 0;

  clear_cp2k_aics_api_while_python_alive();
  if (add_cp2k_aics_package_directory() != 0)
    goto fail;
  module = PyImport_ImportModule("cp2k_aics");
  if (module == NULL)
    goto fail;
  prepare = get_callable(module, "prepare_initialize");
  if (prepare == NULL)
    goto fail;
  initialize = get_callable(module, "initialize_prepared");
  if (initialize == NULL)
    goto fail;
  solve = get_callable(module, "solve");
  if (solve == NULL)
    goto fail;
  finalize = get_callable(module, "finalize");
  if (finalize == NULL)
    goto fail;

  cp2k_aics_bridge.module = module;
  cp2k_aics_bridge.prepare = prepare;
  cp2k_aics_bridge.initialize = initialize;
  cp2k_aics_bridge.solve = solve;
  cp2k_aics_bridge.finalize = finalize;
  if (runtime_diagnostics_enabled() &&
      PyRun_SimpleString(
          "import sys, cp2k_aics, numpy, mpi4py, petsc4py\n"
          "print('CP2K_AICS_DIAGNOSTIC sys.executable=' + sys.executable, "
          "file=sys.stderr)\n"
          "print('CP2K_AICS_DIAGNOSTIC sys.prefix=' + sys.prefix, "
          "file=sys.stderr)\n"
          "print('CP2K_AICS_DIAGNOSTIC cp2k_aics=' + cp2k_aics.__file__, "
          "file=sys.stderr)\n"
          "print('CP2K_AICS_DIAGNOSTIC numpy=' + numpy.__file__, "
          "file=sys.stderr)\n"
          "print('CP2K_AICS_DIAGNOSTIC mpi4py=' + mpi4py.__file__, "
          "file=sys.stderr)\n"
          "print('CP2K_AICS_DIAGNOSTIC petsc4py=' + petsc4py.__file__, "
          "file=sys.stderr)\n") != 0) {
    clear_cp2k_aics_api_while_python_alive();
    return -1;
  }
  return 0;

fail:
  Py_XDECREF(finalize);
  Py_XDECREF(solve);
  Py_XDECREF(initialize);
  Py_XDECREF(prepare);
  Py_XDECREF(module);
  return -1;
}

static int call_python_no_args(PyObject *callable) {
  PyObject *result = PyObject_CallObject(callable, NULL);

  if (result == NULL)
    return -1;
  Py_DECREF(result);
  return 0;
}

static void report_python_error(const char *operation) {
  PyObject *exception_type = NULL;
  PyObject *exception_value = NULL;
  PyObject *exception_traceback = NULL;
  PyObject *type_name_object = NULL;
  PyObject *value_text_object = NULL;
  const char *type_name = NULL;
  const char *value_text = NULL;

  if (!PyErr_Occurred()) {
    set_last_error_if_empty("CP2K-AICS Python operation failed without an "
                            "available exception.");
    return;
  }

  PyErr_Fetch(&exception_type, &exception_value, &exception_traceback);
  PyErr_NormalizeException(&exception_type, &exception_value,
                           &exception_traceback);
  if (exception_type != NULL) {
    type_name_object = PyObject_GetAttrString(exception_type, "__name__");
    if (type_name_object != NULL) {
      type_name = PyUnicode_AsUTF8(type_name_object);
      if (type_name == NULL)
        PyErr_Clear();
    } else {
      PyErr_Clear();
    }
  }
  if (exception_value != NULL) {
    value_text_object = PyObject_Str(exception_value);
    if (value_text_object != NULL) {
      value_text = PyUnicode_AsUTF8(value_text_object);
      if (value_text == NULL)
        PyErr_Clear();
    } else {
      PyErr_Clear();
    }
  }
  set_concise_python_error(type_name, value_text, operation);

  Py_XDECREF(value_text_object);
  Py_XDECREF(type_name_object);
  /* String conversion may itself have failed.  Discard that secondary error
   * before either restoring the original exception or releasing it. */
  PyErr_Clear();
  if (runtime_diagnostics_enabled()) {
    PyErr_Restore(exception_type, exception_value, exception_traceback);
    exception_type = NULL;
    exception_value = NULL;
    exception_traceback = NULL;
    PyErr_Print();
  }
  Py_XDECREF(exception_traceback);
  Py_XDECREF(exception_value);
  Py_XDECREF(exception_type);
  PyErr_Clear();
}

int cp2k_aics_initialize(int comm_f_arg, int normal_output_enabled,
                         int poisson_solver, int normal_axis,
                         double relative_tolerance, double absolute_tolerance,
                         int max_iterations) {
  int mpi_initialized = 0;
  MPI_Fint comm_f = (MPI_Fint)comm_f_arg;
  MPI_Comm comm;
  char *cwd = NULL;
  PyObject *bridge_library_path = NULL;
  PyObject *coupling_dir = NULL;
  PyObject *communicator_handle = NULL;
  PyObject *poisson_solver_name = NULL;
  PyObject *relative_tolerance_arg = NULL;
  PyObject *absolute_tolerance_arg = NULL;
  PyObject *max_iterations_arg = NULL;
  PyObject *normal_axis_arg = NULL;
  PyObject *request = NULL;

  python_normal_output_enabled = normal_output_enabled;
  clear_last_error();
  if (cp2k_python_initialize() != 0 || !python_api_ready ||
      !python_owned_by_cp2k || !Py_IsInitialized())
    return -1;
  cp2k_aics_bridge.python_owned = python_owned_by_cp2k;
  if (!cp2k_aics_scope_state_is_valid())
    return -1;
  if (cp2k_aics_bridge.initialized)
    return 0;
  if (cp2k_aics_bridge.pending_request != NULL) {
    set_last_error("CP2K-AICS initialization preparation is already pending.");
    return -1;
  }
  if (cp2k_aics_bridge.owner_scope_depth != 0) {
    set_last_error("CP2K-AICS Python runtime state is inconsistent during "
                   "initialization.");
    return -1;
  }
  if (MPI_Initialized(&mpi_initialized) != MPI_SUCCESS || !mpi_initialized) {
    set_last_error("CP2K-AICS requires MPI to be initialized before Python "
                   "coupling starts.");
    return -1;
  }
  comm = MPI_Comm_f2c(comm_f);
  if (comm == MPI_COMM_NULL) {
    set_last_error("CP2K-AICS received an invalid MPI communicator during "
                   "initialization.");
    return -1;
  }
  comm_f = MPI_Comm_c2f(comm);
  if (load_cp2k_aics_api() != 0) {
    report_python_error("loading the cp2k_aics API");
    return -1;
  }
  bridge_library_path = get_loaded_cp2k_library_path();
  if (bridge_library_path == NULL) {
    report_python_error("identifying the loaded CP2K library");
    return -1;
  }
  cwd = getcwd(NULL, 0);
  if (cwd == NULL) {
    Py_DECREF(bridge_library_path);
    set_last_errorf(
        "CP2K-AICS could not determine the calculation directory: %s.",
        strerror(errno));
    return -1;
  }
  coupling_dir = PyUnicode_DecodeFSDefault(cwd);
  free(cwd);
  if (coupling_dir == NULL) {
    Py_DECREF(bridge_library_path);
    report_python_error("creating the coupling directory argument");
    return -1;
  }
  communicator_handle = PyLong_FromLong((long)comm_f);
  if (communicator_handle == NULL) {
    Py_DECREF(coupling_dir);
    Py_DECREF(bridge_library_path);
    report_python_error("creating the communicator argument");
    return -1;
  }
  if (poisson_solver == 0)
    poisson_solver_name = PyUnicode_FromString("finite_difference");
  else if (poisson_solver == 1)
    poisson_solver_name = PyUnicode_FromString("finite_element");
  else {
    Py_DECREF(communicator_handle);
    Py_DECREF(coupling_dir);
    Py_DECREF(bridge_library_path);
    set_last_errorf("CP2K-AICS received invalid Poisson solver code %d.",
                    poisson_solver);
    return -1;
  }
  if (poisson_solver_name == NULL) {
    Py_DECREF(communicator_handle);
    Py_DECREF(coupling_dir);
    Py_DECREF(bridge_library_path);
    report_python_error("creating the Poisson solver argument");
    return -1;
  }
  relative_tolerance_arg = PyFloat_FromDouble(relative_tolerance);
  absolute_tolerance_arg = PyFloat_FromDouble(absolute_tolerance);
  max_iterations_arg = PyLong_FromLong((long)max_iterations);
  normal_axis_arg = PyLong_FromLong((long)normal_axis);
  if (relative_tolerance_arg == NULL || absolute_tolerance_arg == NULL ||
      max_iterations_arg == NULL || normal_axis_arg == NULL) {
    Py_XDECREF(normal_axis_arg);
    Py_XDECREF(max_iterations_arg);
    Py_XDECREF(absolute_tolerance_arg);
    Py_XDECREF(relative_tolerance_arg);
    Py_DECREF(poisson_solver_name);
    Py_DECREF(communicator_handle);
    Py_DECREF(coupling_dir);
    Py_DECREF(bridge_library_path);
    report_python_error("creating the linear solver arguments");
    return -1;
  }
  request = PyObject_CallFunctionObjArgs(
      cp2k_aics_bridge.prepare, bridge_library_path, coupling_dir,
      communicator_handle, poisson_solver_name, normal_axis_arg,
      relative_tolerance_arg, absolute_tolerance_arg, max_iterations_arg, NULL);
  Py_DECREF(normal_axis_arg);
  Py_DECREF(max_iterations_arg);
  Py_DECREF(absolute_tolerance_arg);
  Py_DECREF(relative_tolerance_arg);
  Py_DECREF(poisson_solver_name);
  Py_DECREF(communicator_handle);
  Py_DECREF(coupling_dir);
  Py_DECREF(bridge_library_path);
  if (request == NULL) {
    report_python_error("cp2k_aics.prepare_initialize");
    return -1;
  }
  cp2k_aics_bridge.pending_request = request;
  return 0;
}

int cp2k_aics_initialize_prepared(void) {
  PyObject *context;

  clear_last_error();
  if (cp2k_aics_bridge.initialized)
    return 0;
  if (!python_api_ready || !python_owned_by_cp2k || !Py_IsInitialized() ||
      cp2k_aics_bridge.pending_request == NULL ||
      cp2k_aics_bridge.initialize == NULL) {
    set_last_error("CP2K-AICS context initialization was requested without "
                   "successful rank-local preparation.");
    return -1;
  }
  context = PyObject_CallFunctionObjArgs(
      cp2k_aics_bridge.initialize, cp2k_aics_bridge.pending_request, NULL);
  Py_DECREF(cp2k_aics_bridge.pending_request);
  cp2k_aics_bridge.pending_request = NULL;
  if (context == NULL) {
    report_python_error("cp2k_aics.initialize_prepared");
    return -1;
  }
  Py_DECREF(context);
  cp2k_aics_bridge.initialized = 1;
  // Depth zero denotes an externally initiated context; it must be finalized
  // explicitly because no bridge scope owns it.
  cp2k_aics_bridge.owner_scope_depth = cp2k_aics_bridge.current_scope_depth;
  return 0;
}

void cp2k_aics_cancel_initialize(void) {
  if (python_api_ready && Py_IsInitialized())
    Py_XDECREF(cp2k_aics_bridge.pending_request);
  cp2k_aics_bridge.pending_request = NULL;
}

int cp2k_aics_solve(void) {
  clear_last_error();
  if (!cp2k_aics_bridge.initialized) {
    set_last_error("CP2K-AICS solve was requested before successful Python "
                   "initialization.");
    return -1;
  }
  if (load_cp2k_aics_api() != 0) {
    report_python_error("loading the cp2k_aics API");
    return -1;
  }
  if (call_python_no_args(cp2k_aics_bridge.solve) != 0) {
    report_python_error("cp2k_aics.solve");
    return -1;
  }
  return 0;
}

int cp2k_aics_scope_enter(void) {
  clear_last_error();
  if (!cp2k_aics_scope_state_is_valid())
    return -1;
  if (cp2k_aics_bridge.current_scope_depth == INT_MAX) {
    set_last_error("CP2K-AICS Python runtime nesting limit was exceeded.");
    return -1;
  }

  cp2k_aics_bridge.current_scope_depth++;
  return 0;
}

int cp2k_aics_finalize(void) {
  int status = 0;

  clear_last_error();

  if (!python_api_ready || !Py_IsInitialized()) {
    if (!cp2k_aics_api_is_cached() && !cp2k_aics_bridge.initialized)
      return 0;
    set_last_error("CP2K-AICS Python runtime state is inconsistent after "
                   "interpreter shutdown.");
    return -1;
  }

  if (!cp2k_aics_api_is_cached() && !cp2k_aics_bridge.initialized)
    return 0;

  if (cp2k_aics_bridge.initialized) {
    if (cp2k_aics_bridge.finalize == NULL) {
      set_last_error("CP2K-AICS Python runtime cannot be finalized because its "
                     "lifecycle API is incomplete.");
      status = -1;
    } else if (call_python_no_args(cp2k_aics_bridge.finalize) != 0) {
      report_python_error("cp2k_aics.finalize");
      status = -1;
    }
  }

  clear_cp2k_aics_api_while_python_alive();
  return status;
}

int cp2k_aics_scope_leave(void) {
  int status = 0;

  clear_last_error();

  if (cp2k_aics_bridge.current_scope_depth <= 0) {
    set_last_error("CP2K-AICS Python runtime lifecycle is unbalanced.");
    return -1;
  }
  if (python_api_ready && Py_IsInitialized() &&
      !cp2k_aics_scope_state_is_valid())
    return -1;

  if (python_api_ready && Py_IsInitialized() &&
      (cp2k_aics_bridge.owner_scope_depth ==
           cp2k_aics_bridge.current_scope_depth ||
       ((cp2k_aics_bridge.python_owned || python_owned_by_cp2k) &&
        cp2k_aics_bridge.current_scope_depth == 1 &&
        cp2k_aics_api_is_cached())))
    status = cp2k_aics_finalize();

  cp2k_aics_bridge.current_scope_depth--;
  return status;
}

int cp2k_python_finalize(void) {
  int status = 0;

  clear_last_error();
  if (!python_api_ready || !python_owned_by_cp2k)
    return 0;
  if (!Py_IsInitialized()) {
    set_last_error("CP2K-AICS Python runtime state is inconsistent during "
                   "process finalization.");
    return -1;
  }
  if (cp2k_aics_bridge.current_scope_depth != 0 ||
      cp2k_aics_bridge.initialized || cp2k_aics_api_is_cached()) {
    set_last_error(
        "CP2K-AICS process finalization found an active run context.");
    return -1;
  }
  if (Py_FinalizeEx() != 0) {
    set_last_error("CP2K-AICS could not finalize the Python runtime.");
    status = -1;
  }
  cp2k_aics_bridge.python_owned = 0;
  python_owned_by_cp2k = 0;
  python_api_ready = 0;
  memset(&python_api, 0, sizeof(python_api));
  free(python_loaded_library);
  python_loaded_library = NULL;
  return status;
}
