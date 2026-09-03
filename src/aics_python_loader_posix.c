/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include "aics_python_loader.h"

#include <dlfcn.h>
#include <errno.h>
#include <link.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

static void *python_handle;

static int resolve_api(void *handle, Cp2kPythonApi *api, char *error,
                       size_t error_size) {
#define RESOLVE(field, symbol)                                                 \
  do {                                                                         \
    void *address = dlsym(handle, #symbol);                                    \
    if (address == NULL) {                                                     \
      snprintf(error, error_size,                                              \
               "required CPython symbol %s is unavailable: %s", #symbol,       \
               dlerror());                                                     \
      return -1;                                                               \
    }                                                                          \
    memcpy(&api->field, &address, sizeof(address));                            \
  } while (0)
  memset(api, 0, sizeof(*api));
  RESOLVE(callable_check, PyCallable_Check);
  RESOLVE(config_clear, PyConfig_Clear);
  RESOLVE(config_init_python, PyConfig_InitPythonConfig);
  RESOLVE(config_read, PyConfig_Read);
  RESOLVE(config_set_bytes_string, PyConfig_SetBytesString);
  RESOLVE(dec_ref, Py_DecRef);
  RESOLVE(err_clear, PyErr_Clear);
  RESOLVE(err_fetch, PyErr_Fetch);
  RESOLVE(err_format, PyErr_Format);
  RESOLVE(err_no_memory, PyErr_NoMemory);
  RESOLVE(err_normalize_exception, PyErr_NormalizeException);
  RESOLVE(err_occurred, PyErr_Occurred);
  RESOLVE(err_print, PyErr_Print);
  RESOLVE(err_restore, PyErr_Restore);
  RESOLVE(err_set_string, PyErr_SetString);
  RESOLVE(exc_runtime_error, PyExc_RuntimeError);
  RESOLVE(exc_type_error, PyExc_TypeError);
  RESOLVE(float_from_double, PyFloat_FromDouble);
  RESOLVE(import_module, PyImport_ImportModule);
  RESOLVE(list_insert, PyList_Insert);
  RESOLVE(long_from_long, PyLong_FromLong);
  RESOLVE(object_call_function_obj_args, PyObject_CallFunctionObjArgs);
  RESOLVE(object_call_object, PyObject_CallObject);
  RESOLVE(object_get_attr_string, PyObject_GetAttrString);
  RESOLVE(object_str, PyObject_Str);
  RESOLVE(run_simple_string, PyRun_SimpleString);
  RESOLVE(sequence_contains, PySequence_Contains);
  RESOLVE(status_exception, PyStatus_Exception);
  RESOLVE(status_is_exit, PyStatus_IsExit);
  RESOLVE(sys_get_object, PySys_GetObject);
  RESOLVE(unicode_as_utf8, PyUnicode_AsUTF8);
  RESOLVE(unicode_decode_fs_default, PyUnicode_DecodeFSDefault);
  RESOLVE(unicode_from_string, PyUnicode_FromString);
  RESOLVE(finalize_ex, Py_FinalizeEx);
  RESOLVE(get_prefix, Py_GetPrefix);
  RESOLVE(get_program_full_path, Py_GetProgramFullPath);
  RESOLVE(get_version, Py_GetVersion);
  RESOLVE(initialize_from_config, Py_InitializeFromConfig);
  RESOLVE(is_initialized, Py_IsInitialized);
#undef RESOLVE
  return 0;
}

int cp2k_python_api_load_exact(const char *library_path, Cp2kPythonApi *api,
                               char *error, size_t error_size) {
  void *handle;
  if (library_path == NULL || library_path[0] != '/') {
    snprintf(error, error_size,
             "Python library must be an absolute shared-library pathname");
    return -1;
  }
  dlerror();
  /* CPython extension modules commonly resolve Py_* references from the global
   * namespace.  Loading only this exact pathname first with RTLD_GLOBAL makes
   * that work without exposing any other CP2K dependency globally. */
  handle = dlopen(library_path, RTLD_NOW | RTLD_GLOBAL);
  if (handle == NULL) {
    snprintf(error, error_size, "cannot load exact Python library '%s': %s",
             library_path, dlerror());
    return -1;
  }
  if (resolve_api(handle, api, error, error_size) != 0) {
    dlclose(handle);
    return -1;
  }
  python_handle = handle;
  return 0;
}

typedef struct {
  char *path;
  char *error;
  size_t error_size;
  int status;
} PythonInspection;

static int find_python_image(struct dl_phdr_info *info, size_t size,
                             void *data) {
  PythonInspection *inspection = data;
  void *handle;
  __typeof__(&Py_IsInitialized) is_initialized;
  void *address;

  (void)size;
  if (info->dlpi_name == NULL || info->dlpi_name[0] == '\0' ||
      strstr(info->dlpi_name, "libpython") == NULL)
    return 0;

  dlerror();
  handle = dlopen(info->dlpi_name, RTLD_NOW | RTLD_NOLOAD);
  if (handle == NULL) {
    snprintf(inspection->error, inspection->error_size,
             "cannot inspect loaded Python library '%s': %s", info->dlpi_name,
             dlerror());
    inspection->status = -1;
    return 1;
  }
  dlerror();
  address = dlsym(handle, "Py_IsInitialized");
  if (address == NULL) {
    snprintf(inspection->error, inspection->error_size,
             "required CPython symbol Py_IsInitialized is unavailable in "
             "'%s': %s",
             info->dlpi_name, dlerror());
    dlclose(handle);
    inspection->status = -1;
    return 1;
  }
  memcpy(&is_initialized, &address, sizeof(address));
  if (is_initialized()) {
    inspection->path = strdup(info->dlpi_name);
    if (inspection->path == NULL) {
      snprintf(inspection->error, inspection->error_size,
               "cannot record active Python library '%s': %s", info->dlpi_name,
               strerror(errno));
      inspection->status = -1;
    } else {
      inspection->status = 1;
    }
    dlclose(handle);
    return 1;
  }
  dlclose(handle);
  return 0;
}

int cp2k_python_is_preinitialized(char **loaded_library_path, char *error,
                                  size_t error_size) {
  PythonInspection inspection = {NULL, error, error_size, 0};
  __typeof__(&Py_IsInitialized) is_initialized;
  void *address;
  Dl_info symbol_info;

  *loaded_library_path = NULL;

  /* A dynamically linked or symbol-exporting Python executable normally makes
   * its active Py_IsInitialized visible here even when no mapped image pathname
   * contains "libpython". */
  dlerror();
  address = dlsym(RTLD_DEFAULT, "Py_IsInitialized");
  if (address != NULL) {
    memcpy(&is_initialized, &address, sizeof(address));
    if (is_initialized()) {
      if (dladdr(address, &symbol_info) != 0 && symbol_info.dli_fname != NULL &&
          symbol_info.dli_fname[0] != '\0')
        *loaded_library_path = strdup(symbol_info.dli_fname);
      else
        *loaded_library_path = strdup("<process-global Python symbols>");
      if (*loaded_library_path == NULL) {
        snprintf(error, error_size,
                 "cannot record the active process-global Python runtime: %s",
                 strerror(errno));
        return -1;
      }
      return 1;
    }
  }

  /* Also inspect every locally loaded libpython image.  Looking at only the
   * first image can miss an active runtime when another image is inactive. */
  dl_iterate_phdr(find_python_image, &inspection);
  *loaded_library_path = inspection.path;
  return inspection.status;
}

int cp2k_python_validate_executable(const char *executable, char *error,
                                    size_t error_size) {
  struct stat metadata;

  if (executable == NULL || executable[0] != '/') {
    snprintf(error, error_size,
             "Python executable must be an absolute pathname");
    return -1;
  }
  if (stat(executable, &metadata) != 0) {
    snprintf(error, error_size, "cannot access Python executable '%s': %s",
             executable, strerror(errno));
    return -1;
  }
  if (!S_ISREG(metadata.st_mode)) {
    snprintf(error, error_size, "Python executable '%s' is not a regular file",
             executable);
    return -1;
  }
  if (access(executable, X_OK) != 0) {
    snprintf(error, error_size, "Python executable '%s' is not executable: %s",
             executable, strerror(errno));
    return -1;
  }
  return 0;
}
