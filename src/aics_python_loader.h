/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#ifndef AICS_PYTHON_LOADER_H
#define AICS_PYTHON_LOADER_H

#include <Python.h>

typedef struct {
  __typeof__(&PyCallable_Check) callable_check;
  __typeof__(&PyConfig_Clear) config_clear;
  __typeof__(&PyConfig_InitPythonConfig) config_init_python;
  __typeof__(&PyConfig_Read) config_read;
  __typeof__(&PyConfig_SetBytesString) config_set_bytes_string;
  __typeof__(&Py_DecRef) dec_ref;
  __typeof__(&PyErr_Format) err_format;
  __typeof__(&PyErr_Clear) err_clear;
  __typeof__(&PyErr_Fetch) err_fetch;
  __typeof__(&PyErr_NoMemory) err_no_memory;
  __typeof__(&PyErr_NormalizeException) err_normalize_exception;
  __typeof__(&PyErr_Occurred) err_occurred;
  __typeof__(&PyErr_Print) err_print;
  __typeof__(&PyErr_Restore) err_restore;
  __typeof__(&PyErr_SetString) err_set_string;
  PyObject **exc_runtime_error;
  PyObject **exc_type_error;
  __typeof__(&PyFloat_FromDouble) float_from_double;
  __typeof__(&PyImport_ImportModule) import_module;
  __typeof__(&PyList_Insert) list_insert;
  __typeof__(&PyLong_FromLong) long_from_long;
  __typeof__(&PyObject_CallFunctionObjArgs) object_call_function_obj_args;
  __typeof__(&PyObject_CallObject) object_call_object;
  __typeof__(&PyObject_GetAttrString) object_get_attr_string;
  __typeof__(&PyObject_Str) object_str;
  __typeof__(&PyRun_SimpleString) run_simple_string;
  __typeof__(&PySequence_Contains) sequence_contains;
  __typeof__(&PyStatus_Exception) status_exception;
  __typeof__(&PyStatus_IsExit) status_is_exit;
  __typeof__(&PySys_GetObject) sys_get_object;
  __typeof__(&PyUnicode_DecodeFSDefault) unicode_decode_fs_default;
  __typeof__(&PyUnicode_AsUTF8) unicode_as_utf8;
  __typeof__(&PyUnicode_FromString) unicode_from_string;
  __typeof__(&Py_FinalizeEx) finalize_ex;
  __typeof__(&Py_GetPrefix) get_prefix;
  __typeof__(&Py_GetProgramFullPath) get_program_full_path;
  __typeof__(&Py_GetVersion) get_version;
  __typeof__(&Py_InitializeFromConfig) initialize_from_config;
  __typeof__(&Py_IsInitialized) is_initialized;
} Cp2kPythonApi;

int cp2k_python_api_load_exact(const char *library_path, Cp2kPythonApi *api,
                               char *error, size_t error_size);
int cp2k_python_is_preinitialized(char **loaded_library_path, char *error,
                                  size_t error_size);
int cp2k_python_validate_executable(const char *executable, char *error,
                                    size_t error_size);

#endif
