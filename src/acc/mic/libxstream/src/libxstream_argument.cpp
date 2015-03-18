/******************************************************************************
** Copyright (c) 2014-2015, Intel Corporation                                **
** All rights reserved.                                                      **
**                                                                           **
** Redistribution and use in source and binary forms, with or without        **
** modification, are permitted provided that the following conditions        **
** are met:                                                                  **
** 1. Redistributions of source code must retain the above copyright         **
**    notice, this list of conditions and the following disclaimer.          **
** 2. Redistributions in binary form must reproduce the above copyright      **
**    notice, this list of conditions and the following disclaimer in the    **
**    documentation and/or other materials provided with the distribution.   **
** 3. Neither the name of the copyright holder nor the names of its          **
**    contributors may be used to endorse or promote products derived        **
**    from this software without specific prior written permission.          **
**                                                                           **
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       **
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT         **
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR     **
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT      **
** HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,    **
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED  **
** TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR    **
** PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    **
** LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      **
** NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        **
** SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              **
******************************************************************************/
/* Hans Pabst (Intel Corp.)
******************************************************************************/
#if defined(LIBXSTREAM_EXPORTED) || defined(__LIBXSTREAM)
#include "libxstream_argument.hpp"

#include <libxstream_begin.h>
#include <algorithm>
#include <cstring>
#include <cstdio>
#include <libxstream_end.h>


int libxstream_construct(libxstream_argument arguments[], size_t arg, libxstream_argument::kind_type kind, const void* value, libxstream_type type, size_t dims, const size_t shape[])
{
  size_t typesize = 0;
  const bool weak_candidate = LIBXSTREAM_TYPE_VOID == type || (LIBXSTREAM_ERROR_NONE == libxstream_get_typesize(type, &typesize) && 1 == typesize);
  LIBXSTREAM_CHECK_CONDITION((((libxstream_argument::kind_invalid == kind || libxstream_argument::kind_inout == kind) && LIBXSTREAM_TYPE_INVALID == type) || LIBXSTREAM_TYPE_INVALID > type)
    && ((0 == dims && 0 == shape) || (0 == dims && 0 != shape && weak_candidate) || (0 < dims))
    && (LIBXSTREAM_MAX_NDIMS) >= dims);

  LIBXSTREAM_ASSERT((LIBXSTREAM_MAX_NARGS) >= arg);
  libxstream_argument& argument = arguments[arg];

#if defined(LIBXSTREAM_DEBUG)
  memset(argument.data.self, 0, sizeof(libxstream_argument)); // avoid false pos. with mem. analysis
#endif
#if defined(LIBXSTREAM_PRINT)
  static const char *const context[] = { "", "input", "output", "inout" };
#endif
  argument.kind = kind;
  argument.dims = dims;

  if (shape) {
    if (0 < dims || !weak_candidate) {
#if defined(LIBXSTREAM_PRINT)
      if (0 == dims && !weak_candidate) {
        LIBXSTREAM_PRINT_WARN("libxstream_fn_%s: signature=0x%llx arg=%lu is strong-typed (ignored shape)!",
          context[kind], reinterpret_cast<unsigned long long>(arguments), static_cast<unsigned long>(arg));
      }
#endif
      argument.type = type;
    }
    else { // 0 == dims && weak_candidate
      argument.type = LIBXSTREAM_TYPE_VOID;
      LIBXSTREAM_CHECK_CONDITION(sizeof(libxstream_argument::data_union) >= *shape);
      argument.shape[0] = shape[0];
    }

#if defined(__INTEL_COMPILER)
#   pragma loop_count min(0), max(LIBXSTREAM_MAX_NDIMS), avg(2)
#endif
    for (size_t i = 0; i < dims; ++i) argument.shape[i] = shape[i];
  }
  else {
#if defined(LIBXSTREAM_PRINT)
    if (0 < dims && 0 == shape) {
      LIBXSTREAM_PRINT_WARN("libxstream_fn_%s: signature=0x%llx arg=%lu is weak-typed (no shape information)!",
        context[kind], reinterpret_cast<unsigned long long>(arguments), static_cast<unsigned long>(arg));
    }
#endif
    std::fill_n(argument.shape, dims, 0);
    argument.type = type;
  }

  return libxstream_argument::kind_invalid != kind ? libxstream_set_value(argument, value) : LIBXSTREAM_ERROR_NONE;
}


int libxstream_construct(libxstream_argument* signature, size_t nargs)
{
  LIBXSTREAM_CHECK_CONDITION((0 != signature || 0 == nargs) && (LIBXSTREAM_MAX_NARGS) >= nargs);

  if (0 != signature) {
#if defined(__INTEL_COMPILER)
#   pragma loop_count min(0), max(LIBXSTREAM_MAX_NARGS), avg(LIBXSTREAM_MAX_NARGS/2)
#endif
    for (size_t i = 0; i < nargs; ++i) {
      LIBXSTREAM_CHECK_CALL(libxstream_construct(signature, i, libxstream_argument::kind_inout, 0, LIBXSTREAM_TYPE_INVALID, 0, 0));
    }
    LIBXSTREAM_CHECK_CALL(libxstream_construct(signature, nargs, libxstream_argument::kind_invalid, 0, LIBXSTREAM_TYPE_INVALID, 0, 0));
  }

  return LIBXSTREAM_ERROR_NONE;
}


LIBXSTREAM_TARGET(mic) libxstream_argument::value_union libxstream_get_value(const libxstream_argument& arg, bool byvalue)
{
  const void* data = 0;

  if (0 != arg.dims || 0 != (libxstream_argument::kind_output & arg.kind) || (byvalue &&
    (LIBXSTREAM_TYPE_C32 > arg.type || (LIBXSTREAM_TYPE_VOID == arg.type && sizeof(void*) >= *arg.shape))))
  {
    char *const dst = reinterpret_cast<char*>(&data);
    for (size_t i = 0; i < sizeof(void*); ++i) dst[i] = arg.data.self[i];
  }
  else { // by-value
    data = arg.data.self; // pointer to data
  }

  libxstream_argument::value_union value;
  value.const_pointer = data;
  return value;
}


LIBXSTREAM_TARGET(mic) int libxstream_set_value(libxstream_argument& arg, const void* data)
{
  if (0 != arg.dims || 0 != (libxstream_argument::kind_output & arg.kind)) { // by-pointer
    *reinterpret_cast<const void**>(&arg) = data;
    LIBXSTREAM_ASSERT(libxstream_get_value(arg).const_pointer == data);
  }
  else { // by-value: copy from the given pointer
    size_t typesize = 0 != arg.dims ? 1 : *arg.shape;
    if (LIBXSTREAM_TYPE_VOID != arg.type) {
      LIBXSTREAM_CHECK_CALL(libxstream_get_typesize(arg.type, &typesize));
    }

    if (data) {
      const char *const src = static_cast<const char*>(data);
      for (size_t i = 0; i < typesize; ++i) arg.data.self[i] = src[i];
    }
    else {
      for (size_t i = 0; i < typesize; ++i) arg.data.self[i] = 0;
    }

    // allows to promote smaller types to pointer-size
    for (size_t i = typesize; i < sizeof(void*); ++i) arg.data.self[i] = 0;
  }

  return LIBXSTREAM_ERROR_NONE;
}

#endif // defined(LIBXSTREAM_EXPORTED) || defined(__LIBXSTREAM)
