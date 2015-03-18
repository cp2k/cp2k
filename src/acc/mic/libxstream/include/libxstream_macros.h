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
#ifndef LIBXSTREAM_MACROS_H
#define LIBXSTREAM_MACROS_H

#include "libxstream_config.h"

#if !defined(LIBXSTREAM_PREFER_OPENMP) || !defined(_OPENMP)
# if (201103L <= __cplusplus)
#   if !defined(LIBXSTREAM_STDFEATURES)
#     define LIBXSTREAM_STDFEATURES
#   endif
#   if !defined(LIBXSTREAM_STDFEATURES_THREADX) && !defined(__MIC__)
#     define LIBXSTREAM_STDFEATURES_THREADX
#   endif
# elif (1600 < _MSC_VER)
#   if !defined(LIBXSTREAM_STDFEATURES)
#     define LIBXSTREAM_STDFEATURES
#   endif
#   if !defined(LIBXSTREAM_STDFEATURES_THREADX)
#     define LIBXSTREAM_STDFEATURES_THREADX
#   endif
# elif ((4 <= __GNUC__ && 5 <= __GNUC_MINOR__) && (1L == __cplusplus)) || (defined(__INTEL_COMPILER) && defined(__GXX_EXPERIMENTAL_CXX0X__))
#   if !defined(LIBXSTREAM_STDFEATURES)
#     define LIBXSTREAM_STDFEATURES
#   endif
# endif
#endif

#if defined(_WIN32) && !defined(__GNUC__)
# define LIBXSTREAM_ATTRIBUTE(A) __declspec(A)
# define LIBXSTREAM_ALIGNED(DECL, N) LIBXSTREAM_ATTRIBUTE(align(N)) DECL
# define LIBXSTREAM_CDECL __cdecl
#elif defined(__GNUC__)
# define LIBXSTREAM_ATTRIBUTE(A) __attribute__((A))
# define LIBXSTREAM_ALIGNED(DECL, N) DECL LIBXSTREAM_ATTRIBUTE(aligned(N))
# define LIBXSTREAM_CDECL LIBXSTREAM_ATTRIBUTE(cdecl)
#endif
#define LIBXSTREAM_ALIGN(A, ALIGNMENT) ((uintptr_t)(-(-((intptr_t)A) & -((intptr_t)(ALIGNMENT)))))

#if defined(_WIN32) && !defined(__GNUC__)
# define LIBXSTREAM_TLS LIBXSTREAM_ATTRIBUTE(thread)
#elif defined(__GNUC__)
# define LIBXSTREAM_TLS __thread
#elif defined(LIBXSTREAM_STDFEATURES)
# define LIBXSTREAM_TLS thread_local
#endif
#if !defined(LIBXSTREAM_TLS)
# define LIBXSTREAM_TLS
#endif

#define LIBXSTREAM_TOSTRING_AUX(SYMBOL) #SYMBOL
#define LIBXSTREAM_TOSTRING(SYMBOL) LIBXSTREAM_TOSTRING_AUX(SYMBOL)

#if defined(__INTEL_COMPILER)
# define LIBXSTREAM_ASSUME_ALIGNED(A, N) __assume_aligned(A, N)
# define LIBXSTREAM_PRAGMA(DIRECTIVE) __pragma(DIRECTIVE)
#elif (199901L <= __STDC_VERSION__)
# define LIBXSTREAM_ASSUME_ALIGNED(A, N)
# define LIBXSTREAM_PRAGMA(DIRECTIVE) _Pragma(LIBXSTREAM_STRINGIFY(DIRECTIVE))
#elif defined(_MSC_VER)
# define LIBXSTREAM_ASSUME_ALIGNED(A, N)
# define LIBXSTREAM_PRAGMA(DIRECTIVE) __pragma(DIRECTIVE)
#else
# define LIBXSTREAM_ASSUME_ALIGNED(A, N)
# define LIBXSTREAM_PRAGMA(DIRECTIVE)
#endif

#if defined(__INTEL_OFFLOAD) && (!defined(_WIN32) || (1400 <= __INTEL_COMPILER))
# define LIBXSTREAM_OFFLOAD 1
# define LIBXSTREAM_TARGET(A) LIBXSTREAM_ATTRIBUTE(target(A))
#else
/*# define LIBXSTREAM_OFFLOAD 0*/
# define LIBXSTREAM_TARGET(A)
#endif

#define LIBXSTREAM_IMPORT_DLL __declspec(dllimport)

#if defined(_WINDLL) && defined(_WIN32)
# if defined(LIBXSTREAM_EXPORTED)
#   define LIBXSTREAM_EXPORT __declspec(dllexport)
# else
#   define LIBXSTREAM_EXPORT LIBXSTREAM_IMPORT_DLL
# endif
#else
# define LIBXSTREAM_EXPORT
#endif

#if defined(__cplusplus)
# define LIBXSTREAM_EXTERN_C extern "C"
# define LIBXSTREAM_EXPORT_C LIBXSTREAM_EXTERN_C LIBXSTREAM_EXPORT
# define LIBXSTREAM_VARIADIC ...
#else
# define LIBXSTREAM_EXTERN_C
# define LIBXSTREAM_EXPORT_C LIBXSTREAM_EXPORT
# define LIBXSTREAM_VARIADIC
#endif /*__cplusplus*/

#if defined(__GNUC__) && !defined(_WIN32) && !defined(__CYGWIN32__)
# define LIBXSTREAM_RESTRICT __restrict__
#elif defined(_MSC_VER)
# define LIBXSTREAM_RESTRICT __restrict
#else
# define LIBXSTREAM_RESTRICT
#endif

#if (defined(LIBXSTREAM_ERROR_DEBUG) || defined(_DEBUG)) && !defined(NDEBUG) && !defined(LIBXSTREAM_DEBUG)
# define LIBXSTREAM_DEBUG
#endif

#if defined(LIBXSTREAM_ERROR_CHECK) && !defined(LIBXSTREAM_CHECK)
# define LIBXSTREAM_CHECK
#endif

LIBXSTREAM_EXPORT_C LIBXSTREAM_TARGET(mic) void libxstream_use_sink(const void*);
LIBXSTREAM_EXPORT_C LIBXSTREAM_TARGET(mic) int libxstream_not_constant(int value);

#if defined(LIBXSTREAM_DEBUG)
# define LIBXSTREAM_USE_SINK(VAR) libxstream_use_sink(VAR)
# define LIBXSTREAM_ASSERT(A) assert(A)
# include "libxstream_begin.h"
# include <assert.h>
# include "libxstream_end.h"
#else
# define LIBXSTREAM_USE_SINK(VAR)
# define LIBXSTREAM_ASSERT(A)
#endif

#define LIBXSTREAM_TRUE  1
#define LIBXSTREAM_FALSE 0

#define LIBXSTREAM_ERROR_NONE       0
#define LIBXSTREAM_ERROR_RUNTIME   -1
#define LIBXSTREAM_ERROR_CONDITION -2

#if defined(LIBXSTREAM_TRACE) && ((1 == ((2*LIBXSTREAM_TRACE+1)/2) && defined(LIBXSTREAM_DEBUG)) || 1 < ((2*LIBXSTREAM_TRACE+1)/2))
# define LIBXSTREAM_PRINT
# define LIBXSTREAM_PRINT_INFO(MESSAGE, ...) fprintf(stderr, "DBG " MESSAGE "\n", __VA_ARGS__)
# define LIBXSTREAM_PRINT_INFO0(MESSAGE) fprintf(stderr, "DBG " MESSAGE "\n")
# define LIBXSTREAM_PRINT_INFOCTX(MESSAGE, ...) fprintf(stderr, "DBG %s: " MESSAGE "\n", __FUNCTION__, __VA_ARGS__)
# define LIBXSTREAM_PRINT_INFOCTX0(MESSAGE) fprintf(stderr, "DBG %s: " MESSAGE "\n", __FUNCTION__)
# define LIBXSTREAM_PRINT_WARN(MESSAGE, ...) fprintf(stderr, "WRN " MESSAGE "\n", __VA_ARGS__)
# define LIBXSTREAM_PRINT_WARN0(MESSAGE) fprintf(stderr, "WRN " MESSAGE "\n")
# define LIBXSTREAM_PRINT_WARNCTX(MESSAGE, ...) fprintf(stderr, "WRN %s: " MESSAGE "\n", __FUNCTION__, __VA_ARGS__)
# define LIBXSTREAM_PRINT_WARNCTX0(MESSAGE) fprintf(stderr, "WRN %s: " MESSAGE "\n", __FUNCTION__)
#else
# define LIBXSTREAM_PRINT_INFO(MESSAGE, ...)
# define LIBXSTREAM_PRINT_INFO0(MESSAGE)
# define LIBXSTREAM_PRINT_INFOCTX(MESSAGE, ...)
# define LIBXSTREAM_PRINT_INFOCTX0(MESSAGE)
# define LIBXSTREAM_PRINT_WARN(MESSAGE, ...)
# define LIBXSTREAM_PRINT_WARN0(MESSAGE)
# define LIBXSTREAM_PRINT_WARNCTX(MESSAGE, ...)
# define LIBXSTREAM_PRINT_WARNCTX0(MESSAGE)
#endif

#if defined(_MSC_VER)
# define LIBXSTREAM_SNPRINTF(S, N, F, ...) _snprintf_s(S, N, _TRUNCATE, F, __VA_ARGS__)
#else
# define LIBXSTREAM_SNPRINTF(S, N, F, ...) snprintf(S, N, F, __VA_ARGS__)
#endif

#define LIBXSTREAM_MIN(A, B) ((A) < (B) ? (A) : (B))
#define LIBXSTREAM_MAX(A, B) ((A) < (B) ? (B) : (A))

#if defined(LIBXSTREAM_CHECK)
# define LIBXSTREAM_CHECK_ERROR(RETURN_VALUE) if (LIBXSTREAM_ERROR_NONE != (RETURN_VALUE)) return RETURN_VALUE
# define LIBXSTREAM_CHECK_CONDITION(CONDITION) if (!(CONDITION)) return LIBXSTREAM_ERROR_CONDITION
# define LIBXSTREAM_CHECK_CONDITION_RETURN(CONDITION) if (!(CONDITION)) return
# define LIBXSTREAM_CHECK_CALL_RETURN(EXPRESSION) if (LIBXSTREAM_ERROR_NONE != (EXPRESSION)) return
# if defined(__cplusplus)
#   define LIBXSTREAM_CHECK_CALL_THROW(EXPRESSION) if (LIBXSTREAM_ERROR_NONE != (EXPRESSION)) throw std::runtime_error(LIBXSTREAM_TOSTRING(EXPRESSION) " at " __FILE__ ":" LIBXSTREAM_TOSTRING(__LINE__))
#   define LIBXSTREAM_CHECK_CONDITION_THROW(CONDITION) if (!(CONDITION)) throw std::runtime_error(LIBXSTREAM_TOSTRING(CONDITION) " at " __FILE__ ":" LIBXSTREAM_TOSTRING(__LINE__))
# else
#   define LIBXSTREAM_CHECK_CALL_THROW(EXPRESSION) do { int result = (EXPRESSION); if (LIBXSTREAM_ERROR_NONE != result) abort(result); } while(libxstream_not_constant(LIBXSTREAM_FALSE))
#   define LIBXSTREAM_CHECK_CONDITION_THROW(CONDITION) if (!(CONDITION)) abort(1)
# endif
# if defined(_OPENMP)
#   if defined(LIBXSTREAM_DEBUG)
#     define LIBXSTREAM_CHECK_CALL(EXPRESSION) LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == (EXPRESSION))
#   else
#     define LIBXSTREAM_CHECK_CALL(EXPRESSION) (EXPRESSION)
#   endif
# else
#   define LIBXSTREAM_CHECK_CALL(EXPRESSION) do { int result = (EXPRESSION); if (LIBXSTREAM_ERROR_NONE != result) return result; } while(libxstream_not_constant(LIBXSTREAM_FALSE))
# endif
#elif defined(LIBXSTREAM_DEBUG)
# define LIBXSTREAM_CHECK_ERROR(RETURN_VALUE) LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == (RETURN_VALUE))
# define LIBXSTREAM_CHECK_CONDITION(CONDITION) LIBXSTREAM_ASSERT(CONDITION)
# define LIBXSTREAM_CHECK_CONDITION_RETURN(CONDITION) LIBXSTREAM_ASSERT(CONDITION)
# define LIBXSTREAM_CHECK_CONDITION_THROW(CONDITION) LIBXSTREAM_ASSERT(CONDITION)
# define LIBXSTREAM_CHECK_CALL_RETURN(EXPRESSION) LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == (EXPRESSION))
# define LIBXSTREAM_CHECK_CALL_THROW(EXPRESSION) LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == (EXPRESSION))
# define LIBXSTREAM_CHECK_CALL(EXPRESSION) LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == (EXPRESSION))
#else
# define LIBXSTREAM_CHECK_ERROR(RETURN_VALUE) LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == (RETURN_VALUE))
# define LIBXSTREAM_CHECK_CONDITION(CONDITION) LIBXSTREAM_ASSERT(CONDITION)
# define LIBXSTREAM_CHECK_CONDITION_RETURN(CONDITION) LIBXSTREAM_ASSERT(CONDITION)
# define LIBXSTREAM_CHECK_CONDITION_THROW(CONDITION) LIBXSTREAM_ASSERT(CONDITION)
# define LIBXSTREAM_CHECK_CALL_RETURN(EXPRESSION) EXPRESSION
# define LIBXSTREAM_CHECK_CALL_THROW(EXPRESSION) EXPRESSION
# define LIBXSTREAM_CHECK_CALL(EXPRESSION) EXPRESSION
#endif
#if defined(LIBXSTREAM_DEBUG)
# define LIBXSTREAM_CHECK_CALL_ASSERT(EXPRESSION) LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == (EXPRESSION))
#else
# define LIBXSTREAM_CHECK_CALL_ASSERT(EXPRESSION) EXPRESSION
#endif

#if defined(LIBXSTREAM_CALL_BYVALUE)
# define LIBXSTREAM_INVAL(TYPE) TYPE
# define LIBXSTREAM_GETVAL(VALUE) VALUE
# define LIBXSTREAM_SETVAL(VALUE) VALUE
#else /*by-pointer*/
# if defined(__cplusplus)
#   define LIBXSTREAM_INVAL(TYPE) const TYPE&
#   define LIBXSTREAM_GETVAL(VALUE) VALUE
#   define LIBXSTREAM_SETVAL(VALUE) VALUE
# else
#   define LIBXSTREAM_INVAL(TYPE) const TYPE*
#   define LIBXSTREAM_GETVAL(VALUE) *VALUE
#   define LIBXSTREAM_SETVAL(VALUE) &VALUE
# endif
#endif

#endif /*LIBXSTREAM_MACROS_H*/
