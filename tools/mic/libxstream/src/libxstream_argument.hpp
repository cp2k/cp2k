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
#ifndef LIBXSTREAM_ARGUMENT_HPP
#define LIBXSTREAM_ARGUMENT_HPP

#include <libxstream.h>

#if defined(LIBXSTREAM_EXPORTED) || defined(__LIBXSTREAM)


template<size_t N> struct libxstream_argument_value {};
template<> struct libxstream_argument_value<4> { typedef float type; };
template<> struct libxstream_argument_value<8> { typedef double type; };


extern "C" struct LIBXSTREAM_TARGET(mic) libxstream_argument {
  enum kind_type {
    kind_invalid  = 0,
    kind_input    = 1,
    kind_output   = 2,
    kind_inout    = kind_output | kind_input,
  };

  // This data member *must* be the first!
  union data_union {
    char self[sizeof(void*)], c;
    uintptr_t pointer;
    int8_t    i8;
    uint8_t   u8;
    int16_t   i16;
    uint16_t  u16;
    int32_t   i32;
    uint32_t  u32;
    int64_t   i64;
    uint64_t  u64;
    float     f32, c32[2];
    double    f64, c64[2];
  } data;

  union value_union {
    void* pointer;
    const void* const_pointer;
    typedef libxstream_argument_value<sizeof(void*)>::type value_type;
    value_type value;
  };

  size_t shape[LIBXSTREAM_MAX_NDIMS];
  kind_type kind;
  libxstream_type type;
  size_t dims;
};


int libxstream_construct(libxstream_argument arguments[], size_t arg, libxstream_argument::kind_type kind, const void* value, libxstream_type type, size_t dims, const size_t shape[]);
int libxstream_construct(libxstream_argument* signature, size_t nargs);

LIBXSTREAM_TARGET(mic) libxstream_argument::value_union libxstream_get_value(const libxstream_argument& arg,
#if defined(LIBXSTREAM_CALL_BYVALUE)
  bool byvalue = true);
#else
  bool byvalue = false);
#endif
LIBXSTREAM_TARGET(mic) int libxstream_set_value(libxstream_argument& arg, const void* data);

#endif // defined(LIBXSTREAM_EXPORTED) || defined(__LIBXSTREAM)
#endif // LIBXSTREAM_ARGUMENT_HPP
