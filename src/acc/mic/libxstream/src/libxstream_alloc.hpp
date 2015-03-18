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
#ifndef LIBXSTREAM_ALLOC_HPP
#define LIBXSTREAM_ALLOC_HPP

#include <libxstream_macros.h>
#include <libxstream_begin.h>
#include <cstddef>
#include <libxstream_end.h>

#if defined(LIBXSTREAM_EXPORTED) || defined(__LIBXSTREAM)


LIBXSTREAM_TARGET(mic) size_t libxstream_gcd(size_t a, size_t b);
LIBXSTREAM_TARGET(mic) size_t libxstream_lcm(size_t a, size_t b);

LIBXSTREAM_TARGET(mic) size_t libxstream_alignment(size_t size, size_t alignment);
LIBXSTREAM_TARGET(mic) size_t libxstream_align(size_t size, size_t alignment);
LIBXSTREAM_TARGET(mic) const void* libxstream_align(const void* address, size_t alignment);
LIBXSTREAM_TARGET(mic) void* libxstream_align(void* address, size_t alignment);

LIBXSTREAM_TARGET(mic) size_t libxstream_linear_size(size_t dims, const size_t shape[], size_t initial_size = 1);
LIBXSTREAM_TARGET(mic) int libxstream_linear_offset(size_t dims, const int offset[], const size_t shape[]);
LIBXSTREAM_TARGET(mic) size_t libxstream_linear_address(size_t dims, const int offset[], const size_t shape[], const size_t pitch[]);

int libxstream_real_allocate(void** memory, size_t size, size_t alignment);
int libxstream_real_deallocate(const void* memory);

int libxstream_virt_allocate(void** memory, size_t size, size_t alignment, const void* data = 0, size_t data_size = 0);
int libxstream_virt_deallocate(const void* memory);

void* libxstream_virt_data(void* memory);
const void* libxstream_virt_data(const void* memory);

#endif // defined(LIBXSTREAM_EXPORTED) || defined(__LIBXSTREAM)
#endif // LIBXSTREAM_ALLOC_HPP
