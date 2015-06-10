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
#ifndef LIBXSTREAM_HPP
#define LIBXSTREAM_HPP

#include <libxstream.h>

#if defined(LIBXSTREAM_EXPORTED) || defined(__LIBXSTREAM)


/** Data type representing a signal. */
typedef uintptr_t libxstream_signal;

typedef void libxstream_lock;
libxstream_lock* libxstream_lock_create();
void libxstream_lock_destroy(libxstream_lock* lock);
void libxstream_lock_acquire(libxstream_lock* lock);
void libxstream_lock_release(libxstream_lock* lock);
bool libxstream_lock_try(libxstream_lock* lock);

size_t nthreads_active();
int this_thread_id();
void this_thread_yield();
void this_thread_sleep(size_t ms = 1);

enum {
  LIBXSTREAM_CALL_UNLOCK    = (2 * (LIBXSTREAM_CALL_INVALID - 1)),
  LIBXSTREAM_CALL_EXTERNAL  = (4 * (LIBXSTREAM_CALL_INVALID - 1))
};

#endif // defined(LIBXSTREAM_EXPORTED) || defined(__LIBXSTREAM)
#endif // LIBXSTREAM_HPP
