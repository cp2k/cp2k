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
#ifndef LIBXSTREAM_STREAM_HPP
#define LIBXSTREAM_STREAM_HPP

#include "libxstream.hpp"

#if defined(LIBXSTREAM_OFFLOAD) && (0 != LIBXSTREAM_OFFLOAD)
# include <offload.h>
#endif

#if defined(LIBXSTREAM_EXPORTED) || defined(__LIBXSTREAM)


struct libxstream_event;


struct libxstream_stream {
public:
  static int enqueue(libxstream_event& event, const libxstream_stream* exclude = 0);

  static int sync(int device);
  static int sync();

public:
  libxstream_stream(int device,
    /**
     * Controls "demuxing" threads and streams i.e., when multiple threads are queuing into the same stream.
     * demux<0: automatic (LIBXSTREAM guesses locks incl. deadlock resolution using LIBXSTREAM_LOCK_RETRY)
     * demux=0: disabled  (application is supposed to call libxstream_stream_lock/libxstream_stream_unlock)
     * demux>0: enabled   (application is supposed to use correct stream synchronization)
     */
    int demux,
    int priority, const char* name);
  ~libxstream_stream();

public:
  int demux() const       { return m_demux; }
  int device() const      { return m_device; }
  int priority() const    { return m_priority; }

  libxstream_signal signal() const;
  int wait(libxstream_signal signal);

  void pending(int thread, libxstream_signal signal);
  libxstream_signal pending(int thread) const;

  int thread() const;
  void begin();
  void end();

  void lock(bool retry);
  void unlock();

#if defined(LIBXSTREAM_OFFLOAD) && (0 != LIBXSTREAM_OFFLOAD) && defined(LIBXSTREAM_ASYNC) && (2 == (2*LIBXSTREAM_ASYNC+1)/2)
  _Offload_stream handle() const;
#endif

#if defined(LIBXSTREAM_PRINT)
  const char* name() const;
#endif

private:
  libxstream_stream(const libxstream_stream& other);
  libxstream_stream& operator=(const libxstream_stream& other);

private:
#if defined(LIBXSTREAM_THREADLOCAL_SIGNALS) && defined(LIBXSTREAM_ASYNC) && (0 != (2*LIBXSTREAM_ASYNC+1)/2)
  libxstream_signal m_pending[LIBXSTREAM_MAX_NTHREADS];
#endif
#if defined(LIBXSTREAM_PRINT)
  char m_name[128];
#endif
  void* m_thread;
#if !(defined(LIBXSTREAM_THREADLOCAL_SIGNALS) && defined(LIBXSTREAM_ASYNC) && (0 != (2*LIBXSTREAM_ASYNC+1)/2))
  libxstream_signal m_signal, *const m_pending;
#endif
#if defined(LIBXSTREAM_LOCK_RETRY) && (0 < (LIBXSTREAM_LOCK_RETRY))
  size_t m_begin, m_end;
#endif
  int m_demux;
  int m_device;
  int m_priority;

#if defined(LIBXSTREAM_OFFLOAD) && (0 != LIBXSTREAM_OFFLOAD) && defined(LIBXSTREAM_ASYNC) && (2 == (2*LIBXSTREAM_ASYNC+1)/2)
  mutable _Offload_stream m_handle; // lazy creation
  mutable size_t m_npartitions;
#endif
};


const libxstream_stream* cast_to_stream(const void* stream);
libxstream_stream* cast_to_stream(void* stream);

const libxstream_stream* cast_to_stream(const libxstream_stream* stream);
libxstream_stream* cast_to_stream(libxstream_stream* stream);

const libxstream_stream* cast_to_stream(const libxstream_stream& stream);
libxstream_stream* cast_to_stream(libxstream_stream& stream);

template<typename T> libxstream_stream* cast_to_stream(T stream) {
  libxstream_use_sink(&stream);
  LIBXSTREAM_ASSERT(0 == stream);
  return static_cast<libxstream_stream*>(0);
}

#endif // defined(LIBXSTREAM_EXPORTED) || defined(__LIBXSTREAM)
#endif // LIBXSTREAM_STREAM_HPP
