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
#include "libxstream_stream.hpp"
#include "libxstream_capture.hpp"
#include "libxstream_event.hpp"

#include <libxstream_begin.h>
#include <algorithm>
#include <string>
#include <cstdio>
#if defined(LIBXSTREAM_STDFEATURES)
# include <atomic>
#endif
#include <libxstream_end.h>

// allows to wait for an event issued prior to the pending signal
//#define LIBXSTREAM_STREAM_WAIT_PAST
// unlocking the stream is only allowed for the thread that locked
//#define LIBXSTREAM_STREAM_UNLOCK_OWNER


namespace libxstream_stream_internal {

class registry_type {
public:
  registry_type()
    : m_istreams(0)
#if !defined(LIBXSTREAM_STDFEATURES)
    , m_lock(libxstream_lock_create())
#endif
  {
    std::fill_n(m_signals, LIBXSTREAM_MAX_NDEVICES, 0);
    std::fill_n(m_streams, LIBXSTREAM_MAX_NDEVICES * LIBXSTREAM_MAX_NSTREAMS, static_cast<libxstream_stream*>(0));
  }

  ~registry_type() {
    const size_t n = max_nstreams();
    for (size_t i = 0; i < n; ++i) {
#if defined(LIBXSTREAM_DEBUG)
      if (0 != m_streams[i]) {
        LIBXSTREAM_PRINT_WARN("dangling stream \"%s\"!", m_streams[i]->name());
      }
#endif
      libxstream_stream_destroy(m_streams[i]);
    }
#if !defined(LIBXSTREAM_STDFEATURES)
    libxstream_lock_destroy(m_lock);
#endif
  }

public:
  libxstream_stream** allocate() {
#if !defined(LIBXSTREAM_STDFEATURES)
    libxstream_lock_acquire(m_lock);
#endif
    libxstream_stream** i = m_streams + (m_istreams++ % (LIBXSTREAM_MAX_NDEVICES * LIBXSTREAM_MAX_NSTREAMS));
    while (0 != *i) i = m_streams + (m_istreams++ % (LIBXSTREAM_MAX_NDEVICES * LIBXSTREAM_MAX_NSTREAMS));
#if !defined(LIBXSTREAM_STDFEATURES)
    libxstream_lock_release(m_lock);
#endif
    return i;
  }

  size_t max_nstreams() const {
    return std::min<size_t>(m_istreams, LIBXSTREAM_MAX_NDEVICES * LIBXSTREAM_MAX_NSTREAMS);
  }

  size_t nstreams(int device) const {
    const size_t n = max_nstreams();
    size_t result = 0;
    for (size_t i = 0; i < n; ++i) {
      result += (0 != m_streams[i] && m_streams[i]->device() == device) ? 1 : 0;
    }
    return result;
  }

  size_t nstreams() const {
    const size_t n = max_nstreams();
    size_t result = 0;
    for (size_t i = 0; i < n; ++i) {
      result += 0 != m_streams[i] ? 1 : 0;
    }
    return result;
  }

  libxstream_signal& signal(int device) {
    LIBXSTREAM_ASSERT(-1 <= device && device <= LIBXSTREAM_MAX_NDEVICES);
    return m_signals[device+1];
  }

  libxstream_stream** streams() {
    return m_streams;
  }

  int enqueue(libxstream_event& event, const libxstream_stream* exclude) {
    LIBXSTREAM_ASSERT(0 == event.expected());
    int result = LIBXSTREAM_ERROR_NONE;
    const size_t n = max_nstreams();
    bool reset = true;

    for (size_t i = 0; i < n; ++i) {
      libxstream_stream *const stream = m_streams[i];

      if (stream != exclude) {
        result = event.enqueue(*stream, reset);
        LIBXSTREAM_CHECK_ERROR(result);
        reset = false;
      }
    }
    if (reset) {
      result = event.reset();
    }

    return result;
  }

  int sync(int device) {
    const size_t n = max_nstreams();
    for (size_t i = 0; i < n; ++i) {
      if (libxstream_stream *const stream = m_streams[i]) {
        const int stream_device = stream->device();
        if (stream_device == device) {
          const int result = stream->wait(0);
          LIBXSTREAM_CHECK_ERROR(result);
        }
      }
    }
    return LIBXSTREAM_ERROR_NONE;
  }

  int sync() {
    const size_t n = max_nstreams();
    for (size_t i = 0; i < n; ++i) {
      if (libxstream_stream *const stream = m_streams[i]) {
        const int result = stream->wait(0);
        LIBXSTREAM_CHECK_ERROR(result);
      }
    }
    return LIBXSTREAM_ERROR_NONE;
  }

#if !defined(LIBXSTREAM_STDFEATURES)
  libxstream_lock* lock() { return m_lock; }
#endif

private:
  // not necessary to be device-specific due to single-threaded offload
  libxstream_signal m_signals[LIBXSTREAM_MAX_NDEVICES + 1];
  libxstream_stream* m_streams[LIBXSTREAM_MAX_NDEVICES*LIBXSTREAM_MAX_NSTREAMS];
#if defined(LIBXSTREAM_STDFEATURES)
  std::atomic<size_t> m_istreams;
#else
  size_t m_istreams;
  libxstream_lock* m_lock;
#endif
} registry;


template<typename A, typename E, typename D>
bool atomic_compare_exchange(A& atomic, E& expected, D desired)
{
#if defined(LIBXSTREAM_STDFEATURES)
  const bool result = std::atomic_compare_exchange_weak(&atomic, &expected, desired);
#elif defined(_OPENMP)
  bool result = false;
# pragma omp critical
  {
    result = atomic == expected;
    if (result) {
      atomic = desired;
    }
    else {
      expected = atomic;
    }
  }
#else // generic
  bool result = false;
  libxstream_lock_acquire(registry.lock());
  result = atomic == expected;
  if (result) {
    atomic = desired;
  }
  else {
    expected = atomic;
  }
  libxstream_lock_release(registry.lock());
#endif
  return result;
}


template<typename A, typename T>
T atomic_store(A& atomic, T value)
{
  T result = value;
#if defined(LIBXSTREAM_STDFEATURES)
  result = std::atomic_exchange(&atomic, value);
#elif defined(_OPENMP)
# pragma omp critical
  {
    result = atomic;
    atomic = value;
  }
#else // generic
  libxstream_lock_acquire(registry.lock());
  result = atomic;
  atomic = value;
  libxstream_lock_release(registry.lock());
#endif
  return result;
}

} // namespace libxstream_stream_internal


/*static*/int libxstream_stream::enqueue(libxstream_event& event, const libxstream_stream* exclude)
{
  return libxstream_stream_internal::registry.enqueue(event, exclude);
}


/*static*/int libxstream_stream::sync(int device)
{
  return libxstream_stream_internal::registry.sync(device);
}


/*static*/int libxstream_stream::sync()
{
  return libxstream_stream_internal::registry.sync();
}


libxstream_stream::libxstream_stream(int device, int demux, int priority, const char* name)
#if defined(LIBXSTREAM_STDFEATURES)
  : m_thread(new std::atomic<int>(-1))
#else
  : m_thread(new int(-1))
#endif
#if !(defined(LIBXSTREAM_THREADLOCAL_SIGNALS) && defined(LIBXSTREAM_ASYNC) && (0 != (2*LIBXSTREAM_ASYNC+1)/2))
  , m_signal(0), m_pending(&m_signal)
#endif
#if defined(LIBXSTREAM_LOCK_RETRY) && (0 < (LIBXSTREAM_LOCK_RETRY))
  , m_begin(0), m_end(0)
#endif
  , m_demux(demux)
  , m_device(device), m_priority(priority)
#if defined(LIBXSTREAM_OFFLOAD) && defined(LIBXSTREAM_ASYNC) && (2 == (2*LIBXSTREAM_ASYNC+1)/2)
  , m_handle(0) // lazy creation
  , m_npartitions(0)
#endif
{
  libxstream_use_sink(name);
#if defined(LIBXSTREAM_THREADLOCAL_SIGNALS) && defined(LIBXSTREAM_ASYNC) && (0 != (2*LIBXSTREAM_ASYNC+1)/2)
  std::fill_n(m_pending, LIBXSTREAM_MAX_NTHREADS, static_cast<libxstream_signal>(0));
#endif
#if defined(LIBXSTREAM_PRINT)
  if (name && 0 != *name) {
    const size_t length = std::min(std::char_traits<char>::length(name), sizeof(m_name) - 1);
    std::copy(name, name + length, m_name);
    m_name[length] = 0;
  }
  else {
    m_name[0] = 0;
  }
#endif
  using namespace libxstream_stream_internal;
  libxstream_stream* *const slot = libxstream_stream_internal::registry.allocate();
  *slot = this;
}


libxstream_stream::~libxstream_stream()
{
  using namespace libxstream_stream_internal;
  libxstream_stream* *const end = registry.streams() + registry.max_nstreams();
  libxstream_stream* *const stream = std::find(registry.streams(), end, this);
  LIBXSTREAM_ASSERT(stream != end);
  *stream = 0; // unregister stream
#if defined(LIBXSTREAM_STDFEATURES)
  delete static_cast<std::atomic<int>*>(m_thread);
#else
  delete static_cast<int*>(m_thread);
#endif
#if defined(LIBXSTREAM_OFFLOAD) && (0 != LIBXSTREAM_OFFLOAD) && !defined(__MIC__) && defined(LIBXSTREAM_ASYNC) && (2 == (2*LIBXSTREAM_ASYNC+1)/2)
  if (0 != m_handle) {
    _Offload_stream_destroy(m_device, m_handle);
  }
#endif
}


libxstream_signal libxstream_stream::signal() const
{
  return ++libxstream_stream_internal::registry.signal(m_device);
}


int libxstream_stream::wait(libxstream_signal signal)
{
  int result = LIBXSTREAM_ERROR_NONE;

  LIBXSTREAM_ASYNC_BEGIN(this, m_pending, signal)
  {
    libxstream_signal *const pending_signals = ptr<libxstream_signal,0>();
    const libxstream_signal signal = val<const libxstream_signal,1>();

#if defined(LIBXSTREAM_THREADLOCAL_SIGNALS) && defined(LIBXSTREAM_ASYNC) && (0 != (2*LIBXSTREAM_ASYNC+1)/2)
    const int nthreads = static_cast<int>(nthreads_active());
    for (int i = 0; i < nthreads; ++i)
#else
    const int i = thread();
#endif
    {
      const libxstream_signal pending_signal = pending_signals[i];
      if (0 != pending_signal) {
#if defined(LIBXSTREAM_OFFLOAD) && defined(LIBXSTREAM_ASYNC) && (0 != (2*LIBXSTREAM_ASYNC+1)/2)
        if (0 <= LIBXSTREAM_ASYNC_DEVICE) {
# if defined(LIBXSTREAM_STREAM_WAIT_PAST)
          const libxstream_signal wait_pending = 0 != signal ? signal : pending_signal;
# else
          const libxstream_signal wait_pending = pending_signal;
# endif
#         pragma offload_wait LIBXSTREAM_ASYNC_TARGET wait(wait_pending)
        }
#endif
        if (0 == signal) {
          pending_signals[i] = 0;
        }
#if defined(LIBXSTREAM_STREAM_WAIT_PAST)
        else {
          i = nthreads; // break
        }
#endif
      }
    }

    if (0 != signal) {
#if defined(LIBXSTREAM_THREADLOCAL_SIGNALS) && defined(LIBXSTREAM_ASYNC) && (0 != (2*LIBXSTREAM_ASYNC+1)/2)
      for (int i = 0; i < nthreads; ++i)
#else
      const int i = thread();
#endif
      {
        if (signal == pending_signals[i]) {
          pending_signals[i] = 0;
        }
      }
    }
  }
  LIBXSTREAM_ASYNC_END(LIBXSTREAM_CALL_DEFAULT | LIBXSTREAM_CALL_WAIT | LIBXSTREAM_CALL_UNLOCK, result);

  return result;
}


void libxstream_stream::pending(int thread, libxstream_signal signal)
{
#if defined(LIBXSTREAM_THREADLOCAL_SIGNALS) && defined(LIBXSTREAM_ASYNC) && (0 != (2*LIBXSTREAM_ASYNC+1)/2)
  LIBXSTREAM_ASSERT(0 <= thread && thread < LIBXSTREAM_MAX_NTHREADS);
  m_pending[thread] = signal;
#else
  libxstream_use_sink(&thread);
  LIBXSTREAM_ASSERT(0 == thread);
  m_signal = signal;
#endif
}


libxstream_signal libxstream_stream::pending(int thread) const
{
#if defined(LIBXSTREAM_THREADLOCAL_SIGNALS) && defined(LIBXSTREAM_ASYNC) && (0 != (2*LIBXSTREAM_ASYNC+1)/2)
  LIBXSTREAM_ASSERT(0 <= thread && thread < LIBXSTREAM_MAX_NTHREADS);
  const libxstream_signal signal = m_pending[thread];
#else
  libxstream_use_sink(&thread);
  LIBXSTREAM_ASSERT(0 == thread);
  const libxstream_signal signal = m_signal;
#endif
  return signal;
}


int libxstream_stream::thread() const
{
  LIBXSTREAM_ASSERT(m_thread);
#if defined(LIBXSTREAM_STDFEATURES)
  return *static_cast<std::atomic<int>*>(m_thread);
#else
  return *static_cast<int*>(m_thread);
#endif
}


void libxstream_stream::begin()
{
#if defined(LIBXSTREAM_LOCK_RETRY) && (0 < (LIBXSTREAM_LOCK_RETRY))
# if defined(LIBXSTREAM_STREAM_UNLOCK_OWNER)
  if (thread() == this_thread_id())
# endif
  {
    ++m_begin;
  }
#endif
}


void libxstream_stream::end()
{
#if defined(LIBXSTREAM_LOCK_RETRY) && (0 < (LIBXSTREAM_LOCK_RETRY))
# if defined(LIBXSTREAM_STREAM_UNLOCK_OWNER)
  if (thread() == this_thread_id())
# endif
  {
    ++m_end;
  }
#endif
}


void libxstream_stream::lock(bool retry)
{
  const int this_thread = this_thread_id();
#if defined(LIBXSTREAM_STDFEATURES)
  std::atomic<int> *const stream_thread = static_cast<std::atomic<int>*>(m_thread);
#else
  volatile int *const stream_thread = static_cast<volatile int*>(m_thread);
#endif

  int lock_thread = *stream_thread;
  if (this_thread != lock_thread) {
    int unlocked = -1;
#if defined(LIBXSTREAM_LOCK_RETRY) && (0 < (LIBXSTREAM_LOCK_RETRY))
    size_t thread_begin = m_begin, thread_end = m_end, nretry = 0;
# if defined(LIBXSTREAM_PRINT)
    size_t delay = 0;
# endif
    if (retry) {
      while (!libxstream_stream_internal::atomic_compare_exchange(*stream_thread, unlocked, this_thread)) {
        static /*const*/ size_t sleep_ms = (LIBXSTREAM_WAIT_LOCK_MS) / (LIBXSTREAM_LOCK_RETRY);

        if ((LIBXSTREAM_LOCK_RETRY) > nretry || m_begin != m_end) {
          nretry += (thread_begin == m_begin && thread_end == m_end) ? 1 : 0;
          if (0 < sleep_ms) {
            this_thread_sleep(sleep_ms);
          }
          else {
            this_thread_yield();
          }
          thread_begin = m_begin;
          thread_end = m_end;
          unlocked = -1;
        }
        else {
          unlocked = lock_thread != *stream_thread ? -1 : lock_thread;
          lock_thread = *stream_thread;
# if defined(LIBXSTREAM_PRINT)
          delay += nretry * sleep_ms;
# endif
          nretry = 0;
        }
      }

      if (-1 != unlocked) {
# if defined(LIBXSTREAM_PRINT)
        LIBXSTREAM_PRINT_WARN("libxstream_stream_unlock: stream=0x%llx released by thread=%i with delay=%lu ms",
          reinterpret_cast<unsigned long long>(this), this_thread, static_cast<unsigned long>(delay));
# else
        LIBXSTREAM_PRINT_WARN("libxstream_stream_unlock: stream=0x%llx released by thread=%i",
          reinterpret_cast<unsigned long long>(this), this_thread);
# endif
      }

      m_begin = 0;
      m_end = 0;
    }
    else
#endif
    {
      while (!libxstream_stream_internal::atomic_compare_exchange(*stream_thread, unlocked, this_thread)) {
        this_thread_yield();
        unlocked = -1;
      }
    }

    LIBXSTREAM_ASSERT(this_thread == *stream_thread);
    LIBXSTREAM_PRINT_INFO("libxstream_stream_lock: stream=0x%llx acquired by thread=%i",
      reinterpret_cast<unsigned long long>(this), this_thread);
  }
}


void libxstream_stream::unlock()
{
#if defined(LIBXSTREAM_STDFEATURES)
  std::atomic<int> *const stream_thread = static_cast<std::atomic<int>*>(m_thread);
#else
  volatile int *const stream_thread = static_cast<volatile int*>(m_thread);
#endif

#if defined(LIBXSTREAM_STREAM_UNLOCK_OWNER)
  int locked = this_thread_id();
  if (libxstream_stream_internal::atomic_compare_exchange(*stream_thread, locked, -1)) {
#else
  if (libxstream_stream_internal::atomic_store(*stream_thread, -1)) {
#endif
    LIBXSTREAM_PRINT_INFO("libxstream_stream_unlock: stream=0x%llx released by thread=%i",
      reinterpret_cast<unsigned long long>(this), this_thread_id());
  }
}


#if defined(LIBXSTREAM_OFFLOAD) && (0 != LIBXSTREAM_OFFLOAD) && defined(LIBXSTREAM_ASYNC) && (2 == (2*LIBXSTREAM_ASYNC+1)/2)
_Offload_stream libxstream_stream::handle() const
{
  const size_t nstreams = libxstream_stream_internal::registry.nstreams(m_device);

  if (nstreams != m_npartitions) {
    if (0 != m_handle) {
      const_cast<libxstream_stream*>(this)->wait(0); // complete pending operations on old stream
      _Offload_stream_destroy(m_device, m_handle);
    }

    // TODO: implement device discovery (number of threads)
    const size_t nthreads = 224;
    // TODO: implement stream priorities (weighting)
    m_handle = _Offload_stream_create(m_device, static_cast<int>(nthreads / nstreams));
    m_npartitions = nstreams;
  }
  return m_handle;
}
#endif


#if defined(LIBXSTREAM_PRINT)
const char* libxstream_stream::name() const
{
  return m_name;
}
#endif


const libxstream_stream* cast_to_stream(const void* stream)
{
  return static_cast<const libxstream_stream*>(stream);
}


libxstream_stream* cast_to_stream(void* stream)
{
  return static_cast<libxstream_stream*>(stream);
}


const libxstream_stream* cast_to_stream(const libxstream_stream* stream)
{
  return stream;
}


libxstream_stream* cast_to_stream(libxstream_stream* stream)
{
  return stream;
}


const libxstream_stream* cast_to_stream(const libxstream_stream& stream)
{
  return &stream;
}


libxstream_stream* cast_to_stream(libxstream_stream& stream)
{
  return &stream;
}

#endif // defined(LIBXSTREAM_EXPORTED) || defined(__LIBXSTREAM)
