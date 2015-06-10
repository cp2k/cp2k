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
#include "libxstream.hpp"
#include "libxstream_alloc.hpp"
#include "libxstream_capture.hpp"
#include "libxstream_context.hpp"
#include "libxstream_event.hpp"
#include "libxstream_offload.hpp"

#include <libxstream_begin.h>
#include <algorithm>
#include <cstring>
#include <limits>
#include <cstdio>
#if defined(LIBXSTREAM_STDFEATURES)
# include <thread>
# include <atomic>
# define LIBXSTREAM_STDMUTEX
# if defined(LIBXSTREAM_STDMUTEX)
#   include <mutex>
# endif
#endif
#include <libxstream_end.h>

#if defined(_OPENMP)
# include <omp.h>
#endif

#if defined(__GNUC__)
# include <pthread.h>
#endif

#if defined(__MKL)
# include <mkl.h>
#endif

#if defined(LIBXSTREAM_OFFLOAD) && (0 != LIBXSTREAM_OFFLOAD)
# include <offload.h>
#endif

#if defined(_WIN32)
# include <windows.h>
#else
# include <unistd.h>
#endif


namespace libxstream_internal {

class context_type {
public:
#if defined(LIBXSTREAM_STDFEATURES)
  typedef std::atomic<size_t> counter_type;
#else
  typedef size_t counter_type;
#endif

public:
  context_type()
    : m_lock(libxstream_lock_create())
    , m_nthreads_active(0)
    , m_device(-2)
  {}

  ~context_type() {
    libxstream_lock_destroy(m_lock);
  }

public:
  libxstream_lock* lock() { return m_lock; }
  int global_device() const { return m_device; }
  
  counter_type& nthreads_active() {
    return m_nthreads_active;
  }

  void global_device(int device) {
    libxstream_lock_acquire(m_lock);
    if (-1 > m_device) {
      m_device = device;
    }
    libxstream_lock_release(m_lock);
  }

  // store the active device per host-thread
  int& device() {
    static LIBXSTREAM_TLS int instance = -2;
    return instance;
  }

private:
  libxstream_lock* m_lock;
  counter_type m_nthreads_active;
  int m_device;
} context;


LIBXSTREAM_TARGET(mic) void mem_info(uint64_t& memory_physical, uint64_t& memory_not_allocated)
{
#if defined(_WIN32)
  MEMORYSTATUSEX status;
  status.dwLength = sizeof(status);
  const BOOL ok = GlobalMemoryStatusEx(&status);
  if (TRUE == ok) {
    memory_not_allocated = status.ullAvailPhys;
    memory_physical = status.ullTotalPhys;
  }
  else {
    memory_not_allocated = 0;
    memory_physical = 0;
  }
#else
  const long memory_pages_size = sysconf(_SC_PAGE_SIZE);
  const long memory_pages_phys = sysconf(_SC_PHYS_PAGES);
  const long memory_pages_avail = sysconf(_SC_AVPHYS_PAGES);
  memory_not_allocated = memory_pages_size * memory_pages_avail;
  memory_physical = memory_pages_size * memory_pages_phys;
#endif
}

} // namespace libxstream_internal


LIBXSTREAM_EXPORT_C LIBXSTREAM_TARGET(mic) void libxstream_use_sink(const void*) {}
LIBXSTREAM_EXPORT_C LIBXSTREAM_TARGET(mic) int libxstream_not_constant(int value) { return value; }


libxstream_lock* libxstream_lock_create()
{
#if defined(LIBXSTREAM_STDFEATURES)
# if defined(LIBXSTREAM_STDMUTEX)
#   if defined(LIBXSTREAM_LOCK_NONRECURSIVE)
  std::mutex *const typed_lock = new std::mutex;
#   else
  std::recursive_mutex *const typed_lock = new std::recursive_mutex;
#   endif
# else
  std::atomic<int> *const typed_lock = new std::atomic<int>(0);
# endif
#elif defined(_OPENMP)
# if defined(LIBXSTREAM_LOCK_NONRECURSIVE)
  omp_lock_t *const typed_lock = new omp_lock_t;
  omp_init_lock(typed_lock);
# else
  omp_nest_lock_t *const typed_lock = new omp_nest_lock_t;
  omp_init_nest_lock(typed_lock);
# endif
#else //defined(__GNUC__)
  pthread_mutexattr_t attributes;
  pthread_mutexattr_init(&attributes);
# if defined(LIBXSTREAM_LOCK_NONRECURSIVE)
  pthread_mutexattr_settype(&attributes, PTHREAD_MUTEX_NORMAL);
# else
  pthread_mutexattr_settype(&attributes, PTHREAD_MUTEX_RECURSIVE);
# endif
  pthread_mutex_t *const typed_lock = new pthread_mutex_t;
  pthread_mutex_init(typed_lock, &attributes);
#endif
  return typed_lock;
}


void libxstream_lock_destroy(libxstream_lock* lock)
{
#if defined(LIBXSTREAM_STDFEATURES)
# if defined(LIBXSTREAM_STDMUTEX)
#   if defined(LIBXSTREAM_LOCK_NONRECURSIVE)
  std::mutex *const typed_lock = static_cast<std::mutex*>(lock);
#   else
  std::recursive_mutex *const typed_lock = static_cast<std::recursive_mutex*>(lock);
#   endif
# else
  std::atomic<int> *const typed_lock = static_cast<std::atomic<int>*>(lock);
# endif
#elif defined(_OPENMP)
# if defined(LIBXSTREAM_LOCK_NONRECURSIVE)
  omp_lock_t *const typed_lock = static_cast<omp_lock_t*>(lock);
  omp_destroy_lock(typed_lock);
# else
  omp_nest_lock_t *const typed_lock = static_cast<omp_nest_lock_t*>(lock);
  omp_destroy_nest_lock(typed_lock);
# endif
#else //defined(__GNUC__)
  pthread_mutex_t *const typed_lock = static_cast<pthread_mutex_t*>(lock);
  pthread_mutex_destroy(typed_lock);
#endif
  delete typed_lock;
}


void libxstream_lock_acquire(libxstream_lock* lock)
{
  LIBXSTREAM_ASSERT(lock);
#if defined(LIBXSTREAM_STDFEATURES)
# if defined(LIBXSTREAM_STDMUTEX)
#   if defined(LIBXSTREAM_LOCK_NONRECURSIVE)
  std::mutex *const typed_lock = static_cast<std::mutex*>(lock);
#   else
  std::recursive_mutex *const typed_lock = static_cast<std::recursive_mutex*>(lock);
#   endif
  typed_lock->lock();
# else
#   if !defined(LIBXSTREAM_LOCK_NONRECURSIVE)
  LIBXSTREAM_ASSERT(false/*TODO: not implemented!*/);
#   endif
  std::atomic<int>& typed_lock = *static_cast<std::atomic<int>*>(lock);
  if (1 < ++typed_lock) {
    while (1 < typed_lock) {
      this_thread_sleep();
    }
  }
# endif
#elif defined(_OPENMP)
# if defined(LIBXSTREAM_LOCK_NONRECURSIVE)
  omp_lock_t *const typed_lock = static_cast<omp_lock_t*>(lock);
  omp_set_lock(typed_lock);
# else
  omp_nest_lock_t *const typed_lock = static_cast<omp_nest_lock_t*>(lock);
  omp_set_nest_lock(typed_lock);
# endif
#else //defined(__GNUC__)
  pthread_mutex_t *const typed_lock = static_cast<pthread_mutex_t*>(lock);
  pthread_mutex_lock(typed_lock);
#endif
}


void libxstream_lock_release(libxstream_lock* lock)
{
  LIBXSTREAM_ASSERT(lock);
#if defined(LIBXSTREAM_STDFEATURES)
# if defined(LIBXSTREAM_STDMUTEX)
#   if defined(LIBXSTREAM_LOCK_NONRECURSIVE)
  std::mutex *const typed_lock = static_cast<std::mutex*>(lock);
#   else
  std::recursive_mutex *const typed_lock = static_cast<std::recursive_mutex*>(lock);
#   endif
  typed_lock->unlock();
# else
#   if !defined(LIBXSTREAM_LOCK_NONRECURSIVE)
  LIBXSTREAM_ASSERT(false/*TODO: not implemented!*/);
#   endif
  std::atomic<int>& typed_lock = *static_cast<std::atomic<int>*>(lock);
  --typed_lock;
# endif
#elif defined(_OPENMP)
# if defined(LIBXSTREAM_LOCK_NONRECURSIVE)
  omp_lock_t *const typed_lock = static_cast<omp_lock_t*>(lock);
  omp_unset_lock(typed_lock);
# else
  omp_nest_lock_t *const typed_lock = static_cast<omp_nest_lock_t*>(lock);
  omp_unset_nest_lock(typed_lock);
# endif
#else //defined(__GNUC__)
  pthread_mutex_t *const typed_lock = static_cast<pthread_mutex_t*>(lock);
  pthread_mutex_unlock(typed_lock);
#endif
}


bool libxstream_lock_try(libxstream_lock* lock)
{
  LIBXSTREAM_ASSERT(lock);
#if defined(LIBXSTREAM_STDFEATURES)
# if defined(LIBXSTREAM_STDMUTEX)
#   if defined(LIBXSTREAM_LOCK_NONRECURSIVE)
  std::mutex *const typed_lock = static_cast<std::mutex*>(lock);
#   else
  std::recursive_mutex *const typed_lock = static_cast<std::recursive_mutex*>(lock);
#   endif
  const bool result = typed_lock->try_lock();
# else
#   if !defined(LIBXSTREAM_LOCK_NONRECURSIVE)
  LIBXSTREAM_ASSERT(false/*TODO: not implemented!*/);
#   endif
  std::atomic<int>& typed_lock = *static_cast<std::atomic<int>*>(lock);
  const bool result = 1 == ++typed_lock;
  if (!result) --typed_lock;
# endif
#elif defined(_OPENMP)
# if defined(LIBXSTREAM_LOCK_NONRECURSIVE)
  omp_lock_t *const typed_lock = static_cast<omp_lock_t*>(lock);
  const bool result = 0 != omp_test_lock(typed_lock);
# else
  omp_nest_lock_t *const typed_lock = static_cast<omp_nest_lock_t*>(lock);
  const bool result = 0 != omp_test_nest_lock(typed_lock);
# endif
#else //defined(__GNUC__)
  pthread_mutex_t *const typed_lock = static_cast<pthread_mutex_t*>(lock);
  const bool result =  0 == pthread_mutex_trylock(typed_lock);
#endif
  return result;
}


size_t nthreads_active()
{
  const size_t result = libxstream_internal::context.nthreads_active();
  LIBXSTREAM_ASSERT(result <= LIBXSTREAM_MAX_NTHREADS);
  return result;
}


int this_thread_id()
{
  static LIBXSTREAM_TLS int id = -1;
  if (0 > id) {
    libxstream_internal::context_type::counter_type& nthreads_active = libxstream_internal::context.nthreads_active();
#if defined(LIBXSTREAM_STDFEATURES)
    id = static_cast<int>(nthreads_active++);
#elif defined(_OPENMP)
    size_t nthreads = 0;
# if (201107 <= _OPENMP)
#   pragma omp atomic capture
# else
#   pragma omp critical
# endif
    nthreads = ++nthreads_active;
    id = static_cast<int>(nthreads - 1);
#else // generic
    libxstream_lock *const lock = libxstream_internal::context.lock();
    libxstream_lock_acquire(lock);
    id = static_cast<int>(nthreads_active++);
    libxstream_lock_release(lock);
#endif
  }
  return id;
}


void this_thread_yield()
{
#if defined(LIBXSTREAM_STDFEATURES) && defined(LIBXSTREAM_STDFEATURES_THREADX)
  std::this_thread::yield();
#elif defined(__GNUC__)
  pthread_yield();
#endif
}


void this_thread_sleep(size_t ms)
{
#if defined(LIBXSTREAM_STDFEATURES) && defined(LIBXSTREAM_STDFEATURES_THREADX)
  typedef std::chrono::milliseconds milliseconds;
  LIBXSTREAM_ASSERT(ms <= static_cast<size_t>(std::numeric_limits<milliseconds::rep>::max() / 1000));
  const milliseconds interval(static_cast<milliseconds::rep>(ms));
  std::this_thread::sleep_for(interval);
#elif defined(_WIN32)
  if (1 < ms) {
    LIBXSTREAM_ASSERT(ms <= std::numeric_limits<DWORD>::max());
    Sleep(static_cast<DWORD>(ms));
  }
  else {
    SwitchToThread();
  }
#else
  const size_t s = ms / 1000;
  ms -= 1000 * s;
  LIBXSTREAM_ASSERT(ms <= static_cast<size_t>(std::numeric_limits<long>::max() / (1000 * 1000)));
  const timespec pause = {
    static_cast<time_t>(s),
    static_cast<long>(ms * 1000 * 1000)
  };
  nanosleep(&pause, 0);
#endif
}


LIBXSTREAM_EXPORT_C int libxstream_get_ndevices(size_t* ndevices)
{
  LIBXSTREAM_CHECK_CONDITION(ndevices);

#if defined(LIBXSTREAM_OFFLOAD) && (0 != LIBXSTREAM_OFFLOAD) && !defined(__MIC__)
  static const int idevices = std::min(_Offload_number_of_devices(), LIBXSTREAM_MAX_NDEVICES);
  LIBXSTREAM_CHECK_CONDITION(0 <= idevices);
  *ndevices = static_cast<size_t>(idevices);
#else
  *ndevices = 1; // host
#endif

#if defined(LIBXSTREAM_PRINT)
  static LIBXSTREAM_TLS bool print = true;
  if (print) {
    LIBXSTREAM_PRINT_INFOCTX("ndevices=%lu", static_cast<unsigned long>(*ndevices));
    print = false;
  }
#endif

  return LIBXSTREAM_ERROR_NONE;
}


LIBXSTREAM_EXPORT_C int libxstream_get_active_device(int* device)
{
  LIBXSTREAM_CHECK_CONDITION(0 != device);
  int result = LIBXSTREAM_ERROR_NONE, active_device = libxstream_internal::context.device();

  if (-1 > active_device) {
    active_device = libxstream_internal::context.global_device();

    if (-1 > active_device) {
      size_t ndevices = 0;
      result = libxstream_get_ndevices(&ndevices);
      active_device = static_cast<int>(ndevices - 1);
      libxstream_internal::context.global_device(active_device);
      libxstream_internal::context.device() = active_device;
    }

    LIBXSTREAM_PRINT_INFOCTX("device=%i (fallback) thread=%i", active_device, this_thread_id());
  }

  *device = active_device;
  return result;
}


LIBXSTREAM_EXPORT_C int libxstream_set_active_device(int device)
{
  size_t ndevices = LIBXSTREAM_MAX_NDEVICES;
  LIBXSTREAM_CHECK_CONDITION(-1 <= device && ndevices >= static_cast<size_t>(device + 1) && LIBXSTREAM_ERROR_NONE == libxstream_get_ndevices(&ndevices) && ndevices >= static_cast<size_t>(device + 1));

  if (-1 > libxstream_internal::context.global_device()) {
    libxstream_internal::context.global_device(device);
  }

  libxstream_internal::context.device() = device;
  LIBXSTREAM_PRINT_INFOCTX("device=%i thread=%i", device, this_thread_id());

  return LIBXSTREAM_ERROR_NONE;
}


LIBXSTREAM_EXPORT_C int libxstream_mem_info(int device, size_t* allocatable, size_t* physical)
{
  LIBXSTREAM_CHECK_CONDITION(allocatable || physical);
  uint64_t memory_physical = 0, memory_allocatable = 0;
  int result = LIBXSTREAM_ERROR_NONE;

#if defined(LIBXSTREAM_OFFLOAD)
  if (0 <= device) {
    LIBXSTREAM_ASYNC_BEGIN(0, device, &memory_physical, &memory_allocatable)
    {
      uint64_t& memory_physical = *ptr<uint64_t,1>();
      uint64_t& memory_allocatable = *ptr<uint64_t,2>();
#     pragma offload target(mic:LIBXSTREAM_ASYNC_DEVICE) //out(memory_physical, memory_allocatable)
      libxstream_internal::mem_info(memory_physical, memory_allocatable);
    }
    LIBXSTREAM_ASYNC_END(LIBXSTREAM_CALL_DEFAULT | LIBXSTREAM_CALL_WAIT, result);
  }
  else {
#else
  {
    libxstream_use_sink(&device);
#endif
    libxstream_internal::mem_info(memory_physical, memory_allocatable);
  }

  LIBXSTREAM_PRINT_INFOCTX("device=%i allocatable=%lu physical=%lu", device,
    static_cast<unsigned long>(memory_allocatable),
    static_cast<unsigned long>(memory_physical));
  LIBXSTREAM_CHECK_CONDITION(0 < memory_physical && 0 < memory_allocatable);

  if (allocatable) {
    *allocatable = static_cast<size_t>(memory_allocatable);
  }

  if (physical) {
    *physical = static_cast<size_t>(memory_physical);
  }

  return result;
}


LIBXSTREAM_EXPORT_C int libxstream_mem_pointer(int device, const void* memory, const void** real)
{
  LIBXSTREAM_CHECK_CONDITION(0 != real);
  int result = LIBXSTREAM_ERROR_NONE;

  if (memory) {
#if defined(LIBXSTREAM_OFFLOAD)
    if (0 <= device) {
      LIBXSTREAM_ASYNC_BEGIN(0, device, memory, real)
      {
        const char* memory = ptr<const char,1>();
        const void*& real = *ptr<const void*,2>();
#       pragma offload target(mic:LIBXSTREAM_ASYNC_DEVICE) in(memory: LIBXSTREAM_OFFLOAD_REFRESH) //out(real)
        real = memory;
      }
      LIBXSTREAM_ASYNC_END(LIBXSTREAM_CALL_DEFAULT | LIBXSTREAM_CALL_WAIT, result);
    }
    else {
#else
    {
      libxstream_use_sink(&device);
#endif
      *real = memory;
    }
  }
  else {
    *real = 0;
  }

  return result;
}


LIBXSTREAM_EXPORT_C int libxstream_mem_allocate(int device, void** memory, size_t size, size_t alignment)
{
  LIBXSTREAM_CHECK_CONDITION(0 != memory);
  int result = LIBXSTREAM_ERROR_NONE;

#if defined(LIBXSTREAM_OFFLOAD)
  if (0 <= device) {
    void* buffer = 0;
    result = libxstream_virt_allocate(&buffer, size, alignment, &device, sizeof(device));

    if (LIBXSTREAM_ERROR_NONE == result && 0 != buffer) {
      LIBXSTREAM_ASYNC_BEGIN(0, device, buffer, size)
      {
        const char* buffer = ptr<const char,1>();
        const size_t size = val<const size_t,2>();
#       pragma offload_transfer target(mic:LIBXSTREAM_ASYNC_DEVICE) nocopy(buffer: length(size) LIBXSTREAM_OFFLOAD_ALLOC)
      }
      LIBXSTREAM_ASYNC_END(LIBXSTREAM_CALL_DEFAULT | LIBXSTREAM_CALL_WAIT, result);
      LIBXSTREAM_CHECK_ERROR(result);
      *memory = buffer;
    }
  }
  else {
#else
  {
    libxstream_use_sink(&device);
#endif
    void* buffer = 0;
    result = libxstream_real_allocate(&buffer, size, alignment);

    if (LIBXSTREAM_ERROR_NONE == result && 0 != buffer) {
#if defined(LIBXSTREAM_OFFLOAD) && defined(LIBXSTREAM_ALLOC_PINNED)
      LIBXSTREAM_ASYNC_BEGIN(0, device, buffer, size)
      {
        const char* buffer = ptr<const char,1>();
        const size_t size = val<const size_t,2>();
#       pragma offload_transfer target(mic) host_pin(buffer: length(size))
      }
      LIBXSTREAM_ASYNC_END(LIBXSTREAM_CALL_DEFAULT | LIBXSTREAM_CALL_WAIT, result);
      LIBXSTREAM_CHECK_ERROR(result);
#endif
      *memory = buffer;
    }
  }

#if defined(LIBXSTREAM_SYNCMEM)
  libxstream_stream::sync(device);
#endif
#if defined(LIBXSTREAM_DEBUG)
  const void* real = 0;
  if (0 <= device && LIBXSTREAM_ERROR_NONE == libxstream_mem_pointer(device, *memory, &real)) {
    LIBXSTREAM_PRINT_INFOCTX("device=%i buffer=0x%llx real=0x%llx size=%lu", device,
      memory ? reinterpret_cast<unsigned long long>(*memory) : 0, reinterpret_cast<unsigned long long>(real),
      static_cast<unsigned long>(size));
  }
  else
#endif
  {
    LIBXSTREAM_PRINT_INFOCTX("device=%i buffer=0x%llx size=%lu", device,
      reinterpret_cast<unsigned long long>(*memory), static_cast<unsigned long>(size));
  }

  return result;
}


LIBXSTREAM_EXPORT_C int libxstream_mem_deallocate(int device, const void* memory)
{
  int result = LIBXSTREAM_ERROR_NONE;
#if defined(LIBXSTREAM_DEBUG)
  const void* real = 0;
  if (0 <= device && LIBXSTREAM_ERROR_NONE == libxstream_mem_pointer(device, memory, &real)) {
    LIBXSTREAM_PRINT_INFOCTX("device=%i buffer=0x%llx real=0x%llx", device,
      reinterpret_cast<unsigned long long>(memory), reinterpret_cast<unsigned long long>(real));
  }
  else
#endif
  {
    LIBXSTREAM_PRINT_INFOCTX("device=%i buffer=0x%llx", device, reinterpret_cast<unsigned long long>(memory));
  }

  if (memory) {
#if defined(LIBXSTREAM_SYNCMEM)
    libxstream_stream::sync(device);
#endif
#if defined(LIBXSTREAM_OFFLOAD)
    if (0 <= device) {
# if defined(LIBXSTREAM_CHECK)
      const int memory_device = *static_cast<const int*>(libxstream_virt_data(memory));
      if (device != memory_device) {
        LIBXSTREAM_PRINT_WARN("device %i does not match allocating device %i!", device, memory_device);
        LIBXSTREAM_CHECK_CONDITION(0 <= memory_device);
        device = memory_device;
      }
# endif
      LIBXSTREAM_ASYNC_BEGIN(0, device, memory)
      {
        const char *const memory = ptr<const char,1>();
#       pragma offload_transfer target(mic:LIBXSTREAM_ASYNC_DEVICE) nocopy(memory: length(0) LIBXSTREAM_OFFLOAD_FREE)
        LIBXSTREAM_CHECK_CALL_ASSERT(status(libxstream_virt_deallocate(memory)));
      }
      LIBXSTREAM_ASYNC_END(LIBXSTREAM_CALL_DEFAULT, result);
    }
    else {
#else
    {
      libxstream_use_sink(&device);
#endif
#if defined(LIBXSTREAM_OFFLOAD) && defined(LIBXSTREAM_ALLOC_PINNED)
      LIBXSTREAM_ASYNC_BEGIN(0, device, memory)
      {
        const char* memory = ptr<const char,1>();
#       pragma offload_transfer target(mic) host_unpin(memory: length(0))
        LIBXSTREAM_CHECK_CALL_ASSERT(status(libxstream_real_deallocate(memory)));
      }
      LIBXSTREAM_ASYNC_END(LIBXSTREAM_CALL_DEFAULT, result);
#else
      result = libxstream_real_deallocate(memory);
#endif
    }
  }

  return result;
}


LIBXSTREAM_EXPORT_C int libxstream_memset_zero(void* memory, size_t size, libxstream_stream* stream)
{
  LIBXSTREAM_PRINT_INFO("libxstream_memset_zero: buffer=0x%llx size=%lu stream=0x%llx",
    reinterpret_cast<unsigned long long>(memory), static_cast<unsigned long>(size),
    reinterpret_cast<unsigned long long>(stream));
  LIBXSTREAM_CHECK_CONDITION(memory && stream);
  int result = LIBXSTREAM_ERROR_NONE;

  LIBXSTREAM_ASYNC_BEGIN(stream, memory, size)
  {
    char* dst = ptr<char,0>();
    const size_t size = val<const size_t,1>();

#if defined(LIBXSTREAM_OFFLOAD)
    if (0 <= LIBXSTREAM_ASYNC_DEVICE) {
      if (LIBXSTREAM_ASYNC_READY) {
#       pragma offload LIBXSTREAM_ASYNC_TARGET_SIGNAL in(size) out(dst: LIBXSTREAM_OFFLOAD_REFRESH)
        memset(dst, 0, size);
      }
      else {
#       pragma offload LIBXSTREAM_ASYNC_TARGET_WAIT in(size) out(dst: LIBXSTREAM_OFFLOAD_REFRESH)
        memset(dst, 0, size);
      }
    }
    else
#endif
    {
      memset(dst, 0, size);
    }
  }
  LIBXSTREAM_ASYNC_END(LIBXSTREAM_CALL_DEFAULT, result);

  return result;
}


LIBXSTREAM_EXPORT_C int libxstream_memcpy_h2d(const void* host_mem, void* dev_mem, size_t size, libxstream_stream* stream)
{
  LIBXSTREAM_PRINT_INFO("libxstream_memcpy_h2d: 0x%llx->0x%llx size=%lu stream=0x%llx", reinterpret_cast<unsigned long long>(host_mem),
    reinterpret_cast<unsigned long long>(dev_mem), static_cast<unsigned long>(size), reinterpret_cast<unsigned long long>(stream));
  LIBXSTREAM_CHECK_CONDITION(host_mem && dev_mem && host_mem != dev_mem && stream);
  int result = LIBXSTREAM_ERROR_NONE;

  LIBXSTREAM_ASYNC_BEGIN(stream, host_mem, dev_mem, size)
  {
    const char *const src = ptr<const char,0>();
    char *const dst = ptr<char,1>();
    const size_t size = val<const size_t,2>();

#if defined(LIBXSTREAM_OFFLOAD)
    if (0 <= LIBXSTREAM_ASYNC_DEVICE) {
      if (LIBXSTREAM_ASYNC_READY) {
#       pragma offload_transfer LIBXSTREAM_ASYNC_TARGET_SIGNAL in(src: length(size) into(dst) LIBXSTREAM_OFFLOAD_REUSE)
      }
      else {
#       pragma offload_transfer LIBXSTREAM_ASYNC_TARGET_WAIT in(src: length(size) into(dst) LIBXSTREAM_OFFLOAD_REUSE)
      }
    }
    else
#endif
    {
#if defined(LIBXSTREAM_ASYNCHOST)
      if (LIBXSTREAM_ASYNC_READY) {
#       pragma omp task depend(out:capture_region_signal) depend(in:LIBXSTREAM_ASYNC_PENDING)
        std::copy(src, src + size, dst);
      }
      else {
#       pragma omp task depend(out:capture_region_signal)
        std::copy(src, src + size, dst);
        ++capture_region_signal_consumed;
      }
#else
      std::copy(src, src + size, dst);
#endif
    }
  }
  LIBXSTREAM_ASYNC_END(LIBXSTREAM_CALL_DEFAULT, result);

  return result;
}


LIBXSTREAM_EXPORT_C int libxstream_memcpy_d2h(const void* dev_mem, void* host_mem, size_t size, libxstream_stream* stream)
{
  LIBXSTREAM_PRINT_INFO("libxstream_memcpy_d2h: 0x%llx->0x%llx size=%lu stream=0x%llx", reinterpret_cast<unsigned long long>(dev_mem),
    reinterpret_cast<unsigned long long>(host_mem), static_cast<unsigned long>(size), reinterpret_cast<unsigned long long>(stream));
  LIBXSTREAM_CHECK_CONDITION(dev_mem && host_mem && dev_mem != host_mem && stream);
  int result = LIBXSTREAM_ERROR_NONE;

  LIBXSTREAM_ASYNC_BEGIN(stream, dev_mem, host_mem, size)
  {
    const char* src = ptr<const char,0>();
    char *const dst = ptr<char,1>();
    const size_t size = val<const size_t,2>();

#if defined(LIBXSTREAM_OFFLOAD)
    if (0 <= LIBXSTREAM_ASYNC_DEVICE) {
      if (LIBXSTREAM_ASYNC_READY) {
#       pragma offload_transfer LIBXSTREAM_ASYNC_TARGET_SIGNAL out(src: length(size) into(dst) LIBXSTREAM_OFFLOAD_REUSE)
      }
      else {
#       pragma offload_transfer LIBXSTREAM_ASYNC_TARGET_WAIT out(src: length(size) into(dst) LIBXSTREAM_OFFLOAD_REUSE)
      }
    }
    else
#endif
    {
      std::copy(src, src + size, dst);
    }
  }
  LIBXSTREAM_ASYNC_END(LIBXSTREAM_CALL_DEFAULT, result);

  return result;
}


LIBXSTREAM_EXPORT_C int libxstream_memcpy_d2d(const void* src, void* dst, size_t size, libxstream_stream* stream)
{
  LIBXSTREAM_PRINT_INFO("libxstream_memcpy_d2d: 0x%llx->0x%llx size=%lu stream=0x%llx", reinterpret_cast<unsigned long long>(src),
    reinterpret_cast<unsigned long long>(dst), static_cast<unsigned long>(size), reinterpret_cast<unsigned long long>(stream));
  LIBXSTREAM_CHECK_CONDITION(src && dst && stream);
  int result = LIBXSTREAM_ERROR_NONE;

  if (src != dst) {
    LIBXSTREAM_ASYNC_BEGIN(stream, src, dst, size)
    {
      const char *const src = ptr<const char,0>();
      char* dst = ptr<char,1>();
      const size_t size = val<const size_t,2>();

#if defined(LIBXSTREAM_OFFLOAD)
      if (0 <= LIBXSTREAM_ASYNC_DEVICE) {
        // TODO: implement cross-device transfer

        if (LIBXSTREAM_ASYNC_READY) {
#         pragma offload LIBXSTREAM_ASYNC_TARGET_SIGNAL in(size) in(src: LIBXSTREAM_OFFLOAD_REFRESH) out(dst: LIBXSTREAM_OFFLOAD_REFRESH)
          memcpy(dst, src, size);
        }
        else {
#         pragma offload LIBXSTREAM_ASYNC_TARGET_WAIT in(size) in(src: LIBXSTREAM_OFFLOAD_REFRESH) out(dst: LIBXSTREAM_OFFLOAD_REFRESH)
          memcpy(dst, src, size);
        }
      }
      else
#endif
      {
        memcpy(dst, src, size);
      }
    }
    LIBXSTREAM_ASYNC_END(LIBXSTREAM_CALL_DEFAULT, result);
  }

  return result;
}


LIBXSTREAM_EXPORT_C int libxstream_stream_priority_range(int* least, int* greatest)
{
  *least = -1;
  *greatest = -1;
  return LIBXSTREAM_ERROR_NONE;
}


LIBXSTREAM_EXPORT_C int libxstream_stream_create(libxstream_stream** stream, int device, int demux, int priority, const char* name)
{
  LIBXSTREAM_CHECK_CONDITION(stream);
  libxstream_stream *const s = new libxstream_stream(device, demux, priority, name);
  LIBXSTREAM_ASSERT(s);
  *stream = s;

#if defined(LIBXSTREAM_PRINT)
  if (name && *name) {
    LIBXSTREAM_PRINT_INFOCTX("stream=0x%llx device=%i demux=%i priority=%i (%s)",
      reinterpret_cast<unsigned long long>(*stream), device, demux, priority, name);
  }
  else {
    LIBXSTREAM_PRINT_INFOCTX("stream=0x%llx device=%i demux=%i priority=%i",
      reinterpret_cast<unsigned long long>(*stream), device, demux, priority);
  }
#endif

  return LIBXSTREAM_ERROR_NONE;
}


LIBXSTREAM_EXPORT_C int libxstream_stream_destroy(const libxstream_stream* stream)
{
#if defined(LIBXSTREAM_PRINT)
  if (stream) {
    const char *const name = stream->name();
    if (name && *name) {
      LIBXSTREAM_PRINT_INFOCTX("stream=0x%llx name=\"%s\"", reinterpret_cast<unsigned long long>(stream), name);
    }
    else {
      LIBXSTREAM_PRINT_INFOCTX("stream=0x%llx", reinterpret_cast<unsigned long long>(stream));
    }
  }
#endif
  delete stream;
  return LIBXSTREAM_ERROR_NONE;
}


LIBXSTREAM_EXPORT_C int libxstream_stream_sync(libxstream_stream* stream)
{
#if defined(LIBXSTREAM_PRINT)
  if (0 != stream) {
    const char *const name = stream->name();
    if (name && *name) {
      LIBXSTREAM_PRINT_INFOCTX("stream=0x%llx name=\"%s\"", reinterpret_cast<unsigned long long>(stream), name);
    }
    else {
      LIBXSTREAM_PRINT_INFOCTX("stream=0x%llx", reinterpret_cast<unsigned long long>(stream));
    }
  }
  else {
    LIBXSTREAM_PRINT_INFOCTX0("synchronize all streams");
  }
#endif

  return stream ? stream->wait(0) : libxstream_stream::sync();
}


LIBXSTREAM_EXPORT_C int libxstream_stream_wait_event(const libxstream_stream* stream, const libxstream_event* event)
{
  LIBXSTREAM_PRINT_INFOCTX("event=0x%llx stream=0x%llx", reinterpret_cast<unsigned long long>(event), reinterpret_cast<unsigned long long>(stream));
  LIBXSTREAM_CHECK_CONDITION(0 != event);
  return libxstream_event(*event).wait(stream);
}


LIBXSTREAM_EXPORT_C int libxstream_stream_lock(libxstream_stream* stream)
{
  LIBXSTREAM_CHECK_CONDITION(stream && 0 == stream->demux());
  // manual locking is supposed to be correct and hence there is no need to retry
  stream->lock(false);
  return LIBXSTREAM_ERROR_NONE;
}


LIBXSTREAM_EXPORT_C int libxstream_stream_unlock(libxstream_stream* stream)
{
  LIBXSTREAM_CHECK_CONDITION(stream && 0 == stream->demux());
  stream->unlock();
  return LIBXSTREAM_ERROR_NONE;
}


LIBXSTREAM_EXPORT_C int libxstream_stream_device(const libxstream_stream* stream, int* device)
{
  LIBXSTREAM_CHECK_CONDITION(stream && device);
  *device = stream->device();
  LIBXSTREAM_PRINT_INFOCTX("stream=0x%llx device=%i", reinterpret_cast<unsigned long long>(stream), *device);
  return LIBXSTREAM_ERROR_NONE;
}


LIBXSTREAM_EXPORT_C int libxstream_stream_demux(const libxstream_stream* stream, int* demux)
{
  LIBXSTREAM_CHECK_CONDITION(stream && demux);
  *demux = stream->demux();
  LIBXSTREAM_PRINT_INFOCTX("stream=0x%llx demux=%i", reinterpret_cast<unsigned long long>(stream), *demux);
  return LIBXSTREAM_ERROR_NONE;
}


LIBXSTREAM_EXPORT_C int libxstream_event_create(libxstream_event** event)
{
  LIBXSTREAM_CHECK_CONDITION(event);
  *event = new libxstream_event;
  LIBXSTREAM_PRINT_INFOCTX("event=0x%llx", reinterpret_cast<unsigned long long>(*event));
  return LIBXSTREAM_ERROR_NONE;
}


LIBXSTREAM_EXPORT_C int libxstream_event_destroy(const libxstream_event* event)
{
  LIBXSTREAM_PRINT_INFOCTX("event=0x%llx", reinterpret_cast<unsigned long long>(event));
  delete event;
  return LIBXSTREAM_ERROR_NONE;
}


LIBXSTREAM_EXPORT_C int libxstream_event_record(libxstream_event* event, libxstream_stream* stream)
{
  LIBXSTREAM_PRINT_INFOCTX("event=0x%llx stream=0x%llx", reinterpret_cast<unsigned long long>(event), reinterpret_cast<unsigned long long>(stream));
  return stream ? event->enqueue(*stream, true) : libxstream_stream::enqueue(*event);
}


LIBXSTREAM_EXPORT_C int libxstream_event_query(const libxstream_event* event, libxstream_bool* occured)
{
  LIBXSTREAM_PRINT_INFOCTX("event=0x%llx", reinterpret_cast<unsigned long long>(event));
  LIBXSTREAM_CHECK_CONDITION(event && occured);

  bool has_occurred = true;
  const int result = event->query(has_occurred);
  *occured = has_occurred ? LIBXSTREAM_TRUE : LIBXSTREAM_FALSE;

  return result;
}


LIBXSTREAM_EXPORT_C int libxstream_event_synchronize(libxstream_event* event)
{
  LIBXSTREAM_PRINT_INFOCTX("event=0x%llx", reinterpret_cast<unsigned long long>(event));
  return event ? event->wait() : libxstream_stream::sync();
}


LIBXSTREAM_EXPORT_C int libxstream_fn_create_signature(libxstream_argument** signature, size_t nargs)
{
  if (0 < nargs) {
    libxstream_argument *const arguments = new libxstream_argument[nargs+1];
    LIBXSTREAM_CHECK_CALL(libxstream_construct(arguments, nargs));
    *signature = arguments;
  }
  else {
    LIBXSTREAM_CHECK_CONDITION(0 != signature && (LIBXSTREAM_MAX_NARGS) >= nargs);
    *signature = 0;
  }
  LIBXSTREAM_PRINT_INFOCTX("signature=0x%llx nargs=%lu", reinterpret_cast<unsigned long long>(*signature), static_cast<unsigned long>(nargs));
  return LIBXSTREAM_ERROR_NONE;
}


LIBXSTREAM_EXPORT_C int libxstream_fn_destroy_signature(const libxstream_argument* signature)
{
#if defined(LIBXSTREAM_PRINT)
  size_t nargs = 0, arity = 0;
  LIBXSTREAM_CHECK_CALL(libxstream_fn_nargs(signature, &nargs));
  LIBXSTREAM_CHECK_CALL(libxstream_get_arity(signature, &arity));
  LIBXSTREAM_PRINT_INFOCTX("signature=0x%llx nargs=%lu arity=%lu",
    reinterpret_cast<unsigned long long>(signature),
    static_cast<unsigned long>(nargs),
    static_cast<unsigned long>(arity));
  if (arity < nargs && (LIBXSTREAM_MAX_NARGS) != nargs) {
    LIBXSTREAM_PRINT_WARNCTX0("unused argument slots!");
  }
#endif
  delete signature;
  return LIBXSTREAM_ERROR_NONE;
}


LIBXSTREAM_EXPORT_C int libxstream_fn_clear_signature(libxstream_argument* signature)
{
  size_t nargs = 0;
  if (signature) {
    LIBXSTREAM_CHECK_CALL(libxstream_fn_nargs(signature, &nargs));
    LIBXSTREAM_CHECK_CALL(libxstream_construct(signature, nargs));
  }
  LIBXSTREAM_PRINT_INFOCTX("signature=0x%llx nargs=%lu", reinterpret_cast<unsigned long long>(signature), static_cast<unsigned long>(nargs));
  return LIBXSTREAM_ERROR_NONE;
}


LIBXSTREAM_EXPORT_C int libxstream_fn_signature(libxstream_argument** signature)
{
  static LIBXSTREAM_TLS libxstream_argument arguments[(LIBXSTREAM_MAX_NARGS)+1];
  LIBXSTREAM_CHECK_CALL(libxstream_construct(arguments, LIBXSTREAM_MAX_NARGS));
  *signature = arguments;
  return LIBXSTREAM_ERROR_NONE;
}


LIBXSTREAM_EXPORT_C int libxstream_fn_input(libxstream_argument* signature, size_t arg, const void* in, libxstream_type type, size_t dims, const size_t shape[])
{
  LIBXSTREAM_CHECK_CONDITION(0 != signature);
#if defined(LIBXSTREAM_DEBUG)
  size_t nargs = 0;
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == libxstream_fn_nargs(signature, &nargs) && arg < nargs);
#endif
  return libxstream_construct(signature, arg, libxstream_argument::kind_input, in, type, dims, shape);
}


LIBXSTREAM_EXPORT_C int libxstream_fn_output(libxstream_argument* signature, size_t arg, void* out, libxstream_type type, size_t dims, const size_t shape[])
{
  LIBXSTREAM_CHECK_CONDITION(0 != signature);
#if defined(LIBXSTREAM_DEBUG)
  size_t nargs = 0;
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == libxstream_fn_nargs(signature, &nargs) && arg < nargs);
#endif
  return libxstream_construct(signature, arg, libxstream_argument::kind_output, out, type, dims, shape);
}


LIBXSTREAM_EXPORT_C int libxstream_fn_inout(libxstream_argument* signature, size_t arg, void* inout, libxstream_type type, size_t dims, const size_t shape[])
{
  LIBXSTREAM_CHECK_CONDITION(0 != signature);
#if defined(LIBXSTREAM_DEBUG)
  size_t nargs = 0;
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == libxstream_fn_nargs(signature, &nargs) && arg < nargs);
#endif
  return libxstream_construct(signature, arg, libxstream_argument::kind_inout, inout, type, dims, shape);
}


LIBXSTREAM_EXPORT_C int libxstream_fn_nargs(const libxstream_argument* signature, size_t* nargs)
{
  LIBXSTREAM_CHECK_CONDITION(0 != nargs);
  size_t n = 0;
  if (signature) {
    while (libxstream_argument::kind_invalid != signature[n].kind) {
      LIBXSTREAM_ASSERT(n < (LIBXSTREAM_MAX_NARGS));
      ++n;
    }
  }
  *nargs = n;
  return LIBXSTREAM_ERROR_NONE;
}


LIBXSTREAM_EXPORT_C int libxstream_fn_call(libxstream_function function, const libxstream_argument* signature, libxstream_stream* stream, int flags)
{
  LIBXSTREAM_PRINT_INFOCTX("function=0x%llx signature=0x%llx stream=0x%llx flags=%i",
    reinterpret_cast<unsigned long long>(function), reinterpret_cast<unsigned long long>(signature),
    reinterpret_cast<unsigned long long>(stream), flags);
  LIBXSTREAM_CHECK_CONDITION(0 != function && 0 != stream);
  return libxstream_offload(function, signature, stream, flags);
}


LIBXSTREAM_EXPORT_C LIBXSTREAM_TARGET(mic) int libxstream_get_typesize(libxstream_type type, size_t* typesize)
{
  LIBXSTREAM_CHECK_CONDITION(0 != typesize);
  int result = LIBXSTREAM_ERROR_NONE;

  switch(type) {
    case LIBXSTREAM_TYPE_I8:    *typesize = 1;  break;
    case LIBXSTREAM_TYPE_I16:   *typesize = 2;  break;
    case LIBXSTREAM_TYPE_I32:   *typesize = 4;  break;
    case LIBXSTREAM_TYPE_I64:   *typesize = 8;  break;
    case LIBXSTREAM_TYPE_F32:   *typesize = 4;  break;
    case LIBXSTREAM_TYPE_F64:   *typesize = 8;  break;
    case LIBXSTREAM_TYPE_C32:   *typesize = 8;  break;
    case LIBXSTREAM_TYPE_C64:   *typesize = 16; break;
    case LIBXSTREAM_TYPE_U8:    *typesize = 1;  break;
    case LIBXSTREAM_TYPE_U16:   *typesize = 2;  break;
    case LIBXSTREAM_TYPE_U32:   *typesize = 4;  break;
    case LIBXSTREAM_TYPE_U64:   *typesize = 8;  break;
    case LIBXSTREAM_TYPE_CHAR:  *typesize = 1;  break;
    default: // LIBXSTREAM_TYPE_VOID, etc.
      result = LIBXSTREAM_ERROR_CONDITION;
  }
  return result;
}


LIBXSTREAM_EXPORT_C LIBXSTREAM_TARGET(mic) int libxstream_get_autotype(size_t typesize, libxstream_type* autotype)
{
  LIBXSTREAM_CHECK_CONDITION(0 != autotype);
  size_t size = 0;
  int i = 0;

  do {
    *autotype = static_cast<libxstream_type>(i);
    if (LIBXSTREAM_ERROR_NONE == libxstream_get_typesize(*autotype, &size) && typesize != size) {
      ++i;
    }
    else {
      i = LIBXSTREAM_TYPE_INVALID; // break
    }
  }
  while (i <= LIBXSTREAM_TYPE_VOID);

  return LIBXSTREAM_ERROR_NONE;
}


LIBXSTREAM_EXPORT_C LIBXSTREAM_TARGET(mic) int libxstream_get_typename(libxstream_type type, const char** name)
{
  LIBXSTREAM_CHECK_CONDITION(0 != name);
  int result = LIBXSTREAM_ERROR_NONE;

  switch(type) {
    case LIBXSTREAM_TYPE_I8:    *name = "i8";   break;
    case LIBXSTREAM_TYPE_U8:    *name = "u8";   break;
    case LIBXSTREAM_TYPE_I16:   *name = "i16";  break;
    case LIBXSTREAM_TYPE_U16:   *name = "u16";  break;
    case LIBXSTREAM_TYPE_I32:   *name = "i32";  break;
    case LIBXSTREAM_TYPE_U32:   *name = "u32";  break;
    case LIBXSTREAM_TYPE_I64:   *name = "i64";  break;
    case LIBXSTREAM_TYPE_U64:   *name = "u64";  break;
    case LIBXSTREAM_TYPE_F32:   *name = "f32";  break;
    case LIBXSTREAM_TYPE_F64:   *name = "f64";  break;
    case LIBXSTREAM_TYPE_C32:   *name = "c32";  break;
    case LIBXSTREAM_TYPE_C64:   *name = "c64";  break;
    case LIBXSTREAM_TYPE_CHAR:  *name = "char"; break;
    case LIBXSTREAM_TYPE_VOID:  *name = "void"; break;
    default:
      result = LIBXSTREAM_ERROR_CONDITION;
  }
  return result;
}


LIBXSTREAM_EXPORT_C LIBXSTREAM_TARGET(mic) int libxstream_get_argument(const void* variable, size_t* arg)
{
  const libxstream_context& context = libxstream_context::instance();
  LIBXSTREAM_CHECK_CONDITION(0 != arg && 0 == (LIBXSTREAM_CALL_EXTERNAL & context.flags));
  const libxstream_argument *const hit = libxstream_find(context, variable);
  int result = LIBXSTREAM_ERROR_NONE;

  if (0 != hit) {
    LIBXSTREAM_ASSERT((LIBXSTREAM_MAX_NARGS) > (hit - context.signature));
    *arg = hit - context.signature;
  }
  else {
    result = LIBXSTREAM_ERROR_CONDITION;
  }

  return result;
}


LIBXSTREAM_EXPORT_C LIBXSTREAM_TARGET(mic) int libxstream_get_arity(const libxstream_argument* signature, size_t* arity)
{
  LIBXSTREAM_CHECK_CONDITION(0 != arity);
  if (0 == signature) {
    const libxstream_context& context = libxstream_context::instance();
    if (0 == (LIBXSTREAM_CALL_EXTERNAL & context.flags)) {
      signature = context.signature;
    }
  }

  size_t n = 0;
  if (signature) {
    while (LIBXSTREAM_TYPE_INVALID != signature[n].type) {
      LIBXSTREAM_ASSERT(n < (LIBXSTREAM_MAX_NARGS));
      LIBXSTREAM_ASSERT(libxstream_argument::kind_invalid != signature[n].kind);
      ++n;
    }
  }
  *arity = n;

  return LIBXSTREAM_ERROR_NONE;
}


LIBXSTREAM_EXPORT_C LIBXSTREAM_TARGET(mic) int libxstream_get_data(const libxstream_argument* signature, size_t arg, const void** data)
{
  LIBXSTREAM_CHECK_CONDITION(0 != data);
#if defined(LIBXSTREAM_DEBUG) && (!defined(LIBXSTREAM_OFFLOAD) || (0 == LIBXSTREAM_OFFLOAD))
  size_t arity = 0;
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == libxstream_get_arity(signature, &arity) && arg < arity);
#endif
  if (0 == signature) {
    const libxstream_context& context = libxstream_context::instance();
    if (0 == (LIBXSTREAM_CALL_EXTERNAL & context.flags)) {
      signature = context.signature;
    }
  }
  LIBXSTREAM_ASSERT(0 != signature);
  *data = libxstream_get_value(signature[arg], false).const_pointer;
  return LIBXSTREAM_ERROR_NONE;
}


LIBXSTREAM_EXPORT_C LIBXSTREAM_TARGET(mic) int libxstream_get_string(const libxstream_argument* signature, size_t arg, const char** value)
{
  LIBXSTREAM_CHECK_CONDITION(0 != value);
#if defined(LIBXSTREAM_DEBUG) && (!defined(LIBXSTREAM_OFFLOAD) || (0 == LIBXSTREAM_OFFLOAD))
  size_t arity = 0;
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == libxstream_get_arity(signature, &arity) && arg < arity);
#endif
  if (0 == signature) {
    const libxstream_context& context = libxstream_context::instance();
    if (0 == (LIBXSTREAM_CALL_EXTERNAL & context.flags)) {
      signature = context.signature;
    }
  }
  LIBXSTREAM_ASSERT(0 != signature);
  const libxstream_argument& argument = signature[arg];
  const void *const data = libxstream_get_value(argument, false).const_pointer;
  static LIBXSTREAM_TLS char buffer[128];
  int result = LIBXSTREAM_ERROR_NONE;

  if (0 < argument.dims || 0 == data) {
    if (LIBXSTREAM_TYPE_INVALID > argument.type) {
      LIBXSTREAM_SNPRINTF(buffer, sizeof(buffer), "0x%llx", reinterpret_cast<unsigned long long>(data));
      *value = buffer;
    }
    else {
      result = LIBXSTREAM_ERROR_CONDITION;
    }
  }
  else { // 0 == argument.dims && 0 != data
    libxstream_type type = argument.type;
    if (LIBXSTREAM_TYPE_VOID == type) {
      LIBXSTREAM_CHECK_CALL(libxstream_get_autotype(*argument.shape, &type));
    }

    switch(type) {
      case LIBXSTREAM_TYPE_CHAR:  LIBXSTREAM_SNPRINTF(buffer, sizeof(buffer), "%c", *static_cast<const char*>(data)); break;
      case LIBXSTREAM_TYPE_I8:    LIBXSTREAM_SNPRINTF(buffer, sizeof(buffer), "%i", *static_cast<const int8_t*>(data)); break;
      case LIBXSTREAM_TYPE_U8:    LIBXSTREAM_SNPRINTF(buffer, sizeof(buffer), "%u", *static_cast<const uint8_t*>(data)); break;
      case LIBXSTREAM_TYPE_I16:   LIBXSTREAM_SNPRINTF(buffer, sizeof(buffer), "%i", *static_cast<const int16_t*>(data)); break;
      case LIBXSTREAM_TYPE_U16:   LIBXSTREAM_SNPRINTF(buffer, sizeof(buffer), "%u", *static_cast<const uint16_t*>(data)); break;
      case LIBXSTREAM_TYPE_I32:   LIBXSTREAM_SNPRINTF(buffer, sizeof(buffer), "%i", *static_cast<const int32_t*>(data)); break;
      case LIBXSTREAM_TYPE_U32:   LIBXSTREAM_SNPRINTF(buffer, sizeof(buffer), "%u", *static_cast<const uint32_t*>(data)); break;
      case LIBXSTREAM_TYPE_I64:   LIBXSTREAM_SNPRINTF(buffer, sizeof(buffer), "%li", *static_cast<const int64_t*>(data)); break;
      case LIBXSTREAM_TYPE_U64:   LIBXSTREAM_SNPRINTF(buffer, sizeof(buffer), "%lu", *static_cast<const uint64_t*>(data)); break;
      case LIBXSTREAM_TYPE_F32:   LIBXSTREAM_SNPRINTF(buffer, sizeof(buffer), "%f", *static_cast<const float*>(data)); break;
      case LIBXSTREAM_TYPE_F64:   LIBXSTREAM_SNPRINTF(buffer, sizeof(buffer), "%f", *static_cast<const double*>(data)); break;
      case LIBXSTREAM_TYPE_C32: {
        const float *const c = static_cast<const float*>(data);
        LIBXSTREAM_SNPRINTF(buffer, sizeof(buffer), "(%f, %f)", c[0], c[1]);
      } break;
      case LIBXSTREAM_TYPE_C64: {
        const double *const c = static_cast<const double*>(data);
        LIBXSTREAM_SNPRINTF(buffer, sizeof(buffer), "(%f, %f)", c[0], c[1]);
      } break;
      default: {
        result = LIBXSTREAM_ERROR_CONDITION;
      }
    }

    *value = buffer;
  }

  return result;
}


LIBXSTREAM_EXPORT_C LIBXSTREAM_TARGET(mic) int libxstream_get_type(const libxstream_argument* signature, size_t arg, libxstream_type* type)
{
  LIBXSTREAM_CHECK_CONDITION(0 != type);
#if defined(LIBXSTREAM_DEBUG) && (!defined(LIBXSTREAM_OFFLOAD) || (0 == LIBXSTREAM_OFFLOAD))
  size_t arity = 0;
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == libxstream_get_arity(signature, &arity) && arg < arity);
#endif
  if (0 == signature) {
    const libxstream_context& context = libxstream_context::instance();
    if (0 == (LIBXSTREAM_CALL_EXTERNAL & context.flags)) {
      signature = context.signature;
    }
  }
  LIBXSTREAM_ASSERT(0 != signature);
  *type = signature[arg].type;
  return LIBXSTREAM_ERROR_NONE;
}


LIBXSTREAM_EXPORT_C LIBXSTREAM_TARGET(mic) int libxstream_get_dims(const libxstream_argument* signature, size_t arg, size_t* dims)
{
  LIBXSTREAM_CHECK_CONDITION(0 != dims);
#if defined(LIBXSTREAM_DEBUG) && (!defined(LIBXSTREAM_OFFLOAD) || (0 == LIBXSTREAM_OFFLOAD))
  size_t arity = 0;
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == libxstream_get_arity(signature, &arity) && arg < arity);
#endif
  if (0 == signature) {
    const libxstream_context& context = libxstream_context::instance();
    if (0 == (LIBXSTREAM_CALL_EXTERNAL & context.flags)) {
      signature = context.signature;
    }
  }
  LIBXSTREAM_ASSERT(0 != signature);
  *dims = signature[arg].dims;
  return LIBXSTREAM_ERROR_NONE;
}


LIBXSTREAM_EXPORT_C LIBXSTREAM_TARGET(mic) int libxstream_get_shape(const libxstream_argument* signature, size_t arg, size_t shape[])
{
  LIBXSTREAM_CHECK_CONDITION(0 != shape);
#if defined(LIBXSTREAM_DEBUG) && (!defined(LIBXSTREAM_OFFLOAD) || (0 == LIBXSTREAM_OFFLOAD))
  size_t arity = 0;
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == libxstream_get_arity(signature, &arity) && arg < arity);
#endif
  if (0 == signature) {
    const libxstream_context& context = libxstream_context::instance();
    if (0 == (LIBXSTREAM_CALL_EXTERNAL & context.flags)) {
      signature = context.signature;
    }
  }
  LIBXSTREAM_ASSERT(0 != signature);
  const libxstream_argument& argument = signature[arg];
  const size_t dims = argument.dims;

  if (0 < dims) {
    const size_t *const src = argument.shape;
#if defined(__INTEL_COMPILER)
#   pragma loop_count min(0), max(LIBXSTREAM_MAX_NDIMS), avg(2)
#endif
    for (size_t i = 0; i < dims; ++i) shape[i] = src[i];
  }
  else {
    *shape = 0;
  }
  return LIBXSTREAM_ERROR_NONE;
}


LIBXSTREAM_EXPORT_C LIBXSTREAM_TARGET(mic) int libxstream_get_size(const libxstream_argument* signature, size_t arg, size_t* size)
{
  LIBXSTREAM_CHECK_CONDITION(0 != size);
#if defined(LIBXSTREAM_DEBUG) && (!defined(LIBXSTREAM_OFFLOAD) || (0 == LIBXSTREAM_OFFLOAD))
  size_t arity = 0;
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == libxstream_get_arity(signature, &arity) && arg < arity);
#endif
  if (0 == signature) {
    const libxstream_context& context = libxstream_context::instance();
    if (0 == (LIBXSTREAM_CALL_EXTERNAL & context.flags)) {
      signature = context.signature;
    }
  }
  LIBXSTREAM_ASSERT(0 != signature);
  const libxstream_argument& argument = signature[arg];
  *size = libxstream_linear_size(argument.dims, argument.shape, 1);
  return LIBXSTREAM_ERROR_NONE;
}


LIBXSTREAM_EXPORT_C LIBXSTREAM_TARGET(mic) int libxstream_get_elemsize(const libxstream_argument* signature, size_t arg, size_t* size)
{
  LIBXSTREAM_CHECK_CONDITION(0 != size);
#if defined(LIBXSTREAM_DEBUG) && (!defined(LIBXSTREAM_OFFLOAD) || (0 == LIBXSTREAM_OFFLOAD))
  size_t arity = 0;
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == libxstream_get_arity(signature, &arity) && arg < arity);
#endif
  if (0 == signature) {
    const libxstream_context& context = libxstream_context::instance();
    if (0 == (LIBXSTREAM_CALL_EXTERNAL & context.flags)) {
      signature = context.signature;
    }
  }
  LIBXSTREAM_ASSERT(0 != signature);
  const libxstream_argument& argument = signature[arg];
  size_t typesize = 0 != argument.dims ? 1 : *argument.shape;
  if (LIBXSTREAM_TYPE_VOID != argument.type) {
    LIBXSTREAM_CHECK_CALL(libxstream_get_typesize(argument.type, &typesize));
  }
  *size = typesize;
  return LIBXSTREAM_ERROR_NONE;
}


LIBXSTREAM_EXPORT_C LIBXSTREAM_TARGET(mic) int libxstream_get_datasize(const libxstream_argument* signature, size_t arg, size_t* size)
{
  LIBXSTREAM_CHECK_CONDITION(0 != size);
#if defined(LIBXSTREAM_DEBUG) && (!defined(LIBXSTREAM_OFFLOAD) || (0 == LIBXSTREAM_OFFLOAD))
  size_t arity = 0;
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == libxstream_get_arity(signature, &arity) && arg < arity);
#endif
  if (0 == signature) {
    const libxstream_context& context = libxstream_context::instance();
    if (0 == (LIBXSTREAM_CALL_EXTERNAL & context.flags)) {
      signature = context.signature;
    }
  }
  LIBXSTREAM_ASSERT(0 != signature);
  const libxstream_argument& argument = signature[arg];
  size_t typesize = 0 != argument.dims ? 1 : *argument.shape;
  if (LIBXSTREAM_TYPE_VOID != argument.type) {
    LIBXSTREAM_CHECK_CALL(libxstream_get_typesize(argument.type, &typesize));
  }
  *size = libxstream_linear_size(argument.dims, argument.shape, typesize);
  return LIBXSTREAM_ERROR_NONE;
}

#endif // defined(LIBXSTREAM_EXPORTED) || defined(__LIBXSTREAM) || defined(__LIBXSTREAM)
