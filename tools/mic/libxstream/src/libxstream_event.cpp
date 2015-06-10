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
#include "libxstream_event.hpp"
#include "libxstream_capture.hpp"

#include <libxstream_begin.h>
#include <algorithm>
#include <libxstream_end.h>

#if defined(LIBXSTREAM_OFFLOAD) && (0 != LIBXSTREAM_OFFLOAD)
# include <offload.h>
#endif

// allows to wait for an event issued prior to the pending signal
//#define LIBXSTREAM_EVENT_WAIT_PAST


libxstream_event::libxstream_event()
  : m_expected(0)
{}


size_t libxstream_event::expected() const
{
  LIBXSTREAM_ASSERT((LIBXSTREAM_MAX_NDEVICES * LIBXSTREAM_MAX_NSTREAMS) >= m_expected);
  return m_expected;
}


int libxstream_event::reset()
{ 
  int result = LIBXSTREAM_ERROR_NONE;

  LIBXSTREAM_ASYNC_BEGIN(0, -2/*invalid device*/, m_slots, &m_expected)
  {
#if defined(LIBXSTREAM_DEBUG)
    std::fill_n(ptr<slot_type,1>(), LIBXSTREAM_MAX_NDEVICES * LIBXSTREAM_MAX_NSTREAMS, slot_type());
#endif
    *ptr<size_t,2>() = 0;
  }
  LIBXSTREAM_ASYNC_END(LIBXSTREAM_CALL_DEFAULT | LIBXSTREAM_CALL_UNLOCK, result);

  return result;
}


int libxstream_event::enqueue(libxstream_stream& stream, bool reset)
{ 
  int result = LIBXSTREAM_ERROR_NONE;

  LIBXSTREAM_ASYNC_BEGIN(stream, m_slots, &m_expected, reset)
  {
    slot_type *const slots = ptr<slot_type,0>();
    size_t& expected = *ptr<size_t,1>();
    const bool reset = val<bool,2>();

    if (reset) {
#if defined(LIBXSTREAM_DEBUG)
      std::fill_n(slots, LIBXSTREAM_MAX_NDEVICES * LIBXSTREAM_MAX_NSTREAMS, slot_type());
#endif
      expected = 0;
    }

#if defined(LIBXSTREAM_DEBUG)
    LIBXSTREAM_ASSERT((LIBXSTREAM_MAX_NDEVICES * LIBXSTREAM_MAX_NSTREAMS) > expected);
#endif
    slot_type& slot = slots[expected];
    slot = slot_type(thread(), *LIBXSTREAM_ASYNC_STREAM);
    ++expected;
  }
  LIBXSTREAM_ASYNC_END(LIBXSTREAM_CALL_DEFAULT | LIBXSTREAM_CALL_UNLOCK, result);

  return result;
}


int libxstream_event::query(bool& occurred, const libxstream_stream* exclude) const
{
  int result = LIBXSTREAM_ERROR_NONE;

  LIBXSTREAM_ASYNC_BEGIN(0, -2/*invalid device*/, &occurred, exclude, m_slots, &m_expected)
  {
    const libxstream_stream *const exclude = ptr<const libxstream_stream,2>();
    slot_type *const slots = ptr<slot_type,3>();
    const size_t expected = *ptr<const size_t,4>();
    bool occurred = true; // everythig "occurred" if nothing is expected

    for (size_t i = 0; i < expected; ++i) {
      slot_type& slot = slots[i];
      const libxstream_signal pending_slot = slot.pending();
      libxstream_stream *const stream = slot.stream();

      if (exclude != stream && 0 != pending_slot) {
        const libxstream_signal pending_stream = stream ? stream->pending(thread()) : 0;

        if (0 != pending_stream) {
#if defined(LIBXSTREAM_EVENT_WAIT_PAST)
          const libxstream_signal signal = pending_slot;
#else
          const libxstream_signal signal = pending_stream;
#endif
#if defined(LIBXSTREAM_OFFLOAD) && (0 != LIBXSTREAM_OFFLOAD) && !defined(__MIC__) && defined(LIBXSTREAM_ASYNC) && (0 != (2*LIBXSTREAM_ASYNC+1)/2)
          if (0 != _Offload_signaled(stream->device(), reinterpret_cast<void*>(signal)))
#endif
          {
            if (signal == pending_stream) {
              stream->pending(thread(), 0);
            }
            slot.pending(0);
          }
        }
        else {
          slot.pending(0);
        }

        occurred = occurred && 0 == slot.pending();
      }
    }

    *ptr<bool,1>() = occurred;
  }
  LIBXSTREAM_ASYNC_END(LIBXSTREAM_CALL_DEFAULT | LIBXSTREAM_CALL_WAIT, result);

  return result;
}


int libxstream_event::wait(const libxstream_stream* exclude)
{
  int result = LIBXSTREAM_ERROR_NONE;

  LIBXSTREAM_ASYNC_BEGIN(0, -2/*invalid device*/, exclude, m_slots, &m_expected)
  {
    const libxstream_stream *const exclude = ptr<const libxstream_stream,1>();
    slot_type *const slots = ptr<slot_type,2>();
    const size_t expected = *ptr<const size_t,3>();
    size_t completed = 0;

    for (size_t i = 0; i < expected; ++i) {
      slot_type& slot = slots[i];
      libxstream_stream *const stream = slot.stream();
      const libxstream_signal pending_slot = slot.pending();

      if (exclude != stream && 0 != pending_slot) {
        const libxstream_signal pending_stream = stream ? stream->pending(thread()) : 0;

        if (0 != pending_stream) {
  #if defined(LIBXSTREAM_EVENT_WAIT_PAST)
          const libxstream_signal signal = pending_slot;
  #else
          const libxstream_signal signal = pending_stream;
  #endif
  #if defined(LIBXSTREAM_OFFLOAD) && defined(LIBXSTREAM_ASYNC) && (0 != (2*LIBXSTREAM_ASYNC+1)/2)
          if (0 <= stream->device()) {
            LIBXSTREAM_ASYNC_DEVICE_UPDATE(stream->device());
  #         pragma offload_wait LIBXSTREAM_ASYNC_TARGET wait(signal)
          }
  #endif
          if (signal == pending_stream) {
            stream->pending(thread(), 0);
          }
          slot.pending(0);
        }
        else {
          slot.pending(0);
        }

        completed += 0 == slot.pending() ? 1 : 0;
      }
    }

    LIBXSTREAM_ASSERT(completed <= expected);
    // reset only if all slots fired (need to iterate all slots otherwise)
    if (completed == *ptr<const size_t,3>()) {
      *ptr<size_t,3>() = 0;
    }
  }
  LIBXSTREAM_ASYNC_END(LIBXSTREAM_CALL_DEFAULT | LIBXSTREAM_CALL_WAIT | LIBXSTREAM_CALL_UNLOCK, result);

  return result;
}


libxstream_event::slot_type::slot_type(int thread, libxstream_stream& stream)
  : m_stream(&stream) // no need to lock the stream
  , m_pending(stream.pending(thread))
{}

#endif // defined(LIBXSTREAM_EXPORTED) || defined(__LIBXSTREAM)
