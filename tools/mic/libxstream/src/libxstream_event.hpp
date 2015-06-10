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
#ifndef LIBXSTREAM_EVENT_HPP
#define LIBXSTREAM_EVENT_HPP

#include "libxstream.hpp"

#if defined(LIBXSTREAM_EXPORTED) || defined(__LIBXSTREAM)


struct libxstream_event {
public:
  libxstream_event();

public:
  // Number of streams the event was recorded for.
  size_t expected() const;

  // Reset the event.
  int reset();

  // Enqueue this event into the given stream; reset to start over.
  int enqueue(libxstream_stream& stream, bool reset);

  // Query whether the event already happened or not.
  int query(bool& occurred, const libxstream_stream* exclude = 0) const;

  // Wait for the event to happen.
  int wait(const libxstream_stream* exclude = 0);

private:
  class slot_type {
    libxstream_stream* m_stream;
    mutable libxstream_signal m_pending;
  public:
    slot_type(): m_stream(0), m_pending(0) {}
    slot_type(int thread, libxstream_stream& stream);
    const libxstream_stream* stream() const { return m_stream; }
    libxstream_stream* stream() { return m_stream; }
    libxstream_signal pending() const { return m_pending; }
    void pending(libxstream_signal signal) { m_pending = signal; }
  } m_slots[LIBXSTREAM_MAX_NDEVICES*LIBXSTREAM_MAX_NSTREAMS];
  size_t m_expected;
};

#endif // defined(LIBXSTREAM_EXPORTED) || defined(__LIBXSTREAM)
#endif // LIBXSTREAM_EVENT_HPP
