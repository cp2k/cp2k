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
#include "multi-dgemm-type.hpp"
#include <libxstream_begin.h>
#include <stdexcept>
#include <algorithm>
#include <cstdlib>
#include <libxstream_end.h>


multi_dgemm_type::host_data_type::host_data_type(libxstream_function process, size_t size, const size_t split[])
  : m_process(process)
  , m_adata(0), m_bdata(0), m_cdata(0), m_idata(0)
  , m_size(size), m_flops(0)
{
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_mem_allocate(-1, reinterpret_cast<void**>(&m_idata), sizeof(size_t) * (size + 1), 0));

  size_t isize = split[0];
  size_t msize = 0, n = 100, nn = n * n;
  for (size_t i = 0; i < isize; ++i) {
    m_flops += nn * (2 * n + 1);
    m_idata[i] = msize;
    msize += nn;
  }
  isize += split[1];
  n = 600, nn = n * n;
  for (size_t i = split[0]; i < isize; ++i) {
    m_flops += nn * (2 * n + 1);
    m_idata[i] = msize;
    msize += nn;
  }
  n = 1000, nn = n * n;
  for (size_t i = isize; i < size; ++i) {
    m_flops += nn * (2 * n + 1);
    m_idata[i] = msize;
    msize += nn;
  }
  m_idata[size] = msize;

  LIBXSTREAM_CHECK_CALL_THROW(libxstream_mem_allocate(-1, reinterpret_cast<void**>(&m_adata), sizeof(double) * msize, 0));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_mem_allocate(-1, reinterpret_cast<void**>(&m_bdata), sizeof(double) * msize, 0));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_mem_allocate(-1, reinterpret_cast<void**>(&m_cdata), sizeof(double) * msize, 0));

  static const double scale = 1.0 / RAND_MAX;
  for (size_t i = 0; i < msize; ++i) {
    m_adata[i] = scale * (2 * std::rand() - RAND_MAX);
    m_bdata[i] = scale * (2 * std::rand() - RAND_MAX);
    m_cdata[i] = 0;
  }
}


multi_dgemm_type::host_data_type::~host_data_type()
{
  /*LIBXSTREAM_CHECK_CALL_THROW*/(libxstream_mem_deallocate(-1, m_adata));
  /*LIBXSTREAM_CHECK_CALL_THROW*/(libxstream_mem_deallocate(-1, m_bdata));
  /*LIBXSTREAM_CHECK_CALL_THROW*/(libxstream_mem_deallocate(-1, m_cdata));
  /*LIBXSTREAM_CHECK_CALL_THROW*/(libxstream_mem_deallocate(-1, m_idata));
}


size_t multi_dgemm_type::host_data_type::max_matrix_size() const
{
  LIBXSTREAM_ASSERT(0 == m_size || 0 == m_idata[0]);
  size_t result = 0, i0 = 0;
  for (size_t i = 0; i < m_size; ++i) {
    const size_t i1 = m_idata[i+1];
    result = std::max(result, i1 - i0);
    i0 = i1;
  }
  return result;
}


size_t multi_dgemm_type::host_data_type::bytes() const
{
  return sizeof(double) * m_idata[m_size] * 3 + sizeof(size_t) * m_size;
}


bool multi_dgemm_type::host_data_type::ready() const
{
  LIBXSTREAM_ASSERT(0 == m_process || (m_adata && m_bdata && m_cdata && m_idata));
  return 0 != m_process;
}


multi_dgemm_type::multi_dgemm_type()
  : m_host_data(0), m_signature(0), m_stream(0), m_event(0)
  , m_adata(0), m_bdata(0), m_cdata(0)
  , m_idata(0), m_max_batch(0)
{}


multi_dgemm_type::~multi_dgemm_type()
{
  /*LIBXSTREAM_CHECK_CALL_THROW*/(deinit());
}


int multi_dgemm_type::deinit()
{
  if (m_host_data) {
    int device = -1;
    LIBXSTREAM_CHECK_CALL(libxstream_stream_device(m_stream, &device));
    LIBXSTREAM_CHECK_CALL(libxstream_fn_destroy_signature(m_signature));
    LIBXSTREAM_CHECK_CALL(libxstream_stream_destroy(m_stream));
    LIBXSTREAM_CHECK_CALL(libxstream_event_destroy(m_event));
    LIBXSTREAM_CHECK_CALL(libxstream_mem_deallocate(device, m_adata));
    LIBXSTREAM_CHECK_CALL(libxstream_mem_deallocate(device, m_bdata));
    LIBXSTREAM_CHECK_CALL(libxstream_mem_deallocate(device, m_cdata));
    LIBXSTREAM_CHECK_CALL(libxstream_mem_deallocate(device, m_idata));
    m_host_data = 0;
#if defined(LIBXSTREAM_DEBUG)
    m_max_batch = 0;
    m_signature = 0;
    m_stream = 0;
    m_event = 0;
    m_adata = 0;
    m_bdata = 0;
    m_cdata = 0;
    m_idata = 0;
#endif
  }

  return LIBXSTREAM_ERROR_NONE;
}


int multi_dgemm_type::init(const char* name, host_data_type& host_data, int device, int demux, size_t max_batch)
{
  LIBXSTREAM_CHECK_CALL(deinit());
  const size_t max_msize = max_batch * host_data.max_matrix_size();
  m_host_data = &host_data;
  m_max_batch = max_batch;

  LIBXSTREAM_CHECK_CALL(libxstream_stream_create(&m_stream, device, demux, 0, name));
  LIBXSTREAM_CHECK_CALL(libxstream_mem_allocate(device, reinterpret_cast<void**>(&m_adata), sizeof(double) * max_msize, 0));
  LIBXSTREAM_CHECK_CALL(libxstream_mem_allocate(device, reinterpret_cast<void**>(&m_bdata), sizeof(double) * max_msize, 0));
  LIBXSTREAM_CHECK_CALL(libxstream_mem_allocate(device, reinterpret_cast<void**>(&m_cdata), sizeof(double) * max_msize, 0));
  LIBXSTREAM_CHECK_CALL(libxstream_mem_allocate(device, reinterpret_cast<void**>(&m_idata), sizeof(size_t) * max_batch, 0));

  LIBXSTREAM_CHECK_CALL(libxstream_fn_create_signature(&m_signature, 6));
  LIBXSTREAM_CHECK_CALL(libxstream_fn_input (m_signature, 2, m_idata, libxstream_map_to_type(m_idata), 1, &max_msize));
  LIBXSTREAM_CHECK_CALL(libxstream_fn_input (m_signature, 3, m_adata, libxstream_map_to_type(m_adata), 1, &max_msize));
  LIBXSTREAM_CHECK_CALL(libxstream_fn_input (m_signature, 4, m_bdata, libxstream_map_to_type(m_bdata), 1, &max_msize));
  LIBXSTREAM_CHECK_CALL(libxstream_fn_output(m_signature, 5, m_cdata, libxstream_map_to_type(m_cdata), 1, &max_msize));

  return LIBXSTREAM_ERROR_NONE;
}


int multi_dgemm_type::operator()(size_t index, size_t size)
{
  LIBXSTREAM_CHECK_CONDITION(ready() && (index + size) <= m_host_data->size());

  if (0 < size) {
    if (0 == demux()) {
      // This manual synchronization prevents multiple threads from queuing work into the *same* stream (at the same time).
      // This is only needed if the stream was created without demux support in order to rely on manual synchronization.
      LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_stream_lock(m_stream));
    }
    const size_t i0 = m_host_data->idata()[index], i1 = m_host_data->idata()[index+size];
    LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_memcpy_h2d(m_host_data->adata() + i0, m_adata, sizeof(double) * (i1 - i0), m_stream));
    LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_memcpy_h2d(m_host_data->bdata() + i0, m_bdata, sizeof(double) * (i1 - i0), m_stream));
    // transferring cdata is part of the benchmark; since it is all zeros we could do better with libxstream_memset_zero
    LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_memcpy_h2d(m_host_data->cdata() + i0, m_cdata, sizeof(double) * (i1 - i0), m_stream));
    LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_memcpy_h2d(m_host_data->idata() + index, m_idata, sizeof(size_t) * size, m_stream));
#if defined(LIBXSTREAM_DEBUG)
    size_t n = 0;
    LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == libxstream_fn_nargs(m_signature, &n) && 6 == n);
#endif
    const size_t nn = i1 - m_host_data->idata()[index+size-1];
    LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_fn_input(m_signature, 0, &size, libxstream_map_to_type(size), 0, 0));
    LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_fn_input(m_signature, 1,   &nn, libxstream_map_to_type(nn  ), 0, 0));
    LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == libxstream_get_arity(m_signature, &n) && 6 == n);
    LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_fn_call(m_host_data->process(), m_signature, m_stream, LIBXSTREAM_CALL_DEFAULT));
    LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_memcpy_d2h(m_cdata, m_host_data->cdata() + i0, sizeof(double) * (i1 - i0), m_stream));
    if (0 == demux()) {
      LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_stream_unlock(m_stream));
    }
  }

  return LIBXSTREAM_ERROR_NONE;
}


libxstream_event* multi_dgemm_type::event()
{
  if (0 == m_event) {
    libxstream_event_create(&m_event);
    LIBXSTREAM_ASSERT(0 != m_event);
  }
  return m_event;
}


bool multi_dgemm_type::ready() const
{
  LIBXSTREAM_ASSERT(0 == m_host_data || (m_signature && m_stream && m_adata && m_bdata && m_cdata && m_idata));
  return 0 != m_host_data;
}


int multi_dgemm_type::demux() const
{
  LIBXSTREAM_ASSERT(ready());
  int value = 0;
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_stream_demux(m_stream, &value));
  return value;
}


size_t multi_dgemm_type::bytes() const
{
  LIBXSTREAM_ASSERT(ready());
  return m_max_batch * m_host_data->max_matrix_size() * (3 * sizeof(double) + sizeof(size_t));
}
