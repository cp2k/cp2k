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
#ifndef MULTI_DGEMM_TYPE_HPP
#define MULTI_DGEMM_TYPE_HPP

#include <libxstream.h>


class multi_dgemm_type {
public:
  class host_data_type {
  public:
    host_data_type(libxstream_function process, size_t size, const size_t split[]);
    ~host_data_type();
  public:
    libxstream_function process()   { return m_process; }
    const double* adata() const     { return m_adata; }
    const double* bdata() const     { return m_bdata; }
    double* cdata()                 { return m_cdata; }
    const size_t* idata() const     { return m_idata; }
    size_t size() const             { return m_size; }
    size_t flops() const            { return m_flops; }
    size_t max_matrix_size() const;
    size_t bytes() const;
    bool ready() const;
  private:
    libxstream_function m_process;
    double *m_adata, *m_bdata, *m_cdata;
    size_t *m_idata, m_size, m_flops;
  };

public:
  multi_dgemm_type();
  ~multi_dgemm_type();

private:
  int deinit();

public:
  int init(const char* name, host_data_type& host_data, int device, int demux, size_t max_batch);
  int operator()(size_t index, size_t size);

  libxstream_stream* stream() { return m_stream; }
  libxstream_event* event();
  size_t bytes() const;
  bool ready() const;
  int demux() const;

private:
  host_data_type* m_host_data;
  libxstream_argument* m_signature;
  libxstream_stream* m_stream;
  libxstream_event* m_event;

  double *m_adata, *m_bdata, *m_cdata;
  size_t *m_idata, m_max_batch;
};

#endif // MULTI_DGEMM_TYPE_HPP
