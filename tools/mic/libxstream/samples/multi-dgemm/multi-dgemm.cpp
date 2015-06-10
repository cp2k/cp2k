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

#include "../../include/libxstream_begin.h"
#include <stdexcept>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <vector>
#include <cmath>
#if defined(_OPENMP)
# include <omp.h>
#endif
#include "../../include/libxstream_end.h"

//#define MULTI_DGEMM_USE_NESTED
#define MULTI_DGEMM_USE_SYNC 1
#define MULTI_DGEMM_USE_CHECK

#define DGEMM dgemm_


LIBXSTREAM_EXTERN_C LIBXSTREAM_TARGET(mic) void DGEMM(
  const char*, const char*, const int*, const int*, const int*,
  const double*, const double*, const int*, const double*, const int*,
  const double*, double*, const int*);


LIBXSTREAM_TARGET(mic) void process(LIBXSTREAM_INVAL(size_t) size, LIBXSTREAM_INVAL(size_t) nn, const size_t* idata,
  const double* adata, const double* bdata, double* cdata)
{
  if (0 < LIBXSTREAM_GETVAL(size)) {
    static const double alpha = 1, beta = 1;
    static const char trans = 'N';
    const int isize = static_cast<int>(size);
    const size_t base = idata[0];

#if defined(_OPENMP) && defined(MULTI_DGEMM_USE_NESTED)
    const int nthreads = omp_get_max_threads() / LIBXSTREAM_GETVAL(size);
    const int dynamic = omp_get_dynamic(), nested = omp_get_nested();
    omp_set_dynamic(0);
    omp_set_nested(1);
#   pragma omp parallel for schedule(dynamic,1) num_threads(LIBXSTREAM_GETVAL(size))
#endif
    for (int i = 0; i < isize; ++i) {
#if defined(_OPENMP) && defined(MULTI_DGEMM_USE_NESTED)
      omp_set_num_threads(nthreads);
#endif
      LIBXSTREAM_ASSERT(base <= idata[i]);
      const size_t i0 = idata[i], i1 = (i + 1) < isize ? idata[i+1] : (i0 + LIBXSTREAM_GETVAL(nn)), n2 = i1 - i0, offset = i0 - base;
      const int n = static_cast<int>(std::sqrt(static_cast<double>(n2)) + 0.5);
      DGEMM(&trans, &trans, &n, &n, &n, &alpha, adata + offset, &n, bdata + offset, &n, &beta, cdata + offset, &n);
    }

#if defined(_OPENMP) && defined(MULTI_DGEMM_USE_NESTED)
    omp_set_dynamic(dynamic);
    omp_set_nested(nested);
#endif
  }
}


int main(int argc, char* argv[])
{
  try {
    const int nitems = std::max(1 < argc ? std::atoi(argv[1]) : 60, 0);
    const int nbatch = std::max(2 < argc ? std::atoi(argv[2]) : 10, 1);
    const int nstreams = std::min(std::max(3 < argc ? std::atoi(argv[3]) : 2, 1), LIBXSTREAM_MAX_NSTREAMS);
#if defined(MULTI_DGEMM_USE_SYNC)
    const int demux = 4 < argc ? std::atoi(argv[4]) : 1;
#else
    const int demux = -1;
#endif
    size_t ndevices = 0;
    if (LIBXSTREAM_ERROR_NONE != libxstream_get_ndevices(&ndevices) || 0 == ndevices) {
      throw std::runtime_error("no device found!");
    }
#if !defined(_OPENMP)
    fprintf(stderr, "Warning: OpenMP support needed for performance results.\n");
#endif

    fprintf(stdout, "Initializing %i device%s and host data...", static_cast<int>(ndevices), 1 == ndevices ? "" : "s");
    const size_t split[] = { size_t(nitems * 18.0 / 250.0 + 0.5), size_t(nitems * 74.0 / 250.0 + 0.5) };
    multi_dgemm_type::host_data_type host_data(reinterpret_cast<libxstream_function>(&process), nitems, split);
    fprintf(stdout, " %.1f MB\n", host_data.bytes() * 1E-6);

    fprintf(stdout, "Initializing %i stream%s per device...", nstreams, 1 < nstreams ? "s" : "");
    const size_t nstreams_total = ndevices * nstreams;
    multi_dgemm_type multi_dgemm[LIBXSTREAM_MAX_NSTREAMS];
    for (size_t i = 0; i < nstreams_total; ++i) {
      char name[128];
      LIBXSTREAM_SNPRINTF(name, sizeof(name), "Stream %i", static_cast<int>(i + 1));
      LIBXSTREAM_CHECK_CALL_THROW(multi_dgemm[i].init(name, host_data, static_cast<int>(i % ndevices), demux, static_cast<size_t>(nbatch)));
    }
    if (0 < nstreams_total) {
      fprintf(stdout, " %.1f MB\n", nstreams * multi_dgemm[0].bytes() * 1E-6);
    }

    const int nbatches = (nitems + nbatch - 1) / nbatch;
    fprintf(stdout, "Running %i batch%s of %i item%s...\n", nbatches,
      1 < nbatches ? "es" : "", std::min(nbatch, nitems),
      1 < nbatch ? "s" : "");

#if defined(_OPENMP)
# if !defined(LIBXSTREAM_OFFLOAD)
    omp_set_dynamic(0);
    omp_set_nested(0);
# endif
    const double start = omp_get_wtime();
#   pragma omp parallel for schedule(dynamic)
#endif
    for (int i = 0; i < nitems; i += nbatch) {
      const size_t j = i / nbatch, n = j % nstreams_total;
      multi_dgemm_type& call = multi_dgemm[n];
      LIBXSTREAM_CHECK_CALL_ASSERT(call(i, std::min(nbatch, nitems - i)));
#if defined(MULTI_DGEMM_USE_SYNC) && (1 <= MULTI_DGEMM_USE_SYNC)
      LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_event_record(call.event(), call.stream()));
#endif
      // synchronize every Nth iteration with N being the total number of streams
      if (n == (nstreams_total - 1)) {
        for (size_t k = 0; k < nstreams_total; ++k) {
#if defined(MULTI_DGEMM_USE_SYNC)
# if (2 <= (MULTI_DGEMM_USE_SYNC))
          // wait for an event within a stream
          LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_stream_wait_event(multi_dgemm[0].stream(), multi_dgemm[k].event()));
# elif (1 <= (MULTI_DGEMM_USE_SYNC))
          // wait for an event on the host
          LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_event_synchronize(multi_dgemm[k].event()));
# else
          // wait for all work in a stream
          LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_stream_sync(multi_dgemm[k].stream()));
# endif
#endif
        }
      }
    }

    // sync all streams to complete any pending work
    LIBXSTREAM_CHECK_CALL_THROW(libxstream_stream_sync(0));

#if defined(_OPENMP)
    const double duration = omp_get_wtime() - start;
    fprintf(stdout, "Performance: %.1f GFLOPS/s (%s)\n", host_data.flops() * 1E-9 / duration,
      0 == demux ? "manual locking" : (0 < demux ? "synchronization" : "automatic locking"));
    fprintf(stdout, "Duration: %.1f s\n", duration);
#endif

#if defined(MULTI_DGEMM_USE_CHECK)
    std::vector<double> expected(host_data.max_matrix_size());
    const size_t testbatchsize = 1;
    double max_error = 0;
    size_t i0 = 0;
    for (int i = 0; i < nitems; ++i) {
      const size_t i1 = host_data.idata()[i+1];
      const int nn = static_cast<int>(i1 - i0);
      std::fill_n(&expected[0], nn, 0.0);
      process(LIBXSTREAM_SETVAL(testbatchsize), LIBXSTREAM_SETVAL(nn), host_data.idata() + i, host_data.adata() + i0, host_data.bdata() + i0, &expected[0]);
      for (int n = 0; n < nn; ++n) max_error = std::max(max_error, std::abs(expected[n] - host_data.cdata()[i0+n]));
      i0 = i1;
    }
    fprintf(stdout, "Error: %g\n", max_error);
#endif
    fprintf(stdout, "Finished\n");
  }
  catch(const std::exception& e) {
    fprintf(stderr, "Error: %s\n", e.what());
    return EXIT_FAILURE;
  }
  catch(...) {
    fprintf(stderr, "Error: unknown exception caught!\n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
