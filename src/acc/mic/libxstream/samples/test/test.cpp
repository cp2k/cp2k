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
#include "test.hpp"
#include <libxstream_begin.h>
#include <stdexcept>
#include <algorithm>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cstdio>
#if defined(_OPENMP)
# include <omp.h>
#endif
#include <libxstream_end.h>


namespace test_internal {

LIBXSTREAM_TARGET(mic) void check(libxstream_bool& result, LIBXSTREAM_INVAL(char) pattern, const void* buffer, LIBXSTREAM_INVAL(size_t) size)
{
  size_t value = 0;
  bool ok = true;
  // check function is called with using LIBXSTREAM hence introspection may not be available
  if (LIBXSTREAM_ERROR_NONE == libxstream_get_argument(buffer, &value)) {
    ok = ok && 2 == value;
    ok = ok && LIBXSTREAM_ERROR_NONE == libxstream_get_shape(0, value, &value) && LIBXSTREAM_GETVAL(size) == value;
    ok = ok && LIBXSTREAM_ERROR_NONE == libxstream_get_arity(0, &value) && 4 == value;
    ok = ok && LIBXSTREAM_ERROR_NONE == libxstream_get_dims(0, 0, &value) && 0 == value;
    ok = ok && LIBXSTREAM_ERROR_NONE == libxstream_get_dims(0, 1, &value) && 0 == value;
    ok = ok && LIBXSTREAM_ERROR_NONE == libxstream_get_dims(0, 2, &value) && 1 == value;
    ok = ok && LIBXSTREAM_ERROR_NONE == libxstream_get_dims(0, 3, &value) && 0 == value;
    const void* data = 0;
    ok = ok && LIBXSTREAM_ERROR_NONE == libxstream_get_data(0, 0, &data) && result == *static_cast<const libxstream_bool*>(data);
    ok = ok && LIBXSTREAM_ERROR_NONE == libxstream_get_data(0, 1, &data) && LIBXSTREAM_GETVAL(pattern) == *static_cast<const char*>(data);
    ok = ok && LIBXSTREAM_ERROR_NONE == libxstream_get_data(0, 2, &data) && buffer == data;
    ok = ok && LIBXSTREAM_ERROR_NONE == libxstream_get_data(0, 3, &data) && LIBXSTREAM_GETVAL(size) == *static_cast<const size_t*>(data);
  }

  const char *const values = reinterpret_cast<const char*>(buffer);
  for (size_t i = 0; i < LIBXSTREAM_GETVAL(size) && ok; ++i) {
    ok = LIBXSTREAM_GETVAL(pattern) == values[i];
  }

  result = ok ? LIBXSTREAM_TRUE : LIBXSTREAM_FALSE;
}

LIBXSTREAM_TARGET(mic) void complex_c(libxstream_bool* ok,
  const float* c32, LIBXSTREAM_INVAL(float) freal, LIBXSTREAM_INVAL(float) fimag,
  const double* c64, LIBXSTREAM_INVAL(double) dreal, LIBXSTREAM_INVAL(double) dimag)
{
  LIBXSTREAM_ASSERT(ok && c32 && c64);
  libxstream_bool result = LIBXSTREAM_TRUE;
  result = LIBXSTREAM_FALSE != result && LIBXSTREAM_GETVAL(freal) == c32[0];
  result = LIBXSTREAM_FALSE != result && LIBXSTREAM_GETVAL(fimag) == c32[1];
  result = LIBXSTREAM_FALSE != result && LIBXSTREAM_GETVAL(dreal) == c64[0];
  result = LIBXSTREAM_FALSE != result && LIBXSTREAM_GETVAL(dimag) == c64[1];
  *ok = result;
}

LIBXSTREAM_TARGET(mic) void complex_cpp(libxstream_bool& ok,
  const std::complex<float>& c32, LIBXSTREAM_INVAL(float) freal, LIBXSTREAM_INVAL(float) fimag,
  const std::complex<double>& c64, LIBXSTREAM_INVAL(double) dreal, LIBXSTREAM_INVAL(double) dimag)
{
  bool result = true;
  result = result && LIBXSTREAM_GETVAL(freal) == c32.real();
  result = result && LIBXSTREAM_GETVAL(fimag) == c32.imag();
  result = result && LIBXSTREAM_GETVAL(dreal) == c64.real();
  result = result && LIBXSTREAM_GETVAL(dimag) == c64.imag();
  ok = result ? LIBXSTREAM_TRUE : LIBXSTREAM_FALSE;
}

} // namespace test_internal

/* workaround for issue "cannot find address of function"; compile using "make.sh -g" */
const libxstream_function check = reinterpret_cast<libxstream_function>(test_internal::check);
const libxstream_function complex_c = reinterpret_cast<libxstream_function>(test_internal::complex_c);
const libxstream_function complex_cpp = reinterpret_cast<libxstream_function>(test_internal::complex_cpp);


test_type::test_type(int device)
  : m_signature(0), m_stream(0), m_event(0)
  , m_host_mem(0), m_dev_mem1(0), m_dev_mem2(0)
{
  size_t mem_free = 0, mem_avail = 0;
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_mem_info(device, &mem_free, &mem_avail));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_stream_create(&m_stream, device, 0, 0, 0));

  const size_t size = 4711u * 1024u;
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_mem_allocate(-1, &m_host_mem, size, 0));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_mem_allocate(device, &m_dev_mem1, size, 0));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_mem_allocate(device, &m_dev_mem2, size, 0));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_mem_info(device, &mem_free, &mem_avail));

  const void* real = 0;
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_mem_pointer(-1, m_host_mem, &real));
  LIBXSTREAM_CHECK_CONDITION_THROW(m_host_mem == real);
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_mem_pointer(device, m_dev_mem1, &real));
#if defined(LIBXSTREAM_OFFLOAD) && (0 != LIBXSTREAM_OFFLOAD)
  LIBXSTREAM_CHECK_CONDITION_THROW(m_dev_mem1 != real);
#endif
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_mem_pointer(device, m_dev_mem2, &real));
#if defined(LIBXSTREAM_OFFLOAD) && (0 != LIBXSTREAM_OFFLOAD)
  LIBXSTREAM_CHECK_CONDITION_THROW(m_dev_mem2 != real);
#endif

  const char pattern_a = 'a', pattern_b = 'b';
  LIBXSTREAM_ASSERT(pattern_a != pattern_b);
  std::fill_n(reinterpret_cast<char*>(m_host_mem), size, pattern_a);
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_memcpy_h2d(m_host_mem, m_dev_mem1, size, m_stream));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_memcpy_d2d(m_dev_mem1, m_dev_mem2, size, m_stream));

  size_t typesize = 0;
  unsigned char *const puc = 0, auc[16], *const* ppuc = &puc, *apuc[16];
  LIBXSTREAM_CHECK_CONDITION_THROW(LIBXSTREAM_ERROR_NONE == libxstream_get_typesize(libxstream_map_to_type(puc ), &typesize) && 1 == typesize);
  LIBXSTREAM_CHECK_CONDITION_THROW(LIBXSTREAM_ERROR_NONE == libxstream_get_typesize(libxstream_map_to_type(auc ), &typesize) && 1 == typesize);
  LIBXSTREAM_CHECK_CONDITION_THROW(LIBXSTREAM_ERROR_NONE == libxstream_get_typesize(libxstream_map_to_type(ppuc), &typesize) && sizeof(void*) == typesize);
  LIBXSTREAM_CHECK_CONDITION_THROW(LIBXSTREAM_ERROR_NONE == libxstream_get_typesize(libxstream_map_to_type(apuc), &typesize) && sizeof(void*) == typesize);

  size_t nargs = 0, arity = 0;
  libxstream_bool ok = LIBXSTREAM_FALSE;
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_fn_create_signature(&m_signature, 4));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_fn_nargs(m_signature, &nargs));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_get_arity(m_signature, &arity));
  LIBXSTREAM_CHECK_CONDITION_THROW(4 == nargs && 0 == arity);

  LIBXSTREAM_CHECK_CALL_THROW(libxstream_fn_output(m_signature, 0, &ok, libxstream_map_to<libxstream_bool>::type(), 0, 0));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_fn_nargs(m_signature, &nargs));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_get_arity(m_signature, &arity));
  LIBXSTREAM_CHECK_CONDITION_THROW(4 == nargs && 1 == arity);
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_get_elemsize(m_signature, 0, &typesize));
  LIBXSTREAM_CHECK_CONDITION_THROW(sizeof(libxstream_bool) == typesize);
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_get_datasize(m_signature, 0, &typesize));
  LIBXSTREAM_CHECK_CONDITION_THROW(sizeof(libxstream_bool) == typesize);

  LIBXSTREAM_CHECK_CALL_THROW(libxstream_fn_input (m_signature, 1, &pattern_a, libxstream_map_to<char>::type(), 0, 0));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_fn_nargs(m_signature, &nargs));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_get_arity(m_signature, &arity));
  LIBXSTREAM_CHECK_CONDITION_THROW(4 == nargs && 2 == arity);
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_get_elemsize(m_signature, 1, &typesize));
  LIBXSTREAM_CHECK_CONDITION_THROW(1 == typesize);
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_get_datasize(m_signature, 1, &typesize));
  LIBXSTREAM_CHECK_CONDITION_THROW(1 == typesize);

  LIBXSTREAM_CHECK_CALL_THROW(libxstream_fn_input (m_signature, 2, m_dev_mem1, LIBXSTREAM_TYPE_VOID, 1, &size));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_fn_nargs(m_signature, &nargs));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_get_arity(m_signature, &arity));
  LIBXSTREAM_CHECK_CONDITION_THROW(4 == nargs && 3 == arity);
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_get_elemsize(m_signature, 2, &typesize));
  LIBXSTREAM_CHECK_CONDITION_THROW(1 == typesize);
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_get_datasize(m_signature, 2, &typesize));
  LIBXSTREAM_CHECK_CONDITION_THROW(size == typesize);
  // for testing purpose the following argument is weak-typed instead of (..., libxstream_map_to<size_t>::type(), 0, 0)
  typesize = sizeof(size_t);
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_fn_input (m_signature, 3, &size, LIBXSTREAM_TYPE_VOID, 0, &typesize));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_fn_nargs(m_signature, &nargs));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_get_arity(m_signature, &arity));
  LIBXSTREAM_CHECK_CONDITION_THROW(4 == nargs && 4 == arity);
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_get_elemsize(m_signature, 3, &typesize));
  LIBXSTREAM_CHECK_CONDITION_THROW(sizeof(size_t) == typesize);
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_get_datasize(m_signature, 3, &typesize));
  LIBXSTREAM_CHECK_CONDITION_THROW(sizeof(size_t) == typesize);

  //const libxstream_function check = reinterpret_cast<libxstream_function>(test_internal::check);
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_fn_call(check, m_signature, m_stream, LIBXSTREAM_CALL_DEFAULT));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_event_create(&m_event));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_event_record(m_event, m_stream));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_event_synchronize(m_event));
  LIBXSTREAM_CHECK_CONDITION_THROW(LIBXSTREAM_FALSE != ok);

  libxstream_argument* signature = 0;
  const std::complex<float>  c32( 1.05f, 19.81f);
  const std::complex<double> c64(25.07 , 19.75 );
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_fn_signature(&signature));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_fn_output(signature, 0,  &ok, libxstream_map_to_type(ok ), 0, 0));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_fn_input (signature, 1, &c32, libxstream_map_to_type(c32), 0, 0));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_fn_input (signature, 2,  reinterpret_cast<const float*>(&c32) + 0, libxstream_map_to_type(c32.real()), 0, 0));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_fn_input (signature, 3,  reinterpret_cast<const float*>(&c32) + 1, libxstream_map_to_type(c32.imag()), 0, 0));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_fn_input (signature, 4, &c64, libxstream_map_to_type(c64), 0, 0));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_fn_input (signature, 5, reinterpret_cast<const double*>(&c64) + 0, libxstream_map_to_type(c64.real()), 0, 0));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_fn_input (signature, 6, reinterpret_cast<const double*>(&c64) + 1, libxstream_map_to_type(c64.imag()), 0, 0));

  //const libxstream_function complex_c = reinterpret_cast<libxstream_function>(test_internal::complex_c);
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_fn_call(complex_c, signature, m_stream, LIBXSTREAM_CALL_DEFAULT));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_stream_sync(m_stream));
  LIBXSTREAM_CHECK_CONDITION_THROW(LIBXSTREAM_FALSE != ok);

  //const libxstream_function complex_cpp = reinterpret_cast<libxstream_function>(test_internal::complex_cpp);
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_fn_call(complex_cpp, signature, m_stream, LIBXSTREAM_CALL_DEFAULT));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_stream_sync(m_stream));
  LIBXSTREAM_CHECK_CONDITION_THROW(LIBXSTREAM_FALSE != ok);

  std::fill_n(reinterpret_cast<char*>(m_host_mem), size, pattern_b);
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_memcpy_d2h(m_dev_mem2, m_host_mem, size, m_stream));

  const size_t size2 = size / 2;
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_memset_zero(m_dev_mem1, size2, m_stream));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_memset_zero(reinterpret_cast<char*>(m_dev_mem1) + size2, size - size2, m_stream));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_event_record(m_event, m_stream));

  int has_occured = 0;
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_event_query(m_event, &has_occured));
  if (0 == has_occured) {
    LIBXSTREAM_CHECK_CALL_THROW(libxstream_event_synchronize(m_event));
  }

  test_internal::check(ok, LIBXSTREAM_SETVAL(pattern_a), m_host_mem, LIBXSTREAM_SETVAL(size));
  LIBXSTREAM_CHECK_CONDITION_THROW(LIBXSTREAM_FALSE != ok);

  LIBXSTREAM_CHECK_CALL_THROW(libxstream_memcpy_d2h(m_dev_mem1, m_host_mem, size2, m_stream));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_memcpy_d2h(reinterpret_cast<const char*>(m_dev_mem1) + size2, reinterpret_cast<char*>(m_host_mem) + size2, size - size2, m_stream));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_event_record(m_event, m_stream));
  LIBXSTREAM_CHECK_CALL_THROW(libxstream_stream_sync(m_stream));

  LIBXSTREAM_CHECK_CALL_THROW(libxstream_event_query(m_event, &has_occured));
  LIBXSTREAM_CHECK_CONDITION_THROW(0 != has_occured);

  const char zero = 0;
  test_internal::check(ok, LIBXSTREAM_SETVAL(zero), m_host_mem, LIBXSTREAM_SETVAL(size));
  LIBXSTREAM_CHECK_CONDITION_THROW(LIBXSTREAM_FALSE != ok);
}


test_type::~test_type()
{
  int device = -1;
  LIBXSTREAM_CHECK_CALL_RETURN(libxstream_stream_device(m_stream, &device));
  LIBXSTREAM_CHECK_CALL_RETURN(libxstream_mem_deallocate(-1, m_host_mem));
  LIBXSTREAM_CHECK_CALL_RETURN(libxstream_mem_deallocate(device, m_dev_mem1));
  LIBXSTREAM_CHECK_CALL_RETURN(libxstream_mem_deallocate(device, m_dev_mem2));
  LIBXSTREAM_CHECK_CALL_RETURN(libxstream_fn_destroy_signature(m_signature));
  LIBXSTREAM_CHECK_CALL_RETURN(libxstream_stream_destroy(m_stream));
  LIBXSTREAM_CHECK_CALL_RETURN(libxstream_event_destroy(m_event));
  fprintf(stdout, "TST successfully completed.\n");
}


int main(int argc, char* argv[])
{
  try {
#if defined(_OPENMP)
    const int ntasks = std::max(1 < argc ? std::atoi(argv[1]) : omp_get_max_threads(), 1);
#else
    const int ntasks = std::max(1 < argc ? std::atoi(argv[1]) : 1, 1);
#endif

    size_t ndevices = 0;
    if (LIBXSTREAM_ERROR_NONE != libxstream_get_ndevices(&ndevices) || 0 == ndevices) {
      throw std::runtime_error("no device found!");
    }

#if defined(_OPENMP)
    // chunksize: limit memory consumption on high core count systems
    const int chunksize = std::max(ntasks / LIBXSTREAM_MAX_NDEVICES, 1);
#   pragma omp parallel for schedule(dynamic,chunksize)
#endif
    for (int i = 0; i < ntasks; ++i) {
      const test_type test(i % ndevices);
    }
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
