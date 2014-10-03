/******************************************************************************
** Copyright (c) 2013-2014, Intel Corporation                                **
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
/* Christopher Dahnken (Intel Corp.), Hans Pabst (Intel Corp.),
 * Alfio Lazzaro (CRAY Inc.), and Gilles Fourestey (CSCS)
******************************************************************************/
#include <immintrin.h>
#include <stdio.h>

#ifdef __MIC__
inline __m512d _MM512_LOADU_PD(const double* a) {
  __m512d va= _mm512_setzero_pd();
  va=_mm512_loadunpacklo_pd(va, &a[0]);
  va=_mm512_loadunpackhi_pd(va, &a[8]);
  return va;
}

inline void _MM512_STOREU_PD(double* a,__m512d v) {
  _mm512_packstorelo_pd(&a[0], v);
  _mm512_packstorehi_pd(&a[8], v);
}

inline __m512d _MM512_MASK_LOADU_PD(const double* a, char mask) {
  __m512d va= _mm512_setzero_pd();
  va=_mm512_mask_loadunpacklo_pd(va, mask, &a[0]);
  va=_mm512_mask_loadunpackhi_pd(va, mask, &a[8]);
  return va;
}

inline void _MM512_MASK_STOREU_PD(double* a,__m512d v, char mask) {
  _mm512_mask_packstorelo_pd(&a[0], mask, v);
  _mm512_mask_packstorehi_pd(&a[8], mask, v);
}
#endif
