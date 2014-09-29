#include <immintrin.h>
#include <stdio.h>

#ifdef __MIC__
inline __m512d _MM512_LOADU_PD(double* a) {
  __m512d va= _mm512_setzero_pd();
  va=_mm512_loadunpacklo_pd(va, &a[0]);
  va=_mm512_loadunpackhi_pd(va, &a[8]);
  return va;
}

inline void _MM512_STOREU_PD(double* a,__m512d v) {
  _mm512_packstorelo_pd(&a[0], v);
  _mm512_packstorehi_pd(&a[8], v);
}

inline __m512d _MM512_MASK_LOADU_PD(double* a, char mask) {
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
