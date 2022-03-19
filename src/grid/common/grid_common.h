/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#ifndef GRID_COMMON_H
#define GRID_COMMON_H

#define GRID_STRINGIFY(SYMBOL) #SYMBOL

// GCC added the simd pragma with version 6.
#if defined(__GNUC__) && !defined(__INTEL_COMPILER) && __GNUC__ < 6
#define GRID_PRAGMA_SIMD(OBJS, N)
// Intel added the simd pragma with version 19.00.
#elif defined(__INTEL_COMPILER) && __INTEL_COMPILER < 1900
#define GRID_PRAGMA_SIMD(OBJS, N)
// All compilers support the same syntax defined by the OpenMP standard.
#else
#define GRID_PRAGMA_SIMD(OBJS, N)                                              \
  _Pragma(GRID_STRINGIFY(omp simd linear OBJS simdlen(N)))
#endif

// GCC added the unroll pragma with version 8 and...
#if defined(__GNUC__) && !defined(__INTEL_COMPILER) && __GNUC__ < 8
#define GRID_PRAGMA_UNROLL(N)
#define GRID_PRAGMA_UNROLL_UP_TO(N)
// ...chose a custom syntax.
#elif defined(__GNUC__) && !defined(__INTEL_COMPILER) && __GNUC__ >= 8
#define GRID_PRAGMA_UNROLL(N) _Pragma(GRID_STRINGIFY(GCC unroll N))
#define GRID_PRAGMA_UNROLL_UP_TO(N) _Pragma(GRID_STRINGIFY(GCC unroll N))
// Most other compilers support a common syntax.
#else
#define GRID_PRAGMA_UNROLL(N) _Pragma(GRID_STRINGIFY(unroll(N)))
#define GRID_PRAGMA_UNROLL_UP_TO(N) _Pragma("unroll")
#endif

#if defined(__CUDACC__) || defined(__HIPCC__)
#define GRID_HOST_DEVICE __host__ __device__
#else
#define GRID_HOST_DEVICE
#endif

/*******************************************************************************
 * \brief Factorial function, e.g. fac(5) = 5! = 120.
 * \author Ole Schuett
 ******************************************************************************/
GRID_HOST_DEVICE static inline double fac(const int i) {
  static const double table[] = {
      0.10000000000000000000E+01, 0.10000000000000000000E+01,
      0.20000000000000000000E+01, 0.60000000000000000000E+01,
      0.24000000000000000000E+02, 0.12000000000000000000E+03,
      0.72000000000000000000E+03, 0.50400000000000000000E+04,
      0.40320000000000000000E+05, 0.36288000000000000000E+06,
      0.36288000000000000000E+07, 0.39916800000000000000E+08,
      0.47900160000000000000E+09, 0.62270208000000000000E+10,
      0.87178291200000000000E+11, 0.13076743680000000000E+13,
      0.20922789888000000000E+14, 0.35568742809600000000E+15,
      0.64023737057280000000E+16, 0.12164510040883200000E+18,
      0.24329020081766400000E+19, 0.51090942171709440000E+20,
      0.11240007277776076800E+22, 0.25852016738884976640E+23,
      0.62044840173323943936E+24, 0.15511210043330985984E+26,
      0.40329146112660563558E+27, 0.10888869450418352161E+29,
      0.30488834461171386050E+30, 0.88417619937397019545E+31,
      0.26525285981219105864E+33};
  return table[i];
}

/*******************************************************************************
 * \brief Number of Cartesian orbitals up to given angular momentum quantum.
 * \author Ole Schuett
 ******************************************************************************/
GRID_HOST_DEVICE static inline int ncoset(const int l) {
  static const int table[] = {1,  // l=0
                              4,  // l=1
                              10, // l=2 ...
                              20,  35,  56,  84,  120, 165, 220,  286,
                              364, 455, 560, 680, 816, 969, 1140, 1330};
  return table[l];
}

/*******************************************************************************
 * \brief Maps three angular momentum components to a single zero based index.
 * \author Ole Schuett
 ******************************************************************************/
GRID_HOST_DEVICE static inline int coset(int lx, int ly, int lz) {
  const int l = lx + ly + lz;
  if (l == 0) {
    return 0;
  } else {
    return ncoset(l - 1) + ((l - lx) * (l - lx + 1)) / 2 + lz;
  }
}

/*******************************************************************************
 * \brief Returns the smaller of two given integer (missing from the C standard)
 * \author Ole Schuett
 ******************************************************************************/
GRID_HOST_DEVICE static inline int imin(int x, int y) {
  return (x < y ? x : y);
}

/*******************************************************************************
 * \brief Returns the larger of two given integer (missing from the C standard)
 * \author Ole Schuett
 ******************************************************************************/
GRID_HOST_DEVICE static inline int imax(int x, int y) {
  return (x > y ? x : y);
}

/*******************************************************************************
 * \brief Equivalent of Fortran's MODULO, which always return a positive number.
 *        https://gcc.gnu.org/onlinedocs/gfortran/MODULO.html
 * \author Ole Schuett
 ******************************************************************************/
GRID_HOST_DEVICE static inline int modulo(int a, int m) {
  return ((a % m + m) % m);
}

/*******************************************************************************
 * \brief Orbital angular momentum.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  int l[3];
} orbital;

/*******************************************************************************
 * \brief Increase i'th component of given orbital angular momentum.
 * \author Ole Schuett
 ******************************************************************************/
GRID_HOST_DEVICE static inline orbital up(const int i, const orbital a) {
  orbital b = a;
  b.l[i] += 1;
  return b;
}

/*******************************************************************************
 * \brief Decrease i'th component of given orbital angular momentum.
 * \author Ole Schuett
 ******************************************************************************/
GRID_HOST_DEVICE static inline orbital down(const int i, const orbital a) {
  orbital b = a;
  b.l[i] = imax(0, a.l[i] - 1);
  return b;
}

/*******************************************************************************
 * \brief Return coset index of given orbital angular momentum.
 * \author Ole Schuett
 ******************************************************************************/
GRID_HOST_DEVICE static inline int idx(const orbital a) {
  return coset(a.l[0], a.l[1], a.l[2]);
}

#endif // GRID_COMMON_H

// EOF
