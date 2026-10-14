/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#ifndef FFT_UTILS_H
#define FFT_UTILS_H

#include "fft_timer.h"

#include <complex.h>
#include <omp.h>
#include <stdio.h>
#include <string.h>

#define FFT_FORWARD 1
#define FFT_BACKWARD -1

/*******************************************************************************
 * \brief Returns the smaller of two given integer (missing from the C standard)
 * \author Ole Schuett
 ******************************************************************************/
static inline int imin(int x, int y) { return (x < y ? x : y); }

/*******************************************************************************
 * \brief Returns the larger of two given integer (missing from the C standard)
 * \author Ole Schuett
 ******************************************************************************/
static inline int imax(int x, int y) { return (x > y ? x : y); }

/*******************************************************************************
 * \brief Returns the smaller of two given integer (missing from the C standard)
 * \author Frederick Stein
 ******************************************************************************/
static inline int dmin(double x, double y) { return (x < y ? x : y); }

/*******************************************************************************
 * \brief Returns the larger of two given integer (missing from the C standard)
 * \author Frederick Stein
 ******************************************************************************/
static inline int dmax(double x, double y) { return (x > y ? x : y); }

/*******************************************************************************
 * \brief Equivalent of Fortran's MODULO which always returns a positive number.
 *        https://gcc.gnu.org/onlinedocs/gfortran/MODULO.html
 * \author Ole Schuett
 ******************************************************************************/
static inline int modulo(int a, int m) { return ((a % m + m) % m); }

/*******************************************************************************
 * \brief Calculates the product of three numbers.
 * \author Frederick Stein
 ******************************************************************************/
static inline long int product3(const int array3[3]) {
  return (long int)array3[0] * array3[1] * array3[2];
}

#endif /* FFT_UTILS_H */

// EOF
