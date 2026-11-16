!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Double-double arithmetic: real and complex numbers stored as unevaluated sums hi + lo of
!>        two double-precision numbers, about 31 significant digits, following Y. Hida, X. S. Li,
!>        D. H. Bailey, Proc. 15th IEEE Symposium on Computer Arithmetic, 155 (2001).
!> \note
!>      The error-free transformations rely on IEEE double-precision rounding of every single
!>      operation. Contraction to fused multiply-adds (e.g. gfortran with -march=native) leaves
!>      the results correct but reduces the precision of the product to about 77 bits;
!>      -ffp-contract=off restores the full precision.
!> \par History
!>      11.2026 moved out of gw_nevanlinna
! **************************************************************************************************
MODULE double_double
   USE kinds,                           ONLY: dp,&
                                              int_8

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'double_double'

   ! clears the lowest 27 of the 52 stored mantissa bits; splitting by bit masking instead of
   ! Dekker's multiplication keeps the split exact when the compiler contracts to FMA
   INTEGER(KIND=int_8), PARAMETER :: dd_split_mask = -134217728_int_8

   PUBLIC :: dd_type, zdd_type, zdd_zero, zdd_one, zdd_i
   PUBLIC :: dd_add, dd_sub, dd_mul, dd_div
   PUBLIC :: zdd, zdd_to_dp, zdd_add, zdd_sub, zdd_mul, zdd_div, zdd_conjg, zdd_abs2

! **************************************************************************************************
!> \brief Double-double number hi + lo with |lo| <= ulp(hi)/2
! **************************************************************************************************
   TYPE dd_type
      REAL(KIND=dp)                                      :: hi = 0.0_dp, lo = 0.0_dp
   END TYPE dd_type

! **************************************************************************************************
!> \brief Complex double-double number
! **************************************************************************************************
   TYPE zdd_type
      TYPE(dd_type)                                      :: re = dd_type(), im = dd_type()
   END TYPE zdd_type

   TYPE(zdd_type), PARAMETER :: zdd_zero = zdd_type(dd_type(0.0_dp, 0.0_dp), dd_type(0.0_dp, 0.0_dp)), &
                                zdd_one = zdd_type(dd_type(1.0_dp, 0.0_dp), dd_type(0.0_dp, 0.0_dp)), &
                                zdd_i = zdd_type(dd_type(0.0_dp, 0.0_dp), dd_type(1.0_dp, 0.0_dp))

CONTAINS

! **************************************************************************************************
!> \brief Error-free sum, a + b = s + e exactly
!> \param a summand
!> \param b summand
!> \param s rounded sum
!> \param e rounding error
! **************************************************************************************************
   ELEMENTAL SUBROUTINE two_sum(a, b, s, e)
      REAL(KIND=dp), INTENT(IN)                          :: a, b
      REAL(KIND=dp), INTENT(OUT)                         :: s, e

      REAL(KIND=dp)                                      :: bb

      s = a + b
      bb = s - a
      e = (a - (s - bb)) + (b - bb)

   END SUBROUTINE two_sum

! **************************************************************************************************
!> \brief Error-free sum for |a| >= |b|
!> \param a larger summand
!> \param b smaller summand
!> \param s rounded sum
!> \param e rounding error
! **************************************************************************************************
   ELEMENTAL SUBROUTINE quick_two_sum(a, b, s, e)
      REAL(KIND=dp), INTENT(IN)                          :: a, b
      REAL(KIND=dp), INTENT(OUT)                         :: s, e

      s = a + b
      e = b - (s - a)

   END SUBROUTINE quick_two_sum

! **************************************************************************************************
!> \brief Error-free product, a*b = p + e to within a relative 2**(-106)
!> \param a factor
!> \param b factor
!> \param p rounded product
!> \param e rounding error
! **************************************************************************************************
   ELEMENTAL SUBROUTINE two_prod(a, b, p, e)
      REAL(KIND=dp), INTENT(IN)                          :: a, b
      REAL(KIND=dp), INTENT(OUT)                         :: p, e

      REAL(KIND=dp)                                      :: a_hi, a_lo, b_hi, b_lo

      p = a*b
      a_hi = TRANSFER(IAND(TRANSFER(a, 0_int_8), dd_split_mask), 0.0_dp)
      a_lo = a - a_hi
      b_hi = TRANSFER(IAND(TRANSFER(b, 0_int_8), dd_split_mask), 0.0_dp)
      b_lo = b - b_hi
      e = ((a_hi*b_hi - p) + a_hi*b_lo + a_lo*b_hi) + a_lo*b_lo

   END SUBROUTINE two_prod

! **************************************************************************************************
!> \brief Sum of two double-double numbers
!> \param a summand
!> \param b summand
!> \return a + b
! **************************************************************************************************
   ELEMENTAL FUNCTION dd_add(a, b) RESULT(c)
      TYPE(dd_type), INTENT(IN)                          :: a, b
      TYPE(dd_type)                                      :: c

      REAL(KIND=dp)                                      :: e, e2, f, s, s2, t

      CALL two_sum(a%hi, b%hi, s, e)
      CALL two_sum(a%lo, b%lo, t, f)
      e = e + t
      CALL quick_two_sum(s, e, s2, e2)
      e2 = e2 + f
      CALL quick_two_sum(s2, e2, c%hi, c%lo)

   END FUNCTION dd_add

! **************************************************************************************************
!> \brief Difference of two double-double numbers
!> \param a minuend
!> \param b subtrahend
!> \return a - b
! **************************************************************************************************
   ELEMENTAL FUNCTION dd_sub(a, b) RESULT(c)
      TYPE(dd_type), INTENT(IN)                          :: a, b
      TYPE(dd_type)                                      :: c

      c = dd_add(a, dd_type(-b%hi, -b%lo))

   END FUNCTION dd_sub

! **************************************************************************************************
!> \brief Product of two double-double numbers
!> \param a factor
!> \param b factor
!> \return a*b
! **************************************************************************************************
   ELEMENTAL FUNCTION dd_mul(a, b) RESULT(c)
      TYPE(dd_type), INTENT(IN)                          :: a, b
      TYPE(dd_type)                                      :: c

      REAL(KIND=dp)                                      :: e, p

      CALL two_prod(a%hi, b%hi, p, e)
      e = e + (a%hi*b%lo + a%lo*b%hi)
      CALL quick_two_sum(p, e, c%hi, c%lo)

   END FUNCTION dd_mul

! **************************************************************************************************
!> \brief Quotient of two double-double numbers, by three steps of long division
!> \param a dividend
!> \param b divisor
!> \return a/b
! **************************************************************************************************
   ELEMENTAL FUNCTION dd_div(a, b) RESULT(c)
      TYPE(dd_type), INTENT(IN)                          :: a, b
      TYPE(dd_type)                                      :: c

      REAL(KIND=dp)                                      :: q1, q2, q3
      TYPE(dd_type)                                      :: q, r

      q1 = a%hi/b%hi
      r = dd_sub(a, dd_mul(dd_type(q1, 0.0_dp), b))
      q2 = r%hi/b%hi
      r = dd_sub(r, dd_mul(dd_type(q2, 0.0_dp), b))
      q3 = r%hi/b%hi
      CALL quick_two_sum(q1, q2, q%hi, q%lo)
      c = dd_add(q, dd_type(q3, 0.0_dp))

   END FUNCTION dd_div

! **************************************************************************************************
!> \brief Complex double-double number from a double-precision one
!> \param z complex number
!> \return z
! **************************************************************************************************
   ELEMENTAL FUNCTION zdd(z) RESULT(c)
      COMPLEX(KIND=dp), INTENT(IN)                       :: z
      TYPE(zdd_type)                                     :: c

      c%re = dd_type(REAL(z, KIND=dp), 0.0_dp)
      c%im = dd_type(AIMAG(z), 0.0_dp)

   END FUNCTION zdd

! **************************************************************************************************
!> \brief Rounds a complex double-double number to double precision
!> \param a complex double-double number
!> \return a in double precision
! **************************************************************************************************
   ELEMENTAL FUNCTION zdd_to_dp(a) RESULT(z)
      TYPE(zdd_type), INTENT(IN)                         :: a
      COMPLEX(KIND=dp)                                   :: z

      z = CMPLX(a%re%hi + a%re%lo, a%im%hi + a%im%lo, KIND=dp)

   END FUNCTION zdd_to_dp

! **************************************************************************************************
!> \brief Sum of two complex double-double numbers
!> \param a summand
!> \param b summand
!> \return a + b
! **************************************************************************************************
   ELEMENTAL FUNCTION zdd_add(a, b) RESULT(c)
      TYPE(zdd_type), INTENT(IN)                         :: a, b
      TYPE(zdd_type)                                     :: c

      c%re = dd_add(a%re, b%re)
      c%im = dd_add(a%im, b%im)

   END FUNCTION zdd_add

! **************************************************************************************************
!> \brief Difference of two complex double-double numbers
!> \param a minuend
!> \param b subtrahend
!> \return a - b
! **************************************************************************************************
   ELEMENTAL FUNCTION zdd_sub(a, b) RESULT(c)
      TYPE(zdd_type), INTENT(IN)                         :: a, b
      TYPE(zdd_type)                                     :: c

      c%re = dd_sub(a%re, b%re)
      c%im = dd_sub(a%im, b%im)

   END FUNCTION zdd_sub

! **************************************************************************************************
!> \brief Product of two complex double-double numbers
!> \param a factor
!> \param b factor
!> \return a*b
! **************************************************************************************************
   ELEMENTAL FUNCTION zdd_mul(a, b) RESULT(c)
      TYPE(zdd_type), INTENT(IN)                         :: a, b
      TYPE(zdd_type)                                     :: c

      c%re = dd_sub(dd_mul(a%re, b%re), dd_mul(a%im, b%im))
      c%im = dd_add(dd_mul(a%re, b%im), dd_mul(a%im, b%re))

   END FUNCTION zdd_mul

! **************************************************************************************************
!> \brief Quotient of two complex double-double numbers
!> \param a dividend
!> \param b divisor
!> \return a/b
! **************************************************************************************************
   ELEMENTAL FUNCTION zdd_div(a, b) RESULT(c)
      TYPE(zdd_type), INTENT(IN)                         :: a, b
      TYPE(zdd_type)                                     :: c

      TYPE(dd_type)                                      :: d

      d = dd_add(dd_mul(b%re, b%re), dd_mul(b%im, b%im))
      c%re = dd_div(dd_add(dd_mul(a%re, b%re), dd_mul(a%im, b%im)), d)
      c%im = dd_div(dd_sub(dd_mul(a%im, b%re), dd_mul(a%re, b%im)), d)

   END FUNCTION zdd_div

! **************************************************************************************************
!> \brief Complex conjugate of a complex double-double number
!> \param a complex double-double number
!> \return conjugate of a
! **************************************************************************************************
   ELEMENTAL FUNCTION zdd_conjg(a) RESULT(c)
      TYPE(zdd_type), INTENT(IN)                         :: a
      TYPE(zdd_type)                                     :: c

      c%re = a%re
      c%im = dd_type(-a%im%hi, -a%im%lo)

   END FUNCTION zdd_conjg

! **************************************************************************************************
!> \brief Squared modulus of a complex double-double number, rounded to double precision
!> \param a complex double-double number
!> \return |a|**2
! **************************************************************************************************
   ELEMENTAL FUNCTION zdd_abs2(a) RESULT(r)
      TYPE(zdd_type), INTENT(IN)                         :: a
      REAL(KIND=dp)                                      :: r

      TYPE(dd_type)                                      :: s

      s = dd_add(dd_mul(a%re, a%re), dd_mul(a%im, a%im))
      r = s%hi + s%lo

   END FUNCTION zdd_abs2

END MODULE double_double
