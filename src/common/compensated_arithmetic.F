!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!> \brief Two-component FP64 sums and products retaining cancellation before final rounding.
MODULE compensated_arithmetic
   USE iso_c_binding,                   ONLY: c_double
   USE kinds,                           ONLY: dp

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: twofold_real_type, twofold_add, twofold_multiply, twofold_product, twofold_divide, twofold_value
   TYPE twofold_real_type
      REAL(dp) :: hi = 0.0_dp, lo = 0.0_dp
   END TYPE twofold_real_type
   INTERFACE OPERATOR(+)
      MODULE PROCEDURE twofold_add
   END INTERFACE
   INTERFACE OPERATOR(-)
      MODULE PROCEDURE twofold_subtract
   END INTERFACE
   INTERFACE OPERATOR(*)
      MODULE PROCEDURE twofold_multiply
   END INTERFACE
   INTERFACE
      FUNCTION fused_product_sum(a, b, c) RESULT(value) BIND(C, NAME="fma")
         IMPORT :: c_double
         REAL(c_double), VALUE, INTENT(IN)::a, b, c
         REAL(c_double)::value
      END FUNCTION fused_product_sum
   END INTERFACE
CONTAINS
! **************************************************************************************************
!> \brief Error-free sum of two FP64 inputs, protected from reassociation across rounding steps.
!> \param a First input.
!> \param b Second input.
!> \param hi Rounded sum.
!> \param lo Its rounding remainder.
! **************************************************************************************************
   SUBROUTINE two_sum(a, b, hi, lo)
      REAL(dp), INTENT(IN)                               :: a, b
      REAL(dp), INTENT(OUT)                              :: hi, lo

      REAL(dp), VOLATILE                                 :: remainder_a, remainder_b, rounded, &
                                                            virtual_a, virtual_b

      rounded = a + b
      virtual_b = rounded - a
      virtual_a = rounded - virtual_b
      remainder_b = b - virtual_b
      remainder_a = a - virtual_a
      hi = rounded
      lo = remainder_a + remainder_b
   END SUBROUTINE two_sum

! **************************************************************************************************
!> \brief Product and its rounding remainder from the standard C fused multiply-add.
!> \param a First factor.
!> \param b Second factor.
!> \param hi Rounded product.
!> \param lo Its rounding remainder.
! **************************************************************************************************
   SUBROUTINE twofold_product(a, b, hi, lo)
      REAL(dp), INTENT(IN)                               :: a, b
      REAL(dp), INTENT(OUT)                              :: hi, lo

      REAL(dp), VOLATILE                                 :: rounded

      rounded = a*b
      hi = rounded
      lo = fused_product_sum(a, b, -rounded)
   END SUBROUTINE twofold_product

! **************************************************************************************************
!> \brief Add normalized two-component FP64 numbers.
!> \param a First input.
!> \param b Second input.
!> \return Compensated sum.
! **************************************************************************************************
   IMPURE ELEMENTAL FUNCTION twofold_add(a, b) RESULT(value)
      TYPE(twofold_real_type), INTENT(IN)                :: a, b
      TYPE(twofold_real_type)                            :: value

      REAL(dp)                                           :: hi, lo

      CALL two_sum(a%hi, b%hi, hi, lo)
      lo = lo + a%lo + b%lo
      CALL two_sum(hi, lo, value%hi, value%lo)
   END FUNCTION twofold_add

! **************************************************************************************************
!> \brief Multiply normalized two-component FP64 numbers.
!> \param a First input.
!> \param b Second input.
!> \return Compensated product.
! **************************************************************************************************
   IMPURE ELEMENTAL FUNCTION twofold_multiply(a, b) RESULT(value)
      TYPE(twofold_real_type), INTENT(IN)                :: a, b
      TYPE(twofold_real_type)                            :: value

      REAL(dp)                                           :: hi, lo

      CALL twofold_product(a%hi, b%hi, hi, lo)
      lo = lo + a%hi*b%lo + a%lo*b%hi + a%lo*b%lo
      CALL two_sum(hi, lo, value%hi, value%lo)
   END FUNCTION twofold_multiply

! **************************************************************************************************
!> \brief Difference retaining both input components.
!> \param a First input.
!> \param b Second input.
!> \return Compensated difference.
! **************************************************************************************************
   IMPURE ELEMENTAL FUNCTION twofold_subtract(a, b) RESULT(value)
      TYPE(twofold_real_type), INTENT(IN)                :: a, b
      TYPE(twofold_real_type)                            :: value

      value = twofold_add(a, twofold_real_type(-b%hi, -b%lo))
   END FUNCTION twofold_subtract

! **************************************************************************************************
!> \brief Quotient corrected with compensated multiplication remainders.
!> \param a Numerator.
!> \param b Nonzero normalized denominator.
!> \return Compensated quotient.
! **************************************************************************************************
   IMPURE ELEMENTAL FUNCTION twofold_divide(a, b) RESULT(value)
      TYPE(twofold_real_type), INTENT(IN)                :: a, b
      TYPE(twofold_real_type)                            :: value

      TYPE(twofold_real_type)                            :: remainder

      value = twofold_real_type(a%hi/b%hi, 0.0_dp)
      remainder = a - b*value
      value = value + twofold_real_type((remainder%hi + remainder%lo)/b%hi, 0.0_dp)
      remainder = a - b*value
      value = value + twofold_real_type((remainder%hi + remainder%lo)/b%hi, 0.0_dp)
   END FUNCTION twofold_divide

! **************************************************************************************************
!> \brief Round a normalized two-component value only at its consumer.
!> \param a Compensated value.
!> \return FP64 value.
! **************************************************************************************************
   PURE ELEMENTAL FUNCTION twofold_value(a) RESULT(value)
      TYPE(twofold_real_type), INTENT(IN)                :: a
      REAL(dp)                                           :: value

      value = a%hi + a%lo
   END FUNCTION twofold_value

END MODULE compensated_arithmetic
