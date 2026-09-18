!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Exact inversion-subgroup EBR and symmetry-indicator analysis (spinful TR).
!> \note These are symmetry signatures, not a complete topological classification.
! **************************************************************************************************
MODULE topology_tqc
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: inversion_indicators, inversion_ebr_signature
CONTAINS

! **************************************************************************************************
!> \brief Fu--Kane indices and inversion indicator with an explicit sign convention.
!> \param odd Number of odd Kramers pairs at each TRIM, x-bit fastest
!> \param pairs Number of occupied Kramers pairs
!> \param DIMENSION Two (kz=0 plane) or three (full reciprocal torus)
!> \param strong Plane index in 2D, strong index in 3D
!> \param weak Three weak indices in 3D; -1 (unavailable) in 2D
!> \param z4 sum(odd) mod 4 in 3D; -1 (unavailable) in 2D
!> \param status Zero on valid input
! **************************************************************************************************
   SUBROUTINE inversion_indicators(odd, pairs, DIMENSION, strong, weak, z4, status)
      INTEGER, INTENT(IN)                                :: odd(:), pairs, DIMENSION
      INTEGER, INTENT(OUT)                               :: strong, weak(3), z4, status

      INTEGER                                            :: axis, k

      status = -1
      strong = -1
      weak = -1
      z4 = -1
      IF (DIMENSION < 2 .OR. DIMENSION > 3) RETURN
      IF (SIZE(odd) /= 2**DIMENSION .OR. pairs < 1 .OR. pairs > 1000000) RETURN
      IF (ANY(odd < 0) .OR. ANY(odd > pairs)) RETURN
      strong = MOD(SUM(odd), 2)
      IF (DIMENSION == 3) THEN
         z4 = MOD(SUM(odd), 4)
         weak = 0
         DO axis = 1, 3
            DO k = 0, 7
               IF (BTEST(k, axis - 1)) weak(axis) = weak(axis) + odd(k + 1)
            END DO
         END DO
         weak = MOD(weak, 2)
      END IF
      status = 0
   END SUBROUTINE inversion_indicators

! **************************************************************************************************
!> \brief Exact signed/nonnegative EBR-signature test for the inversion subgroup.
!> \param odd Odd Kramers pairs, x-bit fastest, 4 or 8 TRIM
!> \param pairs Selected number of Kramers pairs
!> \param signed_atomic Signature belongs to the integer lattice of atomic signatures
!> \param nonnegative_atomic A physical nonnegative EBR-signature decomposition exists
!> \param coefficients Even/odd local-parity EBR coefficients at the 2^d inversion centers
!> \param status Zero on valid input, coefficients valid only if nonnegative_atomic
!> \note EBR at r=a/2 has inversion character s*(-1)^(a.k). Fourier inversion
!>       determines the differences between even/odd multiplicities exactly.
!>       No database or floating-point integer decision is used.
! **************************************************************************************************
   SUBROUTINE inversion_ebr_signature(odd, pairs, signed_atomic, nonnegative_atomic, coefficients, status)
      INTEGER, INTENT(IN)                                :: odd(:), pairs
      LOGICAL, INTENT(OUT)                               :: signed_atomic, nonnegative_atomic
      INTEGER, INTENT(OUT)                               :: coefficients(:, :), status

      INTEGER                                            :: a, k, n, numerator, spare, &
                                                            difference(SIZE(odd))

      status = -1
      signed_atomic = .FALSE.
      nonnegative_atomic = .FALSE.
      coefficients = 0
      n = SIZE(odd)
      IF (n /= 4 .AND. n /= 8) RETURN
      IF (SIZE(coefficients, 1) /= 2 .OR. SIZE(coefficients, 2) /= n) RETURN
      IF (pairs < 1 .OR. pairs > 1000000 .OR. ANY(odd < 0) .OR. ANY(odd > pairs)) RETURN
      status = 0
      DO a = 0, n - 1
         numerator = 0
         DO k = 0, n - 1
            numerator = numerator + (-1)**POPCNT(IAND(a, k))*(pairs - 2*odd(k + 1))
         END DO
         IF (MOD(numerator, n) /= 0) RETURN
         difference(a + 1) = numerator/n
      END DO
      IF (MOD(pairs - SUM(difference), 2) /= 0) RETURN
      signed_atomic = .TRUE.
      spare = pairs - SUM(ABS(difference))
      IF (spare < 0) RETURN
      nonnegative_atomic = .TRUE.
      coefficients(1, :) = MAX(difference, 0)
      coefficients(2, :) = MAX(-difference, 0)
      coefficients(:, 1) = coefficients(:, 1) + spare/2
   END SUBROUTINE inversion_ebr_signature
END MODULE topology_tqc
