!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
PROGRAM topology_symmetry_unittest
   USE ieee_arithmetic,                 ONLY: ieee_quiet_nan,&
                                              ieee_value
   USE kinds,                           ONLY: dp
   USE topology_symmetry,               ONLY: character_multiplicities,&
                                              inversion_representation
   USE topology_tqc,                    ONLY: inversion_ebr_signature,&
                                              inversion_indicators

   IMPLICIT NONE
   COMPLEX(KIND=dp) :: table(6, 3), group_chars(6), metric(2, 2), coeff(4, 4), phase(2), tmp(4)
   REAL(KIND=dp) :: error, energy_error, energies(4)
   INTEGER :: multiplicities(3), status, counts(2), odd(8), ebr(2, 8), a, k, sign, strong, weak(3), z4
   LOGICAL :: signed_atomic, nonnegative_atomic

   ! Complete C3v/S3 character table per element, including a two-dimensional irrep.
   table(:, 1) = CMPLX([1, 1, 1, 1, 1, 1], 0, dp)
   table(:, 2) = CMPLX([1, 1, 1, -1, -1, -1], 0, dp)
   table(:, 3) = CMPLX([2, -1, -1, 0, 0, 0], 0, dp)
   group_chars = table(:, 1) + 2*table(:, 3)
   CALL character_multiplicities(group_chars, table, multiplicities, 1.e-10_dp, status)
   IF (status /= 0 .OR. ANY(multiplicities /= [1, 0, 2])) ERROR STOP 'C3v decomposition failed'
   group_chars(2) = group_chars(2) + 0.1_dp
   CALL character_multiplicities(group_chars, table, multiplicities, 1.e-10_dp, status)
   IF (status == 0) ERROR STOP 'Invalid character accepted'
   metric = CMPLX(0, 0, dp)
   metric(1, 1) = 2.0_dp
   metric(2, 2) = 0.5_dp
   coeff = CMPLX(0, 0, dp)
   coeff(1, 1) = SQRT(0.5_dp)
   coeff(2, 2) = SQRT(2.0_dp)
   coeff(3, 3) = SQRT(0.5_dp)
   coeff(4, 4) = SQRT(2.0_dp)
   phase = CMPLX([1, -1], 0, dp)
   energies = 0.0_dp
   CALL inversion_representation(metric, coeff, [1, 2], phase, energies, .TRUE., &
                                 1.e-9_dp, 1.e-9_dp, counts, error, energy_error, status)
   IF (status /= 0 .OR. ANY(counts /= [2, 2])) ERROR STOP 'Metric/spinor parity failed'
   tmp = coeff(:, 1)
   coeff(:, 1) = (tmp + CMPLX(0, 1, dp)*coeff(:, 2))/SQRT(2.0_dp)
   coeff(:, 2) = (CMPLX(0, 1, dp)*tmp + coeff(:, 2))/SQRT(2.0_dp)
   CALL inversion_representation(metric, coeff, [1, 2], phase, energies, .TRUE., &
                                 1.e-9_dp, 1.e-9_dp, counts, error, energy_error, status)
   IF (status /= 0 .OR. ANY(counts /= [2, 2])) ERROR STOP 'Degenerate-gauge parity failed'
   CALL inversion_representation(metric, coeff(:, :1), [1, 2], phase, energies(:1), .FALSE., &
                                 1.e-9_dp, 1.e-9_dp, counts, error, energy_error, status)
   IF (status == 0) ERROR STOP 'Nonclosed parity subspace accepted'
   CALL inversion_representation(metric, coeff(:, :2), [1, 2], phase, energies(:2), .TRUE., &
                                 1.e-9_dp, 1.e-9_dp, counts, error, energy_error, status)
   IF (status == 0) ERROR STOP 'Broken physical time reversal accepted'
   energies(2) = 1.0_dp
   CALL inversion_representation(metric, coeff, [1, 2], phase, energies, .TRUE., &
                                 1.e-9_dp, 1.e-9_dp, counts, error, energy_error, status)
   IF (status == 0) ERROR STOP 'Noncommuting eigenspace accepted'
   ! Scalar two-site inversion with a genuinely off-diagonal AO metric.
   metric = CMPLX(0.2_dp, 0, dp)
   metric(1, 1) = 1.0_dp
   metric(2, 2) = 1.0_dp
   coeff(:2, 1) = CMPLX([1.0_dp, 1.0_dp]/SQRT(2.4_dp), 0, dp)
   coeff(:2, 2) = CMPLX([1.0_dp, -1.0_dp]/SQRT(1.6_dp), 0, dp)
   phase = CMPLX(1, 0, dp)
   CALL inversion_representation(metric, coeff(:2, :2), [2, 1], phase, energies(:2), .FALSE., &
                                 1.e-9_dp, 1.e-9_dp, counts, error, energy_error, status)
   IF (status /= 0 .OR. ANY(counts /= [1, 1])) ERROR STOP 'Scalar two-site inversion failed'
   CALL inversion_representation(metric, coeff(:2, :1), [2, 1], -phase, energies(:1), .FALSE., &
                                 1.e-9_dp, 1.e-9_dp, counts, error, energy_error, status)
   IF (status /= 0 .OR. ANY(counts /= [0, 1])) ERROR STOP 'Zone-boundary inversion phase failed'
   metric(1, 2) = CMPLX(0.2_dp, 0.01_dp, dp)
   CALL inversion_representation(metric, coeff(:2, :2), [2, 1], phase, energies(:2), .FALSE., &
                                 1.e-9_dp, 1.e-9_dp, counts, error, energy_error, status)
   IF (status == 0) ERROR STOP 'Invalid metric accepted'
   metric(1, 2) = CMPLX(ieee_value(0.0_dp, ieee_quiet_nan), 0, dp)
   CALL inversion_representation(metric, coeff(:2, :2), [2, 1], phase, energies(:2), .FALSE., &
                                 1.e-9_dp, 1.e-9_dp, counts, error, energy_error, status)
   IF (status == 0) ERROR STOP 'Nonfinite metric accepted'
   DO a = 0, 7
      DO sign = -1, 1, 2
         DO k = 0, 7
            odd(k + 1) = (1 - sign*(-1)**POPCNT(IAND(a, k)))/2
         END DO
         CALL inversion_indicators(odd, 1, 3, strong, weak, z4, status)
         IF (status /= 0 .OR. strong /= 0 .OR. ANY(weak /= 0) .OR. z4 /= 0) ERROR STOP 'Atomic indicator failed'
         CALL inversion_ebr_signature(odd, 1, signed_atomic, nonnegative_atomic, ebr, status)
         IF (status /= 0 .OR. .NOT. nonnegative_atomic .OR. .NOT. signed_atomic) ERROR STOP 'Atomic EBR failed'
         IF (SUM(ebr) /= 1) ERROR STOP 'EBR rank failed'
      END DO
   END DO
   odd = 0
   odd(1) = 1
   CALL inversion_indicators(odd, 1, 3, strong, weak, z4, status)
   IF (status /= 0 .OR. strong /= 1 .OR. z4 /= 1) ERROR STOP 'Strong index failed'
   CALL inversion_ebr_signature(odd, 1, signed_atomic, nonnegative_atomic, ebr, status)
   IF (status /= 0 .OR. signed_atomic .OR. nonnegative_atomic) ERROR STOP 'Stable obstruction missed'
   odd(1) = 2
   CALL inversion_indicators(odd, 2, 3, strong, weak, z4, status)
   IF (status /= 0 .OR. strong /= 0 .OR. z4 /= 2) ERROR STOP 'Inversion Z4=2 missed'
   odd(1) = 4
   CALL inversion_ebr_signature(odd, 4, signed_atomic, nonnegative_atomic, ebr, status)
   IF (status /= 0 .OR. .NOT. signed_atomic .OR. nonnegative_atomic) ERROR STOP 'Signed-only EBR missed'
   odd = -1
   CALL inversion_ebr_signature(odd, 4, signed_atomic, nonnegative_atomic, ebr, status)
   IF (status == 0) ERROR STOP 'Invalid multiplicity accepted'
   CALL exhaustive_signatures(2)
   CALL exhaustive_signatures(3)
   PRINT *, 'Native characters, parity, inversion indicators and EBR tests passed.'
CONTAINS

! **************************************************************************************************
!> \brief Compare the exact solver against independent enumeration of rank-two atomic sums.
!> \param DIMENSION Dimension of the inversion subgroup
! **************************************************************************************************
   SUBROUTINE exhaustive_signatures(DIMENSION)
      INTEGER, INTENT(IN)                                :: dimension

      INTEGER :: code, digit, i, ierr, indicator, indices(3), j, n, nu, x, y, &
         values(2**DIMENSION), signatures(2**DIMENSION, 2**(DIMENSION + 1)), &
         reconstructed(2**DIMENSION), coefficients(2, 2**DIMENSION)
      LOGICAL :: indicators_zero, nonnegative_ok, possible(0:3**(2**DIMENSION) - 1), signed_ok

      n = 2**DIMENSION
      DO x = 0, n - 1
         DO i = 0, n - 1
            signatures(i + 1, 2*x + 1) = (1 - (-1)**POPCNT(IAND(x, i)))/2
            signatures(i + 1, 2*x + 2) = 1 - signatures(i + 1, 2*x + 1)
         END DO
      END DO
      possible = .FALSE.
      DO x = 1, 2*n
         DO y = 1, 2*n
            code = SUM((signatures(:, x) + signatures(:, y))*[(3**i, i=0, n - 1)])
            possible(code) = .TRUE.
         END DO
      END DO
      DO code = 0, 3**n - 1
         digit = code
         DO i = 1, n
            values(i) = MOD(digit, 3)
            digit = digit/3
         END DO
         CALL inversion_ebr_signature(values, 2, signed_ok, nonnegative_ok, coefficients, ierr)
         IF (ierr /= 0 .OR. (nonnegative_ok .NEQV. possible(code))) ERROR STOP 'Enumerated EBR mismatch'
         CALL inversion_indicators(values, 2, DIMENSION, nu, indices, indicator, ierr)
         IF (ierr /= 0) ERROR STOP 'Enumerated indicator failure'
         indicators_zero = nu == 0
         IF (DIMENSION == 3) indicators_zero = indicator == 0 .AND. ALL(indices == 0)
         IF (signed_ok .NEQV. indicators_zero) ERROR STOP 'Indicator/atomic-lattice mismatch'
         IF (.NOT. nonnegative_ok) CYCLE
         reconstructed = 0
         DO j = 1, n
            reconstructed = reconstructed + coefficients(1, j)*signatures(:, 2*j - 1) + &
                            coefficients(2, j)*signatures(:, 2*j)
         END DO
         IF (ANY(reconstructed /= values) .OR. SUM(coefficients) /= 2) ERROR STOP 'EBR reconstruction failed'
      END DO
   END SUBROUTINE exhaustive_signatures
END PROGRAM topology_symmetry_unittest
