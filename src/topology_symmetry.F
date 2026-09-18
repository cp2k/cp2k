!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Native character decomposition and metric-aware inversion representations.
! **************************************************************************************************
MODULE topology_symmetry
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: character_multiplicities, inversion_representation
CONTAINS

! **************************************************************************************************
!> \brief Decompose a character against a complete unitary character table.
!> \param CHARACTER One character per group element (not per conjugacy class)
!> \param table Characters, group element first and irrep second; identity first
!> \param multiplicities Nonnegative integer multiplicities, invalid unless status=0
!> \param tolerance Absolute numerical tolerance
!> \param status Zero on success, negative for invalid/incomplete/nonintegral data
!> \note All projective irreps must use the same factor system. This kernel does not
!>       infer a space group or certify the supplied character table's provenance.
! **************************************************************************************************
   SUBROUTINE character_multiplicities(CHARACTER, table, multiplicities, tolerance, status)
      COMPLEX(KIND=dp), INTENT(IN)                       :: character(:), table(:, :)
      INTEGER, INTENT(OUT)                               :: multiplicities(:)
      REAL(KIND=dp), INTENT(IN)                          :: tolerance
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(KIND=dp)                                   :: value
      INTEGER                                            :: dimension, i, j, nirrep, order, total

      status = -1
      multiplicities = -1
      order = SIZE(CHARACTER)
      nirrep = SIZE(table, 2)
      IF (order < 1 .OR. order > 1024 .OR. nirrep < 1 .OR. SIZE(table, 1) /= order) RETURN
      IF (SIZE(multiplicities) /= nirrep .OR. tolerance <= 0.0_dp) RETURN
      IF (.NOT. ieee_is_finite(tolerance)) RETURN
      IF (.NOT. ALL(ieee_is_finite(REAL(CHARACTER, dp)))) RETURN
      IF (.NOT. ALL(ieee_is_finite(AIMAG(CHARACTER)))) RETURN
      IF (.NOT. ALL(ieee_is_finite(REAL(table, dp)))) RETURN
      IF (.NOT. ALL(ieee_is_finite(AIMAG(table)))) RETURN
      IF (MAXVAL(ABS(table)) > REAL(order, dp) + tolerance) RETURN
      IF (MAXVAL(ABS(CHARACTER)) > 1.e6_dp) RETURN
      total = 0
      DO i = 1, nirrep
         DIMENSION = NINT(REAL(table(1, i), dp))
         IF (DIMENSION < 1 .OR. DIMENSION > order) RETURN
         IF (ABS(table(1, i) - DIMENSION) > tolerance) RETURN
         total = total + DIMENSION**2
         IF (total > order) RETURN
         DO j = 1, nirrep
            value = SUM(CONJG(table(:, i))*table(:, j))/REAL(order, dp)
            IF (i == j) value = value - 1.0_dp
            IF (ABS(value) > tolerance) RETURN
         END DO
      END DO
      IF (total /= order) RETURN
      status = -2
      DO i = 1, nirrep
         value = SUM(CONJG(table(:, i))*CHARACTER)/REAL(order, dp)
         j = NINT(REAL(value, dp))
         IF (j < 0 .OR. ABS(value - j) > tolerance) RETURN
         multiplicities(i) = j
      END DO
      IF (MAXVAL(ABS(MATMUL(table, CMPLX(multiplicities, 0, dp)) - CHARACTER)) > tolerance) RETURN
      status = 0
   END SUBROUTINE character_multiplicities

! **************************************************************************************************
!> \brief Inversion irreps of an isolated scalar or spinor eigenspace in an AO metric.
!> \param metric Scalar AO overlap at a TRIM
!> \param coeff Selected coefficients, spinor components stacked by AO
!> \param mapping Target AO index under inversion
!> \param phase Source-AO multiplier, including orbital parity and lattice phase
!> \param energies Selected eigenvalues in Hartree
!> \param check_tr Check physical spinful time reversal as well as inversion
!> \param tolerance Dimensionless metric/subspace/character tolerance
!> \param energy_tolerance Eigenvalue-commutator tolerance in Hartree
!> \param counts Multiplicities of even/odd inversion irreps, counting states not pairs
!> \param error Largest dimensionless residual
!> \param energy_error Largest projected symmetry/eigenvalue commutator
!> \param status Zero on success; invalid output for nonzero status
!> \param diagnostics Optional metric, normalization, inversion and time-reversal residuals
! **************************************************************************************************
   SUBROUTINE inversion_representation(metric, coeff, mapping, phase, energies, check_tr, &
                                       tolerance, energy_tolerance, counts, error, energy_error, status, diagnostics)
      COMPLEX(KIND=dp), INTENT(IN)                       :: metric(:, :), coeff(:, :)
      INTEGER, INTENT(IN)                                :: mapping(:)
      COMPLEX(KIND=dp), INTENT(IN)                       :: phase(:)
      REAL(KIND=dp), INTENT(IN)                          :: energies(:)
      LOGICAL, INTENT(IN)                                :: check_tr
      REAL(KIND=dp), INTENT(IN)                          :: tolerance, energy_tolerance
      INTEGER, INTENT(OUT)                               :: counts(2)
      REAL(KIND=dp), INTENT(OUT)                         :: error, energy_error
      INTEGER, INTENT(OUT)                               :: status
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: diagnostics(4)

      COMPLEX(KIND=dp)                                   :: characters(2), table(2, 2)
      COMPLEX(KIND=dp), ALLOCATABLE                      :: pc(:, :), projected(:, :), sc(:, :), &
                                                            trc(:, :), work(:, :)
      INTEGER                                            :: first, i, j, last, nao, nspin, rank, s
      REAL(KIND=dp)                                      :: metric_scale

      status = -1
      counts = -1
      error = HUGE(1.0_dp)
      energy_error = HUGE(1.0_dp)
      IF (PRESENT(diagnostics)) diagnostics = HUGE(1.0_dp)
      nao = SIZE(metric, 1)
      rank = SIZE(coeff, 2)
      IF (nao < 1 .OR. rank < 1 .OR. SIZE(metric, 2) /= nao) RETURN
      IF (MOD(SIZE(coeff, 1), nao) /= 0) RETURN
      nspin = SIZE(coeff, 1)/nao
      IF (nspin < 1 .OR. nspin > 2) RETURN
      IF (check_tr .AND. (nspin /= 2 .OR. MOD(rank, 2) /= 0)) RETURN
      IF (SIZE(mapping) /= nao .OR. SIZE(phase) /= nao .OR. SIZE(energies) /= rank) RETURN
      IF (tolerance <= 0.0_dp .OR. energy_tolerance <= 0.0_dp) RETURN
      IF (.NOT. ieee_is_finite(tolerance) .OR. .NOT. ieee_is_finite(energy_tolerance)) RETURN
      IF (.NOT. ALL(ieee_is_finite(REAL(metric, dp))) .OR. .NOT. ALL(ieee_is_finite(AIMAG(metric)))) RETURN
      IF (.NOT. ALL(ieee_is_finite(REAL(coeff, dp))) .OR. .NOT. ALL(ieee_is_finite(AIMAG(coeff)))) RETURN
      IF (.NOT. ALL(ieee_is_finite(REAL(phase, dp))) .OR. .NOT. ALL(ieee_is_finite(AIMAG(phase)))) RETURN
      IF (.NOT. ALL(ieee_is_finite(energies))) RETURN
      IF (ANY(mapping < 1) .OR. ANY(mapping > nao)) RETURN
      IF (ANY(ABS(ABS(phase) - 1.0_dp) > tolerance)) RETURN
      DO i = 1, nao
         IF (mapping(mapping(i)) /= i) RETURN
         IF (ABS(phase(i)*phase(mapping(i)) - 1.0_dp) > tolerance) RETURN
      END DO
      metric_scale = MAX(1.0_dp, MAXVAL(ABS(metric)))
      error = MAXVAL(ABS(metric - CONJG(TRANSPOSE(metric))))/metric_scale
      DO j = 1, nao
         DO i = 1, nao
            error = MAX(error, ABS(CONJG(phase(i))*metric(mapping(i), mapping(j))*phase(j) - metric(i, j))/ &
                        metric_scale)
         END DO
      END DO
      IF (PRESENT(diagnostics)) diagnostics(1) = error
      ALLOCATE (sc(nao*nspin, rank), pc(nao*nspin, rank), projected(rank, rank), work(rank, rank))
      DO s = 1, nspin
         first = (s - 1)*nao
         last = s*nao
         sc(first + 1:last, :) = MATMUL(metric, coeff(first + 1:last, :))
         DO i = 1, nao
            pc(first + mapping(i), :) = phase(i)*coeff(first + i, :)
         END DO
      END DO
      work(:, :) = MATMUL(CONJG(TRANSPOSE(coeff)), sc)
      DO i = 1, rank
         work(i, i) = work(i, i) - 1.0_dp
      END DO
      error = MAX(error, MAXVAL(ABS(work)))
      IF (PRESENT(diagnostics)) diagnostics(2) = MAXVAL(ABS(work))
      projected(:, :) = MATMUL(CONJG(TRANSPOSE(sc)), pc)
      IF (.NOT. ALL(ieee_is_finite(REAL(projected, dp)))) RETURN
      IF (.NOT. ALL(ieee_is_finite(AIMAG(projected)))) RETURN
      work(:, :) = MATMUL(CONJG(TRANSPOSE(projected)), projected)
      DO i = 1, rank
         work(i, i) = work(i, i) - 1.0_dp
      END DO
      error = MAX(error, MAXVAL(ABS(work)), MAXVAL(ABS(projected - CONJG(TRANSPOSE(projected)))))
      IF (PRESENT(diagnostics)) diagnostics(3) = MAX(MAXVAL(ABS(work)), &
                                                     MAXVAL(ABS(projected - CONJG(TRANSPOSE(projected)))))
      IF (PRESENT(diagnostics)) THEN
         IF (.NOT. check_tr) diagnostics(4) = 0.0_dp
      END IF
      energy_error = 0.0_dp
      DO j = 1, rank
         DO i = 1, rank
            energy_error = MAX(energy_error, ABS(projected(i, j)*(energies(i) - energies(j))))
         END DO
      END DO
      characters(1) = CMPLX(rank, 0, dp)
      characters(2) = CMPLX(0, 0, dp)
      DO i = 1, rank
         characters(2) = characters(2) + projected(i, i)
      END DO
      table(:, 1) = CMPLX([1, 1], 0, dp)
      table(:, 2) = CMPLX([1, -1], 0, dp)
      CALL character_multiplicities(characters, table, counts, tolerance*rank, status)
      IF (status /= 0) RETURN
      status = -3
      IF (check_tr) THEN
         ! In a real Gaussian Bloch basis at a TRIM, Theta=(i sigma_y) K.
         error = MAX(error, MAXVAL(ABS(AIMAG(metric)))/MAX(1.0_dp, MAXVAL(ABS(metric))))
         ALLOCATE (trc(2*nao, rank))
         trc(:nao, :) = CONJG(coeff(nao + 1:, :))
         trc(nao + 1:, :) = -CONJG(coeff(:nao, :))
         projected(:, :) = MATMUL(CONJG(TRANSPOSE(sc)), trc)
         work(:, :) = MATMUL(CONJG(TRANSPOSE(projected)), projected)
         DO i = 1, rank
            work(i, i) = work(i, i) - 1.0_dp
         END DO
         error = MAX(error, MAXVAL(ABS(work)), MAXVAL(ABS(projected + TRANSPOSE(projected))))
         IF (PRESENT(diagnostics)) diagnostics(4) = MAX(MAXVAL(ABS(work)), &
                                                        MAXVAL(ABS(projected + TRANSPOSE(projected))))
         DO j = 1, rank
            DO i = 1, rank
               energy_error = MAX(energy_error, ABS(projected(i, j)*(energies(i) - energies(j))))
            END DO
         END DO
         IF (ANY(MOD(counts, 2) /= 0)) RETURN
      END IF
      IF (.NOT. ieee_is_finite(error) .OR. .NOT. ieee_is_finite(energy_error)) RETURN
      IF (error > tolerance .OR. energy_error > energy_tolerance) RETURN
      status = 0
   END SUBROUTINE inversion_representation
END MODULE topology_symmetry
