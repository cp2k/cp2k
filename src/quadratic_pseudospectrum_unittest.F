!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
PROGRAM quadratic_pseudospectrum_unittest
#if defined(__HAS_IEEE_EXCEPTIONS)
   USE ieee_exceptions, ONLY: IEEE_ALL, ieee_get_halting_mode, ieee_set_halting_mode
#endif
   USE kinds, ONLY: dp
   USE quadratic_pseudospectrum, ONLY: quadratic_dense, &
                                       quadratic_result_type

   IMPLICIT NONE
   COMPLEX(dp) :: op(2, 2, 4), s(2, 2), change(2, 2), rotated(2, 2, 4), metric(2, 2)
   COMPLEX(dp), ALLOCATABLE :: defect(:, :, :), overlap(:, :), projected(:, :)
   REAL(dp) :: shift(4), weights(4), expected, center
   TYPE(quadratic_result_type) :: dense, iterative, transformed
   INTEGER :: status, i, n

   s = CMPLX(0.0_dp, 0.0_dp, dp)
   s(1, 1) = 1.0_dp
   s(2, 2) = 1.0_dp
   op = CMPLX(0.0_dp, 0.0_dp, dp)
   op(1, 1, 1) = 1.0_dp
   op(2, 2, 1) = -1.0_dp
   op(1, 2, 2) = CMPLX(0.0_dp, 1.0_dp, dp)
   op(2, 1, 2) = CMPLX(0.0_dp, -1.0_dp, dp)
   shift = [0.3_dp, 0.2_dp, 0.0_dp, 0.0_dp]
   weights = [1.0_dp, 0.7_dp, 0.7_dp, 0.7_dp]
   expected = 1.0_dp + shift(1)**2 + weights(2)**2*(1.0_dp + shift(2)**2) &
              - 2.0_dp*SQRT(shift(1)**2 + weights(2)**4*shift(2)**2)
   CALL quadratic_dense(op, s, shift, weights, 1, .FALSE., 2, 200, 1.e-10_dp, 1.e-10_dp, dense, status)
   IF (status /= 0) ERROR STOP 'Quadratic Pauli oracle failed'
   IF (ABS(dense%eigenvalues(1) - expected) > 1.e-12_dp) ERROR STOP 'Analytic quadratic eigenvalue mismatch'
   CALL quadratic_dense(op, s, shift, weights, 1, .TRUE., 2, 200, 1.e-10_dp, 1.e-10_dp, iterative, status)
   IF (status /= 0) ERROR STOP 'Quadratic Pauli iteration failed'
   IF (MAXVAL(ABS(dense%eigenvalues - iterative%eigenvalues)) > 1.e-12_dp) ERROR STOP 'Pauli Ritz mismatch'
   change = RESHAPE([CMPLX(1.2_dp, 0.1_dp, dp), CMPLX(0.15_dp, -0.1_dp, dp), &
                     CMPLX(0.1_dp, 0.2_dp, dp), CMPLX(0.8_dp, 0.0_dp, dp)], [2, 2])
   metric = MATMUL(CONJG(TRANSPOSE(change)), change)
   DO i = 1, 4
      rotated(:, :, i) = MATMUL(CONJG(TRANSPOSE(change)), MATMUL(op(:, :, i), change))
   END DO
   CALL quadratic_dense(rotated, metric, shift, weights, 1, .FALSE., 2, 200, 1.e-10_dp, 1.e-10_dp, &
                        transformed, status)
   IF (status /= 0) ERROR STOP 'Nonorthogonal complex basis failed'
   IF (MAXVAL(ABS(dense%eigenvalues - transformed%eigenvalues)) > 1.e-12_dp) ERROR STOP 'Basis covariance failed'
   IF (MAXVAL(ABS(dense%errors - transformed%errors)) > 1.e-12_dp) ERROR STOP 'Observable residual covariance failed'
   IF (MAXVAL(ABS(dense%means - transformed%means)) > 1.e-12_dp) ERROR STOP 'Observable mean covariance failed'
   DO i = 1, 4
      op(:, :, i) = op(:, :, i) + REAL(i, dp)*s
      shift(i) = shift(i) + REAL(i, dp)
   END DO
   CALL quadratic_dense(op, s, shift, weights, 1, .FALSE., 2, 200, 1.e-10_dp, 1.e-10_dp, transformed, status)
   IF (status /= 0) ERROR STOP 'Shifted quadratic oracle failed'
   IF (MAXVAL(ABS(dense%eigenvalues - transformed%eigenvalues)) > 1.e-12_dp) ERROR STOP 'Origin invariance failed'
   IF (MAXVAL(ABS(dense%errors - transformed%errors)) > 1.e-12_dp) ERROR STOP 'Residual origin invariance failed'

   n = 40
   center = 20.0_dp
   ALLOCATE (defect(n, n, 4), overlap(n, n), projected(3, 3))
   defect = CMPLX(0.0_dp, 0.0_dp, dp)
   overlap = CMPLX(0.0_dp, 0.0_dp, dp)
   DO i = 1, n
      overlap(i, i) = 1.0_dp
      defect(i, i, 1) = 5.0_dp
      defect(i, i, 2) = REAL(i, dp) - center
      IF (i < n) THEN
         defect(i, i + 1, 1) = CMPLX(-0.1_dp, 0.02_dp, dp)
         defect(i + 1, i, 1) = CONJG(defect(i, i + 1, 1))
      END IF
   END DO
   defect(20, 20, 1) = 0.0_dp
   shift = 0.0_dp
   weights = [1.0_dp, 0.2_dp, 0.2_dp, 0.2_dp]
   CALL quadratic_dense(defect, overlap, shift, weights, 3, .FALSE., n, 500, 1.e-9_dp, 1.e-10_dp, dense, status)
   IF (status /= 0) ERROR STOP 'Defect oracle failed'
   IF (ABS(dense%states(20, 1))**2 < 0.95_dp) ERROR STOP 'Defect state not localized'
   CALL quadratic_dense(defect, overlap, shift, weights, 3, .TRUE., 20, 500, 1.e-9_dp, 1.e-10_dp, iterative, status)
   IF (status /= 0) ERROR STOP 'Defect iteration failed'
   IF (MAXVAL(ABS(dense%eigenvalues - iterative%eigenvalues)) > 1.e-8_dp) ERROR STOP 'Defect Ritz eigenvalue mismatch'
   projected(:, :) = MATMUL(CONJG(TRANSPOSE(dense%states)), iterative%states)
   IF (SUM(ABS(projected)**2) < 3.0_dp - 1.e-7_dp) ERROR STOP 'Defect invariant subspace mismatch'
   CALL quadratic_dense(defect, overlap, shift, weights, 3, .TRUE., 12, 1, 1.e-12_dp, 1.e-10_dp, iterative, status)
   IF (status /= 1) ERROR STOP 'Unconverged quadratic solve accepted'
   ! Exactly degenerate joint states: compare the projector, not individual eigenvectors.
   defect = CMPLX(0.0_dp, 0.0_dp, dp)
   DO i = 1, n
      defect(i, i, 1) = 2.0_dp
   END DO
   defect(19, 19, 1) = 0.0_dp
   defect(20, 20, 1) = 0.0_dp
   CALL quadratic_dense(defect, overlap, shift, weights, 2, .TRUE., 8, 200, 1.e-10_dp, 1.e-10_dp, iterative, status)
   IF (status /= 0) ERROR STOP 'Degenerate quadratic iteration failed'
   IF (MAXVAL(ABS(iterative%eigenvalues)) > 1.e-12_dp) ERROR STOP 'Degenerate eigenvalues failed'
   IF (ABS(SUM(ABS(iterative%states(19:20, :))**2) - 2.0_dp) > 1.e-12_dp) ERROR STOP 'Degenerate projector failed'
   overlap(1, 1) = -1.0_dp
   CALL quadratic_dense(defect, overlap, shift, weights, 3, .FALSE., n, 200, 1.e-9_dp, 1.e-10_dp, dense, status)
   IF (status /= -2) ERROR STOP 'Indefinite metric accepted'
   overlap(1, 1) = 1.e-14_dp
   CALL quadratic_dense(defect, overlap, shift, weights, 3, .FALSE., n, 200, 1.e-9_dp, 1.e-10_dp, dense, status)
   IF (status /= -2) ERROR STOP 'Ill-conditioned metric accepted'
   CALL periodic_tests()
   CALL translation_tests()
   CALL folded_translation_tests()
   CALL compressed_translation_tests()
CONTAINS
! **************************************************************************************************
!> \brief Compression of a three-state unitary to two AOs retains the missing normal Gram exactly.
! **************************************************************************************************
   SUBROUTINE compressed_translation_tests()
      COMPLEX(dp)                                        :: change(2, 2), changed(2, 2, 2), &
                                                            metric(2, 2), op(2, 2, 2), query(2), &
                                                            transformed(2, 2)
      INTEGER                                            :: i, status
      REAL(dp)                                           :: center, diagonal(2), err, expected(2), &
                                                            split, weights(2)
      TYPE(quadratic_result_type)                        :: full, projected, rotated

      metric = 0
      metric(1, 1) = 1
      metric(2, 2) = 1
      op = 0
      op(1, 1, 1) = 0.2_dp
      op(2, 2, 1) = -0.4_dp
      op(1, 2, 2) = 1
      query = [CMPLX(0.1_dp, 0.0_dp, dp), EXP(CMPLX(0.0_dp, 0.73_dp, dp))]
      weights = [1.0_dp, 0.7_dp]
      diagonal = [0.1_dp**2, 0.5_dp**2] + 2*weights(2)**2
      center = SUM(diagonal)/2
      split = SQRT(((diagonal(1) - diagonal(2))/2)**2 + weights(2)**4)
      expected = [center - split, center + split]
      CALL quadratic_dense(op, metric, query, weights, 2, .FALSE., 2, 100, 1.e-11_dp, 1.e-10_dp, &
                           full, status, unitary=[.FALSE., .TRUE.])
      IF (status /= 0) ERROR STOP 'Compressed unitary solve failed'
      err = MAXVAL(ABS(full%eigenvalues - expected))
      IF (err > 1.e-12_dp) ERROR STOP 'Full translation normal Gram differs from analytic result'
      IF (MAXVAL(ABS(full%errors**2 - full%projected_errors**2 - full%leakage**2)) > 1.e-12_dp) THEN
         ERROR STOP 'Continuum/projected residual decomposition failed'
      END IF
      IF (MAXVAL(ABS(full%leakage(2, :)**2 - ABS(full%states(1, :))**2)) > 1.e-12_dp) THEN
         ERROR STOP 'Translation leakage differs from discarded ambient component'
      END IF
      CALL quadratic_dense(op, metric, query, weights, 2, .FALSE., 2, 100, 1.e-11_dp, 1.e-10_dp, projected, status)
      IF (status /= 0) ERROR STOP 'Projected unitary comparison failed'
      IF (MAXVAL(ABS(projected%eigenvalues - full%eigenvalues)) < 0.1_dp) ERROR STOP 'Projection negative control failed'
      change = RESHAPE([CMPLX(1.2_dp, 0.1_dp, dp), CMPLX(0.15_dp, -0.1_dp, dp), &
                        CMPLX(0.1_dp, 0.2_dp, dp), CMPLX(0.8_dp, 0.0_dp, dp)], [2, 2])
      transformed = MATMUL(CONJG(TRANSPOSE(change)), change)
      DO i = 1, 2
         changed(:, :, i) = MATMUL(CONJG(TRANSPOSE(change)), MATMUL(op(:, :, i), change))
      END DO
      CALL quadratic_dense(changed, transformed, query, weights, 2, .TRUE., 2, 100, 1.e-11_dp, 1.e-10_dp, &
                           rotated, status, unitary=[.FALSE., .TRUE.])
      IF (status /= 0) ERROR STOP 'Compressed unitary metric iteration failed'
      IF (MAXVAL(ABS(rotated%eigenvalues - full%eigenvalues)) > 1.e-12_dp) ERROR STOP 'Unitary metric covariance failed'
      IF (MAXVAL(ABS(rotated%leakage - full%leakage)) > 1.e-12_dp) ERROR STOP 'Leakage metric covariance failed'
      op(:, :, 2) = 2*metric
      CALL quadratic_dense(op, metric, query, weights, 2, .FALSE., 2, 100, 1.e-11_dp, 1.e-10_dp, &
                           rotated, status, unitary=[.FALSE., .TRUE.])
      IF (status /= -4) ERROR STOP 'Invalid expanding unitary compression accepted'
      PRINT '(A,ES15.6)', 'Compressed unitary full-normal analytic error: ', err
   END SUBROUTINE compressed_translation_tests

! **************************************************************************************************
!> \brief A trimerized periodic chain must equal its independently assembled three-band Bloch blocks.
! **************************************************************************************************
   SUBROUTINE folded_translation_tests()
      INTEGER, PARAMETER                                 :: cells = 6, n = 3*cells

      COMPLEX(dp)                                        :: bloch(3, 3, 2), identity(3, 3), &
                                                            metric(n, n), operators(n, n, 2), &
                                                            phase, query(2)
      INTEGER                                            :: first, i, ie, ik, iq, j, loc, status
      REAL(dp)                                           :: energy(5), expected(n), hopping, &
                                                            joined(n), k, maximum_error, q, &
                                                            weights(2)
      TYPE(quadratic_result_type)                        :: full, sector

      operators = 0
      metric = 0
      identity = 0
      DO i = 1, 3
         identity(i, i) = 1
      END DO
      DO i = 1, n
         metric(i, i) = 1
         j = MODULO(i, n) + 1
         operators(i, j, 2) = 1
         hopping = 0.8_dp
         IF (MODULO(i, 3) == 0) hopping = 1.2_dp
         operators(i, j, 1) = -hopping
         operators(j, i, 1) = -hopping
      END DO
      energy = [-2.2_dp, -1.0_dp, 0.0_dp, 0.7_dp, 2.0_dp]
      weights = [1.0_dp, 0.7_dp]
      maximum_error = 0
      DO ik = 0, 12
         k = -ACOS(-1.0_dp) + 2*ACOS(-1.0_dp)*REAL(ik, dp)/13
         DO ie = 1, SIZE(energy)
            query = [CMPLX(energy(ie), 0.0_dp, dp), EXP(CMPLX(0.0_dp, k, dp))]
            CALL quadratic_dense(operators, metric, query, weights, n, .FALSE., n, 100, &
                                 1.e-10_dp, 1.e-10_dp, full, status)
            IF (status /= 0) ERROR STOP 'Full trimerized quadratic scan failed'
            DO iq = 0, cells - 1
               q = 2*ACOS(-1.0_dp)*REAL(iq, dp)/cells
               phase = EXP(CMPLX(0.0_dp, q, dp))
               bloch = 0
               bloch(1, 2, 2) = 1
               bloch(2, 3, 2) = 1
               bloch(3, 1, 2) = phase
               bloch(1, 2, 1) = -0.8_dp
               bloch(2, 3, 1) = -0.8_dp
               bloch(3, 1, 1) = -1.2_dp*phase
               bloch(:, :, 1) = bloch(:, :, 1) + CONJG(TRANSPOSE(bloch(:, :, 1)))
               CALL quadratic_dense(bloch, identity, query, weights, 3, .FALSE., 3, 100, &
                                    1.e-10_dp, 1.e-10_dp, sector, status)
               IF (status /= 0) ERROR STOP 'Bloch-block quadratic scan failed'
               first = 3*iq + 1
               joined(first:first + 2) = sector%eigenvalues
            END DO
            DO i = 1, n
               loc = MINLOC(joined, DIM=1)
               expected(i) = joined(loc)
               joined(loc) = HUGE(1.0_dp)
            END DO
            maximum_error = MAX(maximum_error, MAXVAL(ABS(full%eigenvalues - expected)))
         END DO
      END DO
      IF (maximum_error > 1.e-11_dp) ERROR STOP 'Full/Bloch trimerized quadratic spectra differ'
      PRINT '(A,I0,A,ES15.6)', 'Trimerized full/Bloch spectra: ', 13*SIZE(energy)*n, &
         ' eigenvalues; maximum error ', maximum_error
   END SUBROUTINE folded_translation_tests

! **************************************************************************************************
!> \brief Exact Bloch phases, nonnormal open translations and an independent generalized Gram oracle.
! **************************************************************************************************
   SUBROUTINE translation_tests()
      INTEGER, PARAMETER                                 :: n = 18

      COMPLEX(dp) :: changed(n, n, 2), identity(n, n), metric(n, n), offset(2), &
         operators(n, n, 2), projector(n, n), query(2), split(n, n, 3), transform(n, n)
      INTEGER                                            :: i, j, status
      REAL(dp)                                           :: covariance_error, expected(n), &
                                                            gap_error, hopping, mean_error, &
                                                            nonnormal_error, phase, weights(2)
      TYPE(quadratic_result_type)                        :: altered, incorrect, iterative, original

      identity = 0
      operators = 0
      DO i = 1, n
         identity(i, i) = 1
         operators(i, MODULO(i, n) + 1, 2) = 1
      END DO
      operators(:, :, 1) = -operators(:, :, 2) - CONJG(TRANSPOSE(operators(:, :, 2)))
      phase = 2*ACOS(-1.0_dp)*5/REAL(n, dp)
      query = [CMPLX(-2*COS(phase) + 0.031_dp, 0.0_dp, dp), EXP(CMPLX(0.0_dp, phase, dp))]
      weights = [1.0_dp, 1.3_dp]
      CALL quadratic_dense(operators, identity, query, weights, 1, .FALSE., n, 200, &
                           1.e-10_dp, 1.e-10_dp, original, status)
      IF (status /= 0) ERROR STOP 'Unitary translation quadratic solve failed'
      IF (ABS(original%eigenvalues(1) - 0.031_dp**2) > 1.e-12_dp) ERROR STOP 'Exact Bloch energy/phase failed'
      IF (ABS(original%complex_means(2, 1) - query(2)) > 1.e-12_dp) ERROR STOP 'Bloch translation mean lost phase'
      IF (original%errors(2, 1) > 1.e-12_dp) ERROR STOP 'Bloch translation residual failed'
      CALL quadratic_dense(operators, identity, REAL(query, dp), weights, 1, .FALSE., n, 200, &
                           1.e-10_dp, 1.e-10_dp, altered, status)
      IF (status /= -1) ERROR STOP 'Hermitian interface accepted non-Hermitian translation'

      ! The truncated translation is nonnormal. Its adjoint cannot be replaced by itself
      ! or by an independent sum of squares of its Hermitian real/imaginary parts.
      operators(n, 1, 2) = 0
      operators(:, :, 1) = 0
      DO i = 1, n - 1
         hopping = 0.8_dp
         IF (MODULO(i, 3) == 0) hopping = 1.2_dp
         operators(i, i + 1, 1) = -hopping
         operators(i + 1, i, 1) = -hopping
      END DO
      operators(4, 4, 1) = 0.17_dp
      query = [CMPLX(0.33_dp, 0.0_dp, dp), EXP(CMPLX(0.0_dp, 0.73_dp, dp))]
      CALL quadratic_dense(operators, identity, query, weights, n, .FALSE., n, 200, &
                           1.e-10_dp, 1.e-10_dp, original, status)
      IF (status /= 0) ERROR STOP 'Open trimerized translation solve failed'
      CALL gram_reference(operators, identity, query, weights, expected)
      gap_error = MAXVAL(ABS(original%eigenvalues - expected))
      IF (gap_error > 1.e-11_dp) ERROR STOP 'General quadratic Gram oracle mismatch'
      CALL quadratic_dense(operators, identity, query, weights, 2, .TRUE., 10, 2000, &
                           1.e-10_dp, 1.e-10_dp, iterative, status)
      IF (status /= 0) ERROR STOP 'Nonnormal restarted iteration failed'
      IF (MAXVAL(ABS(iterative%eigenvalues - expected(:2))) > 1.e-10_dp) ERROR STOP 'Nonnormal Ritz spectrum mismatch'
      projector = MATMUL(iterative%states, CONJG(TRANSPOSE(iterative%states))) - &
                  MATMUL(original%states(:, :2), CONJG(TRANSPOSE(original%states(:, :2))))
      IF (MAXVAL(ABS(projector)) > 1.e-8_dp) ERROR STOP 'Nonnormal invariant-subspace mismatch'
      split(:, :, 1) = operators(:, :, 1)
      split(:, :, 2) = 0.5_dp*(operators(:, :, 2) + CONJG(TRANSPOSE(operators(:, :, 2))))
      split(:, :, 3) = CMPLX(0.0_dp, -0.5_dp, dp)*(operators(:, :, 2) - CONJG(TRANSPOSE(operators(:, :, 2))))
      CALL quadratic_dense(split, identity, [REAL(query(1), dp), REAL(query(2), dp), AIMAG(query(2))], &
                           [weights(1), weights(2), weights(2)], n, .FALSE., n, 200, &
                           1.e-10_dp, 1.e-10_dp, incorrect, status)
      IF (status /= 0) ERROR STOP 'Hermitian-split comparison solve failed'
      nonnormal_error = MAXVAL(ABS(original%eigenvalues - incorrect%eigenvalues))
      IF (nonnormal_error < 0.01_dp) ERROR STOP 'Open-translation counterexample lost nonnormality'

      transform = identity
      DO i = 1, n - 1
         transform(i, i) = 1 + 0.02_dp*i
         transform(i, i + 1) = CMPLX(0.21_dp, -0.12_dp, dp)
      END DO
      metric = MATMUL(CONJG(TRANSPOSE(transform)), transform)
      DO j = 1, 2
         changed(:, :, j) = MATMUL(CONJG(TRANSPOSE(transform)), MATMUL(operators(:, :, j), transform))
      END DO
      CALL quadratic_dense(changed, metric, query, weights, n, .FALSE., n, 200, &
                           1.e-10_dp, 1.e-10_dp, altered, status)
      IF (status /= 0) ERROR STOP 'Nonnormal complex AO basis failed'
      CALL gram_reference(changed, metric, query, weights, expected)
      covariance_error = MAXVAL(ABS(original%eigenvalues - altered%eigenvalues))
      IF (covariance_error > 1.e-11_dp) ERROR STOP 'Nonnormal metric covariance failed'
      IF (MAXVAL(ABS(altered%eigenvalues - expected)) > 1.e-11_dp) ERROR STOP 'Nonorthogonal Gram oracle mismatch'
      mean_error = MAXVAL(ABS(original%complex_means - altered%complex_means))
      IF (mean_error > 1.e-10_dp) ERROR STOP 'Complex expectation covariance failed'
      IF (MAXVAL(ABS(original%errors - altered%errors)) > 1.e-10_dp) ERROR STOP 'Complex residual covariance failed'
      offset = [CMPLX(0.21_dp, -0.12_dp, dp), CMPLX(-0.43_dp, 0.31_dp, dp)]
      DO j = 1, 2
         changed(:, :, j) = changed(:, :, j) + offset(j)*metric
      END DO
      CALL quadratic_dense(changed, metric, query + offset, weights, n, .FALSE., n, 200, &
                           1.e-10_dp, 1.e-10_dp, altered, status)
      IF (status /= 0) ERROR STOP 'Complex operator/query origin solve failed'
      IF (MAXVAL(ABS(original%eigenvalues - altered%eigenvalues)) > 1.e-11_dp) ERROR STOP 'Complex origin covariance failed'
      DO j = 1, 2
         IF (MAXVAL(ABS(altered%complex_means(j, :) - original%complex_means(j, :) - offset(j))) > 1.e-10_dp) THEN
            ERROR STOP 'Complex expectation origin covariance failed'
         END IF
      END DO
      PRINT '(A,4ES15.6)', 'Translation Gram/metric/complex-mean/split-error: ', &
         gap_error, covariance_error, mean_error, nonnormal_error
      changed(:, :, 1) = MATMUL(CONJG(TRANSPOSE(transform)), MATMUL(operators(:, :, 1), transform))
      changed(:, :, 2) = metric
      query(2) = 1
      CALL quadratic_dense(changed, metric, query, weights, 2, .TRUE., 10, 2000, &
                           1.e-10_dp, 1.e-10_dp, iterative, status, unitary=[.FALSE., .TRUE.])
      IF (status /= 0) ERROR STOP 'Restarted identity translation failed'
      IF (iterative%iterations < 2) ERROR STOP 'Identity translation did not exercise restarts'
      IF (MAXVAL(iterative%errors(2, :)) > 5.e-8_dp) ERROR STOP 'Identity full residual exceeds roundoff'
      IF (MAXVAL(iterative%projected_errors(2, :)) > 1.e-12_dp) ERROR STOP 'Identity projected residual drift'
      IF (MAXVAL(iterative%leakage(2, :)) > 5.e-8_dp) ERROR STOP 'Identity projection leakage exceeds roundoff'
   END SUBROUTINE translation_tests

! **************************************************************************************************
!> \brief Independent explicit normal-product reference using LU solves and a generalized eigensolve.
!> \param operators Covariant operators
!> \param metric Positive AO overlap
!> \param query Complex operator queries
!> \param weights Nonnegative residual scales
!> \param values Complete ascending generalized spectrum
! **************************************************************************************************
   SUBROUTINE gram_reference(operators, metric, query, weights, values)
      COMPLEX(dp), INTENT(IN) :: operators(:, :, :), metric(:, :), query(:)
      REAL(dp), INTENT(IN) :: weights(:)
      REAL(dp), INTENT(OUT) :: values(:)
      COMPLEX(dp), ALLOCATABLE :: gram(:, :), shifted(:, :), solved(:, :), lu(:, :), work(:)
      REAL(dp), ALLOCATABLE :: rwork(:)
      INTEGER, ALLOCATABLE :: pivots(:)
      INTEGER :: n, i, status
#if defined(__HAS_IEEE_EXCEPTIONS)
      LOGICAL :: halt(SIZE(IEEE_ALL))
#endif

      n = SIZE(metric, 1)
      ALLOCATE (gram(n, n), shifted(n, n), solved(n, n), lu(n, n), pivots(n), work(4*n), rwork(4*n))
      gram = 0
      DO i = 1, SIZE(query)
         lu(:, :) = metric
         shifted(:, :) = operators(:, :, i) - query(i)*metric
         solved(:, :) = shifted
#if defined(__HAS_IEEE_EXCEPTIONS)
         CALL ieee_get_halting_mode(IEEE_ALL, halt)
         CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
         CALL zgesv(n, n, lu, n, pivots, solved, n, status)
#if defined(__HAS_IEEE_EXCEPTIONS)
         CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif
         IF (status /= 0) ERROR STOP 'Explicit Gram metric solve failed'
         gram(:, :) = gram + weights(i)**2*MATMUL(CONJG(TRANSPOSE(shifted)), solved)
      END DO
      IF (MAXVAL(ABS(gram - CONJG(TRANSPOSE(gram)))) > 1.e-11_dp) ERROR STOP 'Explicit Gram is not Hermitian'
      lu(:, :) = metric
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_get_halting_mode(IEEE_ALL, halt)
      CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
      CALL zhegv(1, 'N', 'L', n, gram, n, lu, n, values, work, SIZE(work), rwork, status)
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif
      IF (status /= 0) ERROR STOP 'Explicit generalized Gram spectrum failed'
   END SUBROUTINE gram_reference

! **************************************************************************************************
!> \brief Three-dimensional trigonometric coordinates, basis covariance and torus-query invariance.
! **************************************************************************************************
   SUBROUTINE periodic_tests()
      INTEGER, PARAMETER                                 :: n = 8

      COMPLEX(dp)                                        :: changed(n, n, 7), identity(n, n), &
                                                            metric(n, n), operators(n, n, 7), &
                                                            transform(n, n)
      INTEGER                                            :: i, j, k, status
      REAL(dp)                                           :: diagonal(n), expected(n), move(3), &
                                                            radius(3), shift(7), theta(3), value, &
                                                            weights(7)
      TYPE(quadratic_result_type)                        :: altered, iterative, original

      identity = 0
      operators = 0
      theta = [0.21_dp, -0.37_dp, 1.43_dp]
      radius = [0.7_dp, 1.2_dp, 0.9_dp]
      shift(1) = 0.73_dp
      weights(1) = 1
      DO j = 1, 3
         shift(2*j:2*j + 1) = [COS(theta(j)), SIN(theta(j))]
         weights(2*j:2*j + 1) = radius(j)
      END DO
      DO i = 1, n
         identity(i, i) = 1
         operators(i, i, 1) = 0.2_dp*i
         DO j = 1, 3
            value = 2*ACOS(-1.0_dp)*MODULO(REAL(i*(2*j - 1), dp)/REAL(n + j, dp), 1.0_dp)
            operators(i, i, 2*j) = COS(value)
            operators(i, i, 2*j + 1) = SIN(value)
         END DO
         diagonal(i) = SUM((weights*ABS(operators(i, i, :) - shift))**2)
      END DO
      DO i = 1, n
         k = MINLOC(diagonal, DIM=1)
         expected(i) = diagonal(k)
         diagonal(k) = HUGE(1.0_dp)
      END DO
      CALL quadratic_dense(operators, identity, shift, weights, n, .FALSE., n, 100, 1.e-10_dp, 1.e-10_dp, original, status)
      IF (status /= 0) ERROR STOP 'Periodic quadratic diagonal oracle failed'
      IF (MAXVAL(ABS(original%eigenvalues - expected)) > 1.e-12_dp) ERROR STOP 'Seven-observable analytic spectrum failed'
      DO i = 1, n - 1
         operators(i, i + 1, 1) = CMPLX(-0.12_dp, 0.04_dp, dp)
         operators(i + 1, i, 1) = CONJG(operators(i, i + 1, 1))
      END DO
      ! The six initial vectors leave only one or two slots for the three-vector residual block.
      DO k = n - 1, n
         CALL quadratic_dense(operators(:k, :k, :), identity(:k, :k), shift, weights, 3, .FALSE., k, 100, &
                              1.e-10_dp, 1.e-10_dp, original, status)
         IF (status /= 0) ERROR STOP 'Partial residual block reference failed'
         CALL quadratic_dense(operators(:k, :k, :), identity(:k, :k), shift, weights, 3, .TRUE., k, 2, &
                              1.e-10_dp, 1.e-10_dp, iterative, status)
         IF (status /= 0) ERROR STOP 'Partial residual block did not fill the available space'
         IF (iterative%iterations /= 2) ERROR STOP 'Partial residual block iteration count mismatch'
         IF (MAXVAL(ABS(original%eigenvalues - iterative%eigenvalues)) > 1.e-12_dp) THEN
            ERROR STOP 'Partial residual block eigenvalue mismatch'
         END IF
      END DO
      CALL quadratic_dense(operators, identity, shift, weights, 3, .FALSE., n, 100, 1.e-10_dp, 1.e-10_dp, original, status)
      IF (status /= 0) ERROR STOP 'Periodic quadratic hopping reference failed'
      CALL quadratic_dense(operators, identity, shift, weights, 3, .TRUE., n, 100, 1.e-10_dp, 1.e-10_dp, iterative, status)
      IF (status /= 0) ERROR STOP 'Periodic quadratic iteration failed'
      IF (MAXVAL(ABS(original%eigenvalues - iterative%eigenvalues)) > 1.e-12_dp) ERROR STOP 'Periodic iteration mismatch'
      transform = identity
      DO i = 1, n - 1
         transform(i, i + 1) = CMPLX(0.13_dp, -0.08_dp, dp)
      END DO
      metric = MATMUL(CONJG(TRANSPOSE(transform)), transform)
      DO j = 1, 7
         changed(:, :, j) = MATMUL(CONJG(TRANSPOSE(transform)), MATMUL(operators(:, :, j), transform))
      END DO
      CALL quadratic_dense(changed, metric, shift, weights, 3, .FALSE., n, 100, 1.e-10_dp, 1.e-10_dp, altered, status)
      IF (status /= 0) ERROR STOP 'Periodic complex metric failed'
      IF (MAXVAL(ABS(original%eigenvalues - altered%eigenvalues)) > 1.e-12_dp) ERROR STOP 'Periodic metric covariance failed'
      CALL quadratic_dense(changed, metric, shift, weights, 3, .TRUE., n, 2, 1.e-10_dp, 1.e-10_dp, iterative, status)
      IF (status /= 0) ERROR STOP 'Partial residual block with complex metric failed'
      IF (MAXVAL(ABS(altered%eigenvalues - iterative%eigenvalues)) > 1.e-12_dp) THEN
         ERROR STOP 'Partial residual block metric covariance failed'
      END IF
      move = [0.43_dp, -0.67_dp, 0.31_dp]
      changed = operators
      DO j = 1, 3
         changed(:, :, 2*j) = COS(move(j))*operators(:, :, 2*j) - SIN(move(j))*operators(:, :, 2*j + 1)
         changed(:, :, 2*j + 1) = SIN(move(j))*operators(:, :, 2*j) + COS(move(j))*operators(:, :, 2*j + 1)
         shift(2*j:2*j + 1) = [COS(theta(j) + move(j)), SIN(theta(j) + move(j))]
      END DO
      CALL quadratic_dense(changed, identity, shift, weights, 3, .FALSE., n, 100, 1.e-10_dp, 1.e-10_dp, altered, status)
      IF (status /= 0) ERROR STOP 'Shifted periodic coordinates failed'
  IF (MAXVAL(ABS(original%eigenvalues - altered%eigenvalues)) > 1.e-12_dp) ERROR STOP 'Periodic coordinate-origin covariance failed'
      DO j = 1, 3
         value = theta(j) + 2*ACOS(-1.0_dp)*(j - 2)
         shift(2*j:2*j + 1) = [COS(value), SIN(value)]
      END DO
      CALL quadratic_dense(operators, identity, shift, weights, 3, .FALSE., n, 100, 1.e-10_dp, 1.e-10_dp, altered, status)
      IF (status /= 0) ERROR STOP 'Periodic query wrapping failed'
      IF (MAXVAL(ABS(original%eigenvalues - altered%eigenvalues)) > 1.e-12_dp) ERROR STOP 'Periodic query invariance failed'
      CALL quadratic_dense(operators, identity, shift, weights(:6), 1, .FALSE., n, 100, 1.e-10_dp, 1.e-10_dp, altered, status)
      IF (status /= -1) ERROR STOP 'Mismatched observable weights accepted'
      CALL quadratic_dense(operators(:, :, :0), identity, shift(:0), weights(:0), 1, .FALSE., n, 100, &
                           1.e-10_dp, 1.e-10_dp, altered, status)
      IF (status /= -1) ERROR STOP 'Empty observable tuple accepted'
      PRINT *, 'Periodic quadratic tuple, metric, iteration and origin checks passed.'
   END SUBROUTINE periodic_tests
END PROGRAM quadratic_pseudospectrum_unittest
