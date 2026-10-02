!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Dense reference spectral localizers in a nonorthogonal finite basis.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE spectral_localizer
   USE ieee_arithmetic, ONLY: ieee_is_finite
#if defined(__HAS_IEEE_EXCEPTIONS)
   USE ieee_exceptions, ONLY: IEEE_ALL, ieee_get_halting_mode, ieee_set_halting_mode
#endif
   USE kinds, ONLY: dp

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: localizer_evaluate, localizer_matrix_evaluate, localizer_inertia, localizer_pfaffian_sign
   PUBLIC :: localizer_periodic_evaluate, localizer_aii_skew
   PUBLIC :: localizer_3d_evaluate, localizer_3d_matrix_evaluate

   INTERFACE localizer_inertia
      MODULE PROCEDURE localizer_inertia_complex, localizer_inertia_real
   END INTERFACE localizer_inertia

CONTAINS

! **************************************************************************************************
!> \brief Finite three-dimensional class AII localizer, using all three Cartesian moments.
!> \param h Covariant physical-spin Hamiltonian, all up AOs followed by all down AOs
!> \param s Positive physical-spin AO metric
!> \param r Three covariant Hermitian position operators
!> \param point Cartesian query coordinates
!> \param energy Absolute energy query
!> \param kappa Positive energy/length scale
!> \param matrix_tol Relative matrix-structure tolerance
!> \param metric_tol Minimum overlap eigenvalue ratio
!> \param gap_tol Absolute unresolved-gap threshold
!> \param invariant Strong local Z2 index, valid only for status zero
!> \param gap Smallest metric-normalized singular value of the chiral block
!> \param metric_ratio Smallest/largest physical metric eigenvalue
!> \param status Same convention as localizer_evaluate
! **************************************************************************************************
   SUBROUTINE localizer_3d_evaluate(h, s, r, point, energy, kappa, matrix_tol, metric_tol, gap_tol, &
                                    invariant, gap, metric_ratio, status)
      COMPLEX(dp), INTENT(IN)                            :: h(:, :), s(:, :), r(:, :, :)
      REAL(dp), INTENT(IN)                               :: point(3), energy, kappa, matrix_tol, &
                                                            metric_tol, gap_tol
      INTEGER, INTENT(OUT)                               :: invariant
      REAL(dp), INTENT(OUT)                              :: gap, metric_ratio
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(dp), PARAMETER                             :: ii = CMPLX(0.0_dp, 1.0_dp, dp)

      COMPLEX(dp), ALLOCATABLE                           :: block(:, :), shifted(:, :, :)
      INTEGER                                            :: j, n

      invariant = 0
      gap = 0.0_dp
      metric_ratio = 0.0_dp
      status = -1
      n = SIZE(h, 1)
      IF (n < 2 .OR. MOD(n, 2) /= 0 .OR. SIZE(h, 2) /= n) RETURN
      IF (ANY(SHAPE(s) /= [n, n]) .OR. ANY(SHAPE(r) /= [n, n, 3])) RETURN
      IF (.NOT. ALL(ieee_is_finite([point, energy, kappa, matrix_tol, metric_tol, gap_tol]))) RETURN
      IF (kappa <= 0.0_dp .OR. matrix_tol <= 0.0_dp .OR. matrix_tol >= 1.0_dp) RETURN
      IF (metric_tol <= 0.0_dp .OR. metric_tol >= 1.0_dp .OR. gap_tol < 0.0_dp) RETURN
      IF (.NOT. hermitian(h, matrix_tol) .OR. .NOT. hermitian(s, matrix_tol)) RETURN
      DO j = 1, 3
         IF (.NOT. hermitian(r(:, :, j), matrix_tol)) RETURN
      END DO
      status = -3
      IF (.NOT. self_dual(h, matrix_tol) .OR. .NOT. self_dual(s, matrix_tol)) RETURN
      DO j = 1, 3
         IF (.NOT. self_dual(r(:, :, j), matrix_tol)) RETURN
      END DO
      ALLOCATE (block(2*n, 2*n), shifted(n, n, 3))
      DO j = 1, 3
         shifted(:, :, j) = kappa*(r(:, :, j) - point(j)*s)
      END DO
      ! A = I_aux (H-E S) - i sum_j sigma_j (kappa (X_j-x_j S)).
      block(:n, :n) = h - energy*s - ii*shifted(:, :, 3)
      block(n + 1:, n + 1:) = h - energy*s + ii*shifted(:, :, 3)
      block(:n, n + 1:) = -ii*shifted(:, :, 1) - shifted(:, :, 2)
      block(n + 1:, :n) = -ii*shifted(:, :, 1) + shifted(:, :, 2)
      CALL localizer_3d_matrix_evaluate(block, s, matrix_tol, metric_tol, gap_tol, invariant, gap, metric_ratio, status)
   END SUBROUTINE localizer_3d_evaluate

! **************************************************************************************************
!> \brief Metric gap and real determinant sign of a three-dimensional AII chiral block.
!> \param block Covariant A=I_aux(H-E S)-i kappa sigma.(X-x S); auxiliary index outermost
!> \param s Positive physical-spin AO metric
!> \param matrix_tol Relative structure tolerance
!> \param metric_tol Minimum overlap eigenvalue ratio
!> \param gap_tol Absolute unresolved-gap threshold
!> \param invariant Zero for positive determinant, one for negative; valid only on success
!> \param gap Minimum absolute generalized eigenvalue of [0,A;A^dagger,0]
!> \param metric_ratio Smallest/largest physical metric eigenvalue
!> \param status Same convention as localizer_evaluate
!> \note A positive atomic Hamiltonian gives positive determinant in every nonsingular AO basis.
!>       Unlike the 2D Pfaffian, no basis-orientation correction is required for this determinant.
! **************************************************************************************************
   SUBROUTINE localizer_3d_matrix_evaluate(block, s, matrix_tol, metric_tol, gap_tol, &
                                           invariant, gap, metric_ratio, status)
      COMPLEX(dp), INTENT(IN) :: block(:, :), s(:, :)
      REAL(dp), INTENT(IN) :: matrix_tol, metric_tol, gap_tol
      INTEGER, INTENT(OUT) :: invariant, status
      REAL(dp), INTENT(OUT) :: gap, metric_ratio
      COMPLEX(dp), ALLOCATABLE :: pencil(:, :), metric(:, :)
      COMPLEX(dp) :: value
      REAL(dp), ALLOCATABLE :: real_block(:, :)
      REAL(dp) :: residual, scale, metric_scale, error_budget
      INTEGER :: n, m, q, i, j, pi, pj, si, sj, signature, determinant
#if defined(__HAS_IEEE_EXCEPTIONS)
      LOGICAL :: halt(SIZE(IEEE_ALL))
#endif

      status = -1
      invariant = 0
      gap = 0.0_dp
      metric_ratio = 0.0_dp
      n = SIZE(s, 1)
      IF (n < 2 .OR. MOD(n, 2) /= 0 .OR. SIZE(s, 2) /= n) RETURN
      m = 2*n
      IF (ANY(SHAPE(block) /= [m, m])) RETURN
      IF (.NOT. ALL(ieee_is_finite([matrix_tol, metric_tol, gap_tol]))) RETURN
      IF (matrix_tol <= 0.0_dp .OR. matrix_tol >= 1.0_dp) RETURN
      IF (metric_tol <= 0.0_dp .OR. metric_tol >= 1.0_dp .OR. gap_tol < 0.0_dp) RETURN
      IF (.NOT. hermitian(s, matrix_tol)) RETURN
      IF (.NOT. ALL(ieee_is_finite(REAL(block, dp))) .OR. .NOT. ALL(ieee_is_finite(AIMAG(block)))) RETURN
      status = -3
      IF (.NOT. self_dual(s, matrix_tol)) RETURN
      ALLOCATE (real_block(m, m))
      q = n/2
      residual = 0.0_dp
      ! Q=(I-i C)/sqrt(2), C=J_aux tensor J_spin: Q^dagger A Q is REAL, not skew.
      DO j = 1, m
         pj = (3 - (j - 1)/q)*q + MODULO(j - 1, q) + 1
         sj = MERGE(1, -1, j <= q .OR. j > 3*q)
         DO i = 1, m
            pi = (3 - (i - 1)/q)*q + MODULO(i - 1, q) + 1
            si = MERGE(1, -1, i <= q .OR. i > 3*q)
            value = 0.5_dp*(block(i, j) + si*sj*block(pi, pj) + &
                            CMPLX(0.0_dp, 1.0_dp, dp)*(si*block(pi, j) - sj*block(i, pj)))
            residual = MAX(residual, ABS(AIMAG(value)))
            real_block(i, j) = REAL(value, dp)
         END DO
      END DO
      scale = MAX(TINY(1.0_dp), MAXVAL(ABS(block)))
      IF (.NOT. ALL(ieee_is_finite(real_block)) .OR. .NOT. ieee_is_finite(residual)) RETURN
      IF (residual > matrix_tol*scale) RETURN

      ALLOCATE (pencil(2*m, 2*m), metric(m, m))
      pencil = 0.0_dp
      pencil(:m, m + 1:) = block
      pencil(m + 1:, :m) = CONJG(TRANSPOSE(block))
      metric = 0.0_dp
      metric(:n, :n) = s
      metric(n + 1:, n + 1:) = s
      CALL localizer_matrix_evaluate(pencil, metric, .FALSE., matrix_tol, metric_tol, gap_tol, &
                                     signature, gap, metric_ratio, status)
      IF (status /= 0) RETURN
      status = -4
      IF (signature /= 0) RETURN
      ! ||S||_max <= lambda_max(S), hence this lower bound is conservative.
      metric_scale = metric_ratio*MAXVAL(ABS(s))
      status = 1
      IF (gap <= REAL(m, dp)*(residual + 100.0_dp*EPSILON(1.0_dp)*scale)/metric_scale) RETURN
      DEALLOCATE (pencil, metric)
      real_block = real_block/scale
      ! Frobenius residuals bound the spectral perturbation after metric normalization.
      error_budget = gap*metric_scale/scale - REAL(m, dp)*(residual/scale + 100.0_dp*EPSILON(1.0_dp))
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_get_halting_mode(IEEE_ALL, halt)
      CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
      CALL real_determinant_sign(real_block, error_budget, determinant, status)
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif
      IF (status /= 0) RETURN
      invariant = (1 - determinant)/2
      status = 0
   END SUBROUTINE localizer_3d_matrix_evaluate

! **************************************************************************************************
!> \brief LU determinant sign with a QR fallback when pivot growth exhausts the error budget.
!> \param matrix Scaled real nonsingular chiral block
!> \param error_budget Maximum admissible Frobenius backward error from the independent metric gap
!> \param determinant Sign of the determinant, only valid on success
!> \param status Zero on success, one for unresolved numerical error, minus four for LAPACK failure
!> \note The caller disables IEEE halting while LAPACK is active and restores it afterwards.
! **************************************************************************************************
   SUBROUTINE real_determinant_sign(matrix, error_budget, determinant, status)
      REAL(dp), INTENT(IN)                               :: matrix(:, :), error_budget
      INTEGER, INTENT(OUT)                               :: determinant, status

      INTEGER                                            :: i, info, method, n
      INTEGER, ALLOCATABLE                               :: pivots(:)
      REAL(dp)                                           :: error, roundoff
      REAL(dp), ALLOCATABLE                              :: factors(:, :), permuted(:, :), &
                                                            PRODUCT(:, :), row(:), tau(:), work(:)

      n = SIZE(matrix, 1)
      ALLOCATE (factors(n, n), PRODUCT(n, n), permuted(n, n), row(n), pivots(n), tau(n), work(64*n))
      determinant = 0
      status = -4
      DO method = 1, 2
         factors(:, :) = matrix
         permuted(:, :) = matrix
         IF (method == 1) THEN
            CALL dgetrf(n, n, factors, n, pivots, info)
         ELSE
            CALL dgeqrf(n, n, factors, n, tau, work, SIZE(work), info)
         END IF
         IF (info /= 0 .OR. .NOT. ALL(ieee_is_finite(factors))) CYCLE
         PRODUCT(:, :) = factors
         DO i = 1, n
            PRODUCT(i + 1:, i) = 0.0_dp
         END DO
         ! Include rounding of both the triangular product and the residual evaluation.
         roundoff = 8.0_dp*EPSILON(1.0_dp)*REAL(n, dp)*SQRT(SUM(product**2))
         IF (method == 1) THEN
            CALL dtrmm('L', 'L', 'N', 'U', n, n, 1.0_dp, factors, n, product, n)
            DO i = 1, n
               row(:) = permuted(i, :)
               permuted(i, :) = permuted(pivots(i), :)
               permuted(pivots(i), :) = row
            END DO
         ELSE
            CALL dormqr('L', 'N', n, n, n, factors, n, tau, product, n, work, SIZE(work), info)
            IF (info /= 0) CYCLE
         END IF
         IF (.NOT. ALL(ieee_is_finite(product))) CYCLE
         error = SQRT(SUM((product - permuted)**2))
         IF (method == 1) THEN
            ! Bound cancellation in the reconstructed product by |L| |U|, not by |L U|.
            permuted(:, :) = ABS(factors)
            PRODUCT(:, :) = permuted
            DO i = 1, n
               PRODUCT(i + 1:, i) = 0.0_dp
            END DO
            CALL dtrmm('L', 'L', 'N', 'U', n, n, 1.0_dp, permuted, n, product, n)
            roundoff = 8.0_dp*EPSILON(1.0_dp)*REAL(n, dp)*SQRT(SUM(product**2))
         END IF
         error = error + roundoff
         status = 1
         IF (.NOT. ieee_is_finite(error) .OR. error >= error_budget) CYCLE
         determinant = 1
         DO i = 1, n
            IF (factors(i, i) == 0.0_dp) EXIT
            IF (factors(i, i) < 0.0_dp) determinant = -determinant
            IF (method == 1) THEN
               IF (pivots(i) /= i) determinant = -determinant
            ELSE
               IF (tau(i) /= 0.0_dp) determinant = -determinant
            END IF
         END DO
         IF (i <= n) CYCLE
         status = 0
         RETURN
      END DO
      determinant = 0
   END SUBROUTINE real_determinant_sign

! **************************************************************************************************
!> \brief Two-dimensional periodic localizer from covariant trigonometric position operators.
!> \param h Periodic Hamiltonian, including all boundary-crossing matrix elements
!> \param s Positive definite overlap
!> \param cosine Covariant cos(G_j.r) matrices, not cos of projected AO positions
!> \param sine Covariant sin(G_j.r) matrices
!> \param origin Two dimensionless phase offsets G_j.r0
!> \param energy Absolute reference energy
!> \param eta Positive energy scale multiplying the torus coordinates
!> \param aii Use physical spinful time reversal
!> \param matrix_tol Relative structure tolerance
!> \param metric_tol Minimum overlap eigenvalue ratio
!> \param gap_tol Absolute unresolved-gap threshold in energy units
!> \param invariant Index in the same orientation convention as localizer_evaluate
!> \param gap Minimum absolute generalized eigenvalue, in energy units
!> \param metric_ratio Smallest/largest overlap eigenvalue ratio
!> \param status Same status convention as localizer_evaluate
!> \note Adapted from Doll, Loring, Schulz-Baldes (2025), Eq. (1): minus eta times
!>       their localizer, conjugated by diag(I,-I), to retain our (H,X,Y) orientation.
!>       A resolved finite-matrix signature is not a certificate of bulk convergence.
! **************************************************************************************************
   SUBROUTINE localizer_periodic_evaluate(h, s, cosine, sine, origin, energy, eta, aii, &
                                          matrix_tol, metric_tol, gap_tol, invariant, gap, metric_ratio, status)
      COMPLEX(dp), INTENT(IN)                            :: h(:, :), s(:, :), cosine(:, :, :), &
                                                            sine(:, :, :)
      REAL(dp), INTENT(IN)                               :: origin(2), energy, eta
      LOGICAL, INTENT(IN)                                :: aii
      REAL(dp), INTENT(IN)                               :: matrix_tol, metric_tol, gap_tol
      INTEGER, INTENT(OUT)                               :: invariant
      REAL(dp), INTENT(OUT)                              :: gap, metric_ratio
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(dp), ALLOCATABLE                           :: l(:, :), mass(:, :), &
                                                            shifted_sine(:, :, :)
      INTEGER                                            :: j, n

      status = -1
      invariant = 0
      gap = 0.0_dp
      metric_ratio = 0.0_dp
      n = SIZE(h, 1)
      IF (n < 1 .OR. SIZE(h, 2) /= n) RETURN
      IF (ANY(SHAPE(s) /= [n, n])) RETURN
      IF (ANY(SHAPE(cosine) /= [n, n, 2]) .OR. ANY(SHAPE(sine) /= [n, n, 2])) RETURN
      IF (.NOT. ALL(ieee_is_finite([origin, energy, eta, matrix_tol]))) RETURN
      IF (eta <= 0.0_dp .OR. matrix_tol <= 0.0_dp .OR. matrix_tol >= 1.0_dp) RETURN
      IF (.NOT. hermitian(h, matrix_tol) .OR. .NOT. hermitian(s, matrix_tol)) RETURN
      DO j = 1, 2
         IF (.NOT. hermitian(cosine(:, :, j), matrix_tol)) RETURN
         IF (.NOT. hermitian(sine(:, :, j), matrix_tol)) RETURN
      END DO
      IF (aii) THEN
         status = -3
         IF (MOD(n, 2) /= 0) RETURN
         IF (.NOT. self_dual(h, matrix_tol) .OR. .NOT. self_dual(s, matrix_tol)) RETURN
         DO j = 1, 2
            IF (.NOT. self_dual(cosine(:, :, j), matrix_tol)) RETURN
            IF (.NOT. self_dual(sine(:, :, j), matrix_tol)) RETURN
         END DO
      END IF
      ALLOCATE (mass(n, n), shifted_sine(n, n, 2), l(2*n, 2*n))
      mass(:, :) = h - energy*s
      DO j = 1, 2
         mass(:, :) = mass - eta*(s - COS(origin(j))*cosine(:, :, j) - SIN(origin(j))*sine(:, :, j))
         shifted_sine(:, :, j) = eta*(COS(origin(j))*sine(:, :, j) - SIN(origin(j))*cosine(:, :, j))
      END DO
      l(:n, :n) = mass
      l(n + 1:, n + 1:) = -mass
      l(:n, n + 1:) = shifted_sine(:, :, 1) - CMPLX(0.0_dp, 1.0_dp, dp)*shifted_sine(:, :, 2)
      l(n + 1:, :n) = CONJG(TRANSPOSE(l(:n, n + 1:)))
      CALL localizer_matrix_evaluate(l, s, aii, matrix_tol, metric_tol, gap_tol, &
                                     invariant, gap, metric_ratio, status)
   END SUBROUTINE localizer_periodic_evaluate

! **************************************************************************************************
!> \brief Evaluate the 2D class A or AII localizer without projecting onto occupied states.
!> \param h Covariant Hamiltonian matrix
!> \param s Positive definite AO overlap
!> \param r Covariant Cartesian position matrices in the two oriented directions
!> \param point Position in the same coordinates and length units as r
!> \param energy Absolute reference energy in the units of h
!> \param kappa Positive energy/length scale
!> \param aii Use spinful time reversal, with all up AOs followed by all down AOs
!> \param matrix_tol Relative Hermiticity and time-reversal tolerance
!> \param metric_tol Minimum allowed smallest/largest overlap eigenvalue ratio
!> \param gap_tol Absolute unresolved-gap threshold in energy units
!> \param invariant Half-signature (A), or 0 trivial / 1 nontrivial (AII); valid only if status=0
!> \param gap Minimum absolute generalized localizer eigenvalue
!> \param metric_ratio Smallest/largest overlap eigenvalue ratio
!> \param status 0 resolved; 1 unresolved gap; -1 invalid data; -2 bad metric; -3 broken TR;
!>               -4 factorization failure; positive values above 1 indicate LAPACK failure
!> \note This is a finite-basis diagnostic, not a certificate of basis or volume convergence.
! **************************************************************************************************
   SUBROUTINE localizer_evaluate(h, s, r, point, energy, kappa, aii, matrix_tol, metric_tol, &
                                 gap_tol, invariant, gap, metric_ratio, status)
      COMPLEX(KIND=dp), INTENT(IN)                       :: h(:, :), s(:, :), r(:, :, :)
      REAL(KIND=dp), INTENT(IN)                          :: point(2), energy, kappa
      LOGICAL, INTENT(IN)                                :: aii
      REAL(KIND=dp), INTENT(IN)                          :: matrix_tol, metric_tol, gap_tol
      INTEGER, INTENT(OUT)                               :: invariant
      REAL(KIND=dp), INTENT(OUT)                         :: gap, metric_ratio
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(KIND=dp), ALLOCATABLE                      :: localizer(:, :)
      INTEGER                                            :: i, n

      status = -1
      invariant = 0
      gap = 0.0_dp
      metric_ratio = 0.0_dp
      n = SIZE(h, 1)
      IF (n < 1 .OR. SIZE(h, 2) /= n) RETURN
      IF (ANY(SHAPE(s) /= [n, n]) .OR. ANY(SHAPE(r) /= [n, n, 2])) RETURN
      IF (.NOT. ALL(ieee_is_finite([point, energy, kappa, matrix_tol, metric_tol, gap_tol]))) RETURN
      IF (kappa <= 0.0_dp .OR. matrix_tol <= 0.0_dp .OR. matrix_tol >= 1.0_dp) RETURN
      IF (metric_tol <= 0.0_dp .OR. metric_tol >= 1.0_dp .OR. gap_tol < 0.0_dp) RETURN
      IF (.NOT. hermitian(h, matrix_tol) .OR. .NOT. hermitian(s, matrix_tol)) RETURN
      DO i = 1, 2
         IF (.NOT. hermitian(r(:, :, i), matrix_tol)) RETURN
      END DO
      IF (aii) THEN
         status = -3
         IF (MOD(n, 2) /= 0) RETURN
         IF (.NOT. self_dual(h, matrix_tol) .OR. .NOT. self_dual(s, matrix_tol)) RETURN
         DO i = 1, 2
            IF (.NOT. self_dual(r(:, :, i), matrix_tol)) RETURN
         END DO
      END IF

      ALLOCATE (localizer(2*n, 2*n))
      localizer(:n, :n) = h - energy*s
      localizer(n + 1:, n + 1:) = -localizer(:n, :n)
      localizer(:n, n + 1:) = kappa*((r(:, :, 1) - point(1)*s) - &
                                     CMPLX(0.0_dp, 1.0_dp, dp)*(r(:, :, 2) - point(2)*s))
      localizer(n + 1:, :n) = CONJG(TRANSPOSE(localizer(:n, n + 1:)))
      CALL localizer_matrix_evaluate(localizer, s, aii, matrix_tol, metric_tol, gap_tol, &
                                     invariant, gap, metric_ratio, status)
   END SUBROUTINE localizer_evaluate

! **************************************************************************************************
!> \brief Dense reference evaluation of an assembled localizer with metric diag(S,S).
!> \param localizer Covariant Hermitian localizer
!> \param s Positive definite physical AO/spinor overlap
!> \param aii Whether to use the physical spinful time-reversal structure
!> \param matrix_tol Relative structure tolerance
!> \param metric_tol Minimum overlap eigenvalue ratio
!> \param gap_tol Absolute unresolved-gap tolerance
!> \param invariant Resolved class A or AII index
!> \param gap Minimum absolute generalized eigenvalue
!> \param metric_ratio Smallest/largest overlap eigenvalue ratio
!> \param status Same status convention as localizer_evaluate
! **************************************************************************************************
   SUBROUTINE localizer_matrix_evaluate(localizer, s, aii, matrix_tol, metric_tol, gap_tol, &
                                        invariant, gap, metric_ratio, status)
      COMPLEX(KIND=dp), INTENT(IN)                       :: localizer(:, :), s(:, :)
      LOGICAL, INTENT(IN)                                :: aii
      REAL(KIND=dp), INTENT(IN)                          :: matrix_tol, metric_tol, gap_tol
      INTEGER, INTENT(OUT)                               :: invariant, status
      REAL(KIND=dp), INTENT(OUT)                         :: gap, metric_ratio

      COMPLEX(KIND=dp), ALLOCATABLE                      :: a(:, :), b(:, :), work(:)
      INTEGER                                            :: info, n, negative, pf, pf_reference, &
                                                            positive, ldl_positive, ldl_negative
#if defined(__HAS_IEEE_EXCEPTIONS)
      LOGICAL :: halt(SIZE(IEEE_ALL))
#endif
      REAL(KIND=dp)                                      :: metric_min, residual, threshold
      REAL(KIND=dp), ALLOCATABLE                         :: eigenvalues(:), rwork(:), skew(:, :)

      status = -1
      invariant = 0
      gap = 0.0_dp
      metric_ratio = 0.0_dp
      n = SIZE(s, 1)
      IF (n < 1 .OR. SIZE(s, 2) /= n) RETURN
      IF (ANY(SHAPE(localizer) /= [2*n, 2*n])) RETURN
      IF (.NOT. ALL(ieee_is_finite([matrix_tol, metric_tol, gap_tol]))) RETURN
      IF (matrix_tol <= 0.0_dp .OR. matrix_tol >= 1.0_dp) RETURN
      IF (metric_tol <= 0.0_dp .OR. metric_tol >= 1.0_dp .OR. gap_tol < 0.0_dp) RETURN
      IF (.NOT. hermitian(s, matrix_tol) .OR. .NOT. hermitian(localizer, matrix_tol)) RETURN
      IF (aii) THEN
         status = -3
         IF (MOD(n, 2) /= 0) RETURN
         IF (.NOT. self_dual(s, matrix_tol)) RETURN
      END IF

      ALLOCATE (a(n, n), eigenvalues(n), work(MAX(1, 2*n)), rwork(MAX(1, 3*n)))
      a(:, :) = s
      ! Follow CP2K's eigensolvers: LAPACK may use exceptional IEEE intermediates internally.
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_get_halting_mode(IEEE_ALL, halt)
      CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
      CALL zheev('N', 'U', n, a, n, eigenvalues, work, SIZE(work), rwork, info)
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif
      status = ABS(info) + 2
      IF (info /= 0 .OR. .NOT. ALL(ieee_is_finite(eigenvalues))) RETURN
      status = -2
      IF (MAXVAL(eigenvalues) <= 0.0_dp) RETURN
      metric_ratio = MINVAL(eigenvalues)/MAXVAL(eigenvalues)
      IF (metric_ratio <= metric_tol) RETURN
      metric_min = MINVAL(eigenvalues)
      DEALLOCATE (a, eigenvalues, work, rwork)

      ALLOCATE (a(2*n, 2*n), b(2*n, 2*n))
      ALLOCATE (eigenvalues(2*n), work(4*n), rwork(6*n))
      status = -1
      IF (.NOT. ALL(ieee_is_finite(REAL(localizer, dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(localizer)))) RETURN
      a(:, :) = localizer
      b(:, :) = CMPLX(0.0_dp, 0.0_dp, dp)
      b(:n, :n) = s
      b(n + 1:, n + 1:) = s
      ! Generalized eigenvalues, unlike eigenvalues of raw AO matrices, are basis invariant.
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_get_halting_mode(IEEE_ALL, halt)
      CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
      CALL zhegv(1, 'N', 'U', 2*n, a, 2*n, b, 2*n, eigenvalues, work, SIZE(work), rwork, info)
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif
      status = ABS(info) + 2
      IF (info /= 0 .OR. .NOT. ALL(ieee_is_finite(eigenvalues))) RETURN
      gap = MINVAL(ABS(eigenvalues))
      ! Bound roundoff and triangle inconsistencies after conversion to the orthonormal metric.
      residual = REAL(2*n, dp)*MAXVAL(ABS(localizer - CONJG(TRANSPOSE(localizer))))/metric_min + &
                 REAL(n, dp)*MAXVAL(ABS(eigenvalues))*MAXVAL(ABS(s - CONJG(TRANSPOSE(s))))/metric_min
      threshold = MAX(gap_tol, residual + 100.0_dp*EPSILON(1.0_dp)*REAL(2*n, dp)* &
                      MAX(MAXVAL(ABS(eigenvalues)), MAXVAL(ABS(localizer))/metric_min))
      status = 1
      IF (gap <= threshold) RETURN
      positive = COUNT(eigenvalues > threshold)
      negative = COUNT(eigenvalues < -threshold)
      CALL localizer_inertia(localizer, ldl_positive, ldl_negative, info)
      status = -4
      IF (info /= 0) RETURN
      IF (positive /= ldl_positive .OR. negative /= ldl_negative) RETURN
      IF (.NOT. aii) THEN
         invariant = (ldl_positive - ldl_negative)/2
      ELSE
         status = -3
         IF (positive /= negative) RETURN
         ALLOCATE (skew(2*n, 2*n))
         CALL localizer_aii_skew(localizer, skew, matrix_tol, residual, info)
         IF (info /= 0) RETURN
         status = 1
         IF (gap <= threshold + REAL(2*n, dp)*residual/metric_min) RETURN
         CALL localizer_pfaffian_sign(skew, pf, info)
         status = -4
         IF (info /= 0) RETURN
         ! Fix the orientation against the positive atomic Hamiltonian in the SAME AO metric.
         a(:, :) = CMPLX(0.0_dp, 0.0_dp, dp)
         a(:n, :n) = s
         a(n + 1:, n + 1:) = -s
         CALL localizer_aii_skew(a, skew, matrix_tol, residual, info)
         IF (info /= 0) RETURN
         IF (REAL(2*n, dp)*residual >= metric_min) RETURN
         CALL localizer_pfaffian_sign(skew, pf_reference, info)
         IF (info /= 0) RETURN
         invariant = (1 - pf*pf_reference)/2
      END IF
      status = 0
   END SUBROUTINE localizer_matrix_evaluate

! **************************************************************************************************
!> \brief Count inertia by pivoted Hermitian LDL^dagger, without orthogonalizing the AO basis.
!> \param matrix Hermitian matrix; input is not modified
!> \param positive Number of positive eigenvalues
!> \param negative Number of negative eigenvalues
!> \param status 0 resolved, 1 singular, -1 invalid data, -4 factorization failure
!> \note Pivot magnitudes are not eigenvalues or a metric-aware protection gap.
! **************************************************************************************************
   SUBROUTINE localizer_inertia_complex(matrix, positive, negative, status)
      COMPLEX(KIND=dp), INTENT(IN)                       :: matrix(:, :)
      INTEGER, INTENT(OUT)                               :: positive, negative, status

      COMPLEX(KIND=dp), ALLOCATABLE                      :: a(:, :), work(:)
      INTEGER                                            :: i, info, n
      INTEGER, ALLOCATABLE                               :: pivots(:)
      REAL(KIND=dp)                                      :: scale
#if defined(__HAS_IEEE_EXCEPTIONS)
      LOGICAL :: halt(SIZE(IEEE_ALL))
#endif

      positive = 0
      negative = 0
      status = -1
      n = SIZE(matrix, 1)
      IF (n < 1 .OR. SIZE(matrix, 2) /= n) RETURN
      IF (.NOT. hermitian(matrix, 1.e-10_dp)) RETURN
      ! Never discard a small imaginary term: only exactly real matrices use LDL^T.
      IF (ALL(AIMAG(matrix) == 0.0_dp)) THEN
         CALL localizer_inertia_real(REAL(matrix, dp), positive, negative, status)
         RETURN
      END IF
      scale = MAXVAL(ABS(matrix))
      status = 1
      IF (scale <= TINY(1.0_dp)) RETURN
      ALLOCATE (a(n, n), pivots(n), work(MAX(1, 64*n)))
      a(:, :) = matrix/scale
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_get_halting_mode(IEEE_ALL, halt)
      CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
      CALL zhetrf('L', n, a, n, pivots, work, SIZE(work), info)
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif
      IF (info > 0) RETURN
      status = -4
      IF (info < 0) RETURN
      IF (.NOT. ALL(ieee_is_finite(REAL(a, dp))) .OR. .NOT. ALL(ieee_is_finite(AIMAG(a)))) RETURN
      CALL block_inertia([(REAL(a(i, i), dp), i=1, n)], [(ABS(a(i + 1, i)), i=1, n - 1)], &
                         pivots, positive, negative, status)
   END SUBROUTINE localizer_inertia_complex

! **************************************************************************************************
!> \brief Native real symmetric LDL^T inertia, without converting to complex storage.
!> \param matrix Real symmetric matrix; input is not modified
!> \param positive Number of positive eigenvalues
!> \param negative Number of negative eigenvalues
!> \param status 0 resolved, 1 singular, -1 invalid data, -4 factorization failure
! **************************************************************************************************
   SUBROUTINE localizer_inertia_real(matrix, positive, negative, status)
      REAL(dp), INTENT(IN) :: matrix(:, :)
      INTEGER, INTENT(OUT) :: positive, negative, status
      REAL(dp), ALLOCATABLE :: a(:, :), work(:)
      INTEGER, ALLOCATABLE :: pivots(:)
      INTEGER :: i, info, n
      REAL(dp) :: scale
#if defined(__HAS_IEEE_EXCEPTIONS)
      LOGICAL :: halt(SIZE(IEEE_ALL))
#endif

      positive = 0
      negative = 0
      status = -1
      n = SIZE(matrix, 1)
      IF (n < 1 .OR. SIZE(matrix, 2) /= n) RETURN
      IF (.NOT. ALL(ieee_is_finite(matrix))) RETURN
      scale = MAXVAL(ABS(matrix))
      IF (MAXVAL(ABS(matrix - TRANSPOSE(matrix))) > 1.e-10_dp*MAX(TINY(1.0_dp), scale)) RETURN
      status = 1
      IF (scale <= TINY(1.0_dp)) RETURN
      ALLOCATE (a(n, n), pivots(n), work(MAX(1, 64*n)))
      a(:, :) = matrix/scale
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_get_halting_mode(IEEE_ALL, halt)
      CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
      CALL dsytrf('L', n, a, n, pivots, work, SIZE(work), info)
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif
      IF (info > 0) RETURN
      status = -4
      IF (info < 0 .OR. .NOT. ALL(ieee_is_finite(a))) RETURN
      CALL block_inertia([(a(i, i), i=1, n)], [(ABS(a(i + 1, i)), i=1, n - 1)], &
                         pivots, positive, negative, status)
   END SUBROUTINE localizer_inertia_real

! **************************************************************************************************
!> \brief Count signs of lower Bunch-Kaufman D blocks, shared by real and complex factors.
!> \param diagonal Real diagonal of the factored matrix
!> \param coupling Magnitudes of first subdiagonal entries
!> \param pivots LAPACK lower-triangle pivot encoding
!> \param positive Positive inertia count
!> \param negative Negative inertia count
!> \param status 0 resolved, 1 singular, -4 invalid block structure
!> \note A Hermitian 2x2 block has the same eigenvalues when its offdiagonal is replaced by its modulus.
! **************************************************************************************************
   SUBROUTINE block_inertia(diagonal, coupling, pivots, positive, negative, status)
      REAL(dp), INTENT(IN)                               :: diagonal(:), coupling(:)
      INTEGER, INTENT(IN)                                :: pivots(:)
      INTEGER, INTENT(OUT)                               :: positive, negative, status

      INTEGER                                            :: i, n
      REAL(dp)                                           :: values(2)
      INTERFACE
         SUBROUTINE dlae2(a, b, c, rt1, rt2)
            IMPORT :: dp
            REAL(dp), INTENT(IN) :: a, b, c
            REAL(dp), INTENT(OUT) :: rt1, rt2
         END SUBROUTINE dlae2
      END INTERFACE

      positive = 0
      negative = 0
      status = -4
      n = SIZE(pivots)
      i = 1
      DO WHILE (i <= n)
         IF (pivots(i) == 0) RETURN
         IF (pivots(i) > 0) THEN
            values = [diagonal(i), 0.0_dp]
            IF (values(1) == 0.0_dp) THEN
               status = 1
               RETURN
            END IF
            i = i + 1
         ELSE
            IF (i == n) RETURN
            IF (pivots(i + 1) /= pivots(i)) RETURN
            ! A 2x2 pivot can have zero diagonal and still contribute one sign of each kind.
            CALL dlae2(diagonal(i), coupling(i), diagonal(i + 1), values(1), values(2))
            IF (ANY(values == 0.0_dp)) THEN
               status = 1
               RETURN
            END IF
            i = i + 2
         END IF
         positive = positive + COUNT(values > 0.0_dp)
         negative = negative + COUNT(values < 0.0_dp)
      END DO
      status = 0
   END SUBROUTINE block_inertia

! **************************************************************************************************
!> \brief Test finiteness and Hermiticity before passing only one triangle to LAPACK.
!> \param a Square matrix
!> \param tolerance Relative tolerance
!> \return Whether the matrix is Hermitian and finite
! **************************************************************************************************
   LOGICAL FUNCTION hermitian(a, tolerance)
      COMPLEX(KIND=dp), INTENT(IN)                       :: a(:, :)
      REAL(KIND=dp), INTENT(IN)                          :: tolerance

      hermitian = .FALSE.
      IF (.NOT. ALL(ieee_is_finite(REAL(a, dp))) .OR. .NOT. ALL(ieee_is_finite(AIMAG(a)))) RETURN
      hermitian = MAXVAL(ABS(a - CONJG(TRANSPOSE(a)))) <= tolerance*MAX(TINY(1.0_dp), MAXVAL(ABS(a)))
   END FUNCTION hermitian

! **************************************************************************************************
!> \brief Check J conj(A) J^dagger=A for J=(0,I;-I,0), without assuming spin conservation.
!> \param a Even-dimensional matrix in the fixed up/down ordering
!> \param tolerance Relative tolerance
!> \return Whether physical fermionic time reversal is respected
! **************************************************************************************************
   LOGICAL FUNCTION self_dual(a, tolerance)
      COMPLEX(KIND=dp), INTENT(IN)                       :: a(:, :)
      REAL(KIND=dp), INTENT(IN)                          :: tolerance

      INTEGER                                            :: n
      REAL(KIND=dp)                                      :: scale

      n = SIZE(a, 1)/2
      scale = MAX(TINY(1.0_dp), MAXVAL(ABS(a)))
      self_dual = MAXVAL(ABS(a(:n, :n) - CONJG(a(n + 1:, n + 1:)))) <= tolerance*scale .AND. &
                  MAXVAL(ABS(a(:n, n + 1:) + CONJG(a(n + 1:, :n)))) <= tolerance*scale
   END FUNCTION self_dual

! **************************************************************************************************
!> \brief Convert the class AII localizer to a real skew matrix using Q=(I-i C)/sqrt(2).
!> \param localizer Hermitian localizer; auxiliary index outside physical up/down index
!> \param skew Real part of i Q^dagger L Q
!> \param tolerance Relative imaginary/skew residual tolerance
!> \param residual Entrywise bound for the discarded real and symmetric components
!> \param status Zero on success; negative if the AII structure is violated
! **************************************************************************************************
   SUBROUTINE localizer_aii_skew(localizer, skew, tolerance, residual, status)
      COMPLEX(KIND=dp), INTENT(IN)                       :: localizer(:, :)
      REAL(KIND=dp), INTENT(OUT)                         :: skew(:, :)
      REAL(KIND=dp), INTENT(IN)                          :: tolerance
      REAL(KIND=dp), INTENT(OUT)                         :: residual
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(KIND=dp)                                   :: value
      INTEGER                                            :: i, j, m, n, pi, pj, si, sj
      REAL(KIND=dp)                                      :: antisymmetric, scale, skew_residual

      m = SIZE(localizer, 1)
      status = -1
      residual = HUGE(1.0_dp)
      IF (m < 4 .OR. MOD(m, 4) /= 0) RETURN
      IF (SIZE(localizer, 2) /= m .OR. ANY(SHAPE(skew) /= [m, m])) RETURN
      IF (.NOT. ieee_is_finite(tolerance)) RETURN
      IF (tolerance <= 0.0_dp .OR. tolerance >= 1.0_dp) RETURN
      IF (.NOT. ALL(ieee_is_finite(REAL(localizer, dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(localizer)))) RETURN
      n = m/4
      residual = 0.0_dp
      ! C is a signed permutation: Q^dagger L Q = (L + CLC + i(CL-LC))/2.
      DO j = 1, m
         pj = (3 - (j - 1)/n)*n + MODULO(j - 1, n) + 1
         sj = MERGE(1, -1, j <= n .OR. j > 3*n)
         DO i = 1, m
            pi = (3 - (i - 1)/n)*n + MODULO(i - 1, n) + 1
            si = MERGE(1, -1, i <= n .OR. i > 3*n)
            value = 0.5_dp*(localizer(i, j) + si*sj*localizer(pi, pj) + &
                            CMPLX(0.0_dp, 1.0_dp, dp)*(si*localizer(pi, j) - sj*localizer(i, pj)))
            IF (.NOT. ieee_is_finite(REAL(value, dp)) .OR. .NOT. ieee_is_finite(AIMAG(value))) THEN
               residual = HUGE(1.0_dp)
               RETURN
            END IF
            residual = MAX(residual, ABS(REAL(value, dp)))
            skew(i, j) = -AIMAG(value)
         END DO
      END DO
      scale = MAX(TINY(1.0_dp), MAXVAL(ABS(localizer)))
      IF (residual > tolerance*scale) RETURN
      skew_residual = 0.0_dp
      DO j = 1, m
         DO i = 1, j
            skew_residual = MAX(skew_residual, ABS(0.5_dp*skew(i, j) + 0.5_dp*skew(j, i)))
         END DO
      END DO
      IF (skew_residual > 0.5_dp*tolerance*scale) RETURN
      ! Project only validated roundoff; callers include this bound in the gap test.
      residual = residual + skew_residual
      DO j = 1, m
         DO i = 1, j
            antisymmetric = 0.5_dp*skew(i, j) - 0.5_dp*skew(j, i)
            skew(i, j) = antisymmetric
            skew(j, i) = -antisymmetric
         END DO
      END DO
      status = 0
   END SUBROUTINE localizer_aii_skew

! **************************************************************************************************
!> \brief Dense pivoted skew elimination, accumulating only the sign of the Pfaffian.
!> \param matrix Real skew-symmetric matrix; input is not modified
!> \param pf_sign Sign, or zero if singular/unresolved
!> \param status Zero on success, -1 invalid data, -2 numerically singular factorization
!> \note Each simultaneous row/column interchange changes the Pfaffian sign, not its square.
! **************************************************************************************************
   SUBROUTINE localizer_pfaffian_sign(matrix, pf_sign, status)
      REAL(KIND=dp), INTENT(IN)                          :: matrix(:, :)
      INTEGER, INTENT(OUT)                               :: pf_sign, status

      INTEGER                                            :: i, j, k, n, p, sign_accum
      REAL(KIND=dp)                                      :: pivot, scale
      REAL(KIND=dp), ALLOCATABLE                         :: a(:, :), row(:)

      n = SIZE(matrix, 1)
      pf_sign = 0
      status = -1
      IF (n < 2 .OR. MOD(n, 2) /= 0 .OR. SIZE(matrix, 2) /= n) RETURN
      IF (.NOT. ALL(ieee_is_finite(matrix))) RETURN
      scale = MAXVAL(ABS(matrix))
      IF (MAXVAL(ABS(matrix + TRANSPOSE(matrix))) > 100.0_dp*EPSILON(1.0_dp)*scale) RETURN
      status = -2
      IF (scale <= TINY(1.0_dp)) RETURN
      ALLOCATE (a(n, n), row(n))
      a(:, :) = matrix/scale
      sign_accum = 1
      DO k = 1, n - 1, 2
         p = k + MAXLOC(ABS(a(k, k + 1:)), DIM=1)
         IF (p /= k + 1) THEN
            row(:) = a(k + 1, :)
            a(k + 1, :) = a(p, :)
            a(p, :) = row
            row(:) = a(:, k + 1)
            a(:, k + 1) = a(:, p)
            a(:, p) = row
            sign_accum = -sign_accum
         END IF
         pivot = a(k, k + 1)
         IF (ABS(pivot) <= 100.0_dp*EPSILON(1.0_dp)*MAXVAL(ABS(a(k:, k:)))) RETURN
         IF (pivot < 0.0_dp) sign_accum = -sign_accum
         DO j = k + 2, n
            DO i = k + 2, j - 1
               a(i, j) = a(i, j) + (a(k + 1, i)*a(k, j) - a(k, i)*a(k + 1, j))/pivot
               a(j, i) = -a(i, j)
            END DO
         END DO
         IF (.NOT. ALL(ieee_is_finite(a(k:, k:)))) RETURN
      END DO
      pf_sign = sign_accum
      status = 0
   END SUBROUTINE localizer_pfaffian_sign

END MODULE spectral_localizer
