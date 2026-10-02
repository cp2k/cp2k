!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Metric-aware quadratic pseudospectra for finite AO operators.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE quadratic_pseudospectrum
   USE ieee_arithmetic, ONLY: ieee_is_finite
#if defined(__HAS_IEEE_EXCEPTIONS)
   USE ieee_exceptions, ONLY: IEEE_ALL, ieee_get_halting_mode, ieee_set_halting_mode
#endif
   USE kinds, ONLY: dp
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: quadratic_dense, quadratic_solve, quadratic_result_type

   TYPE quadratic_result_type
      REAL(dp), ALLOCATABLE :: eigenvalues(:), residuals(:), errors(:, :), means(:, :)
      REAL(dp), ALLOCATABLE :: projected_errors(:, :), leakage(:, :)
      COMPLEX(dp), ALLOCATABLE :: states(:, :), complex_means(:, :)
      INTEGER :: iterations = 0
   END TYPE quadratic_result_type

   INTERFACE quadratic_dense
      MODULE PROCEDURE quadratic_dense_real, quadratic_dense_complex
   END INTERFACE
   INTERFACE quadratic_solve
      MODULE PROCEDURE quadratic_solve_real, quadratic_solve_complex
   END INTERFACE

   ABSTRACT INTERFACE
! **************************************************************************************************
!> \brief Apply a covariant AO observable or the overlap to a block of coefficients.
!> \param OPERATOR Zero for overlap, otherwise the requested observable index
!> \param v Input coefficients
!> \param y Covariant output vectors
!> \param status Zero on success
! **************************************************************************************************
      SUBROUTINE operator_apply(OPERATOR, v, y, status)
         IMPORT dp
         INTEGER, INTENT(IN) :: operator
         COMPLEX(dp), INTENT(IN) :: v(:, :)
         COMPLEX(dp), INTENT(OUT) :: y(:, :)
         INTEGER, INTENT(OUT) :: status
      END SUBROUTINE operator_apply
! **************************************************************************************************
!> \brief Solve the positive AO metric for a block of covariant right hand sides.
!> \param v Right hand sides
!> \param y Contravariant solutions
!> \param status Zero on success
! **************************************************************************************************
      SUBROUTINE metric_solve(v, y, status)
         IMPORT dp
         COMPLEX(dp), INTENT(IN) :: v(:, :)
         COMPLEX(dp), INTENT(OUT) :: y(:, :)
         INTEGER, INTENT(OUT) :: status
      END SUBROUTINE metric_solve
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Dense AO oracle, also providing a dense-metric reference for the iterative solver.
!> \param operators Hermitian observables in a covariant AO basis
!> \param metric Positive AO overlap
!> \param shift Query value for each observable
!> \param weights Nonnegative observable scales
!> \param nstates Number of lowest eigenpairs
!> \param iterative Use the matrix-free Ritz iteration rather than full diagonalization
!> \param max_subspace Maximum trial-space dimension
!> \param max_iter Maximum Ritz iterations
!> \param tolerance Relative eigenpair residual threshold
!> \param metric_tolerance Minimum overlap eigenvalue ratio
!> \param RESULT Eigenpairs, separate observable residuals and expectations
!> \param status 0 success, 1 unconverged, -1 input, -2 metric, -4 numerical failure
! **************************************************************************************************
   SUBROUTINE quadratic_dense_real(operators, metric, shift, weights, nstates, iterative, max_subspace, &
                                   max_iter, tolerance, metric_tolerance, RESULT, status)
      COMPLEX(dp), INTENT(IN)                            :: operators(:, :, :), metric(:, :)
      REAL(dp), INTENT(IN)                               :: shift(:), weights(:)
      INTEGER, INTENT(IN)                                :: nstates
      LOGICAL, INTENT(IN)                                :: iterative
      INTEGER, INTENT(IN)                                :: max_subspace, max_iter
      REAL(dp), INTENT(IN)                               :: tolerance, metric_tolerance
      TYPE(quadratic_result_type), INTENT(OUT)           :: result
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: j

      status = -1
      IF (SIZE(operators, 1) < 1 .OR. SIZE(operators, 1) /= SIZE(operators, 2)) RETURN
      DO j = 1, SIZE(operators, 3)
         IF (.NOT. hermitian(operators(:, :, j))) RETURN
      END DO
      CALL quadratic_dense_complex(operators, metric, CMPLX(shift, 0.0_dp, dp), weights, nstates, iterative, &
                                   max_subspace, max_iter, tolerance, metric_tolerance, RESULT, status)
   END SUBROUTINE quadratic_dense_real

! **************************************************************************************************
!> \brief Dense oracle for general complex queries and possibly non-Hermitian AO operators.
!> \param operators Covariant AO operators, not necessarily normal or Hermitian
!> \param metric Positive AO overlap
!> \param shift Complex query value for each operator
!> \param weights Nonnegative residual scales
!> \param nstates Number of lowest eigenpairs
!> \param iterative Whether to use the bounded Ritz iteration
!> \param max_subspace Maximum trial-space dimension
!> \param max_iter Maximum Ritz iterations
!> \param tolerance Relative eigenpair residual threshold
!> \param metric_tolerance Minimum overlap eigenvalue ratio
!> \param RESULT Eigenpairs, residuals and complex expectations; means holds their real parts
!> \param status 0 success, 1 unconverged, -1 input, -2 metric, -4 numerical failure
!> \param unitary Optional mask for compressions of continuum unitary operators, whose full normal Gram is S
! **************************************************************************************************
   SUBROUTINE quadratic_dense_complex(operators, metric, shift, weights, nstates, iterative, max_subspace, &
                                      max_iter, tolerance, metric_tolerance, RESULT, status, unitary)
      COMPLEX(dp), INTENT(IN) :: operators(:, :, :), metric(:, :), shift(:)
      REAL(dp), INTENT(IN) :: weights(:), tolerance, metric_tolerance
      INTEGER, INTENT(IN) :: nstates, max_subspace, max_iter
      LOGICAL, INTENT(IN) :: iterative
      TYPE(quadratic_result_type), INTENT(OUT) :: result
      INTEGER, INTENT(OUT) :: status
      LOGICAL, INTENT(IN), OPTIONAL :: unitary(:)
      COMPLEX(dp), ALLOCATABLE :: factor(:, :)
      REAL(dp), ALLOCATABLE :: eigenvalues(:)
      INTEGER :: n
#if defined(__HAS_IEEE_EXCEPTIONS)
      LOGICAL :: halt(SIZE(IEEE_ALL))
#endif

      status = -1
      n = SIZE(metric, 1)
      IF (n < 1 .OR. SIZE(metric, 2) /= n) RETURN
      IF (SIZE(shift) < 1 .OR. SIZE(weights) /= SIZE(shift)) RETURN
      IF (ANY(SHAPE(operators) /= [n, n, SIZE(shift)])) RETURN
      IF (.NOT. ieee_is_finite(metric_tolerance)) RETURN
      IF (metric_tolerance <= 0.0_dp .OR. metric_tolerance >= 1.0_dp) RETURN
      IF (.NOT. hermitian(metric)) RETURN
      IF (.NOT. ALL(ieee_is_finite(REAL(operators, dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(operators)))) RETURN
      ALLOCATE (factor(n, n), eigenvalues(n))
      factor(:, :) = metric
      CALL diagonalize(factor, eigenvalues, status)
      IF (status /= 0) RETURN
      status = -2
      IF (eigenvalues(n) <= 0.0_dp) RETURN
      IF (eigenvalues(1) <= metric_tolerance*eigenvalues(n)) RETURN
      factor(:, :) = metric
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_get_halting_mode(IEEE_ALL, halt)
      CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
      CALL zpotrf('L', n, factor, n, status)
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif
      IF (status /= 0) THEN
         status = -2
         RETURN
      END IF
      CALL quadratic_solve(n, shift, weights, nstates, iterative, max_subspace, max_iter, tolerance, &
                           apply, solve, RESULT, status, adjoint, unitary)
   CONTAINS
! **************************************************************************************************
!> \brief Apply a dense covariant operator or the metric.
!> \param OPERATOR Zero for metric, otherwise an observable index
!> \param v Input block
!> \param y Output block
!> \param status Return status
! **************************************************************************************************
      SUBROUTINE apply(OPERATOR, v, y, status)
      INTEGER, INTENT(IN)                                :: operator
      COMPLEX(dp), INTENT(IN)                            :: v(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: y(:, :)
      INTEGER, INTENT(OUT)                               :: status

         status = 0
         IF (OPERATOR == 0) THEN
            y = MATMUL(metric, v)
         ELSE
            y = MATMUL(operators(:, :, OPERATOR), v)
         END IF
      END SUBROUTINE apply
! **************************************************************************************************
!> \brief Apply the actual covariant adjoint, without assuming a normal operator.
!> \param OPERATOR Observable index
!> \param v Input block
!> \param y Adjoint output block
!> \param status Return status
! **************************************************************************************************
      SUBROUTINE adjoint(OPERATOR, v, y, status)
      INTEGER, INTENT(IN)                                :: operator
      COMPLEX(dp), INTENT(IN)                            :: v(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: y(:, :)
      INTEGER, INTENT(OUT)                               :: status

         y = MATMUL(CONJG(TRANSPOSE(operators(:, :, OPERATOR))), v)
         status = 0
      END SUBROUTINE adjoint
! **************************************************************************************************
!> \brief Reuse one Cholesky factorization for AO-metric solutions.
!> \param v Right hand sides
!> \param y Solutions
!> \param status Return status
! **************************************************************************************************
      SUBROUTINE solve(v, y, status)
      COMPLEX(dp), INTENT(IN)                            :: v(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: y(:, :)
      INTEGER, INTENT(OUT)                               :: status

         y = v
#if defined(__HAS_IEEE_EXCEPTIONS)
         CALL ieee_get_halting_mode(IEEE_ALL, halt)
         CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
         CALL zpotrs('L', n, SIZE(v, 2), factor, n, y, n, status)
#if defined(__HAS_IEEE_EXCEPTIONS)
         CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif
         IF (status /= 0) status = -4
      END SUBROUTINE solve
   END SUBROUTINE quadratic_dense_complex

! **************************************************************************************************
!> \brief Solve sum_j w_j^2 A_j S^-1 A_j c = lambda S c without squaring sparse operators.
!> \param n AO/spinor dimension
!> \param shift Query values for the observables
!> \param weights Nonnegative observable scales
!> \param nstates Number of lowest eigenpairs
!> \param iterative Whether to use a bounded restarted block Ritz space
!> \param max_subspace Maximum trial-space size
!> \param max_iter Maximum Ritz iterations
!> \param tolerance Relative S^-1 eigenpair residual threshold
!> \param apply Apply covariant observables or metric (operator zero)
!> \param solve Apply inverse metric using already validated factors
!> \param RESULT Normalized eigenstates and separated physical residuals
!> \param status 0 success, 1 unconverged, -1 invalid input, -4 numerical failure
! **************************************************************************************************
   SUBROUTINE quadratic_solve_real(n, shift, weights, nstates, iterative, max_subspace, max_iter, tolerance, &
                                   apply, solve, RESULT, status)
      INTEGER, INTENT(IN) :: n, nstates, max_subspace, max_iter
      REAL(dp), INTENT(IN) :: shift(:), weights(:), tolerance
      LOGICAL, INTENT(IN) :: iterative
      PROCEDURE(operator_apply) :: apply
      PROCEDURE(metric_solve) :: solve
      TYPE(quadratic_result_type), INTENT(OUT) :: result
      INTEGER, INTENT(OUT) :: status

      CALL quadratic_solve_complex(n, CMPLX(shift, 0.0_dp, dp), weights, nstates, iterative, max_subspace, &
                                   max_iter, tolerance, apply, solve, RESULT, status, apply)
   END SUBROUTINE quadratic_solve_real

! **************************************************************************************************
!> \brief Solve sum_j w_j^2 (A_j-z_j S)^dagger S^-1 (A_j-z_j S) c = lambda S c.
!> \param n AO/spinor dimension
!> \param shift Complex queries for arbitrary, possibly nonnormal operators
!> \param weights Nonnegative residual scales
!> \param nstates Number of lowest eigenpairs
!> \param iterative Whether to use a bounded restarted block Ritz space
!> \param max_subspace Maximum trial-space size
!> \param max_iter Maximum Ritz iterations
!> \param tolerance Relative S^-1 eigenpair residual threshold
!> \param apply Apply covariant operators or the positive metric (operator zero)
!> \param solve Solve with the already validated positive metric
!> \param RESULT Eigenpairs, residuals and complex expectations; means holds their real parts
!> \param status 0 success, 1 unconverged, -1 input, -4 numerical failure
!> \param adjoint Apply A_j^dagger, not the coefficient-space S-adjoint; positive operator indices
!> \param unitary True for compressions of continuum unitary operators; retain full norms and projection leakage
! **************************************************************************************************
   SUBROUTINE quadratic_solve_complex(n, shift, weights, nstates, iterative, max_subspace, max_iter, tolerance, &
                                      apply, solve, RESULT, status, adjoint, unitary)
      INTEGER, INTENT(IN) :: n, nstates, max_subspace, max_iter
      COMPLEX(dp), INTENT(IN) :: shift(:)
      REAL(dp), INTENT(IN) :: weights(:), tolerance
      LOGICAL, INTENT(IN) :: iterative
      PROCEDURE(operator_apply) :: apply, adjoint
      PROCEDURE(metric_solve) :: solve
      TYPE(quadratic_result_type), INTENT(OUT) :: result
      INTEGER, INTENT(OUT) :: status
      LOGICAL, INTENT(IN), OPTIONAL :: unitary(:)
      COMPLEX(dp), ALLOCATABLE :: v(:, :), sv(:, :), qv(:, :), projected(:, :), sc(:, :), qc(:, :), &
                                  r(:, :), z(:, :), trial(:, :), strial(:, :), work(:, :), invwork(:, :), &
                                  restarted(:, :)
      REAL(dp), ALLOCATABLE :: eigenvalues(:)
      REAL(dp) :: square, scale, norm_before, norm_after, full_square, lost_square
      COMPLEX(dp) :: coefficient
      INTEGER :: capacity, m, old, i, j, p, iteration, seed, initial, pass, keep
      LOGICAL :: accepted
      LOGICAL, ALLOCATABLE :: continuum(:)

      status = -1
      IF (n < 1 .OR. nstates < 1 .OR. nstates > n .OR. max_iter < 1) RETURN
      IF (SIZE(shift) < 1 .OR. SIZE(weights) /= SIZE(shift)) RETURN
      ALLOCATE (continuum(SIZE(shift)))
      continuum = .FALSE.
      IF (PRESENT(unitary)) THEN
         IF (SIZE(unitary) /= SIZE(shift)) RETURN
         continuum(:) = unitary
      END IF
      IF (.NOT. ALL(ieee_is_finite([REAL(shift, dp), AIMAG(shift), weights, tolerance]))) RETURN
      IF (ANY(weights < 0.0_dp) .OR. MAXVAL(weights) <= 0.0_dp .OR. tolerance <= 0.0_dp) RETURN
      capacity = n
      IF (iterative) capacity = MIN(n, max_subspace)
      IF (capacity < nstates .OR. (capacity < n .AND. capacity < 2*nstates)) RETURN
      ALLOCATE (v(n, capacity), sv(n, capacity), qv(n, capacity), sc(n, nstates), qc(n, nstates), &
                r(n, nstates), z(n, nstates), trial(n, 1), strial(n, 1), work(n, nstates), invwork(n, nstates))
      ALLOCATE (result%states(n, nstates), result%eigenvalues(nstates), result%residuals(nstates), &
                result%errors(SIZE(shift), nstates), result%means(SIZE(shift), nstates), &
                result%complex_means(SIZE(shift), nstates), &
                result%projected_errors(SIZE(shift), nstates), result%leakage(SIZE(shift), nstates))
      result%leakage = 0.0_dp
      m = 0
      seed = 0
      initial = n
      IF (iterative) initial = MIN(capacity, MAX(nstates + 2, 2*nstates))
      DO j = 1, initial
         trial = CMPLX(0.0_dp, 0.0_dp, dp)
         IF (.NOT. iterative) THEN
            trial(j, 1) = 1.0_dp
         ELSE
            seed = seed + 1
            CALL seed_vector()
         END IF
         CALL append_vector(accepted)
         IF (status /= 0) RETURN
         IF (.NOT. accepted) THEN
            status = -4
            RETURN
         END IF
      END DO
      CALL apply_q(v(:, :m), qv(:, :m), status)
      IF (status /= 0) RETURN
      DO iteration = 1, max_iter
         result%iterations = iteration
         ALLOCATE (projected(m, m), eigenvalues(m))
         projected(:, :) = MATMUL(CONJG(TRANSPOSE(v(:, :m))), qv(:, :m))
         scale = MAX(1.0_dp, MAXVAL(ABS(projected)))
         IF (MAXVAL(ABS(projected - CONJG(TRANSPOSE(projected)))) > 1.e-8_dp*scale) THEN
            status = -4
            RETURN
         END IF
         projected(:, :) = 0.5_dp*(projected + CONJG(TRANSPOSE(projected)))
         CALL diagonalize(projected, eigenvalues, status)
         IF (status /= 0) RETURN
         result%eigenvalues(:) = eigenvalues(:nstates)
         result%states(:, :) = MATMUL(v(:, :m), projected(:, :nstates))
         sc(:, :) = MATMUL(sv(:, :m), projected(:, :nstates))
         qc(:, :) = MATMUL(qv(:, :m), projected(:, :nstates))
         DO pass = 1, 2
            IF (pass == 2) THEN
               ! Cached metric/operator images can drift through repeated Ritz restarts.
               CALL apply(0, result%states, sc, status)
               IF (status /= 0) RETURN
               DO j = 1, nstates
                  square = REAL(DOT_PRODUCT(result%states(:, j), sc(:, j)), dp)
                  IF (.NOT. ieee_is_finite(square) .OR. square <= 0.0_dp) THEN
                     status = -4
                     RETURN
                  END IF
                  result%states(:, j) = result%states(:, j)/SQRT(square)
               END DO
               CALL apply(0, result%states, sc, status)
               IF (status /= 0) RETURN
               CALL apply_q(result%states, qc, status)
               IF (status /= 0) RETURN
            END IF
            DO j = 1, nstates
               r(:, j) = qc(:, j) - result%eigenvalues(j)*sc(:, j)
            END DO
            CALL solve(r, z, status)
            IF (status /= 0) RETURN
            DO j = 1, nstates
               square = REAL(DOT_PRODUCT(r(:, j), z(:, j)), dp)
               IF (.NOT. ieee_is_finite(square) .OR. square < -EPSILON(square)*scale) THEN
                  status = -4
                  RETURN
               END IF
               result%residuals(j) = SQRT(MAX(0.0_dp, square))
            END DO
            IF (ANY(result%residuals > tolerance*MAX(1.0_dp, ABS(result%eigenvalues)))) EXIT
         END DO
         IF (ALL(result%residuals <= tolerance*MAX(1.0_dp, ABS(result%eigenvalues)))) THEN
            DEALLOCATE (projected, eigenvalues)
            EXIT
         END IF
         IF (.NOT. iterative .OR. m == n .OR. iteration == max_iter) THEN
            status = 1
            RETURN
         END IF
         ! Fill remaining slots with a partial residual block before restarting.
         IF (m == capacity) THEN
            ! Retain nearby Ritz vectors across restarts, important for clustered spectra.
            keep = MIN(m, MAX(nstates, capacity/2))
            ALLOCATE (restarted(n, keep))
            restarted(:, :) = MATMUL(v(:, :m), projected(:, :keep))
            v(:, :keep) = restarted
            restarted(:, :) = MATMUL(sv(:, :m), projected(:, :keep))
            sv(:, :keep) = restarted
            restarted(:, :) = MATMUL(qv(:, :m), projected(:, :keep))
            qv(:, :keep) = restarted
            DEALLOCATE (restarted)
            m = keep
         END IF
         DEALLOCATE (projected, eigenvalues)
         old = m
         DO j = 1, nstates
            IF (m == capacity) EXIT
            trial(:, 1) = z(:, j)
            CALL append_vector(accepted)
            IF (status /= 0) RETURN
         END DO
         IF (m == old) THEN
            seed = seed + 1
            CALL seed_vector()
            CALL append_vector(accepted)
            IF (status /= 0) RETURN
            IF (.NOT. accepted) THEN
               status = 1
               RETURN
            END IF
         END IF
         CALL apply_q(v(:, old + 1:m), qv(:, old + 1:m), status)
         IF (status /= 0) RETURN
      END DO
      IF (ANY(result%eigenvalues < -tolerance*MAX(1.0_dp, MAXVAL(ABS(result%eigenvalues))))) THEN
         status = -4
         RETURN
      END IF
      result%eigenvalues = MAX(0.0_dp, result%eigenvalues)
      DO p = 1, SIZE(shift)
         CALL apply(p, result%states, work, status)
         IF (status /= 0) RETURN
         DO j = 1, nstates
            result%complex_means(p, j) = DOT_PRODUCT(result%states(:, j), work(:, j))
            result%means(p, j) = REAL(result%complex_means(p, j), dp)
         END DO
         work(:, :) = work - shift(p)*sc
         CALL solve(work, invwork, status)
         IF (status /= 0) RETURN
         DO j = 1, nstates
            square = REAL(DOT_PRODUCT(work(:, j), invwork(:, j)), dp)
            IF (.NOT. ieee_is_finite(square) .OR. square < -tolerance) THEN
               status = -4
               RETURN
            END IF
            result%projected_errors(p, j) = SQRT(MAX(0.0_dp, square))
            IF (continuum(p)) THEN
               full_square = (1.0_dp + ABS(shift(p))**2)*REAL(DOT_PRODUCT(result%states(:, j), sc(:, j)), dp) &
                             - 2.0_dp*REAL(CONJG(shift(p))*result%complex_means(p, j), dp)
               lost_square = full_square - square
               scale = MAX(1.0_dp, ABS(full_square), ABS(square))
               IF (.NOT. ieee_is_finite(full_square) .OR. full_square < -tolerance*scale .OR. &
                   lost_square < -tolerance*scale) THEN
                  status = -4
                  RETURN
               END IF
               result%leakage(p, j) = SQRT(MAX(0.0_dp, lost_square))
               square = full_square
            END IF
            result%errors(p, j) = SQRT(MAX(0.0_dp, square))
         END DO
      END DO
      DO j = 1, nstates
         square = SUM((weights*result%errors(:, j))**2)
         IF (ABS(square - result%eigenvalues(j)) > 10.0_dp*tolerance*MAX(1.0_dp, square)) THEN
            status = -4
            RETURN
         END IF
      END DO
      status = 0
   CONTAINS
! **************************************************************************************************
!> \brief Deterministic complex trial block, independent of compiler random-number generators.
! **************************************************************************************************
      SUBROUTINE seed_vector()
      INTEGER                                            :: k

         DO k = 1, n
            trial(k, 1) = CMPLX(SIN(REAL(k, dp)*REAL(seed + 1, dp)*1.4142135623730951_dp), &
                                COS(REAL(k, dp)*REAL(seed + 2, dp)*1.7320508075688772_dp), dp)
         END DO
      END SUBROUTINE seed_vector
! **************************************************************************************************
!> \brief Twice-reorthogonalized S-metric expansion, rejecting linearly dependent vectors.
!> \param accepted Whether a new basis vector was appended
! **************************************************************************************************
      SUBROUTINE append_vector(accepted)
      LOGICAL, INTENT(OUT)                               :: accepted

         accepted = .FALSE.
         CALL apply(0, trial, strial, status)
         IF (status /= 0) RETURN
         norm_before = REAL(DOT_PRODUCT(trial(:, 1), strial(:, 1)), dp)
         IF (.NOT. ieee_is_finite(norm_before) .OR. norm_before < 0.0_dp) THEN
            status = -4
            RETURN
         END IF
         DO pass = 1, 2
            DO i = 1, m
               coefficient = DOT_PRODUCT(v(:, i), strial(:, 1))
               trial(:, 1) = trial(:, 1) - coefficient*v(:, i)
               strial(:, 1) = strial(:, 1) - coefficient*sv(:, i)
            END DO
         END DO
         norm_after = REAL(DOT_PRODUCT(trial(:, 1), strial(:, 1)), dp)
         IF (.NOT. ieee_is_finite(norm_after)) THEN
            status = -4
            RETURN
         END IF
         IF (norm_after <= 1.e-24_dp*MAX(TINY(1.0_dp), norm_before)) RETURN
         m = m + 1
         v(:, m) = trial(:, 1)/SQRT(norm_after)
         sv(:, m) = strial(:, 1)/SQRT(norm_after)
         accepted = .TRUE.
      END SUBROUTINE append_vector
! **************************************************************************************************
!> \brief Matrix-free covariant normal products, reusing the positive metric factorization.
!> \param block Input AO coefficients
!> \param image Covariant Q times the input
!> \param status Return status
! **************************************************************************************************
      SUBROUTINE apply_q(block, image, status)
      COMPLEX(dp), INTENT(IN)                            :: block(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: image(:, :)
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(dp), ALLOCATABLE                           :: ab(:, :), sb(:, :), second(:, :), &
                                                            solved(:, :), ss(:, :)
      INTEGER                                            :: nrhs, op

         nrhs = SIZE(block, 2)
         ALLOCATE (sb(n, nrhs), ab(n, nrhs), solved(n, nrhs), second(n, nrhs), ss(n, nrhs))
         CALL apply(0, block, sb, status)
         IF (status /= 0) RETURN
         image = CMPLX(0.0_dp, 0.0_dp, dp)
         DO op = 1, SIZE(shift)
            IF (weights(op) == 0.0_dp) CYCLE
            CALL apply(op, block, ab, status)
            IF (status /= 0) RETURN
            IF (continuum(op)) THEN
               CALL adjoint(op, block, second, status)
               IF (status /= 0) RETURN
               image = image + weights(op)**2*((1.0_dp + ABS(shift(op))**2)*sb &
                                               - CONJG(shift(op))*ab - shift(op)*second)
               CYCLE
            END IF
            ab(:, :) = ab - shift(op)*sb
            CALL solve(ab, solved, status)
            IF (status /= 0) RETURN
            CALL adjoint(op, solved, second, status)
            IF (status /= 0) RETURN
            CALL apply(0, solved, ss, status)
            IF (status /= 0) RETURN
            image = image + weights(op)**2*(second - CONJG(shift(op))*ss)
         END DO
      END SUBROUTINE apply_q
   END SUBROUTINE quadratic_solve_complex

! **************************************************************************************************
!> \brief Hermiticity and finiteness check for the dense reference.
!> \param matrix Input square matrix
!> \return Whether the matrix is finite Hermitian
! **************************************************************************************************
   FUNCTION hermitian(matrix) RESULT(valid)
      COMPLEX(dp), INTENT(IN) :: matrix(:, :)
      LOGICAL :: valid
      valid = .FALSE.
      IF (.NOT. ALL(ieee_is_finite(REAL(matrix, dp))) .OR. .NOT. ALL(ieee_is_finite(AIMAG(matrix)))) RETURN
      valid = MAXVAL(ABS(matrix - CONJG(TRANSPOSE(matrix)))) <= 1.e-10_dp*MAX(TINY(1.0_dp), MAXVAL(ABS(matrix)))
   END FUNCTION hermitian

! **************************************************************************************************
!> \brief Small dense Hermitian eigensolve shared by the oracle and Ritz projection.
!> \param matrix Overwritten with eigenvectors
!> \param values Ascending eigenvalues
!> \param status Return status
! **************************************************************************************************
   SUBROUTINE diagonalize(matrix, values, status)
      COMPLEX(dp), INTENT(INOUT) :: matrix(:, :)
      REAL(dp), INTENT(OUT) :: values(:)
      INTEGER, INTENT(OUT) :: status
      COMPLEX(dp), ALLOCATABLE :: workspace(:)
      REAL(dp), ALLOCATABLE :: rwork(:)
      INTEGER :: n
#if defined(__HAS_IEEE_EXCEPTIONS)
      LOGICAL :: halt(SIZE(IEEE_ALL))
#endif
      n = SIZE(values)
      ALLOCATE (workspace(MAX(1, 2*n)), rwork(MAX(1, 3*n - 2)))
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_get_halting_mode(IEEE_ALL, halt)
      CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
      CALL zheev('V', 'L', n, matrix, n, values, workspace, SIZE(workspace), rwork, status)
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif
      IF (status /= 0 .OR. .NOT. ALL(ieee_is_finite(values))) status = -4
   END SUBROUTINE diagonalize
END MODULE quadratic_pseudospectrum
