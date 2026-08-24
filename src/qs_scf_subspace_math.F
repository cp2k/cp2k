!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Small, matrix-free mathematical kernels used by ADIIS.
!>
!>      The quadratic model uses the convention
!>
!>        f(c) = 1/2 c^T hessian c + linear^T c,
!>
!>      with c on the probability simplex.  Only the symmetric part of the
!>      supplied Hessian contributes to the objective.
! **************************************************************************************************
MODULE qs_scf_subspace_math

   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp

   IMPLICIT NONE

   PRIVATE

   INTEGER, PARAMETER, PUBLIC :: simplex_qp_success = 0
   INTEGER, PARAMETER, PUBLIC :: simplex_qp_invalid_shape = 1
   INTEGER, PARAMETER, PUBLIC :: simplex_qp_nonfinite_input = 2
   INTEGER, PARAMETER, PUBLIC :: simplex_qp_iteration_limit = 3
   INTEGER, PARAMETER, PUBLIC :: simplex_qp_pairwise_stationary = 4

   INTEGER, PARAMETER :: max_exact_dimension = 12

   PUBLIC :: qs_scf_subspace_build_adiis_model, &
             simplex_quadratic_minimize

CONTAINS

! **************************************************************************************************
!> \brief Minimize a quadratic model over the probability simplex.
!> \param hessian Hessian in f(c) = 1/2 c^T hessian c + linear^T c.
!> \param linear Linear part of the quadratic model.
!> \param coeff Resulting non-negative coefficients, normalized to sum to one.
!> \param objective Objective value at coeff.
!> \param status Completion status. Invalid/nonfinite inputs use the preferred-point fallback;
!>      iteration limits return the best legal point found.
!> \param preferred_index History entry preferred for fallbacks and numerical ties; defaults to n.
!>
!>      SCF histories are normally very small.  Up to max_exact_dimension, all
!>      active faces are enumerated.  A global minimizer of a quadratic over a
!>      simplex is a stationary point in the relative interior of one of these
!>      faces (a singular stationary face always contains an equivalent boundary
!>      solution).  Larger problems use multistart, simplex-preserving pairwise
!>      line searches and return simplex_qp_pairwise_stationary.  This only proves
!>      that no improving two-coordinate transfer was found; for an indefinite
!>      quadratic it is a heuristic result, not a general local/global optimum.
! **************************************************************************************************
   SUBROUTINE simplex_quadratic_minimize(hessian, linear, coeff, objective, status, preferred_index)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: hessian
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: linear
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: coeff
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: objective
      INTEGER, INTENT(OUT), OPTIONAL                     :: status
      INTEGER, INTENT(IN), OPTIONAL                      :: preferred_index

      INTEGER                                            :: local_status, n, preferred
      REAL(KIND=dp)                                      :: best_value, model_scale
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: working_linear
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: symmetric_hessian, working_hessian

      n = SIZE(linear)
      preferred = n
      local_status = simplex_qp_success
      best_value = HUGE(0.0_dp)

      coeff = 0.0_dp
      IF (PRESENT(preferred_index)) THEN
         IF (preferred_index >= 1 .AND. preferred_index <= n) THEN
            preferred = preferred_index
         ELSE
            local_status = simplex_qp_invalid_shape
            IF (SIZE(coeff) > 0) coeff(SIZE(coeff)) = 1.0_dp
            IF (PRESENT(objective)) objective = best_value
            IF (PRESENT(status)) status = local_status
            RETURN
         END IF
      END IF
      IF (preferred >= 1 .AND. preferred <= SIZE(coeff)) coeff(preferred) = 1.0_dp

      IF (n < 1 .OR. SIZE(coeff) /= n .OR. SIZE(hessian, 1) /= n .OR. SIZE(hessian, 2) /= n) THEN
         local_status = simplex_qp_invalid_shape
         IF (PRESENT(objective)) objective = best_value
         IF (PRESENT(status)) status = local_status
         RETURN
      END IF

      IF (.NOT. ALL(ieee_is_finite(hessian)) .OR. .NOT. ALL(ieee_is_finite(linear))) THEN
         local_status = simplex_qp_nonfinite_input
         IF (PRESENT(objective)) objective = best_value
         IF (PRESENT(status)) status = local_status
         RETURN
      END IF

      ALLOCATE (symmetric_hessian(n, n), working_hessian(n, n), working_linear(n))
      symmetric_hessian(:, :) = 0.5_dp*hessian + 0.5_dp*TRANSPOSE(hessian)

      ! Remove terms that are constant on the simplex, then normalize the remaining model.
      ! Besides improving the KKT conditioning, this makes coefficient selection invariant to
      working_hessian(:, :) = symmetric_hessian
      working_linear(:) = linear
      model_scale = MAX(MAXVAL(ABS(working_hessian)), MAXVAL(ABS(working_linear)))
      ! Scale first only at extreme magnitudes where subtracting opposite-signed finite
      ! values could overflow.  At ordinary scales, removing the gauge first preserves
      IF (model_scale > 0.25_dp*HUGE(0.0_dp)) THEN
         working_hessian(:, :) = working_hessian/model_scale
         working_linear(:) = working_linear/model_scale
      END IF
      ! Subtracting one scalar from every Hessian element is a H_nn*1*1^T gauge
      ! shift and is therefore constant when sum(coeff)=1.
      working_hessian(:, :) = working_hessian - working_hessian(n, n)
      working_linear(:) = working_linear - working_linear(n)
      model_scale = MAX(MAXVAL(ABS(working_hessian)), MAXVAL(ABS(working_linear)))
      IF (model_scale > 0.0_dp) THEN
         working_hessian(:, :) = working_hessian/model_scale
         working_linear(:) = working_linear/model_scale
      END IF

      CALL best_vertex(working_hessian, working_linear, preferred, coeff, best_value)

      IF (n <= max_exact_dimension) THEN
         CALL enumerate_active_faces(working_hessian, working_linear, preferred, coeff, best_value)
      ELSE
         CALL pairwise_minimize(working_hessian, working_linear, preferred, coeff, best_value, local_status)
      END IF

      ! Remove harmless round-off at active constraints before returning.
      coeff = MAX(coeff, 0.0_dp)
      coeff = coeff/SUM(coeff)
      best_value = quadratic_value(symmetric_hessian, linear, coeff)

      IF (PRESENT(objective)) objective = best_value
      IF (PRESENT(status)) status = local_status

   END SUBROUTINE simplex_quadratic_minimize

! **************************************************************************************************
!> \brief Construct the canonical ADIIS model relative to the newest history entry.
!> \param pf_metric Matrix T_ij=Tr(P_i F_j), already summed over spin/k-points as needed.
!> \param newest Index of the reference history entry.
!> \param hessian Hessian of the canonical simplex quadratic model.
!> \param linear Linear term of the canonical simplex quadratic model.
!> \param valid Whether all dimensions and input values were valid.
! **************************************************************************************************
   SUBROUTINE qs_scf_subspace_build_adiis_model(pf_metric, newest, hessian, linear, valid)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: pf_metric
      INTEGER, INTENT(IN)                                :: newest
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: hessian
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: linear
      LOGICAL, INTENT(OUT), OPTIONAL                     :: valid

      INTEGER                                            :: i, j, n
      LOGICAL                                            :: local_valid
      REAL(KIND=dp)                                      :: model_ij, model_ji, tnn

      n = SIZE(linear)
      hessian = 0.0_dp
      linear = 0.0_dp

      local_valid = n > 0 .AND. newest >= 1 .AND. newest <= n .AND. &
                    SIZE(pf_metric, 1) == n .AND. SIZE(pf_metric, 2) == n .AND. &
                    SIZE(hessian, 1) == n .AND. SIZE(hessian, 2) == n
      IF (local_valid) local_valid = ALL(ieee_is_finite(pf_metric))
      IF (.NOT. local_valid) THEN
         IF (PRESENT(valid)) valid = .FALSE.
         RETURN
      END IF

      tnn = pf_metric(newest, newest)
      DO i = 1, n
         linear(i) = 2.0_dp*(pf_metric(i, newest) - tnn)
      END DO

      DO j = 1, n
         DO i = 1, n
            model_ij = (pf_metric(i, j) - pf_metric(i, newest)) + &
                       (tnn - pf_metric(newest, j))
            model_ji = (pf_metric(j, i) - pf_metric(j, newest)) + &
                       (tnn - pf_metric(newest, i))
            hessian(i, j) = model_ij + model_ji
         END DO
      END DO

      IF (PRESENT(valid)) valid = .TRUE.

   END SUBROUTINE qs_scf_subspace_build_adiis_model

! **************************************************************************************************
!> \brief Select the best simplex vertex, preferring the newest one on ties.
!> \param hessian Symmetric quadratic-model Hessian.
!> \param linear Quadratic-model linear term.
!> \param preferred History index favored when vertex values are tied.
!> \param coeff Coefficients of the selected vertex.
!> \param best_value Objective value at the selected vertex.
! **************************************************************************************************
   SUBROUTINE best_vertex(hessian, linear, preferred, coeff, best_value)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: hessian
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: linear
      INTEGER, INTENT(IN)                                :: preferred
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: coeff
      REAL(KIND=dp), INTENT(OUT)                         :: best_value

      INTEGER                                            :: i, n
      REAL(KIND=dp)                                      :: value

      n = SIZE(linear)
      coeff = 0.0_dp
      coeff(preferred) = 1.0_dp
      best_value = 0.5_dp*hessian(preferred, preferred) + linear(preferred)

      DO i = n, 1, -1
         IF (i == preferred) CYCLE
         value = 0.5_dp*hessian(i, i) + linear(i)
         IF (strictly_better(value, best_value)) THEN
            coeff = 0.0_dp
            coeff(i) = 1.0_dp
            best_value = value
         END IF
      END DO

   END SUBROUTINE best_vertex

! **************************************************************************************************
!> \brief Enumerate stationary points on every active face of a small simplex.
!> \param hessian Symmetric quadratic-model Hessian.
!> \param linear Quadratic-model linear term.
!> \param preferred History index favored when objective values are tied.
!> \param coeff Best feasible coefficients found so far.
!> \param best_value Objective value corresponding to coeff.
! **************************************************************************************************
   SUBROUTINE enumerate_active_faces(hessian, linear, preferred, coeff, best_value)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: hessian
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: linear
      INTEGER, INTENT(IN)                                :: preferred
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: coeff
      REAL(KIND=dp), INTENT(INOUT)                       :: best_value

      INTEGER                                            :: i, j, k, last_mask, mask, n
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: active
      LOGICAL                                            :: solved
      REAL(KIND=dp)                                      :: feasibility_tolerance, scale, value
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: candidate, rhs, solution
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: kkt

      n = SIZE(linear)
      last_mask = ISHFT(1, n) - 1
      ALLOCATE (active(n), candidate(n))

      DO mask = 1, last_mask
         k = POPCNT(mask)
         IF (k < 2) CYCLE

         j = 0
         DO i = 1, n
            IF (BTEST(mask, i - 1)) THEN
               j = j + 1
               active(j) = i
            END IF
         END DO

         ALLOCATE (kkt(k + 1, k + 1), rhs(k + 1), solution(k + 1))
         kkt = 0.0_dp
         rhs = 0.0_dp

         DO j = 1, k
            rhs(j) = -linear(active(j))
            kkt(j, k + 1) = 1.0_dp
            kkt(k + 1, j) = 1.0_dp
            DO i = 1, k
               kkt(i, j) = hessian(active(i), active(j))
            END DO
         END DO
         rhs(k + 1) = 1.0_dp

         CALL solve_dense_linear(kkt, rhs, solution, solved)
         IF (solved) THEN
            scale = MAX(1.0_dp, MAXVAL(ABS(solution(1:k))))
            feasibility_tolerance = 2048.0_dp*EPSILON(1.0_dp)*REAL(n, KIND=dp)*scale
            IF (MINVAL(solution(1:k)) >= -feasibility_tolerance) THEN
               candidate = 0.0_dp
               DO i = 1, k
                  candidate(active(i)) = MAX(solution(i), 0.0_dp)
               END DO
               IF (SUM(candidate) > 0.0_dp) THEN
                  candidate = candidate/SUM(candidate)
                  value = quadratic_value(hessian, linear, candidate)
                  CALL consider_candidate(candidate, value, preferred, coeff, best_value)
               END IF
            END IF
         END IF

         DEALLOCATE (kkt, rhs, solution)
      END DO

   END SUBROUTINE enumerate_active_faces

! **************************************************************************************************
!> \brief Simplex-preserving pairwise coordinate descent for larger histories.
!> \param hessian Symmetric quadratic-model Hessian.
!> \param linear Quadratic-model linear term.
!> \param preferred History index favored when objective values are tied.
!> \param coeff Best multistart coefficients found.
!> \param best_value Objective value corresponding to coeff.
!> \param status Local-convergence or iteration-limit status.
! **************************************************************************************************
   SUBROUTINE pairwise_minimize(hessian, linear, preferred, coeff, best_value, status)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: hessian
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: linear
      INTEGER, INTENT(IN)                                :: preferred
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: coeff
      REAL(KIND=dp), INTENT(INOUT)                       :: best_value
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: n, start, trial_status
      LOGICAL                                            :: select_trial
      REAL(KIND=dp)                                      :: trial_value
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: trial

      n = SIZE(linear)
      ALLOCATE (trial(n))

      ! Starting from every vertex avoids simple multi-coordinate traps while retaining
      ! deterministic newest-first tie breaking.  A uniform start also covers broad interiors.
      best_value = HUGE(0.0_dp)
      status = simplex_qp_iteration_limit
      DO start = n, 1, -1
         trial = 0.0_dp
         trial(start) = 1.0_dp
         trial_value = quadratic_value(hessian, linear, trial)
         CALL pairwise_descent(hessian, linear, trial, trial_value, trial_status)

         select_trial = strictly_better(trial_value, best_value)
         IF (.NOT. select_trial .AND. numerically_equal(trial_value, best_value)) THEN
            select_trial = prefer_newer(trial, coeff, preferred)
         END IF
         IF (select_trial) THEN
            coeff = trial
            best_value = trial_value
            IF (trial_status == simplex_qp_success) THEN
               status = simplex_qp_pairwise_stationary
            ELSE
               status = trial_status
            END IF
         END IF
      END DO

      trial = 1.0_dp/REAL(n, KIND=dp)
      trial_value = quadratic_value(hessian, linear, trial)
      CALL pairwise_descent(hessian, linear, trial, trial_value, trial_status)
      IF (strictly_better(trial_value, best_value)) THEN
         coeff = trial
         best_value = trial_value
         IF (trial_status == simplex_qp_success) THEN
            status = simplex_qp_pairwise_stationary
         ELSE
            status = trial_status
         END IF
      END IF

   END SUBROUTINE pairwise_minimize

! **************************************************************************************************
!> \brief Run pairwise coordinate descent from one legal simplex point.
!> \param hessian Symmetric quadratic-model Hessian.
!> \param linear Quadratic-model linear term.
!> \param coeff Initial and final legal simplex coefficients.
!> \param value Initial and final objective value.
!> \param status Convergence or iteration-limit status.
! **************************************************************************************************
   SUBROUTINE pairwise_descent(hessian, linear, coeff, value, status)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: hessian
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: linear
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: coeff
      REAL(KIND=dp), INTENT(INOUT)                       :: value
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: best_from, best_to, from, iter, &
                                                            max_iter, n, to
      REAL(KIND=dp)                                      :: alpha, alpha_max, best_alpha, &
                                                            best_delta, curvature, delta, &
                                                            derivative, objective_tolerance
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: gradient

      n = SIZE(linear)
      max_iter = MAX(500, 50*n*n)
      ALLOCATE (gradient(n))

      status = simplex_qp_iteration_limit
      DO iter = 1, max_iter
         gradient(:) = MATMUL(hessian, coeff) + linear
         best_delta = 0.0_dp
         best_alpha = 0.0_dp
         best_from = 0
         best_to = 0

         DO from = 1, n
            alpha_max = coeff(from)
            IF (alpha_max <= 64.0_dp*EPSILON(1.0_dp)) CYCLE
            DO to = 1, n
               IF (to == from) CYCLE

               derivative = gradient(to) - gradient(from)
               curvature = hessian(to, to) + hessian(from, from) - 2.0_dp*hessian(to, from)

               IF (curvature > 0.0_dp) THEN
                  ! Compare before dividing so a tiny positive curvature cannot overflow
                  ! -derivative/curvature under CP2K's Debug FPE traps.
                  IF (derivative >= 0.0_dp) THEN
                     alpha = 0.0_dp
                  ELSE IF (-derivative >= alpha_max*curvature) THEN
                     alpha = alpha_max
                  ELSE
                     alpha = -derivative/curvature
                  END IF
               ELSE
                  alpha = alpha_max
               END IF

               delta = alpha*derivative + 0.5_dp*alpha*alpha*curvature
               IF (delta < best_delta) THEN
                  best_delta = delta
                  best_alpha = alpha
                  best_from = from
                  best_to = to
               END IF
            END DO
         END DO

         objective_tolerance = 16.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, ABS(value))
         IF (best_from == 0 .OR. best_delta >= -objective_tolerance) THEN
            status = simplex_qp_success
            EXIT
         END IF

         coeff(best_from) = coeff(best_from) - best_alpha
         coeff(best_to) = coeff(best_to) + best_alpha
         coeff = MAX(coeff, 0.0_dp)
         coeff = coeff/SUM(coeff)
         value = quadratic_value(hessian, linear, coeff)
      END DO

   END SUBROUTINE pairwise_descent

! **************************************************************************************************
!> \brief Dense Gaussian elimination with partial pivoting for tiny KKT systems.
!> \param matrix Square coefficient matrix.
!> \param rhs Right-hand side vector.
!> \param solution Computed solution, or zeros when the system is rejected.
!> \param solved Whether a finite nonsingular solution was obtained.
! **************************************************************************************************
   SUBROUTINE solve_dense_linear(matrix, rhs, solution, solved)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: matrix
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rhs
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: solution
      LOGICAL, INTENT(OUT)                               :: solved

      INTEGER                                            :: i, k, n, pivot
      REAL(KIND=dp)                                      :: factor, pivot_tolerance, scale, temp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: work_rhs, work_row
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: work_matrix

      n = SIZE(rhs)
      solved = .FALSE.
      solution = 0.0_dp
      IF (SIZE(matrix, 1) /= n .OR. SIZE(matrix, 2) /= n .OR. SIZE(solution) /= n) RETURN

      ALLOCATE (work_matrix(n, n), work_rhs(n), work_row(n))
      work_matrix(:, :) = matrix
      work_rhs(:) = rhs
      scale = MAX(1.0_dp, MAXVAL(ABS(work_matrix)))
      pivot_tolerance = 1024.0_dp*EPSILON(1.0_dp)*REAL(n, KIND=dp)*scale

      DO k = 1, n
         pivot = k - 1 + MAXLOC(ABS(work_matrix(k:n, k)), DIM=1)
         IF (ABS(work_matrix(pivot, k)) <= pivot_tolerance) RETURN

         IF (pivot /= k) THEN
            work_row(:) = work_matrix(k, :)
            work_matrix(k, :) = work_matrix(pivot, :)
            work_matrix(pivot, :) = work_row
            temp = work_rhs(k)
            work_rhs(k) = work_rhs(pivot)
            work_rhs(pivot) = temp
         END IF

         DO i = k + 1, n
            factor = work_matrix(i, k)/work_matrix(k, k)
            work_matrix(i, k) = 0.0_dp
            work_matrix(i, k + 1:n) = work_matrix(i, k + 1:n) - &
                                      factor*work_matrix(k, k + 1:n)
            work_rhs(i) = work_rhs(i) - factor*work_rhs(k)
         END DO
      END DO

      DO i = n, 1, -1
         solution(i) = (work_rhs(i) - DOT_PRODUCT(work_matrix(i, i + 1:n), solution(i + 1:n)))/ &
                       work_matrix(i, i)
      END DO

      solved = ALL(ieee_is_finite(solution))

   END SUBROUTINE solve_dense_linear

! **************************************************************************************************
!> \brief Compare and, when appropriate, retain a finite simplex candidate.
!> \param candidate Feasible simplex candidate.
!> \param value Objective value at candidate.
!> \param preferred History index favored when objective values are tied.
!> \param coeff Best coefficients found so far.
!> \param best_value Objective value corresponding to coeff.
! **************************************************************************************************
   SUBROUTINE consider_candidate(candidate, value, preferred, coeff, best_value)

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: candidate
      REAL(KIND=dp), INTENT(IN)                          :: value
      INTEGER, INTENT(IN)                                :: preferred
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: coeff
      REAL(KIND=dp), INTENT(INOUT)                       :: best_value

      IF (.NOT. ieee_is_finite(value)) RETURN
      IF (strictly_better(value, best_value)) THEN
         coeff = candidate
         best_value = value
      ELSE IF (numerically_equal(value, best_value) .AND. prefer_newer(candidate, coeff, preferred)) THEN
         coeff = candidate
         best_value = value
      END IF

   END SUBROUTINE consider_candidate

! **************************************************************************************************
!> \brief Evaluate the canonical quadratic objective.
!> \param hessian Symmetric quadratic-model Hessian.
!> \param linear Quadratic-model linear term.
!> \param coeff Simplex coefficient vector.
!> \return Canonical quadratic objective value.
! **************************************************************************************************
   PURE FUNCTION quadratic_value(hessian, linear, coeff) RESULT(value)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: hessian
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: linear, coeff
      REAL(KIND=dp)                                      :: value

      value = 0.5_dp*DOT_PRODUCT(coeff, MATMUL(hessian, coeff)) + DOT_PRODUCT(linear, coeff)

   END FUNCTION quadratic_value

! **************************************************************************************************
!> \brief Scale-aware strict comparison used for deterministic candidate selection.
!> \param value Candidate objective value.
!> \param reference Current best objective value.
!> \return Whether value improves reference beyond round-off tolerance.
! **************************************************************************************************
   PURE FUNCTION strictly_better(value, reference) RESULT(better)

      REAL(KIND=dp), INTENT(IN)                          :: value, reference
      LOGICAL                                            :: better

      better = value < reference - 512.0_dp*EPSILON(1.0_dp)* &
               MAX(1.0_dp, ABS(value), ABS(reference))

   END FUNCTION strictly_better

! **************************************************************************************************
!> \brief Scale-aware equality comparison.
!> \param value Candidate objective value.
!> \param reference Current best objective value.
!> \return Whether both values agree within round-off tolerance.
! **************************************************************************************************
   PURE FUNCTION numerically_equal(value, reference) RESULT(equal)

      REAL(KIND=dp), INTENT(IN)                          :: value, reference
      LOGICAL                                            :: equal

      equal = ABS(value - reference) <= 512.0_dp*EPSILON(1.0_dp)* &
              MAX(1.0_dp, ABS(value), ABS(reference))

   END FUNCTION numerically_equal

! **************************************************************************************************
!> \brief Deterministic tie break that first favors the caller's preferred history entry.
!> \param candidate Candidate coefficient vector.
!> \param reference Current coefficient vector.
!> \param preferred History index favored first.
!> \return Whether candidate wins the deterministic tie break.
! **************************************************************************************************
   PURE FUNCTION prefer_newer(candidate, reference, preferred) RESULT(prefer)

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: candidate, reference
      INTEGER, INTENT(IN)                                :: preferred
      LOGICAL                                            :: prefer

      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: tolerance

      prefer = .FALSE.
      tolerance = 512.0_dp*EPSILON(1.0_dp)

      IF (candidate(preferred) > reference(preferred) + tolerance) THEN
         prefer = .TRUE.
         RETURN
      ELSE IF (candidate(preferred) < reference(preferred) - tolerance) THEN
         RETURN
      END IF

      DO i = SIZE(candidate), 1, -1
         IF (i == preferred) CYCLE
         IF (candidate(i) > reference(i) + tolerance) THEN
            prefer = .TRUE.
            RETURN
         ELSE IF (candidate(i) < reference(i) - tolerance) THEN
            RETURN
         END IF
      END DO

   END FUNCTION prefer_newer

END MODULE qs_scf_subspace_math
