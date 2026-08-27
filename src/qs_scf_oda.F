!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Safeguarded optimal damping of raw Roothaan SCF steps.
! **************************************************************************************************
MODULE qs_scf_oda

   USE bibliography,                    ONLY: Cances2000,&
                                              Herbst2022,&
                                              cite_reference
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_get_info,&
                                              dbcsr_p_type
   USE cp_dbcsr_contrib,                ONLY: dbcsr_dot
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_rho_types,                    ONLY: qs_rho_type
   USE qs_scf_loop_utils,               ONLY: qs_scf_rho_update
   USE qs_scf_types,                    ONLY: qs_scf_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_oda'
   CHARACTER(len=*), PARAMETER, PRIVATE :: status_no_acceptable = "NO_ACCEPTABLE", &
                                           status_non_descent = "NON_DESCENT", &
                                           status_nonfinite = "NONFINITE", &
                                           status_success = "SUCCESS"

   PUBLIC :: qs_scf_oda_apply

CONTAINS

! **************************************************************************************************
!> \brief Apply a safeguarded ODA line search to an already diagonalized raw SCF endpoint.
!> \param qs_env Quickstep environment.
!> \param scf_env SCF environment; p_mix_new contains the raw endpoint.
!> \param rho Current accepted density at lambda=0.
!> \param ks_env Kohn-Sham environment.
!> \param base_fock Raw F[P0] retained with the accepted ADIIS history state.
!> \param matrix_ks Fock matrix updated during trial evaluations.
!> \param rho_ao Current accepted AO density, indexed by spin and real-space image cell.
!> \param base_energy Evaluated energy components of P0.
!> \param predictor_lambda Initial trial step and prediction for the next consecutive ODA call.
!> \param predictor_valid Whether predictor_lambda can be reused for a consecutive ODA call.
!> \param applied Whether an energy-safeguarded ODA step was accepted.
!> \param state_evaluated Whether the returned density and Fock state were fully evaluated.
!> \param evaluated_energy Energy components of the fully evaluated returned state.
! **************************************************************************************************
   SUBROUTINE qs_scf_oda_apply(qs_env, scf_env, rho, ks_env, base_fock, matrix_ks, rho_ao, base_energy, &
                               predictor_lambda, predictor_valid, applied, state_evaluated, &
                               evaluated_energy)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), INTENT(IN)    :: base_fock
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, rho_ao
      TYPE(qs_energy_type), INTENT(IN)                   :: base_energy
      REAL(KIND=dp), INTENT(INOUT)                       :: predictor_lambda
      LOGICAL, INTENT(INOUT)                             :: predictor_valid
      LOGICAL, INTENT(OUT)                               :: applied, state_evaluated
      TYPE(qs_energy_type), INTENT(OUT)                  :: evaluated_energy

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'qs_scf_oda_apply'
      INTEGER, PARAMETER                                 :: max_backtracking = 9
      REAL(KIND=dp), PARAMETER :: armijo_factor = 1.0E-4_dp, initial_lambda = 0.5_dp, &
         maximum_backtrack_fraction = 0.8_dp, minimum_backtrack_fraction = 0.1_dp, &
         minimum_lambda = 1.0_dp/1024.0_dp, predictor_growth = 1.5_dp

      CHARACTER(len=15)                                  :: status
      CHARACTER(len=default_string_length)               :: name
      INTEGER                                            :: evaluations, handle, i, icell, ispin
      LOGICAL                                            :: model_valid, trial_accepted
      REAL(KIND=dp)                                      :: current_lambda, gradient0, gradient1, &
                                                            lambda, proposed_fraction, &
                                                            proposed_lambda, trial_energy, &
                                                            trial_lambda
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: direction
      TYPE(qs_energy_type), POINTER                      :: energy

      CALL timeset(routineN, handle)
      CALL cite_reference(Cances2000)
      CALL cite_reference(Herbst2022)
      NULLIFY (direction, energy)
      CALL get_qs_env(qs_env, energy=energy)

      applied = .FALSE.
      state_evaluated = .FALSE.
      lambda = 1.0_dp
      trial_energy = base_energy%total
      gradient0 = 0.0_dp
      gradient1 = 0.0_dp
      evaluations = 0
      status = status_nonfinite

      IF (predictor_valid .AND. ieee_is_finite(predictor_lambda) .AND. &
          predictor_lambda >= minimum_lambda .AND. predictor_lambda <= 1.0_dp) THEN
         trial_lambda = predictor_lambda
      ELSE
         trial_lambda = initial_lambda
      END IF
      predictor_lambda = 1.0_dp
      predictor_valid = .FALSE.

      CPASSERT(ASSOCIATED(energy))
      evaluated_energy = base_energy
      CPASSERT(ASSOCIATED(scf_env%p_mix_new))
      CPASSERT(SIZE(scf_env%p_mix_new, 1) == SIZE(rho_ao, 1))
      CPASSERT(SIZE(scf_env%p_mix_new, 2) == SIZE(rho_ao, 2))
      CPASSERT(SIZE(base_fock, 1) == SIZE(rho_ao, 1))
      CPASSERT(SIZE(base_fock, 2) == SIZE(rho_ao, 2))
      CPASSERT(SIZE(matrix_ks, 1) == SIZE(rho_ao, 1))
      CPASSERT(SIZE(matrix_ks, 2) == SIZE(rho_ao, 2))

      oda_search: BLOCK
         CALL create_direction(scf_env%p_mix_new, rho_ao, direction)
         CALL directional_derivative(base_fock, direction, gradient0)
         IF (.NOT. ieee_is_finite(gradient0)) EXIT oda_search
         IF (gradient0 >= 0.0_dp) THEN
            status = status_non_descent
            EXIT oda_search
         END IF

         DO i = 0, max_backtracking
            IF (trial_lambda < minimum_lambda) EXIT
            CALL set_trial_density(rho_ao, scf_env%p_mix_new, direction, trial_lambda)
            CALL evaluate_trial(qs_env, scf_env, rho, ks_env)
            evaluations = evaluations + 1
            current_lambda = trial_lambda
            trial_energy = energy%total
            CALL directional_derivative(matrix_ks, direction, gradient1)
            trial_accepted = .FALSE.
            IF (ieee_is_finite(gradient1)) THEN
               trial_accepted = armijo_satisfied(base_energy%total, gradient0, trial_lambda, &
                                                 trial_energy, armijo_factor)
            END IF
            IF (trial_accepted) THEN
               applied = .TRUE.
               lambda = trial_lambda
               IF (gradient1 < 0.0_dp) THEN
                  predictor_lambda = MIN(1.0_dp, 2.0_dp*lambda)
               ELSE
                  predictor_lambda = MIN(1.0_dp, predictor_growth*lambda)
               END IF
               predictor_valid = .TRUE.
               status = status_success
               EXIT oda_search
            END IF

            ! Scale the interval to [0,1] before using the Hermite model. The model only
            ! proposes a smaller trial; every accepted density is still evaluated explicitly.
            CALL cubic_step(base_energy%total, current_lambda*gradient0, trial_energy, &
                            current_lambda*gradient1, proposed_fraction, model_valid)
            IF (model_valid) THEN
               proposed_lambda = current_lambda*proposed_fraction
               trial_lambda = MIN(maximum_backtrack_fraction*current_lambda, &
                                  MAX(minimum_backtrack_fraction*current_lambda, proposed_lambda))
            ELSE
               trial_lambda = 0.5_dp*current_lambda
            END IF
         END DO

         ! Preserve the pre-ODA raw update as an explicit fallback. Restore its internally
         ! consistent density and Fock state after unsuccessful interior trials.
         IF (ABS(current_lambda - 1.0_dp) > 64.0_dp*EPSILON(1.0_dp)) THEN
            CALL set_trial_density(rho_ao, scf_env%p_mix_new, direction, 1.0_dp)
            CALL evaluate_trial(qs_env, scf_env, rho, ks_env)
            evaluations = evaluations + 1
            CALL directional_derivative(matrix_ks, direction, gradient1)
         END IF
         lambda = 1.0_dp
         status = status_no_acceptable
      END BLOCK oda_search

      IF (evaluations > 0) THEN
         trial_energy = energy%total
         evaluated_energy = energy
         state_evaluated = .TRUE.
      END IF
      IF (applied) THEN
         ! Keep the candidate storage aligned for consumers that inspect it before the next diagonalization.
         DO icell = 1, SIZE(rho_ao, 2)
            DO ispin = 1, SIZE(rho_ao, 1)
               CALL dbcsr_get_info(scf_env%p_mix_new(ispin, icell)%matrix, name=name)
               CALL dbcsr_copy(scf_env%p_mix_new(ispin, icell)%matrix, rho_ao(ispin, icell)%matrix, name=name)
            END DO
         END DO
      END IF
      scf_env%oda_lambda = lambda
      scf_env%oda_energy = trial_energy
      scf_env%oda_gradient0 = gradient0
      scf_env%oda_gradient1 = gradient1
      scf_env%oda_evaluations = evaluations
      scf_env%oda_status = status

      ! Keep the ordinary SCF iteration's energy convention: its printed energy components
      ! belong to the input density P0. The accepted ODA energy is returned separately.
      energy = base_energy
      IF (ASSOCIATED(direction)) CALL dbcsr_deallocate_matrix_set(direction)
      CALL timestop(handle)

   END SUBROUTINE qs_scf_oda_apply

! **************************************************************************************************
!> \brief Propose the minimum of the cubic Hermite energy model on [0,1].
!> \param energy0 Energy at lambda=0.
!> \param gradient0 Directional derivative at lambda=0.
!> \param energy1 Energy at lambda=1.
!> \param gradient1 Directional derivative at lambda=1.
!> \param lambda Proposed damping parameter.
!> \param valid Whether a finite model minimum was obtained.
! **************************************************************************************************
   PURE SUBROUTINE cubic_step(energy0, gradient0, energy1, gradient1, lambda, valid)

      REAL(KIND=dp), INTENT(IN)                          :: energy0, gradient0, energy1, gradient1
      REAL(KIND=dp), INTENT(OUT)                         :: lambda
      LOGICAL, INTENT(OUT)                               :: valid

      INTEGER                                            :: i, nroots
      REAL(KIND=dp) :: a, b, best_value, discriminant, energy_delta, model_value, q, quadratic_a, &
         quadratic_b, quadratic_c, root, scale, sqrt_discriminant, tolerance
      REAL(KIND=dp), DIMENSION(2)                        :: roots

      lambda = 1.0_dp
      valid = .FALSE.
      IF (.NOT. ieee_is_finite(energy0) .OR. .NOT. ieee_is_finite(gradient0) .OR. &
          .NOT. ieee_is_finite(energy1) .OR. .NOT. ieee_is_finite(gradient1)) RETURN

      energy_delta = energy1 - energy0
      scale = MAX(1.0_dp, ABS(energy_delta), ABS(gradient0), ABS(gradient1))
      tolerance = 128.0_dp*EPSILON(1.0_dp)*scale
      IF (gradient0 >= -tolerance) RETURN

      ! p(lambda)-E0 = g0*lambda + a*lambda**2 + b*lambda**3
      a = 3.0_dp*energy_delta - 2.0_dp*gradient0 - gradient1
      b = gradient0 + gradient1 - 2.0_dp*energy_delta
      IF (.NOT. ieee_is_finite(a) .OR. .NOT. ieee_is_finite(b)) RETURN

      quadratic_a = 3.0_dp*b
      quadratic_b = 2.0_dp*a
      quadratic_c = gradient0
      IF (.NOT. ieee_is_finite(quadratic_a) .OR. .NOT. ieee_is_finite(quadratic_b) .OR. &
          .NOT. ieee_is_finite(quadratic_c)) RETURN
      roots = 0.0_dp
      nroots = 0

      IF (ABS(quadratic_a) <= tolerance) THEN
         IF (ABS(quadratic_b) > tolerance) THEN
            nroots = 1
            roots(1) = -quadratic_c/quadratic_b
         END IF
      ELSE
         discriminant = quadratic_b*quadratic_b - 4.0_dp*quadratic_a*quadratic_c
         IF (.NOT. ieee_is_finite(discriminant)) RETURN
         IF (discriminant >= -tolerance*scale) THEN
            discriminant = MAX(0.0_dp, discriminant)
            sqrt_discriminant = SQRT(discriminant)
            q = -0.5_dp*(quadratic_b + SIGN(sqrt_discriminant, quadratic_b))
            IF (ABS(q) > tolerance) THEN
               nroots = 2
               roots(1) = q/quadratic_a
               roots(2) = quadratic_c/q
            ELSE
               nroots = 1
               roots(1) = -quadratic_b/(2.0_dp*quadratic_a)
            END IF
         END IF
      END IF

      best_value = energy_delta
      DO i = 1, nroots
         root = roots(i)
         IF (.NOT. ieee_is_finite(root)) CYCLE
         IF (root <= 0.0_dp .OR. root >= 1.0_dp) CYCLE
         model_value = gradient0*root + a*root*root + b*root*root*root
         IF (.NOT. ieee_is_finite(model_value)) CYCLE
         IF (model_value < best_value - tolerance) THEN
            best_value = model_value
            lambda = root
         END IF
      END DO

      lambda = MIN(1.0_dp, MAX(0.0_dp, lambda))
      valid = .TRUE.

   END SUBROUTINE cubic_step

! **************************************************************************************************
!> \brief Check an Armijo sufficient-decrease condition with roundoff tolerance.
!> \param energy0 Energy at lambda=0.
!> \param gradient0 Directional derivative at lambda=0.
!> \param lambda Trial damping parameter.
!> \param trial_energy Actual energy at the trial density.
!> \param armijo_factor Armijo factor in (0,1).
!> \return Whether the trial satisfies sufficient decrease within roundoff tolerance.
! **************************************************************************************************
   PURE LOGICAL FUNCTION armijo_satisfied(energy0, gradient0, lambda, trial_energy, &
                                          armijo_factor) RESULT(satisfied)

      REAL(KIND=dp), INTENT(IN)                          :: energy0, gradient0, lambda, &
                                                            trial_energy, armijo_factor

      REAL(KIND=dp)                                      :: tolerance

      satisfied = .FALSE.
      IF (.NOT. ieee_is_finite(energy0) .OR. .NOT. ieee_is_finite(gradient0) .OR. &
          .NOT. ieee_is_finite(lambda) .OR. .NOT. ieee_is_finite(trial_energy) .OR. &
          .NOT. ieee_is_finite(armijo_factor)) RETURN
      IF (gradient0 >= 0.0_dp .OR. lambda <= 0.0_dp .OR. lambda > 1.0_dp .OR. &
          armijo_factor <= 0.0_dp .OR. armijo_factor >= 1.0_dp) RETURN

      tolerance = 128.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, ABS(energy0), ABS(trial_energy))
      satisfied = trial_energy <= energy0 + armijo_factor*lambda*gradient0 + tolerance

   END FUNCTION armijo_satisfied

! **************************************************************************************************
!> \brief Allocate and form Delta P = P1-P0.
!> \param endpoint Diagonalized raw endpoint P1.
!> \param base Accepted density P0.
!> \param direction Allocated density direction P1-P0.
! **************************************************************************************************
   SUBROUTINE create_direction(endpoint, base, direction)

      TYPE(dbcsr_p_type), DIMENSION(:, :), INTENT(IN)    :: endpoint, base
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: direction

      INTEGER                                            :: icell, ispin

      CALL dbcsr_allocate_matrix_set(direction, SIZE(endpoint, 1), SIZE(endpoint, 2))
      DO icell = 1, SIZE(endpoint, 2)
         DO ispin = 1, SIZE(endpoint, 1)
            ALLOCATE (direction(ispin, icell)%matrix)
            CALL dbcsr_create(direction(ispin, icell)%matrix, template=endpoint(ispin, icell)%matrix, &
                              name="ODA DENSITY DIRECTION")
            CALL dbcsr_copy(direction(ispin, icell)%matrix, endpoint(ispin, icell)%matrix)
            CALL dbcsr_add(direction(ispin, icell)%matrix, base(ispin, icell)%matrix, &
                           alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
         END DO
      END DO

   END SUBROUTINE create_direction

! **************************************************************************************************
!> \brief Set P(lambda) = P1 + (lambda-1) Delta P.
!> \param density Trial density to overwrite.
!> \param endpoint Diagonalized raw endpoint P1.
!> \param direction Density direction P1-P0.
!> \param lambda Trial damping parameter.
! **************************************************************************************************
   SUBROUTINE set_trial_density(density, endpoint, direction, lambda)

      TYPE(dbcsr_p_type), DIMENSION(:, :), INTENT(INOUT) :: density
      TYPE(dbcsr_p_type), DIMENSION(:, :), INTENT(IN)    :: endpoint, direction
      REAL(KIND=dp), INTENT(IN)                          :: lambda

      CHARACTER(len=default_string_length)               :: name
      INTEGER                                            :: icell, ispin

      DO icell = 1, SIZE(density, 2)
         DO ispin = 1, SIZE(density, 1)
            CALL dbcsr_get_info(density(ispin, icell)%matrix, name=name)
            CALL dbcsr_copy(density(ispin, icell)%matrix, endpoint(ispin, icell)%matrix, name=name)
            IF (lambda /= 1.0_dp) THEN
               CALL dbcsr_add(density(ispin, icell)%matrix, direction(ispin, icell)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=lambda - 1.0_dp)
            END IF
         END DO
      END DO

   END SUBROUTINE set_trial_density

! **************************************************************************************************
!> \brief Sum Tr(F Delta P) over spin channels and real-space image cells.
!> \param fock Raw Fock matrices.
!> \param direction Density direction P1-P0.
!> \param derivative Directional energy derivative.
! **************************************************************************************************
   SUBROUTINE directional_derivative(fock, direction, derivative)

      TYPE(dbcsr_p_type), DIMENSION(:, :), INTENT(IN)    :: fock, direction
      REAL(KIND=dp), INTENT(OUT)                         :: derivative

      INTEGER                                            :: icell, ispin
      REAL(KIND=dp)                                      :: contribution

      derivative = 0.0_dp
      ! The k-point density transform has already folded the k-point weights into
      ! each real-space density cell, matching the contraction used by core energies.
      DO icell = 1, SIZE(fock, 2)
         DO ispin = 1, SIZE(fock, 1)
            CALL dbcsr_dot(fock(ispin, icell)%matrix, direction(ispin, icell)%matrix, contribution)
            derivative = derivative + contribution
         END DO
      END DO

   END SUBROUTINE directional_derivative

! **************************************************************************************************
!> \brief Rebuild the density-dependent energy and raw Fock matrix for one ODA trial.
!> \param qs_env Quickstep environment.
!> \param scf_env SCF environment.
!> \param rho Trial density.
!> \param ks_env Kohn-Sham environment.
! **************************************************************************************************
   SUBROUTINE evaluate_trial(qs_env, scf_env, rho, ks_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL qs_scf_rho_update(rho, qs_env, scf_env, ks_env, mix_rho=.FALSE.)
      CALL qs_ks_update_qs_env(qs_env, just_energy=.FALSE., calculate_forces=.FALSE., &
                               print_active=.FALSE.)

   END SUBROUTINE evaluate_trial

END MODULE qs_scf_oda
