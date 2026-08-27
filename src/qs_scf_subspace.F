!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief History management and matrix combination for ADIIS.
! **************************************************************************************************
MODULE qs_scf_subspace

   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_p_type,&
                                              dbcsr_set
   USE cp_dbcsr_contrib,                ONLY: dbcsr_dot
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE qs_scf_subspace_math,            ONLY: qs_scf_subspace_build_adiis_model,&
                                              qs_scf_subspace_fifo_slot,&
                                              simplex_qp_nonfinite_input,&
                                              simplex_qp_pairwise_stationary,&
                                              simplex_qp_success,&
                                              simplex_quadratic_minimize
   USE qs_scf_subspace_types,           ONLY: qs_scf_subspace_buffer_clear,&
                                              qs_scf_subspace_buffer_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_subspace'

   INTEGER, PARAMETER, PUBLIC :: scf_subspace_invalid_history = 100

   PUBLIC :: qs_scf_subspace_build, qs_scf_subspace_push, qs_scf_subspace_restart

CONTAINS

! **************************************************************************************************
!> \brief Restart the ADIIS history from the current strictly paired P, F[P] state.
!> \param buffer Persistent SCF subspace history.
!> \param fock Current raw F[P], indexed by spin and real-space cell.
!> \param density Current input P that generated fock, with matching indices.
!> \param state_energy ...
!> \param restarted Whether the current state was stored successfully.
! **************************************************************************************************
   SUBROUTINE qs_scf_subspace_restart(buffer, fock, density, state_energy, restarted)

      TYPE(qs_scf_subspace_buffer_type), INTENT(INOUT)   :: buffer
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: fock, density
      REAL(KIND=dp), INTENT(IN)                          :: state_energy
      LOGICAL, INTENT(OUT)                               :: restarted

      restarted = .FALSE.
      CALL qs_scf_subspace_buffer_clear(buffer)
      CALL qs_scf_subspace_push(buffer, fock, density, state_energy, restarted)
      IF (.NOT. restarted) RETURN

      buffer%coefficients = 0.0_dp
      buffer%coefficients(1) = 1.0_dp
      buffer%last_status = simplex_qp_success
      buffer%last_restart = .TRUE.
      buffer%use_combined_fock = .FALSE.
      buffer%last_objective = 0.0_dp
      buffer%last_old_fock_weight = 0.0_dp

   END SUBROUTINE qs_scf_subspace_restart

! **************************************************************************************************
!> \brief Append one accepted, evaluated P,F[P] state to the ADIIS history.
!>
!>      The SCF driver owns state acceptance.  Trial endpoints and arbitrary initial guesses
!>      must not call this routine.
!> \param buffer Persistent SCF subspace history.
!> \param fock Current raw F[P], indexed by spin and real-space cell.
!> \param density Current input P that generated fock, with matching indices.
!> \param state_energy ...
!> \param pushed Whether the state was stored successfully.
! **************************************************************************************************
   SUBROUTINE qs_scf_subspace_push(buffer, fock, density, state_energy, pushed)

      TYPE(qs_scf_subspace_buffer_type), INTENT(INOUT)   :: buffer
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: fock, density
      REAL(KIND=dp), INTENT(IN)                          :: state_energy
      LOGICAL, INTENT(OUT)                               :: pushed

      LOGICAL                                            :: storage_valid

      pushed = .FALSE.
      buffer%use_combined_fock = .FALSE.
      buffer%last_restart = .FALSE.
      buffer%last_objective = 0.0_dp
      buffer%last_old_fock_weight = 0.0_dp
      IF (ALLOCATED(buffer%coefficients)) buffer%coefficients = 0.0_dp

      IF (buffer%nbuffer < 1) THEN
         buffer%last_status = scf_subspace_invalid_history
         RETURN
      END IF
      IF (buffer%ncall < 0 .OR. buffer%nstored < 0 .OR. buffer%nstored > buffer%nbuffer .OR. &
          buffer%nstored /= MIN(buffer%ncall, buffer%nbuffer)) THEN
         buffer%last_status = scf_subspace_invalid_history
         RETURN
      END IF

      CALL ensure_matrix_storage(buffer, fock, density, storage_valid)
      IF (.NOT. storage_valid) THEN
         buffer%last_status = scf_subspace_invalid_history
         RETURN
      END IF

      CALL push_paired_state(buffer, fock, density, state_energy, pushed)
      IF (.NOT. pushed) THEN
         buffer%last_status = simplex_qp_nonfinite_input
         RETURN
      END IF
      buffer%last_status = simplex_qp_success

   END SUBROUTINE qs_scf_subspace_push

! **************************************************************************************************
!> \brief Solve the ADIIS model from previously accepted history and form an effective KS matrix.
!> \param buffer Persistent SCF subspace history.
! **************************************************************************************************
   SUBROUTINE qs_scf_subspace_build(buffer)

      TYPE(qs_scf_subspace_buffer_type), INTENT(INOUT)   :: buffer

      INTEGER                                            :: i, j, m, physical
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: slots
      LOGICAL                                            :: model_valid, solver_usable
      REAL(KIND=dp)                                      :: coefficient_sum, model_delta, &
                                                            model_tolerance
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: difference, linear
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: chronological_pf, hessian

      buffer%use_combined_fock = .FALSE.
      buffer%last_restart = .FALSE.
      buffer%last_objective = 0.0_dp
      buffer%last_old_fock_weight = 0.0_dp
      IF (ALLOCATED(buffer%coefficients)) buffer%coefficients = 0.0_dp

      IF (buffer%nbuffer < 1 .OR. buffer%nstored < 1) THEN
         buffer%last_status = scf_subspace_invalid_history
         RETURN
      END IF
      IF (buffer%ncall < 0 .OR. buffer%nstored < 0 .OR. buffer%nstored > buffer%nbuffer .OR. &
          buffer%nstored /= MIN(buffer%ncall, buffer%nbuffer)) THEN
         buffer%last_status = scf_subspace_invalid_history
         RETURN
      END IF

      IF (.NOT. ASSOCIATED(buffer%fock) .OR. .NOT. ASSOCIATED(buffer%density) .OR. &
          .NOT. ASSOCIATED(buffer%combined_fock)) THEN
         buffer%last_status = scf_subspace_invalid_history
         RETURN
      END IF
      IF (.NOT. ALLOCATED(buffer%generation) .OR. .NOT. ALLOCATED(buffer%pf_metric) .OR. &
          .NOT. ALLOCATED(buffer%coefficients) .OR. .NOT. ALLOCATED(buffer%state_energy)) THEN
         buffer%last_status = scf_subspace_invalid_history
         RETURN
      END IF
      IF (SIZE(buffer%generation) /= buffer%nbuffer .OR. &
          SIZE(buffer%pf_metric, 1) /= buffer%nbuffer .OR. SIZE(buffer%pf_metric, 2) /= buffer%nbuffer .OR. &
          SIZE(buffer%coefficients) /= buffer%nbuffer .OR. SIZE(buffer%state_energy) /= buffer%nbuffer .OR. &
          COUNT(buffer%generation > 0) /= buffer%nstored .OR. ANY(buffer%generation < 0)) THEN
         buffer%last_status = scf_subspace_invalid_history
         RETURN
      END IF

      m = buffer%nstored
      ALLOCATE (slots(m), chronological_pf(m, m), hessian(m, m), linear(m))
      CALL chronological_slots(buffer, slots)

      DO i = 1, m
         DO j = 1, m
            chronological_pf(i, j) = buffer%pf_metric(slots(i), slots(j))
         END DO
      END DO

      CALL qs_scf_subspace_build_adiis_model(chronological_pf, m, hessian, linear, model_valid)

      buffer%coefficients = 0.0_dp
      IF (.NOT. model_valid) THEN
         CALL restart_newest_after_failure(buffer, scf_subspace_invalid_history)
         RETURN
      END IF

      CALL simplex_quadratic_minimize(hessian, linear, buffer%coefficients(1:m), &
                                      buffer%last_objective, buffer%last_status, preferred_index=m)

      solver_usable = buffer%last_status == simplex_qp_success .OR. &
                      buffer%last_status == simplex_qp_pairwise_stationary
      IF (.NOT. solver_usable) THEN
         CALL restart_newest_after_failure(buffer, buffer%last_status)
         RETURN
      END IF
      IF (.NOT. ieee_is_finite(buffer%last_objective) .OR. &
          .NOT. ALL(ieee_is_finite(buffer%coefficients(1:m)))) THEN
         CALL restart_newest_after_failure(buffer, simplex_qp_nonfinite_input)
         RETURN
      END IF

      coefficient_sum = SUM(buffer%coefficients(1:m))
      IF (MINVAL(buffer%coefficients(1:m)) < -100.0_dp*EPSILON(1.0_dp) .OR. &
          ABS(coefficient_sum - 1.0_dp) > 1000.0_dp*EPSILON(1.0_dp)) THEN
         CALL restart_newest_after_failure(buffer, scf_subspace_invalid_history)
         RETURN
      END IF

      ALLOCATE (difference(m))
      difference(:) = buffer%coefficients(1:m)
      difference(m) = difference(m) - 1.0_dp
      model_delta = 0.5_dp*DOT_PRODUCT(difference, MATMUL(hessian, difference)) + &
                    DOT_PRODUCT(difference, hessian(:, m) + linear)
      model_tolerance = 1000.0_dp*EPSILON(1.0_dp)*REAL(m, KIND=dp)* &
                        MAX(1.0_dp, MAXVAL(ABS(hessian)), MAXVAL(ABS(linear)))
      IF (.NOT. ieee_is_finite(model_delta) .OR. model_delta > model_tolerance) THEN
         IF (ieee_is_finite(model_delta)) THEN
            CALL restart_newest_after_failure(buffer, scf_subspace_invalid_history)
         ELSE
            CALL restart_newest_after_failure(buffer, simplex_qp_nonfinite_input)
         END IF
         RETURN
      END IF

      buffer%last_old_fock_weight = MAX(0.0_dp, 1.0_dp - buffer%coefficients(m))
      IF (buffer%last_old_fock_weight <= 100.0_dp*EPSILON(1.0_dp)) RETURN

      CALL zero_fock_combination(buffer)
      DO j = 1, m
         physical = slots(j)
         IF (buffer%coefficients(j) <= 0.0_dp) CYCLE
         CALL add_fock_to_combination(buffer, physical, buffer%coefficients(j))
      END DO
      buffer%use_combined_fock = .TRUE.

   END SUBROUTINE qs_scf_subspace_build

! **************************************************************************************************
!> \brief Retain only the newest accepted state after an unusable ADIIS model.
!> \param buffer Persistent SCF subspace history.
!> \param failure_status Status code that caused the restart.
! **************************************************************************************************
   SUBROUTINE restart_newest_after_failure(buffer, failure_status)

      TYPE(qs_scf_subspace_buffer_type), INTENT(INOUT)   :: buffer
      INTEGER, INTENT(IN)                                :: failure_status

      INTEGER                                            :: icell, ispin, newest
      REAL(KIND=dp)                                      :: contribution, newest_energy, newest_pf

      IF (buffer%nstored < 1) THEN
         buffer%last_status = failure_status
         RETURN
      END IF
      IF (COUNT(buffer%generation > 0) /= buffer%nstored) THEN
         buffer%last_status = failure_status
         RETURN
      END IF
      newest = MAXLOC(buffer%generation, DIM=1)
      newest_energy = buffer%state_energy(newest)
      newest_pf = 0.0_dp
      DO icell = 1, SIZE(buffer%fock, 3)
         DO ispin = 1, SIZE(buffer%fock, 2)
            CALL dbcsr_dot(buffer%density(newest, ispin, icell)%matrix, &
                           buffer%fock(newest, ispin, icell)%matrix, contribution)
            newest_pf = newest_pf + contribution
         END DO
      END DO
      IF (.NOT. ieee_is_finite(newest_energy) .OR. .NOT. ieee_is_finite(newest_pf)) THEN
         buffer%last_status = failure_status
         RETURN
      END IF

      buffer%ncall = 1
      buffer%nstored = 1
      buffer%generation = 0
      buffer%generation(newest) = 1
      buffer%state_energy = HUGE(1.0_dp)
      buffer%state_energy(newest) = newest_energy
      buffer%pf_metric = 0.0_dp
      buffer%pf_metric(newest, newest) = newest_pf
      buffer%coefficients = 0.0_dp
      buffer%coefficients(1) = 1.0_dp
      buffer%last_status = failure_status
      buffer%last_restart = .TRUE.
      buffer%use_combined_fock = .FALSE.
      buffer%last_objective = 0.0_dp
      buffer%last_old_fock_weight = 0.0_dp

   END SUBROUTINE restart_newest_after_failure

! **************************************************************************************************
!> \brief Allocate matrix and scalar storage lazily from current matrix templates.
!> \param buffer Persistent SCF subspace history.
!> \param fock Current Fock-matrix templates.
!> \param density Current density-matrix templates.
!> \param valid Whether dimensions and existing allocation match.
! **************************************************************************************************
   SUBROUTINE ensure_matrix_storage(buffer, fock, density, valid)

      TYPE(qs_scf_subspace_buffer_type), INTENT(INOUT)   :: buffer
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: fock, density
      LOGICAL, INTENT(OUT)                               :: valid

      INTEGER                                            :: icell, islot, ispin, ncell, nspin

      valid = ASSOCIATED(fock) .AND. ASSOCIATED(density)
      IF (.NOT. valid) RETURN
      nspin = SIZE(fock, 1)
      ncell = SIZE(fock, 2)
      valid = nspin > 0 .AND. ncell > 0 .AND. &
              SIZE(density, 1) == nspin .AND. SIZE(density, 2) == ncell
      IF (.NOT. valid) RETURN
      DO icell = 1, ncell
         DO ispin = 1, nspin
            valid = ASSOCIATED(fock(ispin, icell)%matrix) .AND. &
                    ASSOCIATED(density(ispin, icell)%matrix)
            IF (.NOT. valid) RETURN
         END DO
      END DO

      IF (ASSOCIATED(buffer%fock)) THEN
         valid = SIZE(buffer%fock, 1) == buffer%nbuffer .AND. &
                 SIZE(buffer%fock, 2) == nspin .AND. SIZE(buffer%fock, 3) == ncell
         IF (.NOT. valid) RETURN
         valid = ASSOCIATED(buffer%density)
         IF (.NOT. valid) RETURN
         valid = SIZE(buffer%density, 1) == buffer%nbuffer .AND. &
                 SIZE(buffer%density, 2) == nspin .AND. SIZE(buffer%density, 3) == ncell
         IF (.NOT. valid) RETURN
         valid = ASSOCIATED(buffer%combined_fock)
         IF (.NOT. valid) RETURN
         valid = SIZE(buffer%combined_fock, 1) == nspin .AND. &
                 SIZE(buffer%combined_fock, 2) == ncell
         IF (.NOT. valid) RETURN
         valid = ALLOCATED(buffer%generation) .AND. ALLOCATED(buffer%pf_metric) .AND. &
                 ALLOCATED(buffer%coefficients) .AND. ALLOCATED(buffer%state_energy)
         IF (.NOT. valid) RETURN
         valid = SIZE(buffer%generation) == buffer%nbuffer .AND. &
                 SIZE(buffer%pf_metric, 1) == buffer%nbuffer .AND. &
                 SIZE(buffer%pf_metric, 2) == buffer%nbuffer .AND. &
                 SIZE(buffer%coefficients) == buffer%nbuffer .AND. &
                 SIZE(buffer%state_energy) == buffer%nbuffer
         IF (.NOT. valid) RETURN
         DO icell = 1, ncell
            DO ispin = 1, nspin
               valid = ASSOCIATED(buffer%combined_fock(ispin, icell)%matrix)
               IF (.NOT. valid) RETURN
               DO islot = 1, buffer%nbuffer
                  valid = ASSOCIATED(buffer%fock(islot, ispin, icell)%matrix) .AND. &
                          ASSOCIATED(buffer%density(islot, ispin, icell)%matrix)
                  IF (.NOT. valid) RETURN
               END DO
            END DO
         END DO
         RETURN
      END IF

      valid = .NOT. ASSOCIATED(buffer%density) .AND. &
              .NOT. ASSOCIATED(buffer%combined_fock) .AND. .NOT. ALLOCATED(buffer%generation) .AND. &
              .NOT. ALLOCATED(buffer%pf_metric) .AND. .NOT. ALLOCATED(buffer%coefficients) .AND. &
              .NOT. ALLOCATED(buffer%state_energy)
      IF (.NOT. valid) RETURN

      ALLOCATE (buffer%fock(buffer%nbuffer, nspin, ncell))
      ALLOCATE (buffer%density(buffer%nbuffer, nspin, ncell))
      ALLOCATE (buffer%combined_fock(nspin, ncell))
      ALLOCATE (buffer%generation(buffer%nbuffer), buffer%state_energy(buffer%nbuffer))
      ALLOCATE (buffer%pf_metric(buffer%nbuffer, buffer%nbuffer), buffer%coefficients(buffer%nbuffer))
      buffer%generation = 0
      buffer%state_energy = HUGE(1.0_dp)
      buffer%pf_metric = 0.0_dp
      buffer%coefficients = 0.0_dp

      DO icell = 1, ncell
         DO ispin = 1, nspin
            DO islot = 1, buffer%nbuffer
               ALLOCATE (buffer%fock(islot, ispin, icell)%matrix)
               CALL dbcsr_create(buffer%fock(islot, ispin, icell)%matrix, &
                                 template=fock(ispin, icell)%matrix)
               ALLOCATE (buffer%density(islot, ispin, icell)%matrix)
               CALL dbcsr_create(buffer%density(islot, ispin, icell)%matrix, &
                                 template=density(ispin, icell)%matrix)
            END DO
            ALLOCATE (buffer%combined_fock(ispin, icell)%matrix)
            CALL dbcsr_create(buffer%combined_fock(ispin, icell)%matrix, &
                              template=fock(ispin, icell)%matrix)
         END DO
      END DO

   END SUBROUTINE ensure_matrix_storage

! **************************************************************************************************
!> \brief Atomically append one strictly paired P, F[P] state.
!> \param buffer Persistent SCF subspace history.
!> \param fock Current raw F[P].
!> \param density Current P.
!> \param state_energy ...
!> \param pushed Whether all scalar checks passed and the state was committed.
! **************************************************************************************************
   SUBROUTINE push_paired_state(buffer, fock, density, state_energy, pushed)

      TYPE(qs_scf_subspace_buffer_type), INTENT(INOUT)   :: buffer
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: fock, density
      REAL(KIND=dp), INTENT(IN)                          :: state_energy
      LOGICAL, INTENT(OUT)                               :: pushed

      INTEGER                                            :: icell, islot, ispin, target_slot
      LOGICAL                                            :: retained
      REAL(KIND=dp)                                      :: current_current
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: current_old, old_current

      pushed = .FALSE.
      IF (.NOT. ieee_is_finite(state_energy)) RETURN

      ! Keep the baseline deterministic and auditable: fill empty slots in physical
      ! order, then evict the oldest accepted state. Energy remains metadata.
      target_slot = qs_scf_subspace_fifo_slot(buffer%generation)
      IF (target_slot <= 0) RETURN

      ALLOCATE (current_old(buffer%nbuffer), old_current(buffer%nbuffer))
      current_old = 0.0_dp
      old_current = 0.0_dp

      CALL paired_matrix_dot(density, fock, current_current)
      IF (.NOT. ieee_is_finite(current_current)) RETURN

      DO islot = 1, buffer%nbuffer
         retained = buffer%generation(islot) > 0 .AND. islot /= target_slot
         IF (.NOT. retained) CYCLE
         CALL history_current_dot(buffer, islot, fock, density, old_current(islot), current_old(islot))
         IF (.NOT. ieee_is_finite(old_current(islot)) .OR. &
             .NOT. ieee_is_finite(current_old(islot))) RETURN
      END DO

      DO icell = 1, SIZE(fock, 2)
         DO ispin = 1, SIZE(fock, 1)
            CALL dbcsr_copy(buffer%fock(target_slot, ispin, icell)%matrix, fock(ispin, icell)%matrix)
            CALL dbcsr_copy(buffer%density(target_slot, ispin, icell)%matrix, density(ispin, icell)%matrix)
         END DO
      END DO

      buffer%pf_metric(target_slot, :) = 0.0_dp
      buffer%pf_metric(:, target_slot) = 0.0_dp
      buffer%pf_metric(target_slot, target_slot) = current_current
      DO islot = 1, buffer%nbuffer
         retained = buffer%generation(islot) > 0 .AND. islot /= target_slot
         IF (.NOT. retained) CYCLE
         buffer%pf_metric(islot, target_slot) = old_current(islot)
         buffer%pf_metric(target_slot, islot) = current_old(islot)
      END DO

      buffer%ncall = buffer%ncall + 1
      buffer%generation(target_slot) = buffer%ncall
      buffer%state_energy(target_slot) = state_energy
      buffer%nstored = MIN(buffer%nstored + 1, buffer%nbuffer)
      pushed = .TRUE.

   END SUBROUTINE push_paired_state

! **************************************************************************************************
!> \brief Compute Tr(P F), summed across every spin and real-space cell.
!> \param density Density matrices.
!> \param fock Fock matrices.
!> \param value Summed trace contraction.
! **************************************************************************************************
   SUBROUTINE paired_matrix_dot(density, fock, value)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: density, fock
      REAL(KIND=dp), INTENT(OUT)                         :: value

      INTEGER                                            :: icell, ispin
      REAL(KIND=dp)                                      :: contribution

      value = 0.0_dp
      DO icell = 1, SIZE(fock, 2)
         DO ispin = 1, SIZE(fock, 1)
            CALL dbcsr_dot(density(ispin, icell)%matrix, fock(ispin, icell)%matrix, contribution)
            value = value + contribution
         END DO
      END DO

   END SUBROUTINE paired_matrix_dot

! **************************************************************************************************
!> \brief Compute both cross contractions between one old and the current state.
!> \param buffer Persistent SCF subspace history.
!> \param islot Physical old-history slot.
!> \param current_fock Current F[P].
!> \param current_density Current P.
!> \param old_current Tr(P_old F_current).
!> \param current_old Tr(P_current F_old).
! **************************************************************************************************
   SUBROUTINE history_current_dot(buffer, islot, current_fock, current_density, old_current, current_old)

      TYPE(qs_scf_subspace_buffer_type), INTENT(IN)      :: buffer
      INTEGER, INTENT(IN)                                :: islot
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: current_fock, current_density
      REAL(KIND=dp), INTENT(OUT)                         :: old_current, current_old

      INTEGER                                            :: icell, ispin
      REAL(KIND=dp)                                      :: contribution

      old_current = 0.0_dp
      current_old = 0.0_dp
      DO icell = 1, SIZE(current_fock, 2)
         DO ispin = 1, SIZE(current_fock, 1)
            CALL dbcsr_dot(buffer%density(islot, ispin, icell)%matrix, &
                           current_fock(ispin, icell)%matrix, contribution)
            old_current = old_current + contribution
            CALL dbcsr_dot(current_density(ispin, icell)%matrix, &
                           buffer%fock(islot, ispin, icell)%matrix, contribution)
            current_old = current_old + contribution
         END DO
      END DO

   END SUBROUTINE history_current_dot

! **************************************************************************************************
!> \brief Return physical slots ordered from oldest to newest.
!> \param buffer Persistent SCF subspace history.
!> \param slots Chronologically ordered physical slot indices.
! **************************************************************************************************
   PURE SUBROUTINE chronological_slots(buffer, slots)

      TYPE(qs_scf_subspace_buffer_type), INTENT(IN)      :: buffer
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: slots

      INTEGER                                            :: i, islot, j, tmp

      j = 0
      DO islot = 1, buffer%nbuffer
         IF (buffer%generation(islot) <= 0) CYCLE
         j = j + 1
         IF (j <= SIZE(slots)) slots(j) = islot
      END DO

      ! The ADIIS model uses the most recent state as its expansion point, so
      ! return physical slots ordered by insertion generation with the newest
      ! state last. The history size is small, making insertion sort sufficient.
      DO i = 2, SIZE(slots)
         tmp = slots(i)
         j = i - 1
         DO WHILE (j >= 1)
            IF (buffer%generation(slots(j)) <= buffer%generation(tmp)) EXIT
            slots(j + 1) = slots(j)
            j = j - 1
         END DO
         slots(j + 1) = tmp
      END DO

   END SUBROUTINE chronological_slots

! **************************************************************************************************
!> \brief Zero the effective KS-matrix combination.
!> \param buffer Persistent SCF subspace history and combination target.
! **************************************************************************************************
   SUBROUTINE zero_fock_combination(buffer)

      TYPE(qs_scf_subspace_buffer_type), INTENT(INOUT)   :: buffer

      INTEGER                                            :: icell, ispin

      DO icell = 1, SIZE(buffer%combined_fock, 2)
         DO ispin = 1, SIZE(buffer%combined_fock, 1)
            CALL dbcsr_set(buffer%combined_fock(ispin, icell)%matrix, 0.0_dp)
         END DO
      END DO

   END SUBROUTINE zero_fock_combination

! **************************************************************************************************
!> \brief Add one physical history slot to the effective KS-matrix combination.
!> \param buffer Persistent SCF subspace history and combination target.
!> \param physical Physical history slot.
!> \param coefficient Weight of this history matrix.
! **************************************************************************************************
   SUBROUTINE add_fock_to_combination(buffer, physical, coefficient)

      TYPE(qs_scf_subspace_buffer_type), INTENT(INOUT)   :: buffer
      INTEGER, INTENT(IN)                                :: physical
      REAL(KIND=dp), INTENT(IN)                          :: coefficient

      INTEGER                                            :: icell, ispin

      DO icell = 1, SIZE(buffer%combined_fock, 2)
         DO ispin = 1, SIZE(buffer%combined_fock, 1)
            CALL dbcsr_add(buffer%combined_fock(ispin, icell)%matrix, &
                           buffer%fock(physical, ispin, icell)%matrix, &
                           alpha_scalar=1.0_dp, beta_scalar=coefficient)
         END DO
      END DO

   END SUBROUTINE add_fock_to_combination

END MODULE qs_scf_subspace
