!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Data types for the ADIIS SCF subspace accelerator.
! **************************************************************************************************
MODULE qs_scf_subspace_types

   USE cp_dbcsr_api,                    ONLY: dbcsr_deallocate_matrix,&
                                              dbcsr_p_type
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_subspace_types'

   PUBLIC :: qs_scf_subspace_buffer_type, &
             qs_scf_subspace_buffer_clear, &
             qs_scf_subspace_buffer_create, &
             qs_scf_subspace_buffer_release

! **************************************************************************************************
!> \brief History buffer holding strictly paired P and F[P] SCF states.
!>
!>      Matrix history is stored in the real-space cell representation.  This is
!>      the representation of the density that actually built the current KS
!>      matrix, and therefore applies to both Gamma-point and k-point runs.
! **************************************************************************************************
   TYPE qs_scf_subspace_buffer_type
      INTEGER                                           :: nbuffer = 0
      INTEGER                                           :: ncall = 0
      INTEGER                                           :: nstored = 0
      INTEGER                                           :: last_status = 0
      LOGICAL                                           :: last_restart = .FALSE.
      LOGICAL                                           :: use_combined_fock = .FALSE.
      LOGICAL                                           :: diis_state_valid = .FALSE.
      REAL(KIND=dp)                                     :: last_objective = 0.0_dp
      ! Total coefficient weight assigned to all entries older than the current raw Fock.
      REAL(KIND=dp)                                     :: last_old_fock_weight = 0.0_dp
      REAL(KIND=dp)                                     :: diis_weight = 0.0_dp
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER   :: density => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER   :: fock => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER      :: combined_fock => NULL()
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: generation
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE       :: pf_metric
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE          :: coefficients, state_energy
   END TYPE qs_scf_subspace_buffer_type

CONTAINS

! **************************************************************************************************
!> \brief Initialize SCF subspace history metadata.
!> \param buffer Buffer to initialize.
!> \param nbuffer Maximum number of paired history entries.
! **************************************************************************************************
   PURE SUBROUTINE qs_scf_subspace_buffer_create(buffer, nbuffer)

      TYPE(qs_scf_subspace_buffer_type), INTENT(OUT)     :: buffer
      INTEGER, INTENT(IN)                                :: nbuffer

      buffer%nbuffer = MAX(nbuffer, 0)
      buffer%ncall = 0
      buffer%nstored = 0
      buffer%last_status = 0
      buffer%last_restart = .FALSE.
      buffer%use_combined_fock = .FALSE.
      buffer%diis_state_valid = .FALSE.
      buffer%last_objective = 0.0_dp
      buffer%last_old_fock_weight = 0.0_dp
      buffer%diis_weight = 0.0_dp
      NULLIFY (buffer%density, buffer%fock, buffer%combined_fock)

   END SUBROUTINE qs_scf_subspace_buffer_create

! **************************************************************************************************
!> \brief Clear logical history while retaining allocated matrix storage.
!> \param buffer Buffer to clear.
! **************************************************************************************************
   PURE SUBROUTINE qs_scf_subspace_buffer_clear(buffer)

      TYPE(qs_scf_subspace_buffer_type), INTENT(INOUT)   :: buffer

      buffer%ncall = 0
      buffer%nstored = 0
      buffer%last_status = 0
      buffer%last_restart = .FALSE.
      buffer%use_combined_fock = .FALSE.
      buffer%diis_state_valid = .FALSE.
      buffer%last_objective = 0.0_dp
      buffer%last_old_fock_weight = 0.0_dp
      buffer%diis_weight = 0.0_dp
      IF (ALLOCATED(buffer%generation)) buffer%generation = 0
      IF (ALLOCATED(buffer%pf_metric)) buffer%pf_metric = 0.0_dp
      IF (ALLOCATED(buffer%coefficients)) buffer%coefficients = 0.0_dp
      IF (ALLOCATED(buffer%state_energy)) buffer%state_energy = HUGE(1.0_dp)

   END SUBROUTINE qs_scf_subspace_buffer_clear

! **************************************************************************************************
!> \brief Release all matrices and scalar storage owned by a subspace buffer.
!> \param buffer Buffer to release.
! **************************************************************************************************
   SUBROUTINE qs_scf_subspace_buffer_release(buffer)

      TYPE(qs_scf_subspace_buffer_type), INTENT(INOUT)   :: buffer

      INTEGER                                            :: icell, islot, ispin

      IF (ASSOCIATED(buffer%density)) THEN
         DO icell = 1, SIZE(buffer%density, 3)
            DO ispin = 1, SIZE(buffer%density, 2)
               DO islot = 1, SIZE(buffer%density, 1)
                  IF (ASSOCIATED(buffer%density(islot, ispin, icell)%matrix)) THEN
                     CALL dbcsr_deallocate_matrix(buffer%density(islot, ispin, icell)%matrix)
                  END IF
               END DO
            END DO
         END DO
         DEALLOCATE (buffer%density)
      END IF

      IF (ASSOCIATED(buffer%fock)) THEN
         DO icell = 1, SIZE(buffer%fock, 3)
            DO ispin = 1, SIZE(buffer%fock, 2)
               DO islot = 1, SIZE(buffer%fock, 1)
                  IF (ASSOCIATED(buffer%fock(islot, ispin, icell)%matrix)) THEN
                     CALL dbcsr_deallocate_matrix(buffer%fock(islot, ispin, icell)%matrix)
                  END IF
               END DO
            END DO
         END DO
         DEALLOCATE (buffer%fock)
      END IF

      IF (ASSOCIATED(buffer%combined_fock)) THEN
         DO icell = 1, SIZE(buffer%combined_fock, 2)
            DO ispin = 1, SIZE(buffer%combined_fock, 1)
               IF (ASSOCIATED(buffer%combined_fock(ispin, icell)%matrix)) THEN
                  CALL dbcsr_deallocate_matrix(buffer%combined_fock(ispin, icell)%matrix)
               END IF
            END DO
         END DO
         DEALLOCATE (buffer%combined_fock)
      END IF

      IF (ALLOCATED(buffer%generation)) DEALLOCATE (buffer%generation)
      IF (ALLOCATED(buffer%pf_metric)) DEALLOCATE (buffer%pf_metric)
      IF (ALLOCATED(buffer%coefficients)) DEALLOCATE (buffer%coefficients)
      IF (ALLOCATED(buffer%state_energy)) DEALLOCATE (buffer%state_energy)
      CALL qs_scf_subspace_buffer_create(buffer, 0)

   END SUBROUTINE qs_scf_subspace_buffer_release

END MODULE qs_scf_subspace_types
