!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Limited memory BFGS
!> \par History
!>       2019.10 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
MODULE almo_scf_lbfgs_types
   !USE cp_external_control,             ONLY: external_control
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr
   !USE cp_log_handling,                 ONLY: cp_to_string
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_dot,&
                                              dbcsr_release,&
                                              dbcsr_scale,&
                                              dbcsr_type
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'almo_scf_lbfgs_types'

   PUBLIC :: lbfgs_seed, &
             lbfgs_create, &
             lbfgs_release, &
             lbfgs_get_direction, &
             lbfgs_history_type

   TYPE lbfgs_history_type
      INTEGER :: nstore
      ! istore counts the total number of action=2 pushes
      ! istore is designed to become more than nstore eventually
      ! there are two counters: the main variable and gradient
      INTEGER, DIMENSION(2) :: istore
      TYPE(dbcsr_type), DIMENSION(:, :, :), ALLOCATABLE :: matrix
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: rho
   END TYPE lbfgs_history_type

CONTAINS

! **************************************************************************************************
!> \brief interface subroutine to store the first variable/gradient pair
!> \param history ...
!> \param variable ...
!> \param gradient ...
! **************************************************************************************************
   SUBROUTINE lbfgs_seed(history, variable, gradient)

      TYPE(lbfgs_history_type), INTENT(INOUT)            :: history
      TYPE(dbcsr_type), DIMENSION(:), INTENT(IN)         :: variable, gradient

      CALL lbfgs_history_push(history, variable, vartype=1, action=1)
      CALL lbfgs_history_push(history, gradient, vartype=2, action=1)

   END SUBROUTINE lbfgs_seed

! **************************************************************************************************
!> \brief interface subroutine to store a variable/gradient pair
!>        and predict direction
!> \param history ...
!> \param variable ...
!> \param gradient ...
!> \param direction ...
! **************************************************************************************************

   SUBROUTINE lbfgs_get_direction(history, variable, gradient, direction)
      TYPE(lbfgs_history_type), INTENT(INOUT)            :: history
      TYPE(dbcsr_type), DIMENSION(:), INTENT(IN)         :: variable, gradient
      TYPE(dbcsr_type), DIMENSION(:), INTENT(INOUT)      :: direction

      ! action 2 will calculate delta = (new - old)
      ! in the last used storage cell
      CALL lbfgs_history_push(history, variable, vartype=1, action=2)
      CALL lbfgs_history_push(history, gradient, vartype=2, action=2)
      ! compute rho for the last stored value
      CALL lbfgs_history_last_rho(history)

      CALL lbfgs_history_direction(history, gradient, direction)

      ! action 1 will seed the next storage cell
      CALL lbfgs_history_push(history, variable, vartype=1, action=1)
      CALL lbfgs_history_push(history, gradient, vartype=2, action=1)

   END SUBROUTINE lbfgs_get_direction

! **************************************************************************************************
!> \brief create history storage for limited memory bfgs
!> \param history ...
!> \param nspins ...
!> \param nstore ...
! **************************************************************************************************
   SUBROUTINE lbfgs_create(history, nspins, nstore)

      TYPE(lbfgs_history_type), INTENT(INOUT)            :: history
      INTEGER, INTENT(IN)                                :: nspins, nstore

      INTEGER                                            :: nallocate

      nallocate = MAX(1, nstore)
      history%nstore = nallocate
      history%istore(:) = 0  ! total number of action-2 pushes
      ALLOCATE (history%matrix(nspins, nallocate, 2))
      ALLOCATE (history%rho(nspins, nallocate))

   END SUBROUTINE lbfgs_create

! **************************************************************************************************
!> \brief release the bfgs history
!> \param history ...
! **************************************************************************************************
   SUBROUTINE lbfgs_release(history)
      TYPE(lbfgs_history_type), INTENT(INOUT)            :: history

      INTEGER                                            :: ispin, istore, ivartype

      ! delete history
      DO ispin = 1, SIZE(history%matrix, 1)
         DO ivartype = 1, 2
            DO istore = 1, MIN(history%istore(ivartype) + 1, history%nstore)
               !WRITE(*,*) "ZREL: ispin,istore,vartype", ispin, istore, ivartype
               CALL dbcsr_release(history%matrix(ispin, istore, ivartype))
            END DO
         END DO
      END DO
      DEALLOCATE (history%matrix)
      DEALLOCATE (history%rho)

   END SUBROUTINE lbfgs_release

! **************************************************************************************************
!> \brief once all data in the last cell is stored, compute rho
!> \param history ...
! **************************************************************************************************
   SUBROUTINE lbfgs_history_last_rho(history)

      TYPE(lbfgs_history_type), INTENT(INOUT)            :: history

      INTEGER                                            :: ispin, istore

      !logger => cp_get_default_logger()
      !IF (logger%para_env%mepos == logger%para_env%source) THEN
      !   unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      !ELSE
      !   unit_nr = -1
      !ENDIF

      DO ispin = 1, SIZE(history%matrix, 1)

         istore = MOD(history%istore(1) - 1, history%nstore) + 1
         CALL dbcsr_dot(history%matrix(ispin, istore, 1), &
                        history%matrix(ispin, istore, 2), &
                        history%rho(ispin, istore))

         history%rho(ispin, istore) = 1.0_dp/history%rho(ispin, istore)

         !IF (unit_nr > 0) THEN
         !   WRITE (unit_nr, *) "Rho in cell ", istore, " is computed ", history%rho(ispin, istore)
         !ENDIF

      END DO ! ispin

   END SUBROUTINE lbfgs_history_last_rho

! **************************************************************************************************
!> \brief store data in history
!>  vartype - which data piece to store: 1 - variable, 2 - gradient
!>  operation - what to do: 1 - erase existing and store new
!>                          2 - store = new - existing
!> \param history ...
!> \param matrix ...
!> \param vartype ...
!> \param action ...
! **************************************************************************************************
   SUBROUTINE lbfgs_history_push(history, matrix, vartype, action)
      TYPE(lbfgs_history_type), INTENT(INOUT)            :: history
      TYPE(dbcsr_type), DIMENSION(:), INTENT(IN)         :: matrix
      INTEGER, INTENT(IN)                                :: vartype, action

      INTEGER                                            :: ispin, istore

      !logger => cp_get_default_logger()
      !IF (logger%para_env%mepos == logger%para_env%source) THEN
      !   unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      !ELSE
      !   unit_nr = -1
      !ENDIF

      ! increase the counter: it moves the pointer to the next cell
      ! for action==1 this is a "pretend" increase; the pointer will be moved back in the end
      history%istore(vartype) = history%istore(vartype) + 1

      DO ispin = 1, SIZE(history%matrix, 1)

         istore = MOD(history%istore(vartype) - 1, history%nstore) + 1
         !IF (unit_nr > 0) THEN
         !   WRITE (unit_nr, *) "Action ", action, " modifying cell ", istore
         !END IF

         IF (history%istore(vartype) <= history%nstore .AND. &
             action .EQ. 1) THEN
            !WRITE(*,*) "ZCRE: ispin,istore,vartype", ispin, istore, vartype
            CALL dbcsr_create(history%matrix(ispin, istore, vartype), &
                              template=matrix(ispin))
            !IF (unit_nr > 0) THEN
            !   WRITE (unit_nr, *) "Creating new matrix..."
            !END IF
         END IF

         IF (action .EQ. 1) THEN
            CALL dbcsr_copy(history%matrix(ispin, istore, vartype), matrix(ispin))
         ELSE
            CALL dbcsr_add(history%matrix(ispin, istore, vartype), matrix(ispin), -1.0_dp, 1.0_dp)
         END IF

      END DO ! ispin

      ! allow the pointer to move forward only if deltas are stored (action==2)
      ! otherwise return the pointer to the previous value
      IF (action .EQ. 1) THEN
         history%istore(vartype) = history%istore(vartype) - 1
      END IF

   END SUBROUTINE lbfgs_history_push

! **************************************************************************************************
!> \brief use history data to construct dir = -Hinv.grad
!> \param history ...
!> \param gradient ...
!> \param direction ...
! **************************************************************************************************
   SUBROUTINE lbfgs_history_direction(history, gradient, direction)

      TYPE(lbfgs_history_type), INTENT(INOUT)            :: history
      TYPE(dbcsr_type), DIMENSION(:), INTENT(IN)         :: gradient
      TYPE(dbcsr_type), DIMENSION(:), INTENT(INOUT)      :: direction

      INTEGER                                            :: ispin, istore, iterm, nterms
      REAL(KIND=dp)                                      :: beta, gammak
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: alpha
      TYPE(dbcsr_type)                                   :: q

      !logger => cp_get_default_logger()
      !IF (logger%para_env%mepos == logger%para_env%source) THEN
      !   unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      !ELSE
      !   unit_nr = -1
      !ENDIF

      IF (history%istore(1) .NE. history%istore(2)) THEN
         CPABORT("BFGS APIs are not used correctly")
      END IF

      nterms = MIN(history%istore(1), history%nstore)
      !IF (unit_nr > 0) THEN
      !   WRITE (unit_nr, *) "L-BFGS terms used: ", nterms
      !END IF

      ALLOCATE (alpha(nterms))

      DO ispin = 1, SIZE(history%matrix, 1)

         CALL dbcsr_create(q, template=gradient(ispin))

         CALL dbcsr_copy(q, gradient(ispin))

         ! loop over all stored items
         DO iterm = 1, nterms

            ! location: from recent to oldest stored
            istore = MOD(history%istore(1) - iterm, history%nstore) + 1

            !IF (unit_nr > 0) THEN
            !   WRITE (unit_nr, *) "Record locator: ", istore
            !END IF

            CALL dbcsr_dot(history%matrix(ispin, istore, 1), q, alpha(iterm))
            alpha(iterm) = history%rho(ispin, istore)*alpha(iterm)
            CALL dbcsr_add(q, history%matrix(ispin, istore, 2), 1.0_dp, -alpha(iterm))

            ! use the most recent term to
            ! compute gamma_k, Nocedal (7.20) and then get H0
            IF (iterm .EQ. 1) THEN
               CALL dbcsr_dot(history%matrix(ispin, istore, 2), history%matrix(ispin, istore, 2), gammak)
               gammak = 1.0_dp/(gammak*history%rho(ispin, istore))
               !IF (unit_nr > 0) THEN
               !   WRITE (unit_nr, *) "Gamma_k: ", gammak
               !END IF
            END IF

         END DO ! iterm, first loop from recent to oldest

         ! now q stores Nocedal's r = (gamma_k*I).q
         CALL dbcsr_scale(q, gammak)

         ! loop over all stored items
         DO iterm = nterms, 1, -1

            ! location: from oldest to recent stored
            istore = MOD(history%istore(1) - iterm, history%nstore) + 1

            CALL dbcsr_dot(history%matrix(ispin, istore, 2), q, beta)
            beta = history%rho(ispin, istore)*beta
            CALL dbcsr_add(q, history%matrix(ispin, istore, 1), 1.0_dp, alpha(iterm) - beta)

         END DO ! iterm, forst loop from recent to oldest

         !RZK-warning: unclear whether q should be multiplied by minus one
         CALL dbcsr_scale(q, -1.0)
         CALL dbcsr_copy(direction(ispin), q)

         CALL dbcsr_release(q)

      END DO !ispin

      DEALLOCATE (alpha)

   END SUBROUTINE lbfgs_history_direction

END MODULE almo_scf_lbfgs_types

