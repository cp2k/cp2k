!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief types of  preconditioners
!> \par History
!>      Separate types from construction and application
!> \author Joost VandeVondele (09.2002)
! *****************************************************************************
MODULE preconditioner_types
  USE cp_blacs_env,                    ONLY: cp_blacs_env_release,&
                                             cp_blacs_env_retain,&
                                             cp_blacs_env_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_deallocate_matrix,&
                                             cp_dbcsr_p_type,&
                                             cp_dbcsr_release_p,&
                                             cp_dbcsr_type
  USE cp_fm_types,                     ONLY: cp_fm_release,&
                                             cp_fm_type
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_constants,                 ONLY: ot_precond_solver_default
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'preconditioner_types'

  PUBLIC  :: preconditioner_type, preconditioner_p_type
  PUBLIC  :: init_preconditioner
  PUBLIC  :: preconditioner_in_use
  PUBLIC  :: destroy_preconditioner

! *****************************************************************************
  TYPE preconditioner_type
!    PRIVATE
       TYPE(cp_dbcsr_type),    POINTER :: sparse_matrix
       TYPE(cp_fm_type), POINTER :: fm
       TYPE(cp_dbcsr_type), POINTER           :: dbcsr_matrix
       TYPE(cp_dbcsr_type), POINTER           :: max_ev_vector
       TYPE(cp_dbcsr_type), POINTER           :: min_ev_vector
       TYPE(cp_dbcsr_p_type), POINTER, DIMENSION(:) :: inverse_history
       TYPE(cp_para_env_type),  POINTER   :: para_env
       TYPE(cp_blacs_env_type), POINTER   :: ctxt
       INTEGER :: in_use, solver, ihistory
       REAL(KIND = dp), DIMENSION(:), POINTER :: occ_evals,full_evals
       REAL(KIND = dp) :: energy_gap
       REAL(KIND = dp) :: condition_num
  END TYPE preconditioner_type

! *****************************************************************************
  TYPE preconditioner_p_type
     TYPE(preconditioner_type), POINTER :: preconditioner
  END TYPE preconditioner_p_type

CONTAINS

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param preconditioner ...
!> \retval preconditioner_in_use ...
! *****************************************************************************
  FUNCTION preconditioner_in_use(preconditioner)
    TYPE(preconditioner_type)                :: preconditioner
    LOGICAL                                  :: preconditioner_in_use

    preconditioner_in_use = .NOT. (preconditioner%in_use .EQ. 0)
  END FUNCTION

! *****************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param para_env ...
!> \param blacs_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE init_preconditioner(preconditioner_env,para_env,blacs_env, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    NULLIFY(preconditioner_env%sparse_matrix)
    NULLIFY(preconditioner_env%fm)
    NULLIFY(preconditioner_env%dbcsr_matrix)
    NULLIFY(preconditioner_env%occ_evals)
    NULLIFY(preconditioner_env%full_evals)
    NULLIFY(preconditioner_env%inverse_history)
    NULLIFY(preconditioner_env%max_ev_vector)
    NULLIFY(preconditioner_env%min_ev_vector)
    preconditioner_env%solver=ot_precond_solver_default
    preconditioner_env%para_env => para_env
    preconditioner_env%ctxt     => blacs_env
    !inverse is used for filtering in update set it to something huge to
    ! avoid filtering if the information is not available
    preconditioner_env%condition_num=-1.0_dp
    preconditioner_env%ihistory=0

    CALL cp_para_env_retain(preconditioner_env%para_env,error=error)
    CALL cp_blacs_env_retain(preconditioner_env%ctxt,error=error)

  END SUBROUTINE init_preconditioner

! *****************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE destroy_preconditioner(preconditioner_env, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'destroy_preconditioner', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i

    CALL timeset(routineN,handle)

    IF (ASSOCIATED(preconditioner_env%sparse_matrix)) THEN
       CALL cp_dbcsr_deallocate_matrix(preconditioner_env%sparse_matrix,error=error)
       NULLIFY(preconditioner_env%sparse_matrix)
    ENDIF

    IF (ASSOCIATED(preconditioner_env%fm)) THEN
       CALL cp_fm_release(preconditioner_env%fm,error=error)
    ENDIF
    IF (ASSOCIATED(preconditioner_env%dbcsr_matrix)) THEN
       CALL cp_dbcsr_release_p (preconditioner_env%dbcsr_matrix, error=error)
    ENDIF
    IF (ASSOCIATED(preconditioner_env%max_ev_vector)) THEN
       CALL cp_dbcsr_release_p (preconditioner_env%max_ev_vector, error=error)
    ENDIF
    IF (ASSOCIATED(preconditioner_env%min_ev_vector)) THEN
       CALL cp_dbcsr_release_p (preconditioner_env%min_ev_vector, error=error)
    ENDIF
    IF (ASSOCIATED(preconditioner_env%occ_evals)) THEN
       DEALLOCATE(preconditioner_env%occ_evals)
    ENDIF
    IF (ASSOCIATED(preconditioner_env%full_evals)) THEN
       DEALLOCATE(preconditioner_env%full_evals)
    ENDIF
    IF (ASSOCIATED(preconditioner_env%inverse_history)) THEN
       DO i=1,SIZE(preconditioner_env%inverse_history)
          CALL cp_dbcsr_release_p (preconditioner_env%inverse_history(i)%matrix, error=error)
       END DO
       DEALLOCATE(preconditioner_env%inverse_history)
    ENDIF
    CALL cp_para_env_release(preconditioner_env%para_env,error=error)
    CALL cp_blacs_env_release(preconditioner_env%ctxt,error=error)

    preconditioner_env%in_use=0

    CALL timestop(handle)

  END SUBROUTINE destroy_preconditioner

END MODULE preconditioner_types

