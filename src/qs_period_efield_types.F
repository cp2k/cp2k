!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief type for berry phase efield matrices. At the moment only used for
!>        cosmat and sinmat
!> \par History
!>      none
!> \author fschiff (06.2010)
! *****************************************************************************

MODULE qs_period_efield_types

  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_deallocate_matrix_set,&
                                             cp_dbcsr_p_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_period_efield_types'

  PUBLIC:: efield_berry_type, efield_berry_release, init_efield_matrices, &
           set_efield_matrices,deriv_p_type

  TYPE deriv_p_type
     TYPE(cp_dbcsr_p_type),DIMENSION(:), POINTER           :: deriv
  END TYPE

  TYPE efield_berry_type
     TYPE(cp_dbcsr_p_type),DIMENSION(:), POINTER           :: cosmat
     TYPE(cp_dbcsr_p_type),DIMENSION(:), POINTER           :: sinmat
     TYPE(deriv_p_type),DIMENSION(:), POINTER           :: dcosmat
     TYPE(deriv_p_type),DIMENSION(:), POINTER           :: dsinmat
  END TYPE efield_berry_type

CONTAINS


  SUBROUTINE init_efield_matrices(efield,error)
    TYPE(efield_berry_type), POINTER         :: efield
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_efield_matrices', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF(ASSOCIATED(efield))CALL efield_berry_release(efield,error)
    ALLOCATE(efield,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    NULLIFY(efield%dcosmat)
    NULLIFY(efield%dsinmat)
    NULLIFY(efield%cosmat)
    NULLIFY(efield%sinmat)

  END SUBROUTINE init_efield_matrices

  SUBROUTINE set_efield_matrices(efield,sinmat,cosmat,dsinmat,dcosmat,fielddir)

    TYPE(efield_berry_type), POINTER         :: efield
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: sinmat, cosmat, dsinmat, &
                                                dcosmat
    INTEGER, OPTIONAL                        :: fielddir

    CHARACTER(len=*), PARAMETER :: routineN = 'set_efield_matrices', &
      routineP = moduleN//':'//routineN

     IF(PRESENT(cosmat))efield%cosmat=>cosmat
     IF(PRESENT(sinmat))efield%sinmat=>sinmat
     IF(PRESENT(dcosmat))efield%dcosmat(fielddir)%deriv=>dcosmat
     IF(PRESENT(dsinmat))efield%dsinmat(fielddir)%deriv=>dsinmat


  END SUBROUTINE set_efield_matrices

  SUBROUTINE efield_berry_release(efield,error)
    TYPE(efield_berry_type), POINTER         :: efield
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'efield_berry_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure =.FALSE.
    IF(ASSOCIATED(efield))THEN
       IF(ASSOCIATED(efield%sinmat).AND.ASSOCIATED(efield%cosmat))THEN
          CALL cp_dbcsr_deallocate_matrix_set ( efield%cosmat, error )
          CALL cp_dbcsr_deallocate_matrix_set ( efield%sinmat, error )
       END IF
       IF(ASSOCIATED(efield%dsinmat).AND.ASSOCIATED(efield%dcosmat))THEN
          DO i=1,3
             CALL cp_dbcsr_deallocate_matrix_set(efield%dsinmat(i)%deriv,error)
             CALL cp_dbcsr_deallocate_matrix_set(efield%dcosmat(i)%deriv,error)
          END DO
          DEALLOCATE(efield%dcosmat,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(efield%dsinmat,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       DEALLOCATE(efield,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE efield_berry_release

END MODULE qs_period_efield_types
