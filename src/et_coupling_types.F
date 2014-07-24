!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Definition and initialisation of the et_coupling data type.
!> \author Florian Schiffmann (01.2007,fschiff)
! *****************************************************************************
MODULE et_coupling_types

  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_p_type
  USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                             cp_fm_release
  USE f77_blas
  USE kinds,                           ONLY: dp
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'et_coupling_types'

! *** Public data types ***

  PUBLIC :: et_coupling_type

! *** Public subroutines ***

  PUBLIC :: et_coupling_create,&
            et_coupling_release,&
            set_et_coupling_type

! *****************************************************************************
!> \par History
!>      01.2007 created [Florian Schiffmann]
!> \author fschiff
! *****************************************************************************
  TYPE et_coupling_type
    TYPE(cp_fm_p_type),DIMENSION(:), POINTER           :: et_mo_coeff
    TYPE(cp_dbcsr_p_type),DIMENSION(:), POINTER     :: rest_mat
    LOGICAL                                            :: first_run
    LOGICAL                                            :: keep_matrix
    REAL(KIND = dp)                                    :: energy,e1,order_p
  END TYPE

  CONTAINS

! *****************************************************************************
!> \brief ...
!> \param et_coupling ...
!> \param error ...
! *****************************************************************************
    SUBROUTINE et_coupling_create(et_coupling,error)
    TYPE(et_coupling_type), POINTER          :: et_coupling
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'et_coupling_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

      failure=.FALSE.

      ALLOCATE(et_coupling, stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      IF (.NOT. failure) THEN
         NULLIFY(et_coupling%et_mo_coeff)
         NULLIFY(et_coupling%rest_mat)
         et_coupling%first_run=.TRUE.
         et_coupling%keep_matrix=.FALSE.
         ALLOCATE(et_coupling%rest_mat(2))
      END IF

    END SUBROUTINE et_coupling_create

! *****************************************************************************
!> \brief ...
!> \param et_coupling ...
!> \param et_mo_coeff ...
!> \param rest_mat ...
!> \param error ...
! *****************************************************************************
    SUBROUTINE get_et_coupling_type(et_coupling,et_mo_coeff,rest_mat,error)
    TYPE(et_coupling_type), POINTER          :: et_coupling
    TYPE(cp_fm_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: et_mo_coeff
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rest_mat
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_et_coupling_type', &
      routineP = moduleN//':'//routineN

      IF(PRESENT(et_mo_coeff))et_mo_coeff => et_coupling%et_mo_coeff
      IF(PRESENT(rest_mat))rest_mat => et_coupling%rest_mat

    END SUBROUTINE get_et_coupling_type

! *****************************************************************************
!> \brief ...
!> \param et_coupling ...
!> \param et_mo_coeff ...
!> \param rest_mat ...
!> \param error ...
! *****************************************************************************
    SUBROUTINE set_et_coupling_type(et_coupling,et_mo_coeff,rest_mat,error)
    TYPE(et_coupling_type), POINTER          :: et_coupling
    TYPE(cp_fm_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: et_mo_coeff
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rest_mat
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_et_coupling_type', &
      routineP = moduleN//':'//routineN

      IF(PRESENT(et_mo_coeff))  et_coupling%et_mo_coeff  = et_mo_coeff
      IF(PRESENT(rest_mat)) et_coupling%rest_mat => rest_mat

    END SUBROUTINE set_et_coupling_type

! *****************************************************************************
!> \brief ...
!> \param et_coupling ...
!> \param error ...
! *****************************************************************************
    SUBROUTINE et_coupling_release(et_coupling,error)
    TYPE(et_coupling_type), POINTER          :: et_coupling
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'et_coupling_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

      IF(ASSOCIATED(et_coupling%et_mo_coeff))THEN
         DO i=1,SIZE(et_coupling%et_mo_coeff)
            CALL cp_fm_release(et_coupling%et_mo_coeff(i)%matrix,error=error)
         END DO
         DEALLOCATE(et_coupling%et_mo_coeff)
      END IF
      IF(ASSOCIATED(et_coupling%rest_mat))THEN
!         CALL deallocate_matrix_set(et_coupling%rest_mat,error)
         DEALLOCATE(et_coupling%rest_mat)
      END IF

      DEALLOCATE(et_coupling,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END SUBROUTINE et_coupling_release

END MODULE et_coupling_types

