!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
MODULE qs_oce_types

  USE block_p_types,                   ONLY: block_p_type
  USE kinds,                           ONLY: dp
  USE sap_kind_types,                  ONLY: release_sap_int,&
                                             sap_int_type
  USE termination,                     ONLY: stop_memory
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_oce_types'

! *** Define a oce matrix type ***

! *****************************************************************************
  TYPE qlist_type
    TYPE(block_p_type), DIMENSION(:), POINTER :: sac
    REAL(KIND = dp), DIMENSION(:), POINTER    :: r2
    REAL(KIND = dp), DIMENSION(:,:), POINTER  :: r
    INTEGER                                   :: n
    INTEGER, DIMENSION(:), POINTER            :: index_list,list
  END TYPE qlist_type

! *****************************************************************************
  TYPE oce_matrix_type
    TYPE(sap_int_type), DIMENSION(:), POINTER        :: intac
  END TYPE

! *** Public data types ***

  PUBLIC :: oce_matrix_type

! *** Public subroutines ***

  PUBLIC :: allocate_oce_set,&
            create_oce_set,&
            deallocate_oce_set

CONTAINS

! *****************************************************************************
!> \brief   Allocate and initialize the matrix set of oce coefficients.
!> \param oce_set ...
!> \param natom ...
!> \param nkind ...
!> \param ndim ...
!> \param error ...
!> \version 1.0
! *****************************************************************************
  SUBROUTINE allocate_oce_set(oce_set,natom,nkind,ndim,error)
    TYPE(oce_matrix_type), POINTER           :: oce_set
    INTEGER, INTENT(IN)                      :: natom, nkind, ndim
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_oce_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, istat
    LOGICAL                                  :: failure

    failure = .FALSE.

    ALLOCATE(oce_set%intac(nkind*nkind),STAT=istat)
    CPPrecondition(istat==0,cp_warning_level,routineP,error,failure)
    DO i=1,nkind*nkind
      NULLIFY(oce_set%intac(i)%alist)
      NULLIFY(oce_set%intac(i)%asort)
      NULLIFY(oce_set%intac(i)%aindex)
    END DO

  END SUBROUTINE allocate_oce_set

! *****************************************************************************
!> \brief ...
!> \param oce_set ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE create_oce_set(oce_set,error)

    TYPE(oce_matrix_type), POINTER           :: oce_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_oce_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    IF(ASSOCIATED(oce_set)) CALL deallocate_oce_set(oce_set,error=error)

    ALLOCATE (oce_set,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"oce_set",0)

    NULLIFY(oce_set%intac)

  END SUBROUTINE create_oce_set

! *****************************************************************************
!> \brief  Deallocate the matrix set of oce coefficients
!> \param oce_set ...
!> \param error ...
!> \date
!> \author
!> \version 1.0
! *****************************************************************************
  SUBROUTINE deallocate_oce_set(oce_set,error)
    TYPE(oce_matrix_type), POINTER           :: oce_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_oce_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    IF (.NOT.ASSOCIATED(oce_set)) RETURN

    IF(ASSOCIATED(oce_set%intac)) CALL release_sap_int(oce_set%intac,error=error)

    DEALLOCATE (oce_set,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "deallocate_oce_set",0)

  END SUBROUTINE deallocate_oce_set

! *****************************************************************************

END MODULE qs_oce_types
