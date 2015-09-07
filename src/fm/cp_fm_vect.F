!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief routine to handle vectors of full matrixes
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE cp_fm_vect
  USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_retain
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_vect'

  PUBLIC :: cp_fm_vect_dealloc, cp_fm_vect_copy
!***
CONTAINS

! *****************************************************************************
!> \brief deallocate an array of pointers to blacs matrixes
!> \param matrixes the array of matrixes to deallocate
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cp_fm_vect_dealloc(matrixes, error)
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: matrixes
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_vect_dealloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(matrixes)) THEN
       DO i=1,SIZE(matrixes)
          CALL cp_fm_release(matrixes(i)%matrix,error=error)
       END DO
       DEALLOCATE(matrixes,stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    END IF
  END SUBROUTINE cp_fm_vect_dealloc

! *****************************************************************************
!> \brief Does a shallow copy of an array of full matrices (i.e. just retains
!>      the matrices)
!> \param matrixes the matrixes to copy
!> \param copy ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE cp_fm_vect_copy(matrixes, copy, error)
    TYPE(cp_fm_p_type), DIMENSION(:), &
      INTENT(in)                             :: matrixes
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: copy
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_vect_copy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  ALLOCATE(copy(SIZE(matrixes)),stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  DO i=1,SIZE(matrixes)
     copy(i)%matrix => matrixes(i)%matrix
     CALL cp_fm_retain(matrixes(i)%matrix,error=error)
  END DO
END SUBROUTINE cp_fm_vect_copy

END MODULE cp_fm_vect
