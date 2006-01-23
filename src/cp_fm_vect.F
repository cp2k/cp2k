!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****h* cp2k/cp_fm_vect [1.0] *
!!
!!   NAME
!!     cp_fm_vect
!!
!!   FUNCTION
!!     routine to handle vectors of full matrixes
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE cp_fm_vect
  USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_retain,&
                                             cp_fm_set_all
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_fm_vect'

  PUBLIC :: cp_fm_vect_dealloc, cp_fm_vect_copy,&
            cp_fm_vect_set_all
!***
!****************************************************************************
CONTAINS

!!****f* cp_fm_vect/cp_fm_vect_dealloc [1.0] *
!!
!!   NAME
!!     cp_fm_vect_dealloc
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_vect_dealloc(matrixes, error)
!!       Type(cp_fm_p_type), Dimension(:), Pointer:: matrixes
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_fm_vect_dealloc
!!
!!   FUNCTION
!!     deallocate an array of pointers to blacs matrixes
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - matrixes: the array of matrixes to deallocate
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE cp_fm_vect_dealloc(matrixes, error)
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: matrixes
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
!***************************************************************************

!!****f* cp_fm_vect/cp_fm_vect_copy[1.0] *
!!
!!   NAME
!!     cp_fm_vect_copy
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_vect_copy(matrixes, copy, error)
!!       Type(cp_fm_p_type), Dimension(:), Intent (IN):: matrixes
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Type(cp_fm_p_type), Dimension(:), Pointer:: copy
!!     End Subroutine cp_fm_vect_copy
!!
!!   FUNCTION
!!     Does a shallow copy of an array of full matrices (i.e. just retains 
!!     the matrices)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - matrixes: the matrixes to copy
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     09.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_fm_vect_copy(matrixes, copy, error)
    TYPE(cp_fm_p_type), DIMENSION(:), &
      INTENT(in)                             :: matrixes
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: copy
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_vect_copy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  ALLOCATE(copy(SIZE(matrixes)),stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     DO i=1,SIZE(matrixes)
        copy(i)%matrix => matrixes(i)%matrix
        CALL cp_fm_retain(matrixes(i)%matrix,error=error)
     END DO
  END IF
END SUBROUTINE cp_fm_vect_copy
!***************************************************************************

!!****f* cp_fm_vect/cp_fm_vect_set_all [1.0] *
!!
!!   NAME
!!     cp_fm_vect_set_all
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_vect_set_all(matrixes, value, error)
!!       Type(cp_fm_p_type), Dimension(:), Pointer:: matrixes
!!       Real(KIND = dp), Intent (IN), Optional:: value
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_fm_vect_set_all
!!
!!   FUNCTION
!!     sets all the elements of the matrixes to the given value
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - matrixes: the matrixes to change
!!     - value: the value that the element should receive (defaults to 0)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_fm_vect_set_all(matrixes,value,error)
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: matrixes
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: value
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_vect_set_all', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: my_value

  failure=.FALSE.
  my_value=0.0_dp
  
  IF (PRESENT(value)) my_value=value
  
  CPPrecondition(ASSOCIATED(matrixes),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     DO i=1,SIZE(matrixes)
        CPPrecondition(ASSOCIATED(matrixes(i)%matrix),cp_failure_level,routineP,error,failure)
        CALL cp_fm_set_all(matrixes(i)%matrix,my_value,error=error)
     END DO
  END IF
END SUBROUTINE cp_fm_vect_set_all
!***************************************************************************

END MODULE cp_fm_vect
