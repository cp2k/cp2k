!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

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
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_equivalent,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_retain,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_init_random,&
                                             cp_fm_release,&
                                             cp_fm_retain,&
                                             cp_fm_set_all,&
                                             cp_fm_set_element,&
                                             cp_fm_write,&
                                             cp_fm_p_type,&
                                             cp_fm_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_env,                     ONLY: cp_para_env_from_globenv,&
                                             cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_fm_vect'

  PUBLIC :: cp_fm_vect_dealloc, cp_fm_vect_copy,&
       cp_fm_vect_write, cp_fm_vect_set_all
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
    TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: matrixes
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='cp_fm_vect_dealloc',&
         routineP=moduleN//':'//routineN
    INTEGER :: stat, i
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
  TYPE(cp_fm_p_type), DIMENSION(:), INTENT(in) :: matrixes
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: copy
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='cp_fm_vect_copy',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat,i
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

!!****f* cp_fm_vect/cp_fm_vect_write [1.0] *
!!
!!   NAME
!!     cp_fm_vect_write
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_vect_write(matrixes, unit_nr, long_description,&
!!         local, error)
!!       Type(cp_fm_p_type), Dimension(:), Pointer:: matrixes
!!       Integer, Intent (IN):: unit_nr
!!       Logical, Optional, Intent (IN):: long_description, local
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_fm_vect_write
!!
!!   FUNCTION
!!     writes an array of full matrixes to the given output unit
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - matrixes: the array of matrixes to output
!!     - unit_nr: the unit where to output
!!     - long_description: if a long description should be written
!!       (defaults to false)
!!     - local: if the unit is a local unit or a global unit 
!!       (defaults to false, i.e. global)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE cp_fm_vect_write(matrixes, unit_nr, long_description,&
       local, error)
    TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: matrixes
    INTEGER, INTENT(in) :: unit_nr
    LOGICAL, OPTIONAL, INTENT(in) :: long_description, local
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure, my_local, my_long_description, should_w
  CHARACTER(len=*), PARAMETER :: routineN='cp_fm_vect_write',&
        routineP=moduleN//':'//routineN
  INTEGER :: i
  TYPE(cp_logger_type), POINTER :: logger

  failure=.FALSE.; my_local=.FALSE.; my_long_description=.FALSE.
  logger => cp_error_get_logger(error)
  
  IF (PRESENT(local)) my_local=local
  IF (PRESENT(long_description)) my_long_description=long_description
  
  should_w=my_local .OR. logger%para_env%mepos==logger%para_env%source

  IF (ASSOCIATED(matrixes)) THEN
     IF (should_w) WRITE (unit=unit_nr,fmt="(a)") '('
     DO i=1,SIZE(matrixes)
        CALL cp_fm_write(matrix=matrixes(i)%matrix, unit_nr=unit_nr,&
             long_description=my_long_description,&
             local=my_local,error=error)
        IF (i/=SIZE(matrixes).and.should_w) WRITE(unit=unit_nr,fmt="(a)")','
     END DO
     IF (should_w) WRITE (unit=unit_nr,fmt="(a)") ')'
  ELSE
     IF (should_w) THEN
        WRITE (unit=unit_nr,fmt="(a)") ' <cp_fm_p array>:*null*'
     END IF
  END IF
END SUBROUTINE cp_fm_vect_write
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
  TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: matrixes
  REAL(KIND = dp), INTENT(in), OPTIONAL :: value
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='cp_fm_vect_set_all',&
        routineP=moduleN//':'//routineN
  REAL(KIND = dp) :: my_value
  INTEGER :: i

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
