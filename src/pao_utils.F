!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/pao_utils [1.0] *
!!
!!   NAME
!!     pao_utils
!!
!!   FUNCTION
!!     all kinds of utilities
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!*****************************************************************************
MODULE pao_utils
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE sparse_matrix_types
  IMPLICIT NONE
  PRIVATE
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER :: moduleN='pao_utils'

! =========== public objects ==========
! types

! procedures
  PUBLIC :: cp_get_random_integer


! ============== interfaces ===========

!!***
!****************************************************************************


CONTAINS

!!****f* pao_utils/cp_get_random_integer [1.0] *
!!
!!   NAME
!!     cp_get_random_integer
!!
!!   SYNOPSIS
!!     Function cp_get_random_integer(from, to, error) Result(res)
!!       Integer:: res
!!       Integer, Optional, Intent (IN):: from, to
!!       Type(cp_error_type), Optional, Intent (IN):: error
!!     End Function cp_get_random_integer
!!
!!   FUNCTION
!!     This function returns a random integer in the requested range
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     from: minumum random number defaults to 1
!!     to: maximum random number defaults to 100
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling, but does NOT report
!!            the error in this variable (to make it pure)
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!****************************************************************************
  !pure
  FUNCTION cp_get_random_integer(from,to,error) RESULT(res)
    INTEGER :: res
    INTEGER, OPTIONAL, INTENT(in) :: from, to
    TYPE(cp_error_type), OPTIONAL, INTENT(in) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineP='util:cp_get_random_integer'
    REAL :: rand
    INTEGER :: fromi, toi
    TYPE(cp_error_type) :: iError
    failure=.FALSE.
!call cp_init(iError,template_error=error)
    IF (.NOT. PRESENT(from)) THEN
       fromi=1
    ELSE
       fromi=from
    END IF
    IF (.NOT. PRESENT(to)) THEN
       toi=100
    ELSE
       toi=to
    END IF
!CPPrecondition(fromi<=toi,cp_failure_level,routineP,iError,failure)
    IF (.not.fromi<=toi) failure=.TRUE.
    IF (.NOT. failure) THEN
       CALL RANDOM_NUMBER(rand)
       res=fromi+FLOOR(rand*REAL(toi-fromi+1))
!CPPostcondition(res<=toi.and.res>=fromi,cp_warning_level,routineP,iError,failure)
       IF (res>toi) res=toi
       IF (res<fromi) res=fromi
    ELSE
       res=-1
    END IF
!call cp_dealloc_ref(iError,error=error)
  END FUNCTION cp_get_random_integer
!!***
!****************************************************************************



END MODULE pao_utils

