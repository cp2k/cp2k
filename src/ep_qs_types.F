!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2006  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/ep_qs_types [1.0] *
!!
!!   NAME
!!     ep_qs_types
!!
!!   FUNCTION
!!     types for the modifications of the qs calculation to perform an ep
!!     calculation
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     06.2006 created
!!
!!   SOURCE
!****************************************************************************
MODULE ep_qs_types
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ep_qs_types'

  PUBLIC :: ep_qs_type, ep_qs_p_type
  PUBLIC :: ep_qs_retain, ep_qs_release
!!***
!****************************************************************************

!!****s* ep_qs_types/ep_qs_type *
!!
!!   NAME
!!     ep_qs_type
!!
!!   FUNCTION
!!     environement for changes to the normal qs calculation done during an ep
!!     calculation
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - id_nr: identification number for the ep_qs_type
!!     - ref_count: reference count (see cp2k/doc/ReferenceCounting.html)
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
TYPE ep_qs_type
   INTEGER :: id_nr, ref_count
END TYPE ep_qs_type
!***************************************************************************

!!****s* ep_qs_types/ep_qs_p_type *
!!
!!   NAME
!!     ep_qs_p_type
!!
!!   FUNCTION
!!     pointer to an ep_qs_type for arrays of pointers
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - ep_qs_env: the pointer to the ep_qs_type
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
TYPE ep_qs_p_type
   TYPE(ep_qs_type), POINTER :: ep_qs_env
END TYPE ep_qs_p_type
!***************************************************************************

CONTAINS

!!****f* ep_qs_types/ep_qs_retain *
!!
!!   NAME
!!     ep_qs_retain
!!
!!   FUNCTION
!!     retains the given ep_qs_env
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - ep_qs_env: the object to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE ep_qs_retain(ep_qs_env,error)
    TYPE(ep_qs_type), POINTER                :: ep_qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_qs_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(ep_qs_env),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     CPPrecondition(ep_qs_env%ref_count>0,cp_failure_level,routineP,error,failure)
     ep_qs_env%ref_count=ep_qs_env%ref_count+1
  END IF
END SUBROUTINE ep_qs_retain
!***************************************************************************

!!****f* ep_qs_types/ep_qs_release *
!!
!!   NAME
!!     ep_qs_release
!!
!!   FUNCTION
!!     releases the given ep_qs_env
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - ep_qs_env: the object to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE ep_qs_release(ep_qs_env,error)
    TYPE(ep_qs_type), POINTER                :: ep_qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_qs_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  IF (ASSOCIATED(ep_qs_env)) THEN
     CPPrecondition(ep_qs_env%ref_count>0,cp_failure_level,routineP,error,failure)
     ep_qs_env%ref_count=ep_qs_env%ref_count-1
     IF(ep_qs_env%ref_count==0) THEN
        DEALLOCATE(ep_qs_env,stat=stat)
     END IF
  END IF
  NULLIFY(ep_qs_env)
END SUBROUTINE ep_qs_release
!***************************************************************************

END MODULE ep_qs_types
!****************************************************************************
