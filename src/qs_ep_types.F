!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qs_ep_types [1.0] *
!!
!!   NAME
!!     qs_ep_types
!!
!!   FUNCTION
!!     the types needed for the calculation of the enrgy with perturbation
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE qs_ep_types
  USE cp_error_handling,               ONLY: cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: wp=>dp
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_p_type,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_p_sparse_psi,                 ONLY: p_proj_release,&
                                             p_proj_retain,&
                                             qs_p_projection_p_type,&
                                             qs_p_projection_type
  USE qs_p_types,                      ONLY: p_env_release,&
                                             p_env_retain,&
                                             qs_p_env_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_ep_types'

  PUBLIC :: qs_ep_env_type, qs_ep_env_p_type
  PUBLIC :: qs_ep_retain, qs_ep_release
!!***
!****************************************************************************

!!****s* qs_p_ep_types/qs_ep_env_type [1.0] *
!!
!!   NAME
!!     qs_ep_env_type
!!
!!   FUNCTION
!!     environement for an ep calculation
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - id_nr: identification number (unique)
!!     - ref_count: reference count (see doc/ReferenceCounting.html)
!!     - sub_qs_env: qs_env that represent the pieces of the full system
!!     - sub_proj: projection from the full system to the pieces
!!     - main_qs_env: qs_env of the full system
!!     - main_p_env: perturbation environment for the full system
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2002 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE qs_ep_env_type
     INTEGER :: id_nr, ref_count
     TYPE(qs_environment_p_type), DIMENSION(:), POINTER :: sub_qs_env
     TYPE(qs_p_projection_p_type), DIMENSION(:), POINTER :: sub_proj
     TYPE(qs_environment_type), POINTER :: main_qs_env
     TYPE(qs_p_env_type), POINTER :: main_p_env
  END TYPE qs_ep_env_type
!!***
  !****************************************************************************

!!****s* qs_ep_types/qs_ep_env_p_type [1.0] *
!!
!!   NAME
!!     qs_ep_env_p_type
!!
!!   FUNCTION
!!     type to have array of pointers
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - ep_env: the pointer to the ep env
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE qs_ep_env_p_type
     TYPE(qs_ep_env_type), POINTER :: ep_env
  END TYPE qs_ep_env_p_type
!!***
  !****************************************************************************

CONTAINS

!!****f* qs_ep_types/qs_ep_retain [1.0] *
!!
!!   NAME
!!     qs_ep_retain
!!
!!   SYNOPSIS
!!     Subroutine qs_ep_retain(ep_env, error)
!!       Type(qs_ep_env_type), Pointer:: ep_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine qs_ep_retain
!!
!!   FUNCTION
!!     retains an ep environment (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - ep_env: the environment to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE qs_ep_retain(ep_env, error)
  TYPE(qs_ep_env_type), POINTER :: ep_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='qs_ep_retain',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(ep_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(ep_env%ref_count>0,cp_failure_level,routineP,error,failure)
     ep_env%ref_count=ep_env%ref_count+1
  END IF
END SUBROUTINE qs_ep_retain
!***************************************************************************

!!****f* qs_ep_types/qs_ep_release [1.0] *
!!
!!   NAME
!!     qs_ep_release
!!
!!   SYNOPSIS
!!     Subroutine qs_ep_release(ep_env, error)
!!       Type(qs_ep_env_type), Pointer:: ep_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine qs_ep_release
!!
!!   FUNCTION
!!     releases the given ep_env (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - ep_env: the ep_env to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE qs_ep_release(ep_env, error)
  TYPE(qs_ep_env_type), POINTER :: ep_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='qs_ep_release',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat, i
  failure=.FALSE.
  
  IF (ASSOCIATED(ep_env)) THEN
     CPPrecondition(ep_env%ref_count>0,cp_failure_level,routineP,error,failure)
     ep_env%ref_count=ep_env%ref_count-1
     IF (ep_env%ref_count<1) THEN
        IF (ASSOCIATED(ep_env%sub_qs_env)) THEN ! to do
           DEALLOCATE(ep_env%sub_qs_env, stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(ep_env%sub_proj)) THEN
           DO i=1,SIZE(ep_env%sub_proj)
              CALL p_proj_release(ep_env%sub_proj(i)%projection,error=error)
           END DO
           DEALLOCATE(ep_env%sub_proj, stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(ep_env%main_qs_env)) THEN ! to do
           DEALLOCATE(ep_env%main_qs_env,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        CALL p_env_release(ep_env%main_p_env,error=error)
        DEALLOCATE(ep_env, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
END SUBROUTINE qs_ep_release
!***************************************************************************

END MODULE qs_ep_types
