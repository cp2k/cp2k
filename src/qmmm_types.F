!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qmmm_types *
!!
!!   NAME
!!     qmmm_types
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     05.2004 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE qmmm_types
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qmmm_types'
  INTEGER, SAVE, PRIVATE :: last_qmmm_env_id_nr=0

  PUBLIC :: qmmm_env_type
  PUBLIC :: qmmm_env_retain, qmmm_env_release, qmmm_env_create

!***
!****************************************************************************

!!****s* qmmm_types/qmmm_env_type [1.0] *
!!
!!   NAME
!!     qmmm_env_type
!!
!!   FUNCTION
!!     variables needed for QM/MM
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     05.2004 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE qmmm_env_type
     INTEGER :: ref_count, id_nr
     INTEGER, DIMENSION(:), POINTER :: qm_atom_index
     INTEGER, DIMENSION(:), POINTER :: qm_subsys_index
  END TYPE qmmm_env_type
!!***
  !****************************************************************************

CONTAINS

!!****f* qmmm_types/qmmm_env_create [1.0] *
!!
!!   NAME
!!     qmmm_env_create
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!*** **********************************************************************
SUBROUTINE qmmm_env_create(qmmm_env, error)
    TYPE(qmmm_env_type), POINTER             :: qmmm_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  CPPrecondition(.NOT.ASSOCIATED(qmmm_env),cp_failure_level,routineP,error,failure)
  ALLOCATE(qmmm_env, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     qmmm_env%ref_count=1
     last_qmmm_env_id_nr=last_qmmm_env_id_nr+1
     qmmm_env%id_nr=last_qmmm_env_id_nr
     NULLIFY(qmmm_env%qm_atom_index,qmmm_env%qm_subsys_index)
  END IF
END SUBROUTINE qmmm_env_create
!***************************************************************************

!!****f* qmmm_types/qmmm_env_retain *
!!
!!   NAME
!!     qmmm_env_retain
!!
!!   FUNCTION
!!     retains the qmmm_env
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qmmm_env: the qmmm_env to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!*** **********************************************************************
SUBROUTINE qmmm_env_retain(qmmm_env,error)
    TYPE(qmmm_env_type), POINTER             :: qmmm_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(qmmm_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(qmmm_env%ref_count>0,cp_failure_level,routineP,error)
     qmmm_env%ref_count=qmmm_env%ref_count+1
  END IF
END SUBROUTINE qmmm_env_retain
!***************************************************************************


!!****f* qmmm_types/qmmm_env_release [1.0] *
!!
!!   NAME
!!     qmmm_env_release
!!
!!   FUNCTION
!!     releases the given qmmm_env (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qmmm_env: the object to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!*** **********************************************************************
SUBROUTINE qmmm_env_release(qmmm_env,error)
    TYPE(qmmm_env_type), POINTER             :: qmmm_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  IF (ASSOCIATED(qmmm_env)) THEN
     CPPreconditionNoFail(qmmm_env%ref_count>0,cp_failure_level,routineP,error)
     qmmm_env%ref_count=qmmm_env%ref_count-1
     IF (qmmm_env%ref_count==0) THEN
        IF (ASSOCIATED(qmmm_env%qm_atom_index)) THEN
           DEALLOCATE(qmmm_env%qm_atom_index,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(qmmm_env%qm_subsys_index)) THEN
           DEALLOCATE(qmmm_env%qm_subsys_index,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        DEALLOCATE(qmmm_env, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(qmmm_env)
END SUBROUTINE qmmm_env_release
!***************************************************************************


END MODULE qmmm_types
