!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/kg_scf_types [1.0] *
!!
!!   NAME
!!     kg_scf_types
!!
!!   FUNCTION
!!     module that contains the definitions of the scf types
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     fawzi/jgh
!!
!!   MODIFICATION HISTORY
!!
!!***
!****************************************************************************
MODULE kg_scf_types
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

  IMPLICIT NONE
  PRIVATE

  INTEGER, SAVE, PRIVATE :: last_scf_env_id=0

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='kg_scf_types'

  PUBLIC :: kg_scf_env_type, kg_scf_env_p_type,&
            kg_scf_env_create, kg_scf_env_retain, kg_scf_env_release

!****************************************************************************

!!****s* kg_scf_types/kg_scf_env_type [1.0] *
!!
!!   NAME
!!     kg_scf_env_type
!!
!!   FUNCTION
!!     wrapper for temporary and cached objects used in the scf iteration
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     fawzi/jgh
!!
!!   MODIFICATION HISTORY
!!
!!***
!***************************************************************************

TYPE kg_scf_env_type
   INTEGER :: id_nr
   INTEGER :: ref_count
   INTEGER :: iter_count
   INTEGER :: print_count
   LOGICAL :: print_iter_line
   INTEGER :: method
   INTEGER :: mixing_method
   REAL(KIND = dp) :: iter_param
   REAL(KIND = dp) :: iter_delta
   REAL(KIND = dp) :: p_mix
END TYPE kg_scf_env_type

!****************************************************************************
!!****s* kg_scf_types/kg_scf_env_p_type [1.0] *
!!
!!   NAME
!!     kg_scf_env_p_type
!!
!!   FUNCTION
!!     type to have arrays of pointers
!!
!!   AUTHOR
!!     fawzi/jgh
!!
!!   MODIFICATION HISTORY
!!
!!***
!***************************************************************************

  TYPE kg_scf_env_p_type
     TYPE(kg_scf_env_type), POINTER :: scf_env
  END TYPE kg_scf_env_p_type

!****************************************************************************

CONTAINS

!!****f* kg_scf_types/kg_scf_env_create [1.0] *
!!
!!   NAME
!!     kg_scf_env_create
!!
!!   AUTHOR
!!     fawzi/jgh
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************
  SUBROUTINE kg_scf_env_create(scf_env, error)
    TYPE(kg_scf_env_type), POINTER           :: scf_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_scf_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    ALLOCATE(scf_env, stat=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    IF (.NOT. failure) THEN
       scf_env%ref_count=1
       scf_env%print_count=0
       last_scf_env_id=last_scf_env_id+1
       scf_env%id_nr=last_scf_env_id
       scf_env%print_count=0
       scf_env%iter_count=0
!  INTEGER :: method
!  INTEGER :: mixing_method
!  REAL(KIND = dp) :: iter_param
!  REAL(KIND = dp) :: iter_delta
!  REAL(KIND = dp) :: p_mix
    END IF

  END SUBROUTINE kg_scf_env_create
!***************************************************************************

!!****f* kg_scf_types/kg_scf_env_retain [1.0] *
!!
!!   NAME
!!     kg_scf_env_retain
!!
!!   FUNCTION
!!     retains an scf_env (see doc/ReferenceCounting.html)
!!
!!   AUTHOR
!!     fawzi/jgh
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************
SUBROUTINE kg_scf_env_retain(scf_env,error)
    TYPE(kg_scf_env_type), POINTER           :: scf_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_scf_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
     scf_env%ref_count=scf_env%ref_count+1
  END IF
END SUBROUTINE kg_scf_env_retain
!***************************************************************************

!!****f* kg_scf_types/kg_scf_env_release [1.0] *
!!
!!   NAME
!!     kg_scf_env_release
!!
!!   FUNCTION
!!     releases an scf_env (see doc/ReferenceCounting.html)
!!
!!   AUTHOR
!!     fawzi/jgh
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************
SUBROUTINE kg_scf_env_release(scf_env,error)
    TYPE(kg_scf_env_type), POINTER           :: scf_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_scf_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (ASSOCIATED(scf_env)) THEN
     CPPreconditionNoFail(scf_env%ref_count>0,cp_failure_level,routineP,error)
     scf_env%ref_count=scf_env%ref_count-1
     IF (scf_env%ref_count==0) THEN
        DEALLOCATE(scf_env, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(scf_env)
END SUBROUTINE kg_scf_env_release
!***************************************************************************

END MODULE kg_scf_types
