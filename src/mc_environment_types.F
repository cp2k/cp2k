!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2003  CP2K developers group                            !
!-----------------------------------------------------------------------------!
!!****s* cp2k/mc_environment_types [1.0] *
!!
!!   NAME
!!     mc_environment_types
!!
!!   FUNCTION
!!     contains the subroutines for dealing with the mc_env
!!
!!   AUTHOR
!!     MJM Oct. 15-2003
!!
!******************************************************************************

#include "cp_prep_globals.h"

MODULE mc_environment_types
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_error_type, &
                                             cp_assert
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_warning_level,&
                                             cp_to_string
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE force_env_types,                 ONLY: force_env_release,&
                                             force_env_retain,&
                                             force_env_type
  USE mc_types,                        ONLY: mc_simulation_parameters_type

  IMPLICIT NONE

  PRIVATE

 
  TYPE mc_environment_type
     INTEGER :: id_nr, ref_count, in_use
     TYPE (cp_para_env_type), POINTER :: para_env 
     TYPE ( mc_simulation_parameters_type ), POINTER :: mc_par
     TYPE ( force_env_type ), POINTER :: force_env  
  END TYPE mc_environment_type

! *** Public subroutines and data types ***
  PUBLIC :: mc_environment_type, set_mc_env, mc_env_create,&
            get_mc_env, mc_env_release, mc_env_retain

! *** Global parameters ***

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='mc_environment_types'
  INTEGER, SAVE, PRIVATE :: last_mc_env_id=0

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************
!!****s* mc_environment_types.F/mc_env_create
!!
!!   NAME
!!     mc_env_create
!!
!!   SYNOPSIS
!!     Subroutine mc_env_create ( mc_env, para_env, error )
!!       Type(mc_env_type), Pointer :: mc_env
!!       Type(cp_para_env_type), Pointer :: para_env
!!     End Subroutine mc_env_create
!!
!!   FUNCTION
!!     creates and initializes an mc_env
!!
!!   ARGUMENTS
!!     - mc_env: the mc_environment you want to create
!!     - para_env: the cp para environment type you want to associate with
!!                 this mc_env
!!
!!   AUTHOR
!!     MJM
!!
!!**************************************************************************
  SUBROUTINE mc_env_create ( mc_env, para_env, error )

!   ***************************************************************************

    TYPE(mc_environment_type), POINTER       :: mc_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routine_name = "mc_env_create", &
      routineP = moduleN//':'//routine_name

    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------

    ALLOCATE ( mc_env, stat=istat )
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)

    last_mc_env_id=last_mc_env_id+1
    mc_env%id_nr=last_mc_env_id
    mc_env%ref_count=1
    mc_env%in_use=0

    NULLIFY ( mc_env % mc_par )
    NULLIFY ( mc_env % force_env )
    mc_env % para_env => para_env 
    CALL cp_para_env_retain ( mc_env % para_env )

  END SUBROUTINE mc_env_create

! *****************************************************************************

! *****************************************************************************
!!****s* mc_environment_types.F/set_mc_env
!!
!!   NAME
!!     set_mc_env
!!
!!   SYNOPSIS
!!     Subroutine set_mc_env ( mc_env, mc_par, force_env)
!!       Type(mc_env_type), Pointer :: mc_env
!!       Type(mc_simulation_parameters_type), Pointer,Optional:: mc_par
!!       Type(force_env_type), Pointer,Optional:: force_env
!!       Type(cp_para_env_type), Pointer,Optional:: para_env
!!     End Subroutine set_mc_env
!!
!!   FUNCTION
!!     provides a method for attaching various structures to an mc_env
!!
!!   ARGUMENTS
!!     - mc_env: the mc_environment you want to change
!!     - mc_par: the mc parameters you want to associate with this mc_env
!!     - force_env: the force environment type you want to associate
!!                  with this mc_env
!!     - para_env: the cp para environment type you want to associate
!!                 with this mc_env
!!
!!   AUTHOR
!!     MJM
!!
!!**************************************************************************
  SUBROUTINE set_mc_env ( mc_env,mc_par,force_env,para_env)

!   ***************************************************************************
    TYPE(mc_environment_type), POINTER   :: mc_env
    TYPE(cp_para_env_type), POINTER,OPTIONAL :: para_env
    TYPE(mc_simulation_parameters_type), &
      OPTIONAL, POINTER                      :: mc_par
    TYPE(force_env_type), OPTIONAL, POINTER  :: force_env

    IF ( PRESENT ( mc_par ) ) mc_env % mc_par => mc_par
    IF ( PRESENT ( para_env ) ) mc_env % para_env => para_env
    IF ( PRESENT ( force_env )) THEN
       CALL force_env_retain( force_env ) ! accept null pointers?
       IF (ASSOCIATED(mc_env % force_env ))  &
           CALL force_env_release( mc_env % force_env )
       mc_env % force_env => force_env
    END IF
    
!   ---------------------------------------------------------------------------

  END SUBROUTINE set_mc_env

! *****************************************************************************
!!****s* mc_environment_types.F/get_mc_env
!!
!!   NAME
!!     get_mc_env
!!
!!   SYNOPSIS
!!     Subroutine get_mc_env ( mc_env, mc_par, force_env)
!!       Type(mc_env_type), Pointer :: mc_env
!!       Type(mc_simulation_parameters_type), Pointer,Optional:: mc_par
!!       Type(force_env_type), Pointer,Optional:: force_env
!!       Type(cp_para_env_type), Pointer,Optional:: para_env
!!     End Subroutine get_mc_env
!!
!!   FUNCTION
!!     provides a method for getting the various structures attached
!!     to an mc_env
!!
!!   ARGUMENTS
!!     - mc_env: the mc_environment you want to get information on
!!     - mc_par: the mc parameters you want to point to the parameters
!!               associated with this mc_env
!!     - force_env: the force environment type you want to point to the
!!               force environment associated with this mc_env
!!     - para_env: the cp para environment type you want to point to the
!!               para environment associated with this mc_env
!!
!!   AUTHOR
!!     MJM
!!
!!**************************************************************************
  SUBROUTINE get_mc_env ( mc_env, mc_par, force_env, para_env)

    TYPE(mc_environment_type), POINTER    :: mc_env
    TYPE(cp_para_env_type), POINTER,OPTIONAL        :: para_env
    TYPE(mc_simulation_parameters_type), &
      OPTIONAL, POINTER                      :: mc_par
    TYPE(force_env_type), OPTIONAL, POINTER  :: force_env

!   ---------------------------------------------------------------------------

    IF ( PRESENT ( mc_par ) ) mc_par => mc_env % mc_par
    IF ( PRESENT ( force_env ) ) force_env => mc_env % force_env
    IF ( PRESENT ( para_env ) ) para_env => mc_env % para_env

  END SUBROUTINE get_mc_env

! *****************************************************************************
!***************************************************************************

!!****f* mc_environment_types/mc_env_retain [1.0] *
!!
!!   NAME
!!     mc_env_retain
!!
!!   SYNOPSIS
!!     Subroutine mc_env_retain(mc_env, error)
!!       Type(mc_env_type), Pointer:: mc_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine mc_env_retain
!!
!!   FUNCTION
!!     retains the given mc env
!!
!!   NOTES
!!     see doc/ReferenceCounting.html
!!
!!   ARGUMENTS
!!     - mc_env: the force environment to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     MJM
!!
!!*** **********************************************************************
SUBROUTINE mc_env_retain(mc_env, error)
    TYPE(mc_environment_type), POINTER       :: mc_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mc_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(mc_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(mc_env%ref_count>0,cp_failure_level,routineP,error)
     mc_env%ref_count=mc_env%ref_count+1
  END IF
END SUBROUTINE mc_env_retain
!***************************************************************************

!!****f* mc_environment_types/mc_env_release [1.0] *
!!
!!   NAME
!!     mc_env_release
!!
!!   SYNOPSIS
!!     Subroutine mc_env_release(mc_env, error)
!!       Type(mc_env_type), Pointer:: mc_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine mc_env_release
!!
!!   FUNCTION
!!     releases the given mc env
!!
!!   NOTES
!!     see doc/ReferenceCounting.html
!!
!!   ARGUMENTS
!!     - mc_env: the mc environment to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     MJM
!!
!!*** **********************************************************************
SUBROUTINE mc_env_release(mc_env, error)
    TYPE(mc_environment_type), POINTER       :: mc_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mc_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (ASSOCIATED(mc_env)) THEN
     CPPreconditionNoFail(mc_env%ref_count>0,cp_failure_level,routineP,error)
     mc_env%ref_count=mc_env%ref_count-1
     IF (mc_env%ref_count==0) THEN
        mc_env%ref_count=1
        CALL cp_para_env_release ( mc_env % para_env, error = error )
        NULLIFY ( mc_env % mc_par )
        CALL force_env_release( mc_env % force_env) 
        mc_env%ref_count=0
        DEALLOCATE(mc_env,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(mc_env)
END SUBROUTINE mc_env_release

! *****************************************************************************

!***************************************************************************

END MODULE mc_environment_types

!******************************************************************************
