!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2003  CP2K developers group                            !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****s* cp2k/mc_environment_types [1.0] *
!!
!!   NAME
!!     mc_environment_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     MJM Oct. 15-2003
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE mc_environment_types
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_error_type, &
                                             cp_assert
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_warning_level,&
                                             cp_to_string
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type,&
                                             cp_blacs_env_type
  USE force_env_types,                 ONLY: force_env_release,&
                                             force_env_retain,&
                                             force_env_type
  USE mc_control,                      ONLY: mc_simulation_parameters_type

  IMPLICIT NONE

  PRIVATE

 
  TYPE mc_environment_type
     INTEGER :: id_nr, ref_count, in_use
     TYPE (cp_para_env_type), POINTER :: para_env 
     TYPE ( mc_simulation_parameters_type ), POINTER :: mc_par
     TYPE ( force_env_type ), POINTER :: force_env  
     TYPE ( cp_blacs_env_type ), POINTER :: blacs_env  
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

  SUBROUTINE get_mc_env ( mc_env, mc_par, force_env, blacs_env)

!   Purpose: Get the integrator environment...works differently than
!            set_mc_env

!   ***************************************************************************

    TYPE(mc_environment_type), INTENT(IN)    :: mc_env
    TYPE(mc_simulation_parameters_type), &
      OPTIONAL, POINTER                      :: mc_par
    TYPE(force_env_type), OPTIONAL, POINTER  :: force_env
    TYPE(cp_blacs_env_type), OPTIONAL, &
      POINTER                                :: blacs_env

!   ---------------------------------------------------------------------------

    IF ( PRESENT ( mc_par ) ) mc_par => mc_env % mc_par
    IF ( PRESENT ( force_env ) ) force_env => mc_env % force_env
    IF ( PRESENT ( blacs_env ) ) blacs_env => mc_env % blacs_env

  END SUBROUTINE get_mc_env

! *****************************************************************************

  SUBROUTINE mc_env_create ( mc_env, para_env, error )

!   Purpose: Initialise the integrator environment.
!   retain the para_env for this environment (should be used for parallel
!   communications)

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
    NULLIFY ( mc_env % blacs_env )
    mc_env % para_env => para_env 
    CALL cp_para_env_retain ( mc_env % para_env )

  END SUBROUTINE mc_env_create

! *****************************************************************************

  SUBROUTINE set_mc_env ( mc_env, mc_par, force_env, blacs_env )

!   Purpose: Set the integrator environment to the correct program.


!   ***************************************************************************
    TYPE(mc_environment_type), INTENT(INOUT)   :: mc_env
    TYPE(mc_simulation_parameters_type), &
      OPTIONAL, POINTER                      :: mc_par
    TYPE(force_env_type), OPTIONAL, POINTER  :: force_env
    TYPE(cp_blacs_env_type), OPTIONAL, &
      POINTER                                :: blacs_env

    IF ( PRESENT ( mc_par ) ) mc_env % mc_par => mc_par
    IF ( PRESENT ( blacs_env ) ) mc_env % blacs_env => blacs_env
    IF ( PRESENT ( force_env )) THEN
       CALL force_env_retain( force_env ) ! accept null pointers?
       IF (ASSOCIATED(mc_env % force_env ))  &
           CALL force_env_release( mc_env % force_env )
       mc_env % force_env => force_env
    END IF
    
!   ---------------------------------------------------------------------------

  END SUBROUTINE set_mc_env

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
!!   MODIFICATION HISTORY
!!     02.2004 created [MJM]
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
!!   MODIFICATION HISTORY
!!     02.2004 created [MJM]
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
