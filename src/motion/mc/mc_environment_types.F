!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief contains the subroutines for dealing with the mc_env
!> \author MJM Oct. 15-2003
! *****************************************************************************
MODULE mc_environment_types
  
  USE force_env_types,                 ONLY: force_env_type
  USE mc_types,                        ONLY: mc_simpar_type
#include "../../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *****************************************************************************
  TYPE mc_environment_type
     INTEGER :: id_nr, ref_count, in_use
     TYPE ( mc_simpar_type ), POINTER :: mc_par
     TYPE ( force_env_type ), POINTER :: force_env
  END TYPE mc_environment_type

! *****************************************************************************
  TYPE mc_environment_p_type
     TYPE(mc_environment_type),POINTER :: mc_env
  END TYPE mc_environment_p_type

! *** Public subroutines and data types ***
  PUBLIC :: mc_environment_type, mc_environment_p_type,&
            set_mc_env, mc_env_create,&
            get_mc_env, mc_env_release

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mc_environment_types'
  INTEGER, SAVE, PRIVATE :: last_mc_env_id=0

CONTAINS

! *****************************************************************************
!> \brief creates and initializes an mc_env
!> \param mc_env the mc_environment you want to create
!>
!>    Suitable for parallel use.
!> \param error ...
!> \author MJM
! *****************************************************************************
  SUBROUTINE mc_env_create ( mc_env, error )

    TYPE(mc_environment_type), POINTER       :: mc_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mc_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    ALLOCATE ( mc_env, stat=istat )
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)

    last_mc_env_id=last_mc_env_id+1
    mc_env%id_nr=last_mc_env_id
    mc_env%ref_count=1
    mc_env%in_use=0

    NULLIFY ( mc_env % mc_par )
    NULLIFY ( mc_env % force_env )

  END SUBROUTINE mc_env_create

! *****************************************************************************
!> \brief provides a method for attaching various structures to an mc_env
!> \param mc_env the mc_environment you want to change
!> \param mc_par the mc parameters you want to associate with this mc_env
!> \param force_env the force environment type you want to associate
!>                   with this mc_env
!>
!>    Suitable for parallel.
!> \author MJM
! *****************************************************************************
  SUBROUTINE set_mc_env ( mc_env,mc_par,force_env)

    TYPE(mc_environment_type), POINTER       :: mc_env
    TYPE(mc_simpar_type), OPTIONAL, POINTER  :: mc_par
    TYPE(force_env_type), OPTIONAL, POINTER  :: force_env

    IF ( PRESENT ( mc_par ) ) mc_env % mc_par => mc_par
    IF ( PRESENT ( force_env )) THEN
       mc_env % force_env => force_env
    END IF


  END SUBROUTINE set_mc_env

! *****************************************************************************
!> \brief provides a method for getting the various structures attached
!>      to an mc_env
!> \param mc_env the mc_environment you want to get information on
!> \param mc_par the mc parameters you want to point to the parameters
!>                associated with this mc_env
!> \param force_env the force environment type you want to point to the
!>                force environment associated with this mc_env
!>
!>    Suitable for parallel.
!> \author MJM
! *****************************************************************************
  SUBROUTINE get_mc_env ( mc_env, mc_par, force_env)

    TYPE(mc_environment_type), POINTER       :: mc_env
    TYPE(mc_simpar_type), OPTIONAL, POINTER  :: mc_par
    TYPE(force_env_type), OPTIONAL, POINTER  :: force_env

    IF ( PRESENT ( mc_par ) ) mc_par => mc_env % mc_par
    IF ( PRESENT ( force_env ) ) force_env => mc_env % force_env

  END SUBROUTINE get_mc_env

! *****************************************************************************
!> \brief retains the given mc env
!> \param mc_env the force environment to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!>
!>    Suitable for parallel.
!> \author MJM
!> \note
!>      see doc/ReferenceCounting.html
! *****************************************************************************
SUBROUTINE mc_env_retain(mc_env, error)
    TYPE(mc_environment_type), POINTER       :: mc_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mc_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(mc_env),cp_failure_level,routineP,error,failure)
  CPPreconditionNoFail(mc_env%ref_count>0,cp_failure_level,routineP,error)
  mc_env%ref_count=mc_env%ref_count+1
END SUBROUTINE mc_env_retain

! *****************************************************************************
!> \brief releases the given mc env
!> \param mc_env the mc environment to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!>
!>    Suitable for parallel.
!> \author MJM
!> \note
!>      see doc/ReferenceCounting.html
! *****************************************************************************
SUBROUTINE mc_env_release(mc_env, error)
    TYPE(mc_environment_type), POINTER       :: mc_env
    TYPE(cp_error_type), INTENT(inout)       :: error

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
        NULLIFY ( mc_env % mc_par )
        NULLIFY ( mc_env % force_env)
        mc_env%ref_count=0
        DEALLOCATE(mc_env,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(mc_env)
END SUBROUTINE mc_env_release

END MODULE mc_environment_types

