!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!


!!****h* cp2k/gopt_f_types [1.0] *
!!
!!   NAME
!!     gopt_f_types
!!
!!   FUNCTION
!!     contains a functional that calculates the energy and its derivatives
!!     for the geometry optimizer
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     -
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!*****************************************************************************
MODULE gopt_f_types
  USE dimer_types,                     ONLY: dimer_env_create,&
                                             dimer_env_release,&
                                             dimer_env_type
  USE f77_blas
  USE force_env_types,                 ONLY: force_env_get_natom,&
                                             force_env_release,&
                                             force_env_retain,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE gopt_param_types,                ONLY: gopt_param_type
  USE input_constants,                 ONLY: default_dimer_method_id,&
                                             default_ts_method_id
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gopt_f_types'
  INTEGER, PRIVATE, SAVE :: last_gopt_f_id=0, last_gopt_param_id=0

  PUBLIC :: gopt_f_type, gopt_f_create, gopt_f_retain, gopt_f_release

!!****s* gopt_f_types/gopt_f_type [1.0] *
!!
!!   NAME
!!     gopt_f_type
!!
!!   FUNCTION
!!     calculates the potential energy of a system, and its derivatives
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     -
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
  !***************************************************************************
  TYPE gopt_f_type
     INTEGER :: ref_count, id_nr
     INTEGER :: type_id, ts_method_id
     TYPE(force_env_type), POINTER            :: force_env
     TYPE(global_environment_type), POINTER   :: globenv
     TYPE(section_vals_type), POINTER         :: force_env_section,&
                                                 motion_section
     TYPE(dimer_env_type), POINTER            :: dimer_env
  END TYPE gopt_f_type
!!***
  !****************************************************************************

CONTAINS
  
!!****f* gopt_f_types/gopt_f_create [1.0] *
!!
!!   NAME
!!     gopt_f_create
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     gopt_env: the geometry optimization environment to be created
!!     force_env:
!!     globenv:
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE gopt_f_create(gopt_env, gopt_param, force_env, globenv, error)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    TYPE(gopt_param_type), POINTER           :: gopt_param
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gopt_f_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: natom, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    ALLOCATE(gopt_env, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       NULLIFY (gopt_env%dimer_env)
       gopt_env%ref_count=1
       last_gopt_f_id=last_gopt_f_id+1
       gopt_env%id_nr=last_gopt_f_id
       CALL force_env_retain(force_env, error=error)
       gopt_env%force_env => force_env
       gopt_env%globenv => globenv
       gopt_env%force_env_section => force_env%force_env_section
       gopt_env%motion_section =>  section_vals_get_subs_vals(force_env%root_section,"MOTION",error=error)
       gopt_env%type_id = gopt_param%type_id
       IF (gopt_param%type_id== default_ts_method_id) THEN 
          gopt_env%ts_method_id = gopt_param%ts_method_id
          SELECT CASE(gopt_param%ts_method_id)
          CASE(default_dimer_method_id)
             natom = force_env_get_natom(force_env,error)
             CALL dimer_env_create(gopt_env%dimer_env, gopt_param, natom, globenv, error)
          END SELECT
       END IF
    END IF
  END SUBROUTINE gopt_f_create
!***************************************************************************

!!****f* gopt_f_types/gopt_f_retain [1.0] *
!!
!!   NAME
!!     gopt_f_retain
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     gopt_env: the geometry optimization environment to retain
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE gopt_f_retain(gopt_env, error)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gopt_f_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(gopt_env),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CPPreconditionNoFail(gopt_env%ref_count>0,cp_failure_level,routineP,error)
       gopt_env%ref_count=gopt_env%ref_count+1
    END IF
  END SUBROUTINE gopt_f_retain
!***************************************************************************

!!****f* gopt_f_types/gopt_f_release [1.0] *
!!
!!   NAME
!!     gopt_f_release
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     gopt_env: the geometry optimization environment to release
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE gopt_f_release(gopt_env, error)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gopt_f_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(gopt_env)) THEN
       CPPreconditionNoFail(gopt_env%ref_count>0,cp_failure_level,routineP,error)
       gopt_env%ref_count=gopt_env%ref_count-1
       IF (gopt_env%ref_count==0) THEN
          CALL force_env_release(gopt_env%force_env,error=error)
          NULLIFY(gopt_env%force_env,&
                  gopt_env%globenv,&
                  gopt_env%force_env_section,&
                  gopt_env%motion_section)
          CALL dimer_env_release(gopt_env%dimer_env, error)
          DEALLOCATE(gopt_env, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE gopt_f_release
!***************************************************************************

END MODULE gopt_f_types
