!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qs_geoopt [1.0] *
!!
!!   NAME
!!     qs_geoopt
!!
!!   FUNCTION
!!     performs geometry optimization
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
MODULE qs_geoopt
!
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
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
                                             cp_unimplemented_error,&
                                             cp_unimplemented_error_nr
  USE cp_fragment_types,               ONLY: cp_fragment_p_type,&
                                             cp_fragment_type
  USE cp_lbfgs_optimizer_gopt,         ONLY: cp_lbfgs_opt_gopt_type,&
                                             cp_opt_gopt_dealloc_ref,&
                                             cp_opt_gopt_get,&
                                             cp_opt_gopt_init,&
                                             cp_opt_gopt_next
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE force_control,                   ONLY: force_env_calculate_force,&
                                             force_env_create,&
                                             force_env_get,&
                                             force_env_release,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE gopt_f_types,                    ONLY: gopt_f_create,&
                                             gopt_f_create_x0,&
                                             gopt_f_release,&
                                             gopt_f_retain,&
                                             gopt_f_type,&
                                             gopt_param_create,&
                                             gopt_param_read,&
                                             gopt_param_release,&
                                             gopt_param_type
  USE kinds,                           ONLY: wp=>dp
  USE particle_types,                  ONLY: particle_type
  USE qs_energy,                       ONLY: qs_energies
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force,                        ONLY: qs_forces
  USE qs_force_types,                  ONLY: init_qs_force,&
                                             qs_force_type
 IMPLICIT NONE
 PRIVATE

 CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_geoopt'

 PUBLIC :: qs_geo_opt
 CONTAINS

   SUBROUTINE qs_geo_opt(qs_env,globenv, error)

     TYPE(global_environment_type), INTENT(INOUT), TARGET :: globenv
     TYPE(qs_environment_type), POINTER              :: qs_env
     TYPE(cp_error_type), INTENT(inout), OPTIONAL    :: error

!     *** Local variables ***

     TYPE(force_env_type), POINTER                   :: force_env
     TYPE(gopt_param_type), POINTER :: gopt_param
     LOGICAL :: failure, finished
     CHARACTER(len=*), PARAMETER :: routineN='qs_geo_opt',&
          routineP=moduleN//':'//routineN

     NULLIFY (force_env, gopt_param)

     CALL force_env_create(force_env, qs_env=qs_env)
     CALL gopt_param_read(gopt_param, globenv, error=error)
     
     CALL gopt_do_opt(force_env, gopt_param, globenv, error=error)
     
     CALL force_env_release(force_env)
     CALL gopt_param_release(gopt_param, error=error)

  END SUBROUTINE qs_geo_opt


!!****f* m/gopt_do_opt [1.0] *
!!
!!   NAME
!!     gopt_do_opt
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
!!   MODIFICATION HISTORY
!!     08.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE gopt_do_opt(force_env, gopt_param, globenv, error)
  TYPE(force_env_type), POINTER :: force_env
  TYPE(gopt_param_type), POINTER :: gopt_param
  TYPE(global_environment_type), TARGET, INTENT(inout) :: globenv
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure, finished
  CHARACTER(len=*), PARAMETER :: routineN='gopt_do_opt',&
        routineP=moduleN//':'//routineN
  TYPE(gopt_f_type), POINTER :: gopt_env
  TYPE(cp_lbfgs_opt_gopt_type) :: optimizer
  INTEGER :: i, stat
  REAL(kind=wp) :: f
  REAL(kind=wp), DIMENSION(:), POINTER :: x0

  failure=.FALSE.
  NULLIFY(gopt_env, x0)
  
  CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(gopt_param),cp_failure_level,routineP,error,failure)
  CPPrecondition(gopt_param%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN

     CALL gopt_f_create(gopt_env, force_env=force_env, globenv=globenv)

     CALL gopt_f_create_x0(gopt_env, x0, error=error)

     CALL cp_opt_gopt_init(optimizer, global_env=globenv, &
          obj_funct=gopt_env, x0=x0, wanted_relative_f_delta=gopt_param%wanted_rel_f_error,&
          wanted_projected_gradient=gopt_param%wanted_proj_gradient, m=gopt_param%max_h_rank,&
          max_f_per_iter=gopt_param%max_f_per_iter, error=error)

     DEALLOCATE(x0,stat=stat)
     CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)

     DO i=1,gopt_param%max_iter
        IF (.NOT.cp_opt_gopt_next(optimizer)) EXIT
     END DO

     CALL cp_opt_gopt_get(optimizer, f=f, at_end=finished)
     PRINT *,"pippo f=",f,"status=",optimizer%status
     CPPostcondition(finished,cp_warning_level,routineP,error,failure)
     CALL cp_opt_gopt_dealloc_ref(optimizer)
     CALL gopt_f_release(gopt_env)     
     
  END IF
END SUBROUTINE gopt_do_opt
!***************************************************************************

END MODULE qs_geoopt
