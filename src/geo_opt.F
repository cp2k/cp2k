!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/geo_opt [1.0] *
!!
!!   NAME
!!     geo_opt
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
MODULE geo_opt

  USE bfgs_optimizer,                  ONLY: geoopt_bfgs
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE force_env_types,                 ONLY: force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE gopt_f_types,                    ONLY: default_bfgs_method_id,&
                                             gopt_param_read,&
                                             gopt_param_release,&
                                             gopt_param_type
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop

 IMPLICIT NONE
 PRIVATE

 CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='geo_opt'

 PUBLIC :: qs_geo_opt
 CONTAINS

   SUBROUTINE qs_geo_opt(force_env, globenv, error)

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_geo_opt', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(gopt_param_type), POINTER           :: gopt_param

    CALL write_checkpoint_information("entering "//routineN,globenv)

     CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
     CPPrecondition(ASSOCIATED(globenv),cp_failure_level,routineP,error,failure)
     
     NULLIFY (gopt_param)

     CALL gopt_param_read(gopt_param, globenv, error=error)
     
     SELECT CASE (gopt_param%method_id)
     CASE (default_bfgs_method_id)
       CALL geoopt_bfgs(force_env,globenv,gopt_param)
     CASE DEFAULT
       CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
     END SELECT

     CALL gopt_param_release(gopt_param, error=error)

    CALL write_checkpoint_information("leaving "//routineN,globenv)

  END SUBROUTINE qs_geo_opt

END MODULE geo_opt
