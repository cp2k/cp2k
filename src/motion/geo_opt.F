!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief performs geometry optimization
!> \par History
!>      none
! *****************************************************************************
MODULE geo_opt

  USE bfgs_optimizer,                  ONLY: geoopt_bfgs
  USE cg_optimizer,                    ONLY: geoopt_cg
  USE cp_lbfgs_geo,                    ONLY: geoopt_lbfgs
  USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                             cp_iterate,&
                                             cp_rm_iter_level
  USE force_env_types,                 ONLY: force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE gopt_f_methods,                  ONLY: gopt_f_create_x0
  USE gopt_f_types,                    ONLY: gopt_f_create,&
                                             gopt_f_release,&
                                             gopt_f_type
  USE gopt_param_types,                ONLY: gopt_param_read,&
                                             gopt_param_release,&
                                             gopt_param_type
  USE input_constants,                 ONLY: default_bfgs_method_id,&
                                             default_cg_method_id,&
                                             default_lbfgs_method_id
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'geo_opt'

  PUBLIC :: cp_geo_opt, cp_rot_opt

CONTAINS

! *****************************************************************************
!> \brief Main driver to perform geometry optimization
!> \param force_env ...
!> \param globenv ...
!> \param eval_opt_geo ...
!> \param rm_restart_info ...
! *****************************************************************************
  RECURSIVE SUBROUTINE cp_geo_opt(force_env, globenv, eval_opt_geo, rm_restart_info)

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    LOGICAL, INTENT(IN), OPTIONAL            :: eval_opt_geo, rm_restart_info

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_geo_opt', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat, step_start_val
    LOGICAL                                  :: failure, my_rm_restart_info
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x0
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(gopt_f_type), POINTER               :: gopt_env
    TYPE(gopt_param_type), POINTER           :: gopt_param
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                geo_section, root_section

    failure = .FALSE.
    CALL timeset(routineN,handle)
    logger => cp_get_default_logger()
    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,failure)
    CPPrecondition(ASSOCIATED(globenv),cp_failure_level,routineP,failure)
    NULLIFY (gopt_param,force_env_section,gopt_env,x0)
    root_section      => force_env%root_section
    force_env_section => force_env%force_env_section
    geo_section       => section_vals_get_subs_vals(root_section,"MOTION%GEO_OPT")

    CALL gopt_param_read(gopt_param, geo_section)
    CALL gopt_f_create(gopt_env, gopt_param, force_env=force_env, globenv=globenv,&
         geo_opt_section=geo_section, eval_opt_geo=eval_opt_geo)
    CALL gopt_f_create_x0(gopt_env, x0)

    CALL section_vals_val_get(geo_section,"STEP_START_VAL",i_val=step_start_val)
    CALL cp_add_iter_level(logger%iter_info,"GEO_OPT")
    CALL cp_iterate(logger%iter_info,iter_nr=step_start_val)
    CALL cp_geo_opt_low(force_env, globenv, gopt_param, gopt_env,&
         force_env_section, geo_section, x0)
    CALL cp_rm_iter_level(logger%iter_info,"GEO_OPT")

    ! Reset counter for next iteration, unless rm_restart_info==.FALSE.
    my_rm_restart_info = .TRUE.
    IF(PRESENT(rm_restart_info)) my_rm_restart_info = rm_restart_info
    IF(my_rm_restart_info) &
       CALL section_vals_val_set(geo_section,"STEP_START_VAL",i_val=0)

    DEALLOCATE (x0,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,failure)
    CALL gopt_f_release(gopt_env)
    CALL gopt_param_release(gopt_param)
    CALL timestop(handle)

  END SUBROUTINE cp_geo_opt

! *****************************************************************************
!> \brief Main driver to perform rotation optimization for Dimer
!> \param gopt_env ...
!> \param x0 ...
!> \param gopt_param ...
!> \param geo_section ...
! *****************************************************************************
   SUBROUTINE cp_rot_opt(gopt_env, x0, gopt_param, geo_section)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x0
    TYPE(gopt_param_type), POINTER           :: gopt_param
    TYPE(section_vals_type), POINTER         :: geo_section

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_rot_opt', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, step_start_val
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: force_env_section

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY (force_env_section)
    logger => cp_get_default_logger()
    CPPrecondition(ASSOCIATED(gopt_env),cp_failure_level,routineP,failure)
    CPPrecondition(ASSOCIATED(gopt_env%force_env),cp_failure_level,routineP,failure)
    CPPrecondition(ASSOCIATED(gopt_env%globenv),cp_failure_level,routineP,failure)
    force_env_section => gopt_env%force_env%force_env_section

    CALL section_vals_val_get(geo_section,"STEP_START_VAL",i_val=step_start_val)
    CALL cp_add_iter_level(logger%iter_info,"ROT_OPT")
    CALL cp_iterate(logger%iter_info,iter_nr=step_start_val)
    CALL cp_geo_opt_low(gopt_env%force_env, gopt_env%globenv, gopt_param, gopt_env,&
         force_env_section, geo_section, x0)
    CALL cp_rm_iter_level(logger%iter_info,"ROT_OPT")

    ! Reset counter for next iteration
    CALL section_vals_val_set(geo_section,"STEP_START_VAL",i_val=0)
    CALL timestop(handle)

  END SUBROUTINE cp_rot_opt

! *****************************************************************************
!> \brief call to low level geometry optimizers
!> \param force_env ...
!> \param globenv ...
!> \param gopt_param ...
!> \param gopt_env ...
!> \param force_env_section ...
!> \param geo_section ...
!> \param x0 ...
! *****************************************************************************
  RECURSIVE SUBROUTINE cp_geo_opt_low(force_env, globenv, gopt_param, gopt_env, force_env_section,&
       geo_section, x0)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(gopt_param_type), POINTER           :: gopt_param
    TYPE(gopt_f_type), POINTER               :: gopt_env
    TYPE(section_vals_type), POINTER         :: force_env_section, geo_section
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x0

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_geo_opt_low', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,failure)
    CPPrecondition(ASSOCIATED(globenv),cp_failure_level,routineP,failure)
    CPPrecondition(ASSOCIATED(gopt_param),cp_failure_level,routineP,failure)
    CPPrecondition(ASSOCIATED(gopt_env),cp_failure_level,routineP,failure)
    CPPrecondition(ASSOCIATED(x0),cp_failure_level,routineP,failure)
    CPPrecondition(ASSOCIATED(force_env_section),cp_failure_level,routineP,failure)
    CPPrecondition(ASSOCIATED(geo_section),cp_failure_level,routineP,failure)

    SELECT CASE (gopt_param%method_id)
    CASE (default_bfgs_method_id)
       CALL geoopt_bfgs(force_env,gopt_param,globenv,&
            geo_section, gopt_env, x0)
    CASE (default_lbfgs_method_id)
       CALL geoopt_lbfgs(force_env,gopt_param,globenv,&
            geo_section, gopt_env, x0)
    CASE (default_cg_method_id)
       CALL geoopt_cg(force_env,gopt_param,globenv,&
            geo_section, gopt_env, x0)
    CASE DEFAULT
       CPPostcondition(.FALSE.,cp_failure_level,routineP,failure)
    END SELECT

  END SUBROUTINE cp_geo_opt_low

END MODULE geo_opt
