!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief performs CELL optimization
!> \par History
!>      03.2008 - Teodoro Laino [tlaino] - University of Zurich - Cell Optimization
! *****************************************************************************
MODULE cell_opt
  USE bfgs_optimizer,                  ONLY: geoopt_bfgs
  USE cg_optimizer,                    ONLY: geoopt_cg
  USE cp_lbfgs_geo,                    ONLY: geoopt_lbfgs
  USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                             cp_iterate,&
                                             cp_rm_iter_level
  USE f77_blas
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
                                             default_cell_method_id,&
                                             default_cg_method_id,&
                                             default_lbfgs_method_id
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cell_opt'

  PUBLIC :: cp_cell_opt

 CONTAINS

! *****************************************************************************
!> \brief Main driver to perform geometry optimization
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino [tlaino] - University of Zurich - 03.2008
! *****************************************************************************
   SUBROUTINE cp_cell_opt(force_env, globenv, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_cell_opt', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat, step_start_val
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x0
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(gopt_f_type), POINTER               :: gopt_env
    TYPE(gopt_param_type), POINTER           :: gopt_param
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                geo_section, root_section

    failure = .FALSE.
    CALL timeset(routineN,handle)
    logger => cp_error_get_logger(error)
    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(globenv),cp_failure_level,routineP,error,failure)
    NULLIFY (gopt_param,force_env_section,gopt_env,x0)
    root_section      => force_env%root_section
    force_env_section => force_env%force_env_section
    geo_section       => section_vals_get_subs_vals(root_section,"MOTION%CELL_OPT",error=error)

    CALL gopt_param_read(gopt_param, geo_section, type_id=default_cell_method_id, error=error)
    CALL gopt_f_create(gopt_env, gopt_param, force_env=force_env, globenv=globenv,&
         geo_opt_section=geo_section, error=error)
    CALL gopt_f_create_x0(gopt_env, x0, error=error)

    CALL section_vals_val_get(geo_section,"STEP_START_VAL",i_val=step_start_val,error=error)
    CALL cp_add_iter_level(logger%iter_info,"CELL_OPT",error=error)
    CALL cp_iterate(logger%iter_info,iter_nr=step_start_val,error=error)
    CALL cp_cell_opt_low(force_env, globenv, gopt_param, gopt_env,&
         force_env_section, geo_section, x0, error)
    CALL cp_rm_iter_level(logger%iter_info,"CELL_OPT",error=error)

    ! Reset counter for next iteration
    CALL section_vals_val_set(geo_section,"STEP_START_VAL",i_val=0,error=error)
    DEALLOCATE(x0, stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL gopt_f_release(gopt_env, error=error)
    CALL gopt_param_release(gopt_param, error=error)
    CALL timestop(handle)

  END SUBROUTINE cp_cell_opt

! *****************************************************************************
!> \brief call to low level geometry optimizers
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino [tlaino] - University of Zurich - 03.2008
! *****************************************************************************
  SUBROUTINE cp_cell_opt_low(force_env, globenv, gopt_param, gopt_env, force_env_section,&
       geo_section, x0, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(gopt_param_type), POINTER           :: gopt_param
    TYPE(gopt_f_type), POINTER               :: gopt_env
    TYPE(section_vals_type), POINTER         :: force_env_section, geo_section
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x0
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_cell_opt_low', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(globenv),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(gopt_param),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(gopt_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(x0),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(force_env_section),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(geo_section),cp_failure_level,routineP,error,failure)

    SELECT CASE (gopt_param%method_id)
    CASE (default_bfgs_method_id)
       CALL geoopt_bfgs(force_env,gopt_param,globenv,&
            geo_section, gopt_env, x0, error=error)
    CASE (default_lbfgs_method_id)
       CALL geoopt_lbfgs(force_env,gopt_param,globenv,&
            geo_section, gopt_env, x0, error=error)
    CASE (default_cg_method_id)
       CALL geoopt_cg(force_env,gopt_param,globenv,&
            geo_section, gopt_env, x0, error=error)
    CASE DEFAULT
       CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

  END SUBROUTINE cp_cell_opt_low

END MODULE cell_opt
