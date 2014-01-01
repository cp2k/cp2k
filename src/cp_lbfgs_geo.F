!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Main driver for L-BFGS optimizer
!> \par History
!>      none
! *****************************************************************************
MODULE cp_lbfgs_geo
  USE cell_types,                      ONLY: cell_type
  USE cp_external_control,             ONLY: external_control
  USE cp_lbfgs_optimizer_gopt,         ONLY: cp_lbfgs_opt_gopt_type,&
                                             cp_opt_gopt_create,&
                                             cp_opt_gopt_next,&
                                             cp_opt_gopt_release,&
                                             cp_opt_gopt_stop
  USE cp_output_handling,              ONLY: cp_iterate,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE gopt_f_methods,                  ONLY: gopt_f_ii,&
                                             gopt_f_io_finalize,&
                                             print_geo_opt_header,&
                                             print_geo_opt_nc
  USE gopt_f_types,                    ONLY: gopt_f_type
  USE gopt_param_types,                ONLY: gopt_param_type
  USE input_constants,                 ONLY: default_ts_method_id
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_set
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_lbfgs_geo'

  PUBLIC :: geoopt_lbfgs

CONTAINS

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  RECURSIVE SUBROUTINE geoopt_lbfgs(force_env, gopt_param, globenv, geo_section, gopt_env,&
                          x0, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(gopt_param_type), POINTER           :: gopt_param
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: geo_section
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x0
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'geoopt_lbfgs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iter_nr, its, &
                                                output_unit
    LOGICAL                                  :: converged, failure, &
                                                should_stop
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_lbfgs_opt_gopt_type), POINTER    :: optimizer
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: root_section

    CALL timeset(routineN,handle)

    failure = .FALSE.
    NULLIFY (optimizer, para_env)
    logger => cp_error_get_logger(error)
    root_section => force_env%root_section
    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(gopt_param),cp_failure_level,routineP,error,failure)
    CPPrecondition(gopt_param%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
       CALL force_env_get(force_env, para_env=para_env, cell=cell, error=error)

       ! Geometry optimization starts now
       output_unit = cp_print_key_unit_nr(logger,geo_section,"PRINT%PROGRAM_RUN_INFO",&
            extension=".geoLog",error=error)
       CALL print_geo_opt_header(gopt_env, output_unit, "L-BFGS")
       CALL cp_print_key_finished_output(output_unit,logger,geo_section,&
            "PRINT%PROGRAM_RUN_INFO", error=error)

       ! Stop if not implemented
       SELECT CASE (gopt_env%type_id)
       CASE (default_ts_method_id)
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="BFGS method not yet working with DIMER", &
               error=error, error_level=cp_failure_level)
       END SELECT

       CALL cp_opt_gopt_create(optimizer, para_env=para_env, obj_funct=gopt_env,&
            x0=x0, wanted_relative_f_delta=gopt_param%wanted_rel_f_error,&
            wanted_projected_gradient=gopt_param%wanted_proj_gradient, m=gopt_param%max_h_rank,&
            max_f_per_iter=gopt_param%max_f_per_iter, error=error)
       CALL cp_iterate(logger%iter_info,increment=0,iter_nr_out=iter_nr,error=error)
       converged=.FALSE.

       DO its=iter_nr+1,gopt_param%max_iter
          CALL cp_iterate(logger%iter_info,last=(its==gopt_param%max_iter),error=error)
          CALL section_vals_val_set(geo_section,"STEP_START_VAL",i_val=its,error=error)
          CALL gopt_f_ii(its, output_unit)

          ! Real optimization step..
          IF (.NOT.cp_opt_gopt_next(optimizer,geo_section=geo_section, input=root_section,&
               force_env=force_env,gopt_param=gopt_param, converged=converged,error=error)) EXIT

          ! Check for an external exit command
          CALL external_control(should_stop,"GEO",globenv=globenv,error=error)
          IF (should_stop) THEN
             CALL cp_opt_gopt_stop(optimizer, error=error)
             EXIT
          END IF
          IF (its == gopt_param%max_iter) EXIT
       END DO

       IF ((its == gopt_param%max_iter).AND.(.NOT.converged))THEN
          CALL print_geo_opt_nc(gopt_env, output_unit)
       END IF

       ! Write final output information, if converged
       CALL cp_iterate(logger%iter_info,last=.TRUE.,increment=0,error=error)
       CALL gopt_f_io_finalize(gopt_env, force_env, optimizer%x, converged, its, root_section,&
            optimizer%para_env, optimizer%master, error)

       CALL cp_opt_gopt_release(optimizer, error=error)
       CALL cp_print_key_finished_output(output_unit,logger,geo_section,&
            "PRINT%PROGRAM_RUN_INFO", error=error)
    END IF

    CALL timestop(handle)

  END SUBROUTINE geoopt_lbfgs

END MODULE cp_lbfgs_geo
