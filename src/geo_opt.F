!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


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
  USE cg_optimizer,                    ONLY: geoopt_cg
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE cp_lbfgs_optimizer_gopt,         ONLY: cp_lbfgs_opt_gopt_type,&
                                             cp_opt_gopt_create,&
                                             cp_opt_gopt_next,&
                                             cp_opt_gopt_release,&
                                             cp_opt_gopt_stop
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE gopt_f_types,                    ONLY: gopt_f_create,&
                                             gopt_f_create_x0,&
                                             gopt_f_eval_at,&
                                             gopt_f_release,&
                                             gopt_f_type,&
                                             gopt_param_read,&
                                             gopt_param_release,&
                                             gopt_param_type
  USE kinds,                           ONLY: dp
  USE termination,                     ONLY: external_control
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: write_particle_distances
  USE cell_types,                 ONLY: cell_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE input_constants,                 ONLY: default_lbfgs_method_id,&
                                             default_bfgs_method_id,&
                                             default_cg_method_id
  USE input_section_types, ONLY: section_vals_type, section_vals_get_subs_vals,&
       section_vals_get, section_vals_val_get, section_vals_release
  USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                             cp_rm_iter_level
#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

 CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='geo_opt'

 PUBLIC :: qs_geo_opt

 CONTAINS

!!****f* geo_opt/qs_geo_opt [1.0] *
!!
!!   NAME
!!     qs_geo_opt
!!
!!   FUNCTION
!!     Main driver to perform geometry optimization
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** **********************************************************************
   SUBROUTINE qs_geo_opt(force_env, globenv, error)
     IMPLICIT NONE
     ! Arguments
     TYPE(force_env_type), POINTER            :: force_env
     TYPE(global_environment_type), POINTER   :: globenv
     TYPE(cp_error_type), INTENT(inout), &
          OPTIONAL                               :: error
     ! Local Variables
     CHARACTER(len=*), PARAMETER :: routineN = 'qs_geo_opt', &
          routineP = moduleN//':'//routineN
     LOGICAL                                  :: failure
     TYPE(gopt_param_type), POINTER           :: gopt_param
     TYPE(section_vals_type), POINTER :: force_env_section

     CALL write_checkpoint_information("Entering "//routineN,globenv,error)
     CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
     CPPrecondition(ASSOCIATED(globenv),cp_failure_level,routineP,error,failure)
     NULLIFY (gopt_param,force_env_section)
     CALL gopt_param_read(gopt_param, globenv, error=error)

     force_env_section => section_vals_get_subs_vals(globenv%input_file,&
          "FORCE_EVAL",error=error)     

     SELECT CASE (gopt_param%method_id)
     CASE (default_bfgs_method_id)
        CALL geoopt_bfgs(force_env,globenv,gopt_param,force_env_section,error=error)
     CASE (default_lbfgs_method_id)
        CALL gopt_do_opt(force_env,gopt_param,globenv,force_env_section,error=error)
     CASE (default_cg_method_id)
        CALL geoopt_cg(force_env,gopt_param,globenv,force_env_section,error=error)
     CASE DEFAULT
        CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
     END SELECT

     CALL gopt_param_release(gopt_param, error=error)
     CALL write_checkpoint_information("Leaving "//routineN,globenv,error)
   END SUBROUTINE qs_geo_opt

!!****f* geo_opt/gopt_do_opt [1.0] *
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
   SUBROUTINE gopt_do_opt(force_env, gopt_param, globenv, force_env_section,error)
     IMPLICIT NONE
     ! Arguments
     TYPE(force_env_type), POINTER            :: force_env
     TYPE(gopt_param_type), POINTER           :: gopt_param
     TYPE(section_vals_type), POINTER :: force_env_section
     TYPE(global_environment_type), &
          pointer                  :: globenv
     TYPE(cp_error_type), INTENT(inout), &
          OPTIONAL                               :: error
     ! Local Variables
     CHARACTER(len=*), PARAMETER :: routineN = 'gopt_do_opt', &
          routineP = moduleN//':'//routineN
     INTEGER                                  :: handle, i, isubsys, stat
     LOGICAL                                  :: failure, should_stop
     REAL(KIND=dp)                            :: f
     REAL(KIND=dp), DIMENSION(:), POINTER     :: x0
     TYPE(cell_type), POINTER                 :: cell
     TYPE(cp_lbfgs_opt_gopt_type), POINTER    :: optimizer
     TYPE(cp_logger_type), POINTER            :: logger
     TYPE(cp_para_env_type), POINTER          :: para_env
     TYPE(cp_subsystem_p_type), &
          DIMENSION(:), POINTER                  :: subsys
     TYPE(gopt_f_type), POINTER               :: gopt_env
     TYPE(particle_list_type), POINTER        :: particles
     INTEGER :: ihandle
     
     CALL timeset(routineN,"I","",handle)
     failure=.FALSE.
     NULLIFY(gopt_env, optimizer, para_env, cell, particles, subsys, x0)
     logger => cp_error_get_logger(error)
     CALL cp_add_iter_level(logger%iter_info,"GEO_OPT",ihandle,error)
     CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
     CPPrecondition(ASSOCIATED(gopt_param),cp_failure_level,routineP,error,failure)
     CPPrecondition(gopt_param%ref_count>0,cp_failure_level,routineP,error,failure)
     IF (.NOT. failure) THEN

        CALL force_env_get(force_env, para_env=para_env, error=error)

        CALL gopt_f_create(gopt_env, force_env=force_env, globenv=globenv)

        CALL gopt_f_create_x0(gopt_env, x0, error=error)

        CALL cp_opt_gopt_create(optimizer, para_env=para_env, &
             obj_funct=gopt_env, x0=x0, wanted_relative_f_delta=gopt_param%wanted_rel_f_error,&
             wanted_projected_gradient=gopt_param%wanted_proj_gradient, m=gopt_param%max_h_rank,&
             max_f_per_iter=gopt_param%max_f_per_iter, globenv=globenv, error=error)

        DEALLOCATE(x0,stat=stat)
        CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)

        DO i=1,gopt_param%max_iter
           IF (.NOT.cp_opt_gopt_next(optimizer,force_env_section=force_env_section,&
                input=globenv%input_file,force_env=force_env,gopt_param=gopt_param,&
                error=error)) EXIT
           ! check for an external exit command 
           CALL external_control(should_stop,"GEO",globenv,error)
           IF (should_stop) THEN
              CALL cp_opt_gopt_stop(optimizer, error=error)
              EXIT
           END IF
        END DO

        ! write out more data for the minimum
        IF (logger%para_env%mepos==logger%para_env%source) THEN
           WRITE (unit=cp_logger_get_default_unit_nr(logger),&
                fmt="(' Reevaluating wf at minimum')")
        END IF
        CALL gopt_f_eval_at(gopt_env,optimizer%x,f=f,master=optimizer%master,&
             para_env=optimizer%para_env,force_env_section=force_env_section, error=error)
        ! add 'gradient=optimizer%gradient' to revaluate also the force?

        CALL force_env_get(force_env,subsys=subsys,cell=cell,error=error)
        DO isubsys=1,SIZE(subsys)
           CALL cp_subsys_get(subsys(isubsys)%subsys, particles=particles,&
                error=error)
           CALL write_particle_distances(particles%els,&
                cell,force_env_section,error)
        END DO

        CALL cp_opt_gopt_release(optimizer, error=error)
        CALL gopt_f_release(gopt_env, error=error)

     END IF
     CALL cp_rm_iter_level(logger%iter_info,ihandle,error)
     CALL timestop(0.0_dp,handle)
   END SUBROUTINE gopt_do_opt
!***************************************************************************

END MODULE geo_opt
