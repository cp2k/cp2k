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

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE bfgs_optimizer,                  ONLY: geoopt_bfgs
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
  USE cp_subsystem_types,               ONLY: cp_subsystem_p_type,&
                                             cp_subsystem_type,&
                                             cp_subsys_get
  USE cp_lbfgs_optimizer_gopt,         ONLY: cp_lbfgs_opt_gopt_type,&
                                             cp_opt_gopt_create,&
                                             cp_opt_gopt_get,&
                                             cp_opt_gopt_next,&
                                             cp_opt_gopt_release,&
                                             cp_opt_gopt_stop
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE force_env_types,                   ONLY: force_env_calculate_force,&
                                             force_env_create,&
                                             force_env_get,&
                                             force_env_release,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE gopt_f_types,                    ONLY: gopt_f_create,&
                                             gopt_f_create_x0,&
                                             gopt_f_eval_at,&
                                             gopt_f_release,&
                                             gopt_f_retain,&
                                             gopt_f_type,&
                                             gopt_param_create,&
                                             gopt_param_read,&
                                             gopt_param_release,&
                                             gopt_param_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_bcast
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type,&
                                             write_particle_distances
  USE qs_energy,                       ONLY: qs_energies
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force,                        ONLY: qs_forces
  USE qs_force_types,                  ONLY: init_qs_force,&
                                             qs_force_type
  USE simulation_cell,                 ONLY: cell_type
  USE timings,                         ONLY: timeset,&
                                             timestop
 IMPLICIT NONE
 PRIVATE

 CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='geo_opt'

 PUBLIC :: qs_geo_opt
 CONTAINS

   SUBROUTINE qs_geo_opt(force_env, globenv, error)

     TYPE(global_environment_type), POINTER          :: globenv
     TYPE(force_env_type), POINTER                   :: force_env
     TYPE(cp_error_type), INTENT(inout), OPTIONAL    :: error

!     *** Local variables ***

     TYPE(qs_environment_type), POINTER :: qs_env
     TYPE(gopt_param_type), POINTER :: gopt_param
     LOGICAL :: failure
     CHARACTER(len=*), PARAMETER :: routineN='qs_geo_opt',&
          routineP=moduleN//':'//routineN

     CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
     CPPrecondition(ASSOCIATED(globenv),cp_failure_level,routineP,error,failure)
     
     NULLIFY (gopt_param)

     CALL gopt_param_read(gopt_param, globenv, error=error)
     
     IF (gopt_param%method == "BFGS") THEN
       CALL force_env_get(force_env=force_env,qs_env=qs_env)
       CALL geoopt_bfgs(qs_env,globenv,gopt_param)
     ELSE
       CALL gopt_do_opt(force_env,gopt_param,globenv,error=error)
     END IF  

     CALL gopt_param_release(gopt_param, error=error)

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
SUBROUTINE gopt_do_opt(force_env, gopt_param, globenv, error)
  TYPE(force_env_type), POINTER :: force_env
  TYPE(gopt_param_type), POINTER :: gopt_param
  TYPE(global_environment_type), TARGET, INTENT(inout) :: globenv
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure, should_stop
  CHARACTER(len=*), PARAMETER :: routineN='gopt_do_opt',&
        routineP=moduleN//':'//routineN
  TYPE(gopt_f_type), POINTER :: gopt_env
  TYPE(cp_lbfgs_opt_gopt_type), POINTER :: optimizer
  TYPE(cp_para_env_type), POINTER :: para_env
  TYPE(cell_type), POINTER :: cell
  TYPE(cp_logger_type), POINTER :: logger
  TYPE(particle_list_type), POINTER :: particles
  TYPE(cp_subsystem_p_type), DIMENSION(:), POINTER :: subsys 
  INTEGER :: i, stat, isubsys, handle
  REAL(KIND = dp) :: f
  REAL(KIND = dp), DIMENSION(:), POINTER :: x0

  CALL timeset(routineN,"I","",handle)
  failure=.FALSE.
  NULLIFY(gopt_env, optimizer, para_env, cell, particles, subsys, x0)
  logger => cp_error_get_logger(error)
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
          max_f_per_iter=gopt_param%max_f_per_iter, error=error)

     DEALLOCATE(x0,stat=stat)
     CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)

     DO i=1,gopt_param%max_iter
        IF (.NOT.cp_opt_gopt_next(optimizer)) EXIT
        IF (para_env%mepos==para_env%source) &
             INQUIRE(file="EXIT",exist=should_stop)
        CALL mp_bcast(should_stop,para_env%source,para_env%group)
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
          para_env=optimizer%para_env,error=error)
     ! add 'gradient=optimizer%gradient' to revaluate also the force?

     CALL force_env_get(force_env,subsys=subsys,cell=cell,error=error)
     DO isubsys=1,SIZE(subsys)
        CALL cp_subsys_get(subsys(isubsys)%subsys, particles=particles,&
             error=error)
        CALL write_particle_distances(particles%els,&
             cell,&
             globenv)
     END DO

     CALL cp_opt_gopt_release(optimizer, error=error)
     CALL gopt_f_release(gopt_env, error=error)
     
  END IF
  CALL timestop(0.0_dp,handle)
END SUBROUTINE gopt_do_opt
!***************************************************************************

END MODULE geo_opt
