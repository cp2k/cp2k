!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!


!!****** cp2k/cg_optimizer [1.0] *
!!
!!   NAME
!!     cg_optimizer
!!
!!   FUNCTION
!!     Routines for Geometry optimization using  Conjugate Gradients
!!
!!   AUTHOR
!!     Teodoro Laino [teo]
!!     10.2005
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!***************************************************************************

MODULE cg_optimizer
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE cg_utils,                        ONLY: cg_linmin,&
                                             get_conjugate_direction
  USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                             cp_iterate,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_unit_nr,&
                                             cp_rm_iter_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE gopt_f_methods,                  ONLY: check_converg,&
                                             check_rot_conv,&
                                             print_geo_opt_header,&
                                             write_cycle_infos,&
                                             write_final_info,&
                                             write_geo_traj,&
                                             write_rot_cycle_infos
  USE gopt_f_types,                    ONLY: gopt_f_type
  USE gopt_param_types,                ONLY: gopt_param_type
  USE input_cp2k_restarts,             ONLY: write_restart
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE termination,                     ONLY: external_control
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
#include "gopt_f77_methods.h"

  PUBLIC :: geoopt_cg
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cg_optimizer'

CONTAINS

!!****f* cg_optimizer/geoopt_cg [1.0] *
!!
!!   NAME
!!     geoopt_cg
!!
!!   FUNCTION
!!     Driver for conjugate gradient optimization technique
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     10.2005 created [tlaino]
!!
!!*** **********************************************************************
  RECURSIVE SUBROUTINE geoopt_cg(force_env,gopt_param,globenv,force_env_section,&
       geo_section, gopt_env, x0, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(gopt_param_type), POINTER           :: gopt_param
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: force_env_section, geo_section
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x0
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'geoopt_cg', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ihandle, output_unit
    LOGICAL                                  :: failure
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: motion_section, root_section

    CALL timeset(routineN,"I","",handle)
    para_env=>force_env%para_env
    root_section=>force_env%root_section
    logger => cp_error_get_logger(error)
    CALL cp_add_iter_level(logger%iter_info,TRIM(gopt_env%label),ihandle,error)
    failure = .FALSE.
    IF (.NOT.failure) THEN
       motion_section => section_vals_get_subs_vals(root_section,"MOTION",error=error)
       output_unit = cp_print_key_unit_nr(logger,geo_section,"PRINT%PROGRAM_RUN_INFO",&
            extension=".geoLog",error=error)

       CALL force_env_get(force_env,subsys=subsys,error=error)
       CALL cp_subsys_get(subsys=subsys,atomic_kinds=atomic_kinds,particles=particles,error=error)
       CALL print_geo_opt_header(gopt_env, output_unit, "CONJUGATE GRADIENTS")

       particle_set => particles%els
       CALL  cp_cg_main(force_env, x0, gopt_param,  output_unit, globenv,&
            force_env_section, motion_section, gopt_env, error)

       CALL cp_print_key_finished_output(output_unit,logger,geo_section,&
            "PRINT%PROGRAM_RUN_INFO", error=error)
    END IF
    CALL cp_rm_iter_level(logger%iter_info,ihandle,error)
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE geoopt_cg

!!****f* cg_optimizer/cp_cg_main [1.0] *
!!
!!   NAME
!!     cp_cg_main
!!
!!   FUNCTION
!!     This really performs the conjugate gradients optimization
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     10.2005 created [tlaino]
!!
!!*** **********************************************************************
  RECURSIVE SUBROUTINE cp_cg_main(force_env, x0, gopt_param, output_unit, globenv,&
       force_env_section, motion_section, gopt_env,  error)
    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x0
    TYPE(gopt_param_type), POINTER           :: gopt_param
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                motion_section
    TYPE(gopt_f_type), POINTER               :: gopt_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_cg_main', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=5)                         :: wildcard
    INTEGER                                  :: brent_max_iter, handle, its, &
                                                max_steep_steps, maxiter, stat
    LOGICAL                                  :: conv, failure, &
                                                Fletcher_Reeves, should_stop
    REAL(KIND=dp)                            :: brack_limit, brent_tol, emin, &
                                                eold, opt_energy, res_lim, &
                                                step
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: xold
    REAL(KIND=dp), DIMENSION(:), POINTER     :: g, h, xi
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: root_section

    CALL timeset(routineN,handle)

    failure           = .FALSE.
    NULLIFY(logger,g, h, xi)
    root_section => force_env%root_section
    para_env     => force_env%para_env
    logger       => cp_error_get_logger(error)

    IF (gopt_env%dimer_rotation) CALL cp_unimplemented_error(fromWhere=routineP, &
                                 message="DIMER method not yet working with CG.", &
                                 error=error, error_level=cp_failure_level)

    IF (.NOT.failure) THEN
       maxiter            = gopt_param%max_iter
       max_steep_steps    = gopt_param%max_steep_steps
       Fletcher_Reeves    = gopt_param%Fletcher_Reeves
       brent_tol          = gopt_param%brent_tol
       brent_max_iter     = gopt_param%brent_max_iter
       brack_limit        = gopt_param%brack_limit
       res_lim            = gopt_param%restart_limit
       step               = gopt_param%initial_step
       ALLOCATE(g (SIZE(x0)), stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(h (SIZE(x0)), stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(xi(SIZE(x0)), stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(xold(SIZE(x0)), stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

       ! Evaluate energy and forces at the first step
       CALL cp_eval_at(gopt_env, x0, opt_energy, xi, gopt_env%force_env%para_env%mepos,&
            gopt_env%force_env%para_env, error)

       g  = -xi
       h  =   g
       xi =   h
       emin    = HUGE(0.0_dp)
       ! Main Loop
       wildcard = "   SD"
       IF (.NOT.gopt_env%dimer_rotation) THEN
          CALL write_cycle_infos(output_unit,it=0,etot=opt_energy,wildcard=wildcard)
       ELSE
          CALL write_rot_cycle_infos(output_unit,it=0,etot=opt_energy,dimer_env=gopt_env%dimer_env,&
               wildcard=wildcard)
       END IF
       DO its = 1, maxiter
          gopt_env%its = its
          CALL cp_iterate(logger%iter_info,iter_nr=its,error=error)
          xold = x0
          eold = opt_energy
          emin = MIN(emin, opt_energy)

          CALL cg_linmin(gopt_env, x0, xi, opt_energy, step, output_unit, gopt_param, globenv, error)

          ! Check for an external exit command
          CALL external_control(should_stop,"GEO",globenv,error)
          IF(should_stop) EXIT
          
          ! Some IO and Convergence check
          IF (.NOT.gopt_env%dimer_rotation) THEN
             CALL write_restart(force_env=force_env,root_section=root_section,&
                  force_env_section=force_env_section,error=error)
             CALL write_cycle_infos(output_unit,its,opt_energy,opt_energy-eold,emin=emin, wildcard=wildcard)
             CALL check_converg(SIZE(x0),x0-xold,xi,output_unit,conv,gopt_param)
             CALL write_geo_traj(force_env,  root_section, its, opt_energy, error)
          ELSE
             CALL write_rot_cycle_infos(output_unit,its,opt_energy,opt_energy-eold,emin,gopt_env%dimer_env,&
                  wildcard)
             CALL check_rot_conv(gopt_env%dimer_env, output_unit, conv, error)
          END IF

          IF (conv.OR.(its==maxiter))  EXIT
          CALL cp_eval_at(gopt_env, x0, opt_energy, xi, gopt_env%force_env%para_env%mepos,&
               gopt_env%force_env%para_env, error)

          ! Reset Condition or Steepest Descent Requested
          IF ( ABS(DOT_PRODUCT(g,h))/SQRT((DOT_PRODUCT(g,g)*DOT_PRODUCT(h,h))) < res_lim &
               .OR. its+1 <= max_steep_steps) THEN
             ! Steepest Descent
             wildcard = "   SD"
             g  = -xi
             h  =   g
          ELSE
             ! Conjugate Gradients
             wildcard = "   CG"
             CALL get_conjugate_direction(gopt_env%dimer_rotation, Fletcher_Reeves, g, xi, h)
          END IF
          xi  = h
       END DO
       IF(its == maxiter .AND. (.NOT.conv))THEN
          IF(output_unit>0)THEN
             WRITE(UNIT=output_unit,FMT="(/,T2,A)")&
                  "*** MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED ***"
             IF (.NOT.gopt_env%dimer_rotation) THEN
                WRITE(UNIT=output_unit,FMT="(T2,A)")&
                     "***        EXITING GEOMETRY OPTIMIZATION         ***"
             ELSE
                WRITE(UNIT=output_unit,FMT="(T2,A)")&
                     "***        EXITING ROTATION OPTIMIZATION         ***"
             END IF
             CALL m_flush(output_unit)
          END IF
       END IF
       ! Write final particle information, if converged
       CALL cp_iterate(logger%iter_info,last=.TRUE.,iter_nr=its,error=error)
       IF (.NOT.gopt_env%dimer_rotation) THEN
          CALL write_final_info(conv, its, gopt_env, x0, gopt_env%force_env%para_env%mepos,&
               gopt_env%force_env%para_env, force_env, gopt_env%motion_section, root_section,&
               force_env_section, error)
       END IF

       DEALLOCATE(xold, stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(g, stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(h, stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(xi,stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    CALL timestop(handle)

  END SUBROUTINE cp_cg_main

END MODULE cg_optimizer
