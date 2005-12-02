!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

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
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_release,&
                                             force_env_retain,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: write_particle_distances,&
                                             write_qs_particle_coordinates
  USE simulation_cell,                 ONLY: cell_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE input_constants,                 ONLY: default_lbfgs_method_id,&
                                             default_bfgs_method_id,&
                                             default_cg_method_id
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_get,&
                                             section_vals_val_get

  IMPLICIT NONE
  PRIVATE

  REAL(KIND = dp), PARAMETER :: scale = 0.2_dp, scaleI = 1.0_dp/scale
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='gopt_f_types'
  INTEGER, PRIVATE, SAVE :: last_gopt_f_id=0, last_gopt_param_id=0

  PUBLIC :: gopt_f_type
  PUBLIC :: cp_eval_at, gopt_f_eval_at, gopt_f_create, gopt_f_retain,&
            gopt_f_release, gopt_f_create_x0

  PUBLIC :: gopt_param_type
  PUBLIC :: gopt_param_create, gopt_param_read, gopt_param_retain, &
            gopt_param_release

  INTERFACE cp_eval_at
     MODULE PROCEDURE gopt_f_eval_at
  END INTERFACE

!!***
!****************************************************************************
  
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
     TYPE(force_env_type), POINTER :: force_env
     TYPE(global_environment_type), POINTER :: globenv
  END TYPE gopt_f_type
!!***
  !****************************************************************************

!!****s* gopt_f_types/gopt_param_type [1.0] *
!!
!!   NAME
!!     gopt_param_type
!!
!!   FUNCTION
!!     stores the parameters nedeed for a geometry optimization
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - method_id: the optimization method
!!     - wanted_proj_gradient: the wanted accuracy on the projected gradient
!!     - wanted_rel_f_error:
!!     - max_h_rank: maximum rank of the hessian approximation
!!     - max_f_per_iter: maximum number of function evaluations per iteration
!!     - max_iter: maximum number of iterations
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2003 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE gopt_param_type
     INTEGER         :: id_nr, ref_count, method_id, max_h_rank
     INTEGER         :: max_f_per_iter, max_iter
     REAL(KIND = dp) :: wanted_proj_gradient, wanted_rel_f_error,&
                        max_dr,max_force,rms_dr,rms_force
     REAL(KIND = dp) :: brent_tol, brack_limit, restart_limit
     INTEGER         :: brent_max_iter, max_steep_steps
     LOGICAL         :: hesrest, Fletcher_Reeves_ON
     LOGICAL         :: use_fit_line_search
  END TYPE gopt_param_type
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
SUBROUTINE gopt_f_create(gopt_env, force_env, globenv, error)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), &
      INTENT(inout), TARGET                  :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gopt_f_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  ALLOCATE(gopt_env, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     gopt_env%ref_count=1
     last_gopt_f_id=last_gopt_f_id+1
     gopt_env%id_nr=last_gopt_f_id
     CALL force_env_retain(force_env, error=error)
     gopt_env%force_env => force_env
     gopt_env%globenv => globenv
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
        NULLIFY(gopt_env%force_env, gopt_env%globenv)
        DEALLOCATE(gopt_env, stat=stat)
        CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
     END IF
  END IF
END SUBROUTINE gopt_f_release
!***************************************************************************

!!****f* gopt_f_types/gopt_f_eval_at [1.0] *
!!
!!   NAME
!!     gopt_f_eval_at
!!
!!   FUNCTION
!!     evaluete the potential energy and its gradients
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     gopt_env: the geometry optimization environment
!!     x: the position where the function should be evaluated
!!     f: the function value
!!     gradient: the value of its gradient
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 3.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
SUBROUTINE gopt_f_eval_at(gopt_env, x, f, gradient, master, para_env, error)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp), DIMENSION(:), INTENT(in)  :: x
    REAL(KIND=dp), INTENT(out), OPTIONAL     :: f
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: gradient
    INTEGER, INTENT(in)                      :: master
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gopt_f_eval_at', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idg, idir, ip, isubsys
    LOGICAL                                  :: failure
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(particle_list_type), POINTER        :: particles

  failure=.FALSE.
  NULLIFY(particles, subsys, cell)
  logger => cp_error_get_logger(error)
  
  CALL timeset(routineN,"I","",handle)
  
  CALL force_env_get(gopt_env%force_env, subsys=subsys,&
       cell=cell,error=error)
  idg=0
  DO isubsys=1,SIZE(subsys)
     CALL cp_subsys_get(subsys(isubsys)%subsys, particles=particles,&
          error=error)
     DO ip=1,particles%n_els
        DO idir=1,3
           idg=idg+1
           particles%els(ip)%r(idir)=x(idg)*scale
        END DO
     END DO
     IF (logger%print_keys%atomic_coordinates) THEN
          CALL write_qs_particle_coordinates(particles%els,&
               cell,gopt_env%globenv)
       END IF
     IF (logger%print_keys%interatomic_distances) THEN
        CALL write_particle_distances(particles%els,&
             cell,gopt_env%globenv)
     END IF
  END DO

  CALL force_env_calc_energy_force(gopt_env%force_env,&
       calc_force=PRESENT(gradient),error=error)

  IF (PRESENT(f)) THEN
     CALL force_env_get(gopt_env%force_env,potential_energy=f,&
          error=error)
  END IF
  IF (PRESENT(gradient)) THEN
     IF (master==para_env%mepos) THEN ! we are on the master
        CPPrecondition(ASSOCIATED(gradient),cp_failure_level,routineP,error,failure)
        CALL force_env_get(gopt_env%force_env, subsys=subsys,&
             error=error)
        idg=0
        DO isubsys=1,SIZE(subsys)
           CALL cp_subsys_get(subsys(isubsys)%subsys, particles=particles,&
                error=error)
           DO ip=1,particles%n_els
              DO idir=1,3
                 idg=idg+1
                 CPInvariant(SIZE(gradient)>=idg,cp_failure_level,routineP,error,failure)
                 gradient(idg)=-particles%els(ip)%f(idir)*scaleI
              END DO
           END DO
        END DO
     END IF
  END IF
  CALL timestop(0.0_dp,handle)
END SUBROUTINE gopt_f_eval_at
!***************************************************************************

!!****f* gopt_f_types/gopt_f_create_x0 [1.0] *
!!
!!   NAME
!!     gopt_f_create_x0
!!
!!   FUNCTION
!!     returns the value of the parameters for the actual configuration
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     gopt_env: the geometry optimization environment you want the info about
!!     x0: the parameter vector (is allocated by this routine)
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
SUBROUTINE gopt_f_create_x0(gopt_env, x0, error)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x0
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gopt_f_create_x0', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: idg, idir, ip, isubsys, ndg, &
                                                stat
    LOGICAL                                  :: failure
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(particle_list_type), POINTER        :: particles

  failure=.FALSE.
  NULLIFY(particles,subsys)
  
  CALL force_env_get(gopt_env%force_env, subsys=subsys,&
       error=error)
  ndg=0
  DO isubsys=1,SIZE(subsys)
     CALL cp_subsys_get(subsys(isubsys)%subsys, particles=particles,&
          error=error)
     ndg=ndg+particles%n_els*3
  END DO
  
  ALLOCATE(x0(ndg),stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  
  idg=0
  DO isubsys=1,SIZE(subsys)
     CALL cp_subsys_get(subsys(isubsys)%subsys, particles=particles,&
          error=error)
     DO ip=1,particles%n_els
        DO idir=1,3
           idg=idg+1
           CPInvariant(SIZE(x0)>=idg,cp_failure_level,routineP,error,failure)
           x0(idg)=particles%els(ip)%r(idir)*scaleI
        END DO
     END DO
  END DO
END SUBROUTINE gopt_f_create_x0
!***************************************************************************

!!****f* gopt_f_types/gopt_param_create [1.0] *
!!
!!   NAME
!!     gopt_param_create
!!
!!   FUNCTION
!!     creates a new gopt_param object
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - gopt_param: the object to create
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
SUBROUTINE gopt_param_create(gopt_param,error)
    TYPE(gopt_param_type), POINTER           :: gopt_param
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gopt_param_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  ALLOCATE(gopt_param, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  
  IF (.NOT. failure) THEN
     last_gopt_param_id=last_gopt_param_id+1
     gopt_param%id_nr=last_gopt_param_id
     gopt_param%ref_count=1
     gopt_param%method_id=default_bfgs_method_id
     gopt_param%max_h_rank=5
     gopt_param%max_f_per_iter=3
     gopt_param%max_iter=200
     gopt_param%wanted_proj_gradient=1.0e-3_dp
     gopt_param%wanted_rel_f_error=1.0e5_dp
     gopt_param%hesrest=.FALSE.
     gopt_param%Fletcher_Reeves_ON=.FALSE.
     gopt_param%use_fit_line_search=.FALSE.
     gopt_param%max_steep_steps=0
     gopt_param%max_dr = 0.0030_dp
     gopt_param%max_force = 0.00045_dp
     gopt_param%rms_dr = 0.0015_dp
     gopt_param%rms_force = 0.00030_dp
     gopt_param%brent_tol = 1.0E-4_dp
     gopt_param%brent_max_iter = 100
     gopt_param%brack_limit = 100.0_dp
     gopt_param%restart_limit = 0.1_dp
  END IF
END SUBROUTINE gopt_param_create
!***************************************************************************

!!****f* gopt_f_types/gopt_param_read [1.0] *
!!
!!   NAME
!!     gopt_param_read
!!
!!   FUNCTION
!!     reads the parameters of the geopmetry optimizer
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - 
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
SUBROUTINE gopt_param_read(gopt_param, globenv, error)
    TYPE(gopt_param_type), POINTER           :: gopt_param
    TYPE(global_environment_type), &
      INTENT(inout), TARGET                  :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gopt_param_read', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=3)                         :: test_result
    CHARACTER(LEN=60)                        :: keyword, keyword2
    LOGICAL                                  :: failure, section_found
    TYPE(section_vals_type), POINTER         :: gopt_section

    failure=.FALSE.

    IF (.NOT.ASSOCIATED(gopt_param)) THEN
       CALL gopt_param_create(gopt_param, error=error)
    END IF

    NULLIFY(gopt_section)
    gopt_section =>  section_vals_get_subs_vals(globenv%input_file,&
         "MOTION%GEOOPT",error=error)
    failure=.FALSE.
    IF (.NOT.failure) THEN
       CALL section_vals_val_get(gopt_section,"OPTIMIZER",i_val=gopt_param%method_id,error=error)
       CALL section_vals_val_get(gopt_section,"MAX_ITER",i_val=gopt_param%max_iter,error=error)
       CALL section_vals_val_get(gopt_section,"MAX_H_RANK",i_val=gopt_param%max_h_rank,error=error)
       CALL section_vals_val_get(gopt_section,"MAX_F_PER_ITER",i_val=gopt_param%max_f_per_iter,error=error)
       CALL section_vals_val_get(gopt_section,"WANTED_PROJ_GRADIENT",r_val=gopt_param%wanted_proj_gradient,error=error)
       CALL section_vals_val_get(gopt_section,"WANTED_REL_F_ERROR",r_val=gopt_param%wanted_rel_f_error,error=error)
       CALL section_vals_val_get(gopt_section,"FLETCHER_REEVES",l_val=gopt_param%Fletcher_Reeves_ON,error=error)
       CALL section_vals_val_get(gopt_section,"USE_FIT_LINE_SEARCH",l_val=gopt_param%use_fit_line_search,error=error)
       CALL section_vals_val_get(gopt_section,"RESTART",l_val=gopt_param%hesrest,error=error)
       CALL section_vals_val_get(gopt_section,"MAX_STEEP_STEPS",i_val=gopt_param%max_steep_steps,error=error)
       CALL section_vals_val_get(gopt_section,"MAX_DR",r_val=gopt_param%max_dr,error=error)
       CALL section_vals_val_get(gopt_section,"MAX_FORCE",r_val=gopt_param%max_force,error=error)
       CALL section_vals_val_get(gopt_section,"RMS_DR",r_val=gopt_param%rms_dr,error=error)
       CALL section_vals_val_get(gopt_section,"RMS_FORCE",r_val=gopt_param%rms_force,error=error)
       CALL section_vals_val_get(gopt_section,"BRENT_TOL",r_val=gopt_param%brent_tol,error=error)
       CALL section_vals_val_get(gopt_section,"BRENT_MAX_ITER",i_val=gopt_param%brent_max_iter,error=error)
       CALL section_vals_val_get(gopt_section,"BRACK_LIMIT",r_val=gopt_param%brack_limit,error=error)
       CALL section_vals_val_get(gopt_section,"CG_RESTART_LIMIT",r_val=gopt_param%restart_limit,error=error)
    END IF

  END SUBROUTINE gopt_param_read

!!****f* gopt_f_types/gopt_param_retain [1.0] *
!!
!!   NAME
!!     gopt_param_retain
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     gopt_param: the geometry optimization environment to retain
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
SUBROUTINE gopt_param_retain(gopt_param, error)
    TYPE(gopt_param_type), POINTER           :: gopt_param
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gopt_param_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(gopt_param),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     CPPreconditionNoFail(gopt_param%ref_count>0,cp_failure_level,routineP,error)
     gopt_param%ref_count=gopt_param%ref_count+1
  END IF
END SUBROUTINE gopt_param_retain
!***************************************************************************


!!****f* gopt_f_types/gopt_param_release [1.0] *
!!
!!   NAME
!!     gopt_param_release
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     gopt_param: the geometry optimization environment to release
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
SUBROUTINE gopt_param_release(gopt_param, error)
    TYPE(gopt_param_type), POINTER           :: gopt_param
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gopt_param_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (ASSOCIATED(gopt_param)) THEN
     CPPreconditionNoFail(gopt_param%ref_count>0,cp_failure_level,routineP,error)
     gopt_param%ref_count=gopt_param%ref_count-1
     IF (gopt_param%ref_count==0) THEN
        DEALLOCATE(gopt_param, stat=stat)
        CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
     END IF
  END IF
END SUBROUTINE gopt_param_release
!***************************************************************************

END MODULE gopt_f_types
