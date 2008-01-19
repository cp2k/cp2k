!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****h* cp2k/gopt_param_types [1.0] *
!!
!!   NAME
!!     gopt_param_types
!!
!!   FUNCTION
!!     contains typo and related routines to handle parameters controlling the
!!     GEO_OPT module
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
MODULE gopt_param_types
  USE f77_blas
  USE input_constants,                 ONLY: default_bfgs_method_id,&
                                             default_cg_method_id,&
                                             default_lbfgs_method_id,&
                                             default_minimization_method_id,&
                                             default_ts_method_id,&
                                             none_ts_method_id
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gopt_param_types'
  INTEGER, PRIVATE, SAVE :: last_gopt_f_id=0, last_gopt_param_id=0

  PUBLIC :: gopt_param_type
  PUBLIC :: gopt_param_create, gopt_param_read, gopt_param_retain, &
            gopt_param_release

!!****s* gopt_param_types/gopt_param_type [1.0] *
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
     LOGICAL         :: Fletcher_Reeves
     LOGICAL         :: use_fit_line_search
     INTEGER         :: id_nr, ref_count, method_id, type_id
     INTEGER         :: ts_method_id
     INTEGER         :: max_f_per_iter, max_iter, max_h_rank
     INTEGER         :: brent_max_iter, max_steep_steps
     REAL(KIND = dp) :: wanted_proj_gradient, wanted_rel_f_error
     REAL(KIND = dp) :: max_dr, max_force, rms_dr, rms_force
     REAL(KIND = dp) :: brent_tol, brack_limit, restart_limit, initial_step
     REAL(KIND = dp) :: dimer_angle_tol
  END TYPE gopt_param_type
!!***
  !****************************************************************************

CONTAINS
  
!!****f* gopt_param_types/gopt_param_create [1.0] *
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
    TYPE(cp_error_type), INTENT(inout)       :: error

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
    END IF
  END SUBROUTINE gopt_param_create
!***************************************************************************

!!****f* gopt_param_types/gopt_param_read [1.0] *
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
  SUBROUTINE gopt_param_read(gopt_param, gopt_section, type_id, error)
    TYPE(gopt_param_type), POINTER           :: gopt_param
    TYPE(section_vals_type), POINTER         :: gopt_section
    INTEGER, INTENT(IN), OPTIONAL            :: type_id
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gopt_param_read', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(gopt_param),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(gopt_section),cp_failure_level,routineP,error,failure)
    CALL gopt_param_create(gopt_param, error=error)

    failure=.FALSE.
    IF (.NOT.failure) THEN
       IF (PRESENT(type_id)) THEN
          gopt_param%type_id = type_id
       ELSE
          CALL section_vals_val_get(gopt_section,"TYPE",i_val=gopt_param%type_id,error=error)
       END IF
       CALL section_vals_val_get(gopt_section,"OPTIMIZER",i_val=gopt_param%method_id,error=error)

       CALL section_vals_val_get(gopt_section,"MAX_ITER",i_val=gopt_param%max_iter,error=error)
       CALL section_vals_val_get(gopt_section,"MAX_DR",r_val=gopt_param%max_dr,error=error)
       CALL section_vals_val_get(gopt_section,"MAX_FORCE",r_val=gopt_param%max_force,error=error)
       CALL section_vals_val_get(gopt_section,"RMS_DR",r_val=gopt_param%rms_dr,error=error)
       CALL section_vals_val_get(gopt_section,"RMS_FORCE",r_val=gopt_param%rms_force,error=error)

       SELECT CASE(gopt_param%method_id)
       CASE(default_lbfgs_method_id)
          CALL section_vals_val_get(gopt_section,"LBFGS%MAX_H_RANK",i_val=gopt_param%max_h_rank,error=error)
          CALL section_vals_val_get(gopt_section,"LBFGS%MAX_F_PER_ITER",i_val=gopt_param%max_f_per_iter,error=error)
          CALL section_vals_val_get(gopt_section,"LBFGS%WANTED_PROJ_GRADIENT",r_val=gopt_param%wanted_proj_gradient,error=error)
          CALL section_vals_val_get(gopt_section,"LBFGS%WANTED_REL_F_ERROR",r_val=gopt_param%wanted_rel_f_error,error=error)
       CASE(default_bfgs_method_id)
          ! Do nothing
       CASE(default_cg_method_id)
          CALL section_vals_val_get(gopt_section,"CG%FLETCHER_REEVES",l_val=gopt_param%Fletcher_Reeves,error=error)
          CALL section_vals_val_get(gopt_section,"CG%USE_FIT_LINE_SEARCH",l_val=gopt_param%use_fit_line_search,error=error)
          CALL section_vals_val_get(gopt_section,"CG%MAX_STEEP_STEPS",i_val=gopt_param%max_steep_steps,error=error)
          CALL section_vals_val_get(gopt_section,"CG%BRENT_TOL",r_val=gopt_param%brent_tol,error=error)
          CALL section_vals_val_get(gopt_section,"CG%INITIAL_STEP",r_val=gopt_param%initial_step,error=error)
          CALL section_vals_val_get(gopt_section,"CG%BRENT_MAX_ITER",i_val=gopt_param%brent_max_iter,error=error)
          CALL section_vals_val_get(gopt_section,"CG%BRACK_LIMIT",r_val=gopt_param%brack_limit,error=error)
          CALL section_vals_val_get(gopt_section,"CG%RESTART_LIMIT",r_val=gopt_param%restart_limit,error=error)
       END SELECT
       
       SELECT CASE(gopt_param%type_id)
       CASE(default_minimization_method_id)
          ! Do Nothing
          gopt_param%ts_method_id=none_ts_method_id
       CASE(default_ts_method_id)
          CALL section_vals_val_get(gopt_section,"TRANSITION_STATE%METHOD",i_val=gopt_param%ts_method_id,error=error)
       END SELECT
    END IF

  END SUBROUTINE gopt_param_read

!!****f* gopt_param_types/gopt_param_retain [1.0] *
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
    TYPE(cp_error_type), INTENT(inout)       :: error

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

!!****f* gopt_param_types/gopt_param_release [1.0] *
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
    TYPE(cp_error_type), INTENT(inout)       :: error

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

END MODULE gopt_param_types
