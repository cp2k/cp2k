!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief contains typo and related routines to handle parameters controlling the
!>      GEO_OPT module
!> \par History
!>      none
! *****************************************************************************
MODULE gopt_param_types
  
  USE input_constants,                 ONLY: default_bfgs_method_id,&
                                             default_cell_method_id,&
                                             default_cg_method_id,&
                                             default_lbfgs_method_id,&
                                             default_minimization_method_id,&
                                             default_ts_method_id,&
                                             none_ts_method_id
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gopt_param_types'
  INTEGER, PRIVATE, SAVE :: last_gopt_param_id=0

  PUBLIC :: gopt_param_type
  PUBLIC :: gopt_param_read,&
            gopt_param_release

! *****************************************************************************
!> \brief stores the parameters nedeed for a geometry optimization
!> \param method_id the optimization method
!> \param wanted_proj_gradient the wanted accuracy on the projected gradient
!> \param wanted_rel_f_error :
!> \param max_h_rank maximum rank of the hessian approximation
!> \param max_f_per_iter maximum number of function evaluations per iteration
!> \param max_iter maximum number of iterations
!> \par History
!>      08.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE cg_ls_param_type
     INTEGER         :: brent_max_iter, type_id
     REAL(KIND = dp) :: brent_tol, brack_limit, initial_step, max_step
     LOGICAL         :: grad_only
  END TYPE cg_ls_param_type

! *****************************************************************************
  TYPE gopt_param_type
     LOGICAL                  :: Fletcher_Reeves
     INTEGER                  :: id_nr, ref_count
     INTEGER                  :: method_id, type_id
     INTEGER                  :: ts_method_id, cell_method_id, shellcore_method_id
     INTEGER                  :: max_f_per_iter, max_iter, max_h_rank
     INTEGER                  :: max_steep_steps
     REAL(KIND = dp)          :: restart_limit
     REAL(KIND = dp)          :: wanted_proj_gradient, wanted_rel_f_error
     REAL(KIND = dp)          :: max_dr, max_force, rms_dr, rms_force
     REAL(KIND = dp)          :: dimer_angle_tol
     TYPE(cg_ls_param_type)   :: cg_ls
  END TYPE gopt_param_type

CONTAINS

! *****************************************************************************
!> \brief creates a new gopt_param object
!> \param gopt_param the object to create
!> \par History
!>      08.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE gopt_param_create(gopt_param)
    TYPE(gopt_param_type), POINTER           :: gopt_param

    CHARACTER(len=*), PARAMETER :: routineN = 'gopt_param_create', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    ALLOCATE(gopt_param)

    last_gopt_param_id=last_gopt_param_id+1
    gopt_param%id_nr=last_gopt_param_id
    gopt_param%ref_count=1
  END SUBROUTINE gopt_param_create

! *****************************************************************************
!> \brief reads the parameters of the geopmetry optimizer
!> \param gopt_param ...
!> \param gopt_section ...
!> \param type_id ...
!> \par History
!>      08.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE gopt_param_read(gopt_param, gopt_section, type_id)
    TYPE(gopt_param_type), POINTER           :: gopt_param
    TYPE(section_vals_type), POINTER         :: gopt_section
    INTEGER, INTENT(IN), OPTIONAL            :: type_id

    CHARACTER(len=*), PARAMETER :: routineN = 'gopt_param_read', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPASSERT(.NOT.ASSOCIATED(gopt_param))
    CPASSERT(ASSOCIATED(gopt_section))
    CALL gopt_param_create(gopt_param)

    failure=.FALSE.
    IF (PRESENT(type_id)) THEN
       gopt_param%type_id = type_id
    ELSE
       CALL section_vals_val_get(gopt_section,"TYPE",i_val=gopt_param%type_id)
    END IF
    CALL section_vals_val_get(gopt_section,"OPTIMIZER",i_val=gopt_param%method_id)

    CALL section_vals_val_get(gopt_section,"MAX_ITER",i_val=gopt_param%max_iter)
    CALL section_vals_val_get(gopt_section,"MAX_DR",r_val=gopt_param%max_dr)
    CALL section_vals_val_get(gopt_section,"MAX_FORCE",r_val=gopt_param%max_force)
    CALL section_vals_val_get(gopt_section,"RMS_DR",r_val=gopt_param%rms_dr)
    CALL section_vals_val_get(gopt_section,"RMS_FORCE",r_val=gopt_param%rms_force)

    SELECT CASE(gopt_param%method_id)
    CASE(default_lbfgs_method_id)
       CALL section_vals_val_get(gopt_section,"LBFGS%MAX_H_RANK",i_val=gopt_param%max_h_rank)
       CALL section_vals_val_get(gopt_section,"LBFGS%MAX_F_PER_ITER",i_val=gopt_param%max_f_per_iter)
       CALL section_vals_val_get(gopt_section,"LBFGS%WANTED_PROJ_GRADIENT",r_val=gopt_param%wanted_proj_gradient)
       CALL section_vals_val_get(gopt_section,"LBFGS%WANTED_REL_F_ERROR",r_val=gopt_param%wanted_rel_f_error)
    CASE(default_bfgs_method_id)
       ! Do nothing
    CASE(default_cg_method_id)
       CALL section_vals_val_get(gopt_section,"CG%FLETCHER_REEVES",l_val=gopt_param%Fletcher_Reeves)
       CALL section_vals_val_get(gopt_section,"CG%MAX_STEEP_STEPS",i_val=gopt_param%max_steep_steps)
       CALL section_vals_val_get(gopt_section,"CG%RESTART_LIMIT",r_val=gopt_param%restart_limit)
       CALL section_vals_val_get(gopt_section,"CG%LINE_SEARCH%TYPE",i_val=gopt_param%cg_ls%type_id)
       CALL section_vals_val_get(gopt_section,"CG%LINE_SEARCH%GOLD%INITIAL_STEP",&
            r_val=gopt_param%cg_ls%initial_step)
       CALL section_vals_val_get(gopt_section,"CG%LINE_SEARCH%GOLD%BRENT_TOL",&
            r_val=gopt_param%cg_ls%brent_tol)
       CALL section_vals_val_get(gopt_section,"CG%LINE_SEARCH%GOLD%BRENT_MAX_ITER",&
            i_val=gopt_param%cg_ls%brent_max_iter)
       CALL section_vals_val_get(gopt_section,"CG%LINE_SEARCH%GOLD%BRACK_LIMIT",&
            r_val=gopt_param%cg_ls%brack_limit)
       CALL section_vals_val_get(gopt_section,"CG%LINE_SEARCH%2PNT%MAX_ALLOWED_STEP",&
            r_val=gopt_param%cg_ls%max_step)
       CALL section_vals_val_get(gopt_section,"CG%LINE_SEARCH%2PNT%LINMIN_GRAD_ONLY",&
            l_val=gopt_param%cg_ls%grad_only)
    END SELECT

    SELECT CASE(gopt_param%type_id)
    CASE(default_minimization_method_id)
       ! Do Nothing
       gopt_param%ts_method_id=none_ts_method_id
    CASE(default_ts_method_id)
       CALL section_vals_val_get(gopt_section,"TRANSITION_STATE%METHOD",i_val=gopt_param%ts_method_id)
    CASE(default_cell_method_id)
       CALL section_vals_val_get(gopt_section,"TYPE",i_val=gopt_param%cell_method_id)
    END SELECT

  END SUBROUTINE gopt_param_read

! *****************************************************************************
!> \brief ...
!> \param gopt_param the geometry optimization environment to retain
!> \par History
!>      none
! *****************************************************************************
  SUBROUTINE gopt_param_retain(gopt_param)
    TYPE(gopt_param_type), POINTER           :: gopt_param

    CHARACTER(len=*), PARAMETER :: routineN = 'gopt_param_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPASSERT(ASSOCIATED(gopt_param))
    CPASSERT(gopt_param%ref_count>0)
    gopt_param%ref_count=gopt_param%ref_count+1
  END SUBROUTINE gopt_param_retain

! *****************************************************************************
!> \brief ...
!> \param gopt_param the geometry optimization environment to release
!> \par History
!>      none
! *****************************************************************************
  SUBROUTINE gopt_param_release(gopt_param)
    TYPE(gopt_param_type), POINTER           :: gopt_param

    CHARACTER(len=*), PARAMETER :: routineN = 'gopt_param_release', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(gopt_param)) THEN
       CPASSERT(gopt_param%ref_count>0)
       gopt_param%ref_count=gopt_param%ref_count-1
       IF (gopt_param%ref_count==0) THEN
          DEALLOCATE(gopt_param)
       END IF
    END IF
  END SUBROUTINE gopt_param_release

END MODULE gopt_param_types
