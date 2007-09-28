!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!


!!****h* cp2k/free_energy_types *
!!
!!   NAME
!!     free_energy_types
!!
!!   FUNCTION
!!    defines types for metadynamics calculation
!!
!!   NOTES
!!
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!     01.2007 created [tlaino] Teodoro Laino
!!
!!   SOURCE
!!
!****************************************************************************
MODULE free_energy_types
  USE input_constants,                 ONLY: do_fe_ac,&
                                             do_fe_ui
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'free_energy_types'

  PUBLIC :: free_energy_type, fe_env_release, fe_env_create, ui_var_type


  TYPE ui_var_type
     REAL(KIND=dp), DIMENSION(:), POINTER  :: ss 
     INTEGER       :: icolvar
  END TYPE ui_var_type

  TYPE free_energy_type
     INTEGER       :: ncolvar, TYPE
     REAL(KIND=dp) :: eps_conv
     REAL(KIND=dp), DIMENSION(:,:), POINTER    :: covmx
     TYPE(ui_var_type), POINTER, DIMENSION(:)  :: uivar
  END TYPE free_energy_type

CONTAINS

!!****f* free_energy_types/fe_env_create *
!!
!!   NAME
!!     fe_env_create
!!
!!   FUNCTION
!!     creates the fe_env
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Teodoro Laino 01.2007
!!
!!*** **********************************************************************
  SUBROUTINE fe_env_create(fe_env,motion_section,error)
    TYPE(free_energy_type), POINTER          :: fe_env
    TYPE(section_vals_type), POINTER         :: motion_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fe_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, id_method, stat
    LOGICAL                                  :: explicit, failure
    TYPE(section_vals_type), POINTER         :: fe_section, ui_section, &
                                                ui_var_section

    failure=.FALSE.
    CPPreconditionNoFail(.NOT.ASSOCIATED(fe_env),cp_failure_level,routineP,error)
    fe_section => section_vals_get_subs_vals(motion_section,"FREE_ENERGY",error=error)

    CALL section_vals_get(fe_section,explicit=explicit, error=error)
    IF (explicit) THEN
       CALL section_vals_val_get(fe_section,"METHOD",i_val=id_method,error=error)
       SELECT CASE(id_method)
       CASE (do_fe_ui)
          ALLOCATE(fe_env, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          NULLIFY(fe_env%covmx, fe_env%uivar)
          fe_env%type = id_method
          ui_section => section_vals_get_subs_vals(fe_section,"UMBRELLA_INTEGRATION",error=error)
          ui_var_section => section_vals_get_subs_vals(ui_section,"UVAR",error=error)
          CALL section_vals_get(ui_var_section,n_repetition=fe_env%ncolvar,error=error)
          ALLOCATE(fe_env%covmx(fe_env%ncolvar,fe_env%ncolvar), stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          fe_env%covmx = HUGE(0.0_dp)
          ALLOCATE(fe_env%uivar(fe_env%ncolvar), stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          DO i = 1, fe_env%ncolvar
             ! Read Umbrella Integration Variable definition
             CALL section_vals_val_get(ui_var_section,"COLVAR",&
                  i_val=fe_env%uivar(fe_env%ncolvar)%icolvar,i_rep_section=i,error=error)
             NULLIFY(fe_env%uivar(fe_env%ncolvar)%ss)
          END DO
          CALL section_vals_val_get(ui_section,"EPS_CONV",r_val=fe_env%eps_conv,error=error)
       CASE (do_fe_ac)
          ALLOCATE(fe_env, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          NULLIFY(fe_env%covmx, fe_env%uivar)
          ALLOCATE(fe_env%covmx(3,0), stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          fe_env%type = id_method
          CALL section_vals_val_get(fe_section,"ALCHEMICAL_CHANGE%EPS_CONV",r_val=fe_env%eps_conv,error=error)          
       CASE DEFAULT
          ! Do Nothing
       END SELECT
    END IF
  END SUBROUTINE fe_env_create

!!****f* free_energy_types/fe_env_release *
!!
!!   NAME
!!     fe_env_release
!!
!!   FUNCTION
!!     releases the fe_env
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Laino Teodoro 01.2007
!!
!!*** **********************************************************************
  SUBROUTINE fe_env_release(fe_env,error)
    TYPE(free_energy_type), POINTER          :: fe_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fe_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(fe_env)) THEN
       IF (ASSOCIATED(fe_env%covmx)) THEN
          DEALLOCATE(fe_env%covmx,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(fe_env%uivar)) THEN
          DO i = 1, SIZE(fe_env%uivar)
             IF (ASSOCIATED(fe_env%uivar(i)%ss)) THEN
                DEALLOCATE(fe_env%uivar(i)%ss,stat=stat)
                CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
             END IF
          END DO
          DEALLOCATE(fe_env%uivar,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
       DEALLOCATE(fe_env, stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    END IF
  END SUBROUTINE fe_env_release
  !***************************************************************************

END MODULE free_energy_types
