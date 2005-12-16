!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_ddapc_types [1.0] *
!!
!!   NAME
!!     cp_ddapc_types
!!
!!   FUNCTION
!!     contains information regarding the decoupling/recoupling method of Bloechl
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!*****************************************************************************
MODULE cp_ddapc_types
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
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type 
  USE pw_types,                        ONLY: pw_type
  USE cp_ddapc_methods,                ONLY: ddapc_eval_gfunc,&
                                             ddapc_eval_AmI,&
                                             ewald_ddapc_pot
  USE input_section_types,             ONLY: section_type,&
                                             section_create, &
                                             section_add_keyword,&
                                             section_add_subsection,&
                                             section_release,&
                                             section_vals_type,&
                                             section_vals_get_subs_vals,&
                                             section_vals_get,&
                                             section_vals_val_get
  USE simulation_cell,                 ONLY: cell_type

  IMPLICIT NONE
  PRIVATE
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_ddapc_types'  
  INTEGER, PRIVATE, SAVE :: last_cp_ddapc_id=0
  PUBLIC :: cp_ddapc_type, cp_ddapc_create, cp_ddapc_retain, cp_ddapc_release
  
!!****s* cp_ddapc_types/cp_ddapc_type [1.0] *
!!
!!   NAME
!!     cp_ddapc_type
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     -
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     
!!
!!   SOURCE
!!***************************************************************************
  TYPE cp_ddapc_type
     INTEGER :: ref_count, id_nr
     REAL(KIND=dp) :: c0
     REAL(KIND=dp), DIMENSION(:,:), POINTER     :: AmI, Md, Mr, Mt
     REAL(KIND=dp)   , POINTER, DIMENSION(:,:)  :: gfunc     
     REAL(KIND=dp)   , POINTER, DIMENSION(:)    :: w
  END TYPE cp_ddapc_type
  
CONTAINS
  
!!****f* cp_ddapc_types/cp_ddapc_create [1.0] *
!!
!!   NAME
!!     cp_ddapc_create
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!
!!   AUTHOR
!!     Tedoro Laino
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** **********************************************************************
  SUBROUTINE cp_ddapc_create(cp_ddapc_env, qmmm_periodic, particle_set, radii,&
       cell, super_cell, rho_tot_g, gcut, iw, Vol, fitcharge_section, error)
    ! Arguments
    TYPE(cp_ddapc_type), POINTER               :: cp_ddapc_env
    TYPE(particle_type), DIMENSION(:), &
         POINTER                               :: particle_set
    TYPE(section_vals_type), POINTER           :: fitcharge_section
    REAL(kind=dp), DIMENSION(:), POINTER       :: radii
    TYPE(pw_type), POINTER                     :: rho_tot_g
    REAL(KIND=dp), INTENT(IN)                  :: gcut, Vol
    LOGICAL, INTENT(in)                        :: qmmm_periodic
    INTEGER, INTENT(IN)                        :: iw
    TYPE(cell_type), POINTER                   :: cell, super_cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ddapc_create', &
         routineP = moduleN//':'//routineN
    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.    
    IF (ASSOCIATED(cp_ddapc_env)) RETURN

    ALLOCATE(cp_ddapc_env, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       cp_ddapc_env%ref_count=1
       last_cp_ddapc_id=last_cp_ddapc_id+1
       cp_ddapc_env%id_nr=last_cp_ddapc_id
       NULLIFY(cp_ddapc_env%AmI,&
               cp_ddapc_env%Md,&
               cp_ddapc_env%Mt,&
               cp_ddapc_env%Mr,&
               cp_ddapc_env%gfunc,&
               cp_ddapc_env%w)
       !
       ! Evaluates gfunc and AmI 
       !
       CALL ddapc_eval_gfunc(cp_ddapc_env%gfunc, cp_ddapc_env%w, gcut, rho_tot_g, radii, error)
       CALL ddapc_eval_AmI(cp_ddapc_env%AmI, cp_ddapc_env%c0, &
            cp_ddapc_env%gfunc, cp_ddapc_env%w, particle_set, &
            gcut, rho_tot_g, radii, iw, Vol, error)
       !
       ! Evaluate the matrix for the Classical contribution to the coupling/decoupling scheme
       !
       CALL ewald_ddapc_pot(.FALSE.,    cell,       fitcharge_section, particle_set,&
               cp_ddapc_env%Md, radii, error)
       ALLOCATE(cp_ddapc_env%Mt(SIZE(cp_ddapc_env%Md,1),SIZE(cp_ddapc_env%Md,2)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (qmmm_periodic) THEN
          CALL ewald_ddapc_pot(qmmm_periodic, super_cell, fitcharge_section, particle_set,&
               cp_ddapc_env%Mr, radii, error)
          cp_ddapc_env%Mt = cp_ddapc_env%Md + cp_ddapc_env%Mr
       ELSE
          cp_ddapc_env%Mt = cp_ddapc_env%Md 
       END IF
    END IF
  END SUBROUTINE cp_ddapc_create
!***************************************************************************

!!****f* cp_ddapc_types/cp_ddapc_retain [1.0] *
!!
!!   NAME
!!     cp_ddapc_retain
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_ddapc_retain(cp_ddapc_env, error)
    ! Arguments
    TYPE(cp_ddapc_type), POINTER               :: cp_ddapc_env
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ddapc_retain', &
         routineP = moduleN//':'//routineN
    
    LOGICAL                                  :: failure
    
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(cp_ddapc_env),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CPPreconditionNoFail(cp_ddapc_env%ref_count>0,cp_failure_level,routineP,error)
       cp_ddapc_env%ref_count=cp_ddapc_env%ref_count+1
    END IF
  END SUBROUTINE cp_ddapc_retain

!!****f* cp_ddapc_types/cp_ddapc_release [1.0] *
!!
!!   NAME
!!     cp_ddapc_release
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_ddapc_release(cp_ddapc_env, error)
    ! Arguments
    TYPE(cp_ddapc_type), POINTER               :: cp_ddapc_env
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ddapc_release', &
         routineP = moduleN//':'//routineN  
    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    
    failure=.FALSE.    
    IF (ASSOCIATED(cp_ddapc_env)) THEN
       CPPreconditionNoFail(cp_ddapc_env%ref_count>0,cp_failure_level,routineP,error)
       cp_ddapc_env%ref_count=cp_ddapc_env%ref_count-1
       IF (cp_ddapc_env%ref_count==0) THEN
          IF (ASSOCIATED(cp_ddapc_env%AmI)) THEN
             DEALLOCATE(cp_ddapc_env%AmI, stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(cp_ddapc_env%Mt)) THEN
             DEALLOCATE(cp_ddapc_env%Mt, stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(cp_ddapc_env%Md)) THEN
             DEALLOCATE(cp_ddapc_env%Md, stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(cp_ddapc_env%Mr)) THEN
             DEALLOCATE(cp_ddapc_env%Mr, stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(cp_ddapc_env%gfunc)) THEN
             DEALLOCATE(cp_ddapc_env%gfunc, stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(cp_ddapc_env%w)) THEN
             DEALLOCATE(cp_ddapc_env%w, stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          NULLIFY(cp_ddapc_env%AmI,&
                  cp_ddapc_env%Md,&
                  cp_ddapc_env%Mr,&
                  cp_ddapc_env%gfunc,&
                  cp_ddapc_env%w)          
          DEALLOCATE(cp_ddapc_env, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          NULLIFY(cp_ddapc_env)
       END IF
    END IF
  END SUBROUTINE cp_ddapc_release

END MODULE cp_ddapc_types
