!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****h* cp2k/metadynamics_types *
!!
!!   NAME
!!     metadynamics_types
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
!!     01.2005 created [fawzi and ale]
!!
!!   SOURCE
!!
!****************************************************************************
MODULE metadynamics_types
  USE cp_para_env,                     ONLY: cp_para_env_release
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE 

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'metadynamics_types'

  PUBLIC :: hills_env_type, meta_env_type, meta_env_retain, &
            meta_env_release, metavar_type

  TYPE hills_env_type
     LOGICAL                                :: restart
     REAL(KIND=dp)                          :: ww
     INTEGER                                :: n_hills, nt_hills
     REAL(KIND=dp), DIMENSION(:,:), POINTER :: ss_history
     REAL(KIND=dp), DIMENSION(:,:), POINTER :: delta_s_history
     REAL(KIND=dp), DIMENSION(:),   POINTER :: ww_history
  END TYPE hills_env_type

  TYPE metavar_type
     REAL(KIND=dp) :: mass, lambda, ss, ss0, ff_s, vv, vvp, &
                      cv_energy, delta_s, ff_hills
     LOGICAL       :: do_wall_p,do_wall_m
     REAL(KIND=dp) :: wall_p, wall_m
     REAL(KIND=dp) :: TARGET
     INTEGER       :: icolvar,isubsys
  END TYPE metavar_type

  TYPE meta_env_type
     INTEGER       :: ref_count, id_nr
     LOGICAL       :: do_hills
     LOGICAL       :: extended_lagrange
     INTEGER       :: n_colvar
     REAL(KIND=dp) :: ekin_s,cv_energy, dt
     LOGICAL       :: tempcontrol, restart
     REAL(KIND=dp) :: temp_wanted,toll_temp
     INTEGER       :: n_steps
     TYPE(hills_env_type)                      :: hills_env
     TYPE(metavar_type), POINTER, DIMENSION(:) :: metavar
     TYPE(cp_para_env_type), POINTER           :: para_env
     TYPE(section_vals_type), POINTER          :: metadyn_section
  END TYPE meta_env_type

CONTAINS
  
!!****f* metadynamics/meta_env_retain *
!!
!!   NAME
!!     meta_env_retain
!!
!!   FUNCTION
!!     retains the meta_env
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     alessandro laio and fawzi mohamed
!!
!!*** **********************************************************************
SUBROUTINE meta_env_retain(meta_env,error)
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'meta_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPreconditionNoFail(ASSOCIATED(meta_env),cp_failure_level,routineP,error)
  CPPreconditionNoFail(meta_env%ref_count>0,cp_failure_level,routineP,error)
  meta_env%ref_count=meta_env%ref_count+1
END SUBROUTINE meta_env_retain
!***************************************************************************

!!****f* metadynamics/meta_env_release *
!!
!!   NAME
!!     meta_env_release
!!
!!   FUNCTION
!!     releases the meta_env
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     alessandro laio and fawzi mohamed
!!
!!*** **********************************************************************
SUBROUTINE meta_env_release(meta_env,error)
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'meta_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  IF (ASSOCIATED(meta_env)) THEN
     CPPreconditionNoFail(meta_env%ref_count>0,cp_failure_level,routineP,error)
     meta_env%ref_count=meta_env%ref_count-1
     IF (meta_env%ref_count==0) THEN
        CALL cp_para_env_release(meta_env%para_env)
        IF (ASSOCIATED(meta_env%metavar)) THEN
           DEALLOCATE(meta_env%metavar,stat=stat)
           CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
        END IF
        IF (ASSOCIATED(meta_env%hills_env%ss_history)) THEN
           DEALLOCATE(meta_env%hills_env%ss_history,stat=stat)
           CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
        END IF
        IF (ASSOCIATED(meta_env%hills_env%delta_s_history)) THEN
           DEALLOCATE(meta_env%hills_env%delta_s_history,stat=stat)
           CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
        END IF
        IF (ASSOCIATED(meta_env%hills_env%ww_history)) THEN
           DEALLOCATE(meta_env%hills_env%ww_history,stat=stat)
           CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
        END IF
        NULLIFY(meta_env%metadyn_section)
        DEALLOCATE(meta_env, stat=stat)
        CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
     END IF
  END IF
  NULLIFY(meta_env)
END SUBROUTINE meta_env_release
!***************************************************************************

END MODULE metadynamics_types
