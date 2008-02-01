!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief defines types for metadynamics calculation
!> \par History
!>      01.2005 created [fawzi and ale]
! *****************************************************************************
MODULE metadynamics_types
  USE cp_para_env,                     ONLY: cp_para_env_release
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'metadynamics_types'

  PUBLIC :: hills_env_type, meta_env_type, meta_env_retain, &
            meta_env_release, metavar_type, wall_type

! *****************************************************************************
  TYPE hills_env_type
     LOGICAL                                :: restart
     REAL(KIND=dp)                          :: ww
     INTEGER                                :: n_hills, nt_hills
     REAL(KIND=dp), DIMENSION(:,:), POINTER :: ss_history
     REAL(KIND=dp), DIMENSION(:,:), POINTER :: delta_s_history
     REAL(KIND=dp), DIMENSION(:),   POINTER :: ww_history
  END TYPE hills_env_type

! *****************************************************************************
  TYPE wall_type
     INTEGER       :: id_type, id_direction
     REAL(KIND=dp) :: pos
     REAL(KINd=dp) :: k_quadratic, ww_gauss, sigma_gauss
  END TYPE wall_type

! *****************************************************************************
  TYPE metavar_type
     INTEGER       :: icolvar
     LOGICAL       :: do_wall, periodic
     REAL(KIND=dp) :: mass, lambda, vvp
     REAL(KIND=dp) :: cv_energy, delta_s, ff_hills, ss, ss0, ff_s
     TYPE(wall_type), DIMENSION(:), POINTER :: walls
  END TYPE metavar_type

! *****************************************************************************
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

! *****************************************************************************
!> \brief retains the meta_env
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author alessandro laio and fawzi mohamed
! *****************************************************************************
SUBROUTINE meta_env_retain(meta_env,error)
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'meta_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPreconditionNoFail(ASSOCIATED(meta_env),cp_failure_level,routineP,error)
  CPPreconditionNoFail(meta_env%ref_count>0,cp_failure_level,routineP,error)
  meta_env%ref_count=meta_env%ref_count+1
END SUBROUTINE meta_env_retain

! *****************************************************************************
!> \brief releases the meta_env
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author alessandro laio and fawzi mohamed
! *****************************************************************************
SUBROUTINE meta_env_release(meta_env,error)
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'meta_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  IF (ASSOCIATED(meta_env)) THEN
     CPPreconditionNoFail(meta_env%ref_count>0,cp_failure_level,routineP,error)
     meta_env%ref_count=meta_env%ref_count-1
     IF (meta_env%ref_count==0) THEN
        CALL cp_para_env_release(meta_env%para_env,error=error)
        IF (ASSOCIATED(meta_env%metavar)) THEN
           DO i = 1, SIZE(meta_env%metavar)
              IF (ASSOCIATED(meta_env%metavar(i)%walls)) THEN
                 DEALLOCATE(meta_env%metavar(i)%walls,stat=stat)
                 CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)  
              END IF
           END DO
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

END MODULE metadynamics_types
