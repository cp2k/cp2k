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
  USE kinds,                           ONLY: default_path_length,&
                                             dp
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'metadynamics_types'

  PUBLIC :: hills_env_type, meta_env_type, meta_env_retain, set_meta_env, &
            meta_env_release, metavar_type, wall_type, multiple_walkers_type

! *****************************************************************************
!> \brief defines types for HILLS
! *****************************************************************************
  TYPE hills_env_type
     LOGICAL                                   :: restart
     REAL(KIND=dp)                             :: ww
     INTEGER                                   :: n_hills, nt_hills
     REAL(KIND=dp), DIMENSION(:,:), POINTER    :: ss_history
     REAL(KIND=dp), DIMENSION(:,:), POINTER    :: delta_s_history
     REAL(KIND=dp), DIMENSION(:),   POINTER    :: ww_history
  END TYPE hills_env_type

! *****************************************************************************
!> \brief defines types for WALLS
! *****************************************************************************
  TYPE wall_type
     INTEGER                                   :: id_type, id_direction
     REAL(KIND=dp)                             :: pos
     REAL(KIND=dp)                             :: k_quadratic, ww_gauss, sigma_gauss
  END TYPE wall_type

! *****************************************************************************
!> \brief defines types for COLVAR used in the metadynamics
! *****************************************************************************
  TYPE metavar_type
     INTEGER                                   :: icolvar
     LOGICAL                                   :: do_wall, periodic
     REAL(KIND=dp)                             :: mass, lambda, vvp
     REAL(KIND=dp)                             :: cv_energy, delta_s
     REAL(KIND=dp)                             :: ss, ss0, ff_s, ff_hills
     TYPE(wall_type), DIMENSION(:), POINTER    :: walls
  END TYPE metavar_type

! *****************************************************************************
!> \brief defines types for multiple walkers run
! *****************************************************************************
  TYPE multiple_walkers_type
     INTEGER                                   :: n_hills_local
     INTEGER                                   :: walker_id
     INTEGER                                   :: walkers_tot_nr
     INTEGER                                   :: walkers_freq_comm
     INTEGER, DIMENSION(:), POINTER            :: walkers_status
     CHARACTER(LEN=default_path_length),&
          DIMENSION(:), POINTER                :: walkers_file_name
  END TYPE multiple_walkers_type

! *****************************************************************************
!> \brief defines meta_env type
! *****************************************************************************
  TYPE meta_env_type
     INTEGER                                   :: ref_count, id_nr
     LOGICAL                                   :: do_hills, do_multiple_walkers
     LOGICAL                                   :: extended_lagrange
     INTEGER                                   :: n_colvar
     REAL(KIND=dp)                             :: ekin_s,cv_energy, dt
     LOGICAL                                   :: tempcontrol, restart
     REAL(KIND=dp)                             :: temp_wanted,toll_temp
     INTEGER                                   :: n_steps
     ! time pointer should never be allocated itself.. that's a quite bad
     ! choice derived from the md_env.. So whenever the meta_env type is 
     ! used the pointer time should be pointer to the high level time counter
     REAL(KIND=dp), POINTER                    :: time
     TYPE(hills_env_type), POINTER             :: hills_env
     TYPE(metavar_type), POINTER, DIMENSION(:) :: metavar
     TYPE(multiple_walkers_type), POINTER      :: multiple_walkers
     TYPE(cp_para_env_type), POINTER           :: para_env
     TYPE(section_vals_type), POINTER          :: metadyn_section
  END TYPE meta_env_type

CONTAINS

! *****************************************************************************
!> \brief sets the meta_env
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author alessandro laio and fawzi mohamed
! *****************************************************************************
  SUBROUTINE set_meta_env(meta_env, time, error)
    TYPE(meta_env_type), POINTER             :: meta_env
    REAL(KIND=dp), OPTIONAL, POINTER         :: time
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_meta_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    IF (ASSOCIATED(meta_env)) THEN
       IF (PRESENT(time)) THEN
          NULLIFY(meta_env%time)
          meta_env%time => time
       END IF
    END IF
  END SUBROUTINE set_meta_env

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
          ! Hills env
          CALL hills_env_release(meta_env%hills_env,error)
          ! Walkers type
          IF (ASSOCIATED(meta_env%multiple_walkers)) THEN
             IF (ASSOCIATED(meta_env%multiple_walkers%walkers_status)) THEN
                DEALLOCATE(meta_env%multiple_walkers%walkers_status,stat=stat)
                CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
             END IF
             IF (ASSOCIATED(meta_env%multiple_walkers%walkers_file_name)) THEN
                DEALLOCATE(meta_env%multiple_walkers%walkers_file_name,stat=stat)
                CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
             END IF
             DEALLOCATE(meta_env%multiple_walkers,stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          NULLIFY(meta_env%time)
          NULLIFY(meta_env%metadyn_section)
          DEALLOCATE(meta_env, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE meta_env_release

! *****************************************************************************
!> \brief releases the hills_env
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino [tlaino] - University of Zurich 10.2008
! *****************************************************************************
  SUBROUTINE hills_env_release(hills_env, error)
    TYPE(hills_env_type), POINTER            :: hills_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'hills_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(hills_env)) THEN
       IF (ASSOCIATED(hills_env%ss_history)) THEN
          DEALLOCATE(hills_env%ss_history,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(hills_env%delta_s_history)) THEN
          DEALLOCATE(hills_env%delta_s_history,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(hills_env%ww_history)) THEN
          DEALLOCATE(hills_env%ww_history,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
       DEALLOCATE(hills_env,stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    END IF
  END SUBROUTINE hills_env_release

END MODULE metadynamics_types
