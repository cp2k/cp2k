!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief perform classical molecular dynamics and path integral simulations
!> \par History
!>      gt SEPT-23-2002: part is allocated/deallocated/initialized in
!>                       read_coord_vel
!>      CJM rewrite
!> \author CJM-Sept-01-02
! *****************************************************************************
MODULE fist_main
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_type
  USE fist_environment,                ONLY: fist_init
  USE fist_environment_types,          ONLY: fist_env_create,&
                                             fist_env_release,&
                                             fist_env_set,&
                                             fist_environment_type
  USE force_env_methods,               ONLY: force_env_create
  USE force_env_types,                 ONLY: force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_type
  USE qmmm_types,                      ONLY: qmmm_env_mm_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fist_main'

! *** Global variables ***
  PUBLIC :: fist_create_force_env

!!-----------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------!
! FIST FIST FIST FIST FIST FIST FIST FIST FIST FIST FIST FIST FIST FIST FIST  !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Controls program flow for classical MD and path-integrals
!> \param force_env ...
!> \param root_section ...
!> \param para_env ...
!> \param globenv ...
!> \param qmmm ...
!> \param qmmm_env ...
!> \param force_env_section ...
!> \param subsys_section ...
!> \param use_motion_section ...
!> \param prev_subsys ...
!> \param error ...
!> \par Used By
!>      cp2k
!> \author CJM
! *****************************************************************************
  SUBROUTINE fist_create_force_env ( force_env, root_section, para_env, globenv,&
       qmmm, qmmm_env, force_env_section, subsys_section, use_motion_section, prev_subsys, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    LOGICAL, OPTIONAL                        :: qmmm
    TYPE(qmmm_env_mm_type), OPTIONAL, &
      POINTER                                :: qmmm_env
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                subsys_section
    LOGICAL, INTENT(IN)                      :: use_motion_section
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: prev_subsys
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fist_create_force_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: myqmmm
    TYPE(fist_environment_type), POINTER     :: fist_env

    CALL timeset(routineN,handle)
    myqmmm=.FALSE.
    IF (PRESENT(qmmm)) THEN
       myqmmm=qmmm
    END IF

    CALL fist_env_create( fist_env, para_env = para_env ,error=error)
    IF (PRESENT(qmmm_env)) THEN
       CALL fist_env_set (fist_env, qmmm=myqmmm, qmmm_env=qmmm_env, error=error)
    ELSE
       CALL fist_env_set (fist_env, qmmm=myqmmm, error=error)
    END IF
    ! *** Read the input and the database files and perform further  ***
    ! *** initializations for the setup of the FIST environment ***
    CALL fist_init ( fist_env, root_section, para_env, force_env_section,&
         subsys_section, use_motion_section, prev_subsys=prev_subsys, error=error )

    CALL force_env_create ( force_env, root_section, fist_env = fist_env, &
         para_env = para_env, globenv = globenv, &
         force_env_section=force_env_section, error = error )

    CALL fist_env_release ( fist_env, error = error )
    CALL timestop(handle)
  END SUBROUTINE fist_create_force_env

END MODULE fist_main
