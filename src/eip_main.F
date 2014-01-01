!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Main program which drives the EIP calculation
!> \par History
!>      03.2006 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
MODULE eip_main
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE eip_environment,                 ONLY: eip_init
  USE eip_environment_types,           ONLY: eip_env_create,&
                                             eip_env_release,&
                                             eip_environment_type
  USE f77_blas
  USE force_env_methods,               ONLY: force_env_create
  USE force_env_types,                 ONLY: force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'eip_main'

  ! *** Public subroutines ***
  PUBLIC :: eip_create_force_env

!***

CONTAINS

! *****************************************************************************
!> \par History
!>      03.2006 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE eip_create_force_env(force_env,root_section, para_env, globenv,&
       force_env_section, subsys_section, use_motion_section, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                subsys_section
    LOGICAL, INTENT(IN)                      :: use_motion_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'eip_create_force_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(eip_environment_type), POINTER      :: eip_env

    CALL timeset(routineN,handle)
    NULLIFY(eip_env)

    CALL eip_env_create(eip_env, error=error)
    CALL eip_init(eip_env, root_section, para_env, force_env_section=force_env_section,&
         subsys_section=subsys_section, use_motion_section=use_motion_section, error=error)

    CALL force_env_create(force_env,root_section, eip_env=eip_env, para_env=para_env,&
         globenv=globenv, force_env_section=force_env_section, error=error)
    CALL eip_env_release(eip_env, error=error)
    CALL timestop(handle)

  END SUBROUTINE eip_create_force_env

END MODULE eip_main
