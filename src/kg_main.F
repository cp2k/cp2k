!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      jgh: Added new KG_GPW method [01.03.2004]
!> \author gloriolina     (13.11.2003)
! *****************************************************************************
MODULE kg_main
  USE bibliography,                    ONLY: Iannuzzi2006,&
                                             cite_reference
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE force_env_methods,               ONLY: force_env_create
  USE force_env_types,                 ONLY: force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_type
  USE kg_environment,                  ONLY: kg_init
  USE kg_environment_types,            ONLY: kg_env_create,&
                                             kg_env_release,&
                                             kg_environment_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

    IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kg_main'

! *** Global variables ***

  TYPE(kg_environment_type), POINTER :: kg_env

  PUBLIC :: kg_create_force_env

CONTAINS

! *****************************************************************************
  SUBROUTINE kg_create_force_env(force_env,root_section,para_env, globenv,&
       force_env_section, subsys_section, use_motion_section, error)

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                subsys_section
    LOGICAL, INTENT(IN)                      :: use_motion_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'kg_create_force_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    CALL cite_reference(Iannuzzi2006)
    failure=.FALSE.
    CALL kg_env_create(kg_env,para_env=para_env,error=error)

    CALL kg_init( kg_env, root_section, para_env, globenv, force_env_section, &
         subsys_section, use_motion_section, error )
    CALL force_env_create(force_env,root_section, kg_env=kg_env, &
         para_env= para_env, globenv=globenv, force_env_section=force_env_section,&
         error=error)
    CALL kg_env_release(kg_env,error=error)

    CALL timestop(handle)
  END SUBROUTINE kg_create_force_env

END MODULE kg_main

