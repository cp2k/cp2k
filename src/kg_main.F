!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2003  CP2K developers group
!-----------------------------------------------------------------------------!
!!****h* cp2k/kg_main [1.0] *
!!
!!   NAME
!!     kg_main
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     gloriolina     (13.11.2003)
!!
!!   LITERATURE
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!!
MODULE kg_main
! USE qs_energy,                       ONLY: qs_energies
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE force_env_types,                 ONLY: force_env_create,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kg_environment,                  ONLY: kg_init
  USE kg_environment_types,            ONLY: kg_env_create,&
                                             kg_env_release,&
                                             kg_environment_type
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop

    IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: moduleN = "kg_main"

! *** Global variables ***

  TYPE(kg_environment_type), POINTER :: kg_env

  PUBLIC :: kg_create_force_env

! *****************************************************************************

CONTAINS

!!******
!******************************************************************************

!!****f* kg_main/kg_create_force_env [1.0] *
!!
!!   NAME
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

  SUBROUTINE kg_create_force_env(force_env,globenv,error)

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "kg_create_force_env", &
      routineP = moduleN//":"//routineN

    INTEGER                                  :: handle

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    CALL kg_env_create(kg_env,para_env=globenv%para_env)

!   *** Read the input and the database files and perform further  ***
!   *** initializations for the setup of the KG environment ***

    CALL kg_init( kg_env, globenv )

    CALL force_env_create(force_env, kg_env=kg_env, globenv=globenv,&
         error=error)

    CALL kg_env_release(kg_env,error=error)
    
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE kg_create_force_env

! *****************************************************************************

END MODULE kg_main
!!*****
