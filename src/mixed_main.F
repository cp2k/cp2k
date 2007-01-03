!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/mixed_main [1.0] *
!!
!!   NAME
!!     mixed_main
!!
!!   FUNCTION
!!     perform biased molecular dynamics (H= k H1 + (1-k) H2)
!!   AUTHOR
!!     fschiff 11.06 
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE mixed_main
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE mixed_environment,               ONLY: mixed_init
  USE mixed_environment_types,         ONLY: mixed_env_create,&
                                             mixed_env_release,&
                                             mixed_environment_type,&
                                             set_mixed_env
  USE force_env_methods,               ONLY: force_env_create
  USE force_env_types,                 ONLY: force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mixed_main'

! *** Global variables ***

  TYPE(mixed_environment_type), POINTER :: mixed_env

  PUBLIC :: mixed_create_force_env

!!***
!!-----------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------!
!
!-----------------------------------------------------------------------------!
!!****** mixed_main/mixed [1.0] *
!!
!!   NAME
!!     mixed
!!
!!   SYNOPSIS
!!     Subroutine mixed(globenv)
!!       Implicit None
!!       Type(global_environment_type), pointer:: globenv
!!     End Subroutine mixed
!!
!!   FUNCTION
!!     Controls program flow for mixed calculations
!!
!!   AUTHOR
!!     fschiff
!!
!!   MODIFICATION HISTORY
!!
!!   USED BY
!!     cp2k
!!
!!*** *************************************************************************

SUBROUTINE mixed_create_force_env ( force_env, root_section, para_env,globenv,&
                                   force_env_section, method_name_id,error )


! Argument
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: force_env_section
    INTEGER, INTENT(in), OPTIONAL            :: method_name_id
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mixed_create_force_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(mixed_environment_type), POINTER     :: mixed_env

!------------------------------------------------------------------------------

  CALL timeset(routineN,"I","",handle)

  CALL mixed_env_create( mixed_env, para_env = para_env ,error=error)

  CALL set_mixed_env (mixed_env, error=error)



! *** Read the input and the database files and perform further  ***
! *** initializations for the setup of the MIXED environment ***

  CALL mixed_init ( mixed_env, root_section, para_env, globenv, force_env_section, error )

  CALL force_env_create ( force_env, root_section, mixed_env = mixed_env, &
                          para_env = para_env, globenv = globenv, &
                          force_env_section=force_env_section, error = error )

  CALL mixed_env_release ( mixed_env, error = error )

  CALL timestop ( 0.0_dp, handle )

END SUBROUTINE mixed_create_force_env

! *****************************************************************************

END MODULE mixed_main
