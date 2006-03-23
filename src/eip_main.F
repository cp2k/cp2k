!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2006  CP2K developers group                            !
!-----------------------------------------------------------------------------!

!!****h* cp2k/eip_main [1.0] *
!!
!!   NAME
!!     eip_main
!!
!!   FUNCTION
!!     Main program which drives the EIP calculation
!!
!!   NOTES
!!     -
!!     
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     03.2006 initial create [tdk]
!!
!!   SOURCE
!****************************************************************************
MODULE eip_main
  USE eip_environment,                 ONLY: eip_init
  USE eip_environment_types,           ONLY: eip_env_create,&
                                             eip_env_release,&
                                             eip_environment_type
  USE force_env_methods,               ONLY: force_env_create
  USE force_env_types,                 ONLY: force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'eip_main'

  ! *** Public subroutines ***
  PUBLIC :: eip_create_force_env

!***
!****************************************************************************

CONTAINS

!!****f* eip_main/eip_create_force_env [1.0] *
!!
!!   NAME
!!     eip_create_force_env
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     -
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     03.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE eip_create_force_env(force_env, globenv, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'eip_create_force_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(eip_environment_type), POINTER      :: eip_env

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    NULLIFY(eip_env)

    CALL eip_env_create(eip_env, error=error)

    CALL eip_init(eip_env, globenv, error=error)

    CALL force_env_create(force_env, eip_env=eip_env, globenv=globenv, &
                          error=error)

    CALL eip_env_release(eip_env, error=error)

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE eip_create_force_env
!***************************************************************************

END MODULE eip_main
