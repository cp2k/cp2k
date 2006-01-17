!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2006  CP2K developers group                            !
!-----------------------------------------------------------------------------!

!!****h* cp2k/tbmd_main [1.0] *
!!
!!   NAME
!!     tbmd_main
!!
!!   FUNCTION
!!     Main program which drives the tight-binding calculation
!!
!!   NOTES
!!     -
!!     
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!   SOURCE
!****************************************************************************
MODULE tbmd_main
  USE force_env_types,                 ONLY: force_env_type
  USE global_types,                    ONLY: global_environment_type

#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER :: moduleN = 'tbmd_main'

  ! *** Public subroutines ***
  PUBLIC :: tbmd_create_force_env

!***
!****************************************************************************

CONTAINS

!!****f* tbmd_main/tbmd_create_force_env [1.0] *
!!
!!   NAME
!!     tbmd_create_force_env
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
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_create_force_env(force_env, globenv, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

  END SUBROUTINE tbmd_create_force_env
!***************************************************************************

END MODULE tbmd_main
