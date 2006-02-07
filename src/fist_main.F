!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/fist_main [1.0] *
!!
!!   NAME
!!     fist_main
!!
!!   FUNCTION
!!     perform classical molecular dynamics and path integral simulations
!!   AUTHOR
!!     CJM-Sept-01-02
!!
!!   MODIFICATION HISTORY
!!     gt SEPT-23-2002: part is allocated/deallocated/initialized in
!!                      read_coord_vel
!!     CJM rewrite 
!!   SOURCE
!******************************************************************************

MODULE fist_main
  USE fist_environment,                ONLY: fist_init
  USE fist_environment_types,          ONLY: fist_env_create,&
                                             fist_env_release,&
                                             fist_environment_type,&
                                             set_fist_env
  USE force_env_methods,               ONLY: force_env_create
  USE force_env_types,                 ONLY: force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE qmmm_types,                      ONLY: qmmm_env_mm_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fist_main'

! *** Global variables ***

  TYPE(fist_environment_type), POINTER :: fist_env

  PUBLIC :: fist_create_force_env
  
!!***
!!-----------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------!
! FIST FIST FIST FIST FIST FIST FIST FIST FIST FIST FIST FIST FIST FIST FIST  !
!-----------------------------------------------------------------------------!
!!****** fist_main/fist [1.0] *
!!
!!   NAME
!!     fist
!!
!!   SYNOPSIS
!!     Subroutine fist(globenv)
!!       Implicit None
!!       Type(global_environment_type), pointer:: globenv
!!     End Subroutine fist
!!
!!   FUNCTION
!!     Controls program flow for classical MD and path-integrals
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!
!!   USED BY
!!     cp2k
!!
!!*** *************************************************************************

SUBROUTINE fist_create_force_env ( force_env, globenv, qmmm, qmmm_env, error )


! Argument
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    LOGICAL, OPTIONAL                        :: qmmm
    TYPE(qmmm_env_mm_type), OPTIONAL, &
      POINTER                                :: qmmm_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fist_create_force_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: myqmmm
    TYPE(fist_environment_type), POINTER     :: fist_env

!------------------------------------------------------------------------------

  CALL timeset(routineN,"I","",handle)
  myqmmm=.FALSE.
  IF (PRESENT(qmmm)) THEN
     myqmmm=qmmm
  END IF

  CALL fist_env_create( fist_env, para_env = globenv%para_env )
  IF (PRESENT(qmmm_env)) THEN
     CALL set_fist_env (fist_env, qmmm=myqmmm, qmmm_env=qmmm_env, error=error)
  ELSE
     CALL set_fist_env (fist_env, qmmm=myqmmm, error=error)
  END IF


! *** Read the input and the database files and perform further  ***
! *** initializations for the setup of the FIST environment ***

  CALL fist_init ( fist_env, globenv, error )

  CALL force_env_create ( force_env, fist_env = fist_env, &
                          globenv = globenv, &
                          error = error )

  CALL fist_env_release ( fist_env, error = error )
    
  CALL timestop ( 0.0_dp, handle )

  END SUBROUTINE fist_create_force_env

! *****************************************************************************

END MODULE fist_main
