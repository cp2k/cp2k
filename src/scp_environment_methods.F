!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief to complete the setup of the scp_env
!> \author CJM
! *****************************************************************************
MODULE scp_environment_methods
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE distribution_2d_types,           ONLY: distribution_2d_type
  USE f77_blas
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_scp_methods,               ONLY: qs_ks_scp_create
  USE qs_ks_scp_types,                 ONLY: qs_ks_scp_env_type,&
                                             qs_ks_scp_release
  USE scp_environment_types,           ONLY: scp_environment_type,&
                                             set_scp_env
  USE scp_rho_methods,                 ONLY: scp_rho_create
  USE scp_rho_types,                   ONLY: scp_rho_release,&
                                             scp_rho_type
  USE scp_rspw_methods,                ONLY: scp_qs_rspw_rebuild,&
                                             scp_rspw_create
  USE scp_rspw_types,                  ONLY: scp_rspw_release,&
                                             scp_rspw_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: scp_qs_env_setup

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scp_environment_methods'

!-----------------------------------------------------------------------------!

CONTAINS

! *****************************************************************************
  SUBROUTINE scp_qs_env_setup(scp_env,qs_env,error)

    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_qs_env_setup', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_ks_scp_env_type), POINTER        :: ks_scp_env
    TYPE(scp_rho_type), POINTER              :: scp_rho
    TYPE(scp_rspw_type), POINTER             :: rspw

  failure=.FALSE.

  NULLIFY( para_env )
  NULLIFY( blacs_env )
  NULLIFY( scp_rho )
  NULLIFY( rspw )
  NULLIFY( distribution_2d )
  NULLIFY( ks_scp_env )

  CALL get_qs_env ( qs_env = qs_env, blacs_env=blacs_env, & 
                                       para_env=para_env, &
                                         pw_env = pw_env, &
                         distribution_2d = distribution_2d, &
                         error=error)

  CALL pw_env_get ( pw_env, auxbas_pw_pool = auxbas_pw_pool, error=error )

! borrow the parallel, blacs, and 2d-distribution environment from QS
  CALL set_scp_env ( scp_env, distribution_2d = distribution_2d, &
                                            para_env = para_env, &
                                            blacs_env = blacs_env,&
                                            error=error )
  
! borrow the plane-wave environment from QS
  CALL scp_rspw_create ( rspw, error )
  CALL scp_qs_rspw_rebuild ( rspw, qs_env, error )
  CALL set_scp_env ( scp_env, rspw = rspw, error = error )
  CALL scp_rspw_release ( rspw, error = error )

! create rho structures
  CALL pw_env_get ( pw_env, auxbas_pw_pool = auxbas_pw_pool, error=error )
  CALL scp_rho_create ( scp_rho = scp_rho, auxbas_pw_pool=auxbas_pw_pool, error = error)
  CALL set_scp_env ( scp_env, rho_scp = scp_rho, error = error)
  CALL scp_rho_release ( rho_scp = scp_rho, error = error)

! create the SCP KS-env
  CALL qs_ks_scp_create ( ks_scp_env, qs_env, error )
  CALL set_scp_env ( scp_env, ks_scp_env = ks_scp_env, error = error )
  CALL qs_ks_scp_release ( ks_scp_env, error = error )
  
  END SUBROUTINE scp_qs_env_setup

END MODULE scp_environment_methods
