!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C ) 2002 - 2003  CP2K developers group                         !
!-----------------------------------------------------------------------------!
!!****** cp2k/ewald_pw_methods [1.0] *
!!
!!   NAME
!!     ewald_pw_methods
!!
!!   FUNCTION
!!     pw_methods
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE ewald_pw_methods
  USE cell_types,                      ONLY: cell_type
  USE dg_rho0_types,                   ONLY: dg_rho0_get,&
                                             dg_rho0_init,&
                                             dg_rho0_set,&
                                             dg_rho0_type
  USE dg_types,                        ONLY: dg_get,&
                                             dg_set,&
                                             dg_type
  USE dgs,                             ONLY: dg_grid_change
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE ewald_pw_types,                  ONLY: ewald_pw_get,&
                                             ewald_pw_set,&
                                             ewald_pw_type
  USE input_constants,                 ONLY: do_ewald_ewald,&
                                             do_ewald_none,&
                                             do_ewald_pme,&
                                             do_ewald_spme
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_change
  USE pw_poisson_methods,              ONLY: pw_poisson_set
! CJM, pw_poisson_rebuild
  USE pw_poisson_types,                ONLY: pw_poisson_create,&
                                             pw_poisson_release,&
                                             pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_p_type,&
                                             pw_pool_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ewald_pw_methods'
  PUBLIC :: ewald_pw_grid_change

!-----------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------!

!******************************************************************************
!!****** ewald_pw_methods/ewald_pw_grid_change [1.0] *
!!
!!   NAME
!!     ewald_pw_grid_change
!!
!!   FUNCTION
!!     Rescales pw_grids after a box change
!!
!!   AUTHOR
!!     JGH (15-Mar-2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE ewald_pw_grid_change ( ewald_pw, ewald_env, error )
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_pw_grid_change', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ewald_type, o_spline, stat
    LOGICAL                                  :: failure
    REAL(dp)                                 :: alpha
    TYPE(cell_type), POINTER                 :: cell
    TYPE(dg_type), POINTER                   :: dg
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: pw_big_pool, pw_small_pool
    TYPE(section_vals_type), POINTER         :: poisson_section

  failure=.FALSE.
 
  NULLIFY ( pw_big_pool )
  NULLIFY ( pw_small_pool )
  NULLIFY ( cell )
  NULLIFY ( dg, poisson_env, poisson_section )

  CALL ewald_env_get ( ewald_env, ewald_type = ewald_type, &
                       alpha = alpha, o_spline = o_spline,&
                       poisson_section=poisson_section)
  SELECT CASE ( ewald_type )
  CASE ( do_ewald_ewald )
     CALL ewald_pw_get ( ewald_pw, cell = cell, pw_big_pool = pw_big_pool, &
                         dg = dg, poisson_env=poisson_env )
     CALL pw_grid_change ( cell, pw_big_pool % pw_grid )
     CALL ewald_pw_rho0_setup ( ewald_env, pw_big_pool % pw_grid, dg )
     CALL pw_poisson_release(poisson_env,error=error)
     CALL ewald_pw_set ( ewald_pw, pw_big_pool = pw_big_pool, dg = dg, &
          poisson_env=poisson_env )
  CASE ( do_ewald_pme )
     CALL ewald_pw_get ( ewald_pw, cell = cell, pw_big_pool = pw_big_pool, &
                         pw_small_pool = pw_small_pool, dg = dg, &
                         poisson_env=poisson_env )
     IF (.NOT.ASSOCIATED(poisson_env)) THEN
        CALL pw_poisson_create(poisson_env,error=error)
        CALL ewald_pw_set(ewald_pw, poisson_env=poisson_env)
     END IF
     CALL pw_grid_change ( cell, pw_big_pool % pw_grid )
     CALL dg_grid_change ( cell, pw_big_pool % pw_grid, pw_small_pool % pw_grid )
     CALL ewald_pw_rho0_setup ( ewald_env, pw_small_pool % pw_grid, dg )
     CALL ewald_pw_set ( ewald_pw, pw_big_pool = pw_big_pool,  &
                         pw_small_pool = pw_small_pool, dg = dg, &
                         poisson_env=poisson_env )
  CASE ( do_ewald_spme )
     CALL ewald_pw_get ( ewald_pw, cell = cell, pw_big_pool = pw_big_pool,&
          poisson_env=poisson_env)
     IF (.NOT.ASSOCIATED(poisson_env)) THEN
        CALL pw_poisson_create(poisson_env,error=error)
     END IF
     CALL pw_grid_change ( cell, pw_big_pool % pw_grid )
     CALL ewald_pw_set ( ewald_pw, pw_big_pool = pw_big_pool, &
          poisson_env=poisson_env )
  CASE ( do_ewald_none )
  CASE default
     CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
  END SELECT
  IF (ASSOCIATED(poisson_env)) THEN
     ALLOCATE(pw_pools(1),stat=stat)
     CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
     pw_pools(1)%pool => pw_big_pool
     CALL pw_poisson_set(poisson_env,cell=cell,parameters=poisson_section,&
          use_level=1,pw_pools=pw_pools,error=error)
   ! CJM  CALL pw_poisson_rebuild(poisson_env,error=error)
     DEALLOCATE(pw_pools,stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  END IF

END SUBROUTINE ewald_pw_grid_change 
!******************************************************************************
!!****** ewald_pw_methods/ewald_pw_rho0_setup [1.0] *
!!
!!   NAME
!!     ewald_pw_rho0_setup
!!
!!   FUNCTION
!!     Calculates the Fourier transform of the "Ewald function"
!!
!!   AUTHOR
!!     JGH (15-Mar-2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE ewald_pw_rho0_setup ( ewald_env, pw_grid, dg )
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(dg_type), POINTER                   :: dg

    INTEGER                                  :: ewald_type, grid_index, isos
    REAL(dp)                                 :: alpha
    REAL(dp), POINTER                        :: gcc( : ), zet( : )
    TYPE(dg_rho0_type), POINTER              :: dg_rho0

  CALL ewald_env_get ( ewald_env, alpha = alpha, ewald_type = ewald_type )
  CALL dg_get ( dg, dg_rho0 = dg_rho0 )
  CALL dg_rho0_get ( dg_rho0, gcc = gcc, zet = zet )

! This is the first ( and only ) double grid
  grid_index = 1

  IF ( .NOT. ASSOCIATED ( zet ) ) THEN
    ALLOCATE ( zet ( 1 ) , STAT = isos )
    IF ( isos /= 0 ) CALL stop_memory ('ewald_pw_rho0_setup', 'zet', 0 )
  ENDIF

! No contracted Gaussians are used here
  NULLIFY ( gcc )

  zet ( 1 ) =  alpha
  CALL dg_rho0_set (  dg_rho0, TYPE = ewald_type, zet = zet ) 

  CALL dg_rho0_init ( dg_rho0, pw_grid )
  
  CALL dg_set ( dg, dg_rho0 = dg_rho0, grid_index = grid_index ) 

END SUBROUTINE ewald_pw_rho0_setup

!******************************************************************************

END MODULE ewald_pw_methods

!******************************************************************************
