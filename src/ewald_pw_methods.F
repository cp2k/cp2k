!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C ) 2002 - 2003  CP2K developers group                         !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
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
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
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
  USE green_methods,                   ONLY: pw_green_fn_rebuild
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_change
  USE pw_pool_types,                   ONLY: pw_pool_type
  USE simulation_cell,                 ONLY: cell_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE green_types,                     ONLY: greens_function_type
  USE input_constants,                 ONLY: do_ewald_none,&
                                             do_ewald_ewald,&
                                             do_ewald_pme,&
                                             do_ewald_spme  

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='ewald_pw_methods'
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

SUBROUTINE ewald_pw_grid_change ( ewald_pw, ewald_env )
  IMPLICIT NONE
  ! Arguments
  TYPE(ewald_pw_type), POINTER             :: ewald_pw
  TYPE(ewald_environment_type), POINTER    :: ewald_env
  ! Local Variables
  INTEGER                                  :: ewald_type
  INTEGER                                  :: o_spline
  REAL(dp)                                 :: alpha
  TYPE(greens_function_type), POINTER      :: green
  TYPE(cell_type), POINTER                 :: cell
  TYPE(dg_type), POINTER                   :: dg
  TYPE(pw_pool_type), POINTER              :: pw_big_pool, pw_small_pool
  
  NULLIFY ( pw_big_pool )
  NULLIFY ( pw_small_pool )
  NULLIFY ( cell, green )
  NULLIFY ( dg )

  CALL ewald_env_get ( ewald_env, ewald_type = ewald_type, &
                       alpha = alpha, o_spline = o_spline )
  SELECT CASE ( ewald_type )
  CASE ( do_ewald_ewald )
     CALL ewald_pw_get ( ewald_pw, cell = cell, pw_big_pool = pw_big_pool, &
                         dg = dg )
     CALL pw_grid_change ( cell, pw_big_pool % pw_grid )
     CALL ewald_pw_rho0_setup ( ewald_env, pw_big_pool % pw_grid, dg )
     CALL ewald_pw_set ( ewald_pw, pw_big_pool = pw_big_pool, dg = dg )
  CASE ( do_ewald_pme )
     CALL ewald_pw_get ( ewald_pw, cell = cell, green = green, pw_big_pool = pw_big_pool, &
                         pw_small_pool = pw_small_pool, dg = dg )
     CALL pw_grid_change ( cell, pw_big_pool % pw_grid )
     CALL pw_green_fn_rebuild ( green, pw_big_pool )
     CALL dg_grid_change ( cell, pw_big_pool % pw_grid, pw_small_pool % pw_grid )
     CALL ewald_pw_rho0_setup ( ewald_env, pw_small_pool % pw_grid, dg )
     CALL ewald_pw_set ( ewald_pw, pw_big_pool = pw_big_pool,  &
                         pw_small_pool = pw_small_pool, dg = dg )
  CASE ( do_ewald_spme )
     CALL ewald_pw_get ( ewald_pw, cell = cell, green = green, pw_big_pool = pw_big_pool )
     CALL pw_grid_change ( cell, pw_big_pool % pw_grid )
     CALL pw_green_fn_rebuild ( green, pw_big_pool, p3m=o_spline, &
                             alpha = alpha )
     CALL ewald_pw_set ( ewald_pw, pw_big_pool = pw_big_pool )
  END SELECT

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
  IMPLICIT NONE
  ! Arguments
  TYPE(ewald_environment_type), POINTER    :: ewald_env
  TYPE(pw_grid_type), POINTER              :: pw_grid
  TYPE(dg_type), POINTER                   :: dg
  ! Local Variables
  INTEGER                                  :: ewald_type
  INTEGER                                  :: grid_index, isos
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
  CALL dg_rho0_set (  dg_rho0, type = ewald_type, zet = zet ) 

  CALL dg_rho0_init ( dg_rho0, pw_grid )
  
  CALL dg_set ( dg, dg_rho0 = dg_rho0, grid_index = grid_index ) 

END SUBROUTINE ewald_pw_rho0_setup

!******************************************************************************

END MODULE ewald_pw_methods

!******************************************************************************
