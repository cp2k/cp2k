!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001 - 2002  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****s* cp2k/xc [1.0] *
!!
!!   NAME
!!     xc
!!
!!   FUNCTION
!!     Exchange and Correlation functional calculations
!!
!!   AUTHOR
!!     JGH (13-Feb-2001); based on earlier version of apsi
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE xc
  
  USE kinds, ONLY: dbl
  USE termination, ONLY : stop_program
  USE pw_grid_types, ONLY : pw_grid_type
  USE pw_grids, ONLY : pw_grid_compare
  USE pw_types, ONLY : pw_type
  
  PRIVATE
  PUBLIC :: xc_calculate

  INTERFACE xc_calculate
    MODULE PROCEDURE xc_calculate_pw_lda
  END INTERFACE

!!*****
!******************************************************************************
  
CONTAINS

!******************************************************************************
!!****** xc/xc_calculate_pw_lda [1.0] *
!!
!!   NAME
!!     xc_calculate_pw_lda
!!
!!   FUNCTION
!!     Exchange and Correlation functional calculations
!!
!!   AUTHOR
!!     JGH (13-Feb-2001); based on earlier version of apsi
!!
!!   MODIFICATION HISTORY
!!     JGH (13-Jun-2002): adaptation to new functionals
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE xc_calculate_pw_lda ( exc, vxc, xc_pot, xc_stress )
  
  IMPLICIT NONE

! Arguments
  REAL ( dbl ), INTENT ( OUT ) :: exc
  REAL ( dbl ), INTENT ( OUT ) :: vxc
  TYPE ( pw_type ), INTENT ( OUT ), TARGET :: xc_pot
  REAL ( dbl ), DIMENSION ( 3, 3 ), INTENT ( OUT ), OPTIONAL :: xc_stress

! Local
  TYPE ( pw_grid_type ), POINTER :: grid_dens, grid_pot, grid_xc
  TYPE ( pw_type ), POINTER :: pw_dens_a, pw_xcpot_a, pw_dxcpot_a, &
                               pw_dens_b, pw_xcpot_b, pw_dxcpot_b
  REAL ( dbl ) :: pw_exc, pw_vxc
  REAL ( dbl ), DIMENSION ( 3, 3 ) :: xc_stress_local

!------------------------------------------------------------------------------

END SUBROUTINE xc_calculate_pw_lda

!!*****
!******************************************************************************

END MODULE xc

!******************************************************************************
