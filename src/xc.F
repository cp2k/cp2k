!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
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
  USE coefficient_types, ONLY : coeff_type, PLANEWAVES
  USE density_types, ONLY: density_type, D_NOT_USED, D_IN_USE, D_UP_TO_DATE
  USE potential_types, ONLY: potential_type, P_NOT_USED, P_IN_USE, P_UP_TO_DATE
  USE termination, ONLY : stop_program
  USE pw_grid_types, ONLY : pw_grid_type
  USE pw_grids, ONLY : pw_grid_compare
  USE pw_types, ONLY : pw_type
  USE pw_xc, ONLY : pw_calculate_xc
  USE functionals, ONLY : xc_type, init_functionals
  
  TYPE ( xc_type ) :: xc_info

  PRIVATE
  PUBLIC :: set_xc, calculate_xc

!!*****
!******************************************************************************
  
CONTAINS

!******************************************************************************
!!****** xc/set_xc [1.0] *
!!
!!   NAME
!!     set_xc
!!
!!   FUNCTION
!!     Initialize Exchange and Correlation functional calculations
!!
!!   AUTHOR
!!     JGH (14-Feb-2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE set_xc ()
  
  IMPLICIT NONE

! Arguments

!------------------------------------------------------------------------------

  CALL init_functionals()

  xc_info % x_functional = "SLATER"
  xc_info % c_functional = "VWN"
  xc_info % gradient_functional = .FALSE.
  xc_info % smooth_gradient = .FALSE.

END SUBROUTINE set_xc

!!*****
!******************************************************************************
!!****** xc/calculate_xc [1.0] *
!!
!!   NAME
!!     calculate_xc
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

SUBROUTINE calculate_xc ( density, exc, vxc, &
                          alpha_spin, beta_spin, &
                          xc_pot, dxc_pot, xc_stress )
  
  IMPLICIT NONE

! Arguments
  TYPE ( density_type ), INTENT ( IN ) :: density
  REAL ( dbl ), INTENT ( OUT ) :: exc
  REAL ( dbl ), OPTIONAL, INTENT ( OUT ) :: vxc
  INTEGER, OPTIONAL, INTENT ( IN ) :: alpha_spin, beta_spin
  TYPE ( potential_type ), INTENT ( OUT ), OPTIONAL, TARGET :: xc_pot
  TYPE ( potential_type ), INTENT ( OUT ), OPTIONAL, TARGET :: dxc_pot
  REAL ( dbl ), DIMENSION ( 3, 3 ), INTENT ( OUT ), OPTIONAL :: xc_stress

! Local
  LOGICAL :: spin_polarised, calc ( 3 )
  INTEGER :: aspin, bspin, part, den_type
  TYPE ( pw_grid_type ), POINTER :: grid_dens, grid_pot, grid_xc
  TYPE ( pw_type ), POINTER :: pw_dens_a, pw_xcpot_a, pw_dxcpot_a, &
                               pw_dens_b, pw_xcpot_b, pw_dxcpot_b
  REAL ( dbl ) :: pw_exc, pw_vxc
  REAL ( dbl ), DIMENSION ( 3, 3 ) :: xc_stress_local

!------------------------------------------------------------------------------

  calc = .FALSE.

  IF ( PRESENT ( alpha_spin ) ) THEN
    aspin = alpha_spin
  ELSE
    aspin = 1
  END IF

  IF ( PRESENT ( beta_spin ) ) THEN
    bspin = beta_spin
    spin_polarised = .TRUE.
  ELSE
    bspin = 0
    spin_polarised = .FALSE.
  END IF

! check the input
! density
  part = density % active
  IF ( part > density % npart .OR. part < 1 ) THEN
    CALL stop_program ( "calculate_xc", "Invalid density part" )
  END IF
  IF ( aspin > density % nspin .OR. aspin < 1 ) THEN
    CALL stop_program ( "calculate_xc", "Invalid alpha spin component of density" )
  END IF
  IF ( spin_polarised ) THEN
    IF ( bspin > density % nspin .OR. bspin < 1 ) THEN
      CALL stop_program ( "calculate_xc", "Invalid beta spin component of density" )
    END IF
  END IF
  IF ( density % status ( part ) /= D_NOT_USED ) THEN
    CALL stop_program ( "calculate_xc", "Requested density part not in use" )
  END IF
  IF ( density % status ( part ) /= D_UP_TO_DATE ) THEN
    CALL stop_program ( "calculate_xc", "Requested density part not updated" )
  END IF
! xc potential
  IF ( PRESENT ( xc_pot ) ) THEN
    part = xc_pot % active
    IF ( part > xc_pot % npart .OR. part < 1 ) THEN
      CALL stop_program ( "calculate_xc", "Invalid xc_pot part" )
    END IF
    IF ( aspin > xc_pot % nspin .OR. aspin < 1 ) THEN
      CALL stop_program ( "calculate_xc", "Invalid alpha spin component of xc_pot" )
    END IF
    IF ( spin_polarised ) THEN
      IF ( bspin > xc_pot % nspin .OR. bspin < 1 ) THEN
        CALL stop_program ( "calculate_xc", "Invalid beta spin component of xc_pot" )
      END IF
    END IF
    IF ( xc_pot % status ( part ) /= P_NOT_USED ) THEN
      CALL stop_program ( "calculate_xc", "Requested xc_pot part not in use" )
    END IF
    ! mark xc_pot to be calculated
    calc ( 1 ) = .TRUE.
  END IF
! dxc potential
  IF ( PRESENT ( dxc_pot ) ) THEN
    part = dxc_pot % active
    IF ( part > dxc_pot % npart .OR. part < 1 ) THEN
      CALL stop_program ( "calculate_xc", "Invalid dxc_pot part" )
    END IF
    IF ( aspin > dxc_pot % nspin .OR. aspin < 1 ) THEN
      CALL stop_program ( "calculate_xc", "Invalid alpha spin component of dxc_pot" )
    END IF
    IF ( spin_polarised ) THEN
      IF ( bspin > dxc_pot % nspin .OR. bspin < 1 ) THEN
        CALL stop_program ( "calculate_xc", "Invalid beta spin component of dxc_pot" )
      END IF
    END IF
    IF ( dxc_pot % status ( part ) /= P_NOT_USED ) THEN
      CALL stop_program ( "calculate_xc", "Requested dxc_pot part not in use" )
    END IF
    ! mark dxc_pot to be calculated
    calc ( 2 ) = .TRUE.
  END IF
  IF ( PRESENT ( xc_stress ) ) THEN
    ! mark stress to be calculated
    calc ( 3 ) = .TRUE.
  END IF

! check types
  part = density % active
  den_type = density % dpart ( aspin, part ) % in_use
  IF ( PRESENT ( xc_pot ) ) THEN
    part = xc_pot % active
    IF ( den_type /= xc_pot % pot_part ( aspin, part ) % in_use ) THEN
      CALL stop_program ( "calculate_xc", "Density and potential type do not match" )
    END IF
  END IF
  IF ( PRESENT ( dxc_pot ) ) THEN
    part = dxc_pot % active
    IF ( den_type /= dxc_pot % pot_part ( aspin, part ) % in_use ) THEN
      CALL stop_program ( "calculate_xc", "Density and potential type do not match" )
    END IF
  END IF

! now for the actual calculation

! we can have another grid for the XC calculation
!!!!!!  to be changed
  grid_xc => density % dpart ( aspin, density % active ) % pw % pw_grid

  IF ( den_type == PLANEWAVES ) THEN
    
! we have to see if the in_space mode of the plane wave variables is compatible
! we need a real 3d representation in real space for densities and potentials
! However, any other possibility is also allowed and we make the appropriate
! transformations on input or output of the xc routines
! The possible transformation of the density is done later as we might need
! both representations (real and g space) for the gradient calculation

    grid_dens => density % dpart ( aspin, density % active ) % pw % pw_grid
    IF ( PRESENT ( xc_pot ) ) THEN
      grid_pot => xc_pot % pot_part ( aspin, xc_pot % active ) % pw % pw_grid
      IF ( .NOT. pw_grid_compare ( grid_dens , grid_pot ) ) THEN
        CALL stop_program ( "calculate_xc", "PW grids do not match" )
      END IF
    END IF
    IF ( PRESENT ( dxc_pot ) ) THEN
      grid_pot => dxc_pot % pot_part ( aspin, xc_pot % active ) % pw % pw_grid
      IF ( .NOT. pw_grid_compare ( grid_dens , grid_pot ) ) THEN
        CALL stop_program ( "calculate_xc", "PW grids do not match" )
      END IF
    END IF

    IF ( .NOT. spin_polarised ) THEN

      pw_dens_a => density % dpart ( aspin, density % active ) % pw

      IF ( PRESENT ( xc_pot ) ) THEN
        pw_xcpot_a => xc_pot % pot_part ( aspin, xc_pot % active ) % pw
      ELSE
        NULLIFY ( pw_xcpot_a )
      END IF
      IF ( PRESENT ( xc_pot ) ) THEN
        pw_dxcpot_a => dxc_pot % pot_part ( aspin, dxc_pot % active ) % pw
      ELSE
        NULLIFY ( pw_dxcpot_a )
      END IF

      CALL pw_calculate_xc ( xc_info, calc, grid_xc, pw_dens_a, pw_exc, pw_vxc, &
                             pw_xcpot_a, pw_dxcpot_a, xc_stress_local )

    ELSE

      pw_dens_a => density % dpart ( aspin, density % active ) % pw
      pw_dens_b => density % dpart ( bspin, density % active ) % pw

      IF ( PRESENT ( xc_pot ) ) THEN
        pw_xcpot_a => xc_pot % pot_part ( aspin, xc_pot % active ) % pw
        pw_xcpot_b => xc_pot % pot_part ( bspin, xc_pot % active ) % pw
      ELSE
        NULLIFY ( pw_xcpot_a )
        NULLIFY ( pw_xcpot_b )
      END IF
      IF ( PRESENT ( xc_pot ) ) THEN
        pw_dxcpot_a => dxc_pot % pot_part ( aspin, dxc_pot % active ) % pw
        pw_dxcpot_b => dxc_pot % pot_part ( bspin, dxc_pot % active ) % pw
      ELSE
        NULLIFY ( pw_dxcpot_a )
        NULLIFY ( pw_dxcpot_b )
      END IF

      CALL pw_calculate_xc ( xc_info, calc, grid_xc, pw_dens_a, pw_dens_b, &
                             pw_exc, pw_vxc, pw_xcpot_a, pw_xcpot_b, &
                             pw_dxcpot_a, pw_dxcpot_b, xc_stress_local )

    END IF

    exc = pw_exc
    IF ( PRESENT ( vxc ) ) vxc = pw_vxc
    IF ( PRESENT ( xc_stress ) ) xc_stress = xc_stress_local

  ELSE

    CALL stop_program ( "calculate_xc", "Invalid density type" )

  END IF
  
END SUBROUTINE calculate_xc

!!*****
!******************************************************************************

END MODULE xc

!******************************************************************************
