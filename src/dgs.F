!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE dgs
  
  USE coefficient_types, ONLY : coeff_type
  USE fft_tools, ONLY : fft_radix_operations, &
       FFT_RADIX_ALLOWED, FFT_RADIX_DISALLOWED, BWFFT, FWFFT, &
       fft3d => fft_wrap
  USE kinds, ONLY : dbl, sgl
  USE mathconstants, ONLY : twopi, pi
  USE pw_grid_types, ONLY : pw_grid_type
  USE pw_grids, ONLY : pw_grid_setup, pw_find_cutoff
  USE pw_types, ONLY : COMPLEXDATA3D
  USE pws, ONLY : fft_wrap
  USE simulation_cell, ONLY : cell_type, get_hinv, get_cell_param
  USE stop_program, ONLY : stop_prg
  USE structure_factors, ONLY : structure_factor_evaluate
  USE util, ONLY : dotprod_3d
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: dg_grid_setup, dg_find_max_cutoff, dg_sum_patch, &
       dg_sum_patch_force, dg_get_1x_patch, dg_get_2x_patch, &
       dg_get_strucfac, dg_get_rho_tot
  
CONTAINS

!******************************************************************************

SUBROUTINE dg_grid_setup ( box_b, npts_s, epsilon, alpha, grid_s, grid_b, &
     dg_gaussian_type )
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, DIMENSION ( : ), INTENT ( IN ) :: npts_s
  REAL ( dbl ), INTENT ( IN ) :: epsilon, alpha
  TYPE ( cell_type ), INTENT ( IN ) :: box_b
  CHARACTER( LEN = * ), INTENT ( IN ) :: dg_gaussian_type
  
  TYPE ( pw_grid_type ), INTENT ( INOUT ) :: grid_s, grid_b
  
! Locals
  INTEGER :: check
  REAL ( dbl ) :: dr ( 3 ), cutoff
  REAL ( dbl ) :: cell_lengths ( 3 )
  TYPE ( cell_type ) :: unit_box, box_s
  
  INTEGER :: foo ( 3 ) !*apsi
!------------------------------------------------------------------------------
  
  CALL fft_radix_operations ( npts_s ( 1 ), check, &
       operation = FFT_RADIX_ALLOWED )
  IF ( check /= FFT_RADIX_ALLOWED ) THEN
     CALL stop_prg ( "dg_grid_setup", "disallowed small FFT length #1" )
  END IF
  CALL fft_radix_operations ( npts_s ( 2 ), check, &
       operation = FFT_RADIX_ALLOWED )
  IF ( check /= FFT_RADIX_ALLOWED ) THEN
     CALL stop_prg ( "dg_grid_setup", "disallowed small FFT length #2" )
  END IF
  CALL fft_radix_operations ( npts_s ( 3 ), check, &
       operation = FFT_RADIX_ALLOWED )
  IF ( check /= FFT_RADIX_ALLOWED ) THEN
     CALL stop_prg ( "dg_grid_setup", "disallowed small FFT length #3" )
  END IF

  CALL get_cell_param ( box_b, cell_lengths )
  
  CALL dg_get_spacing ( npts_s, epsilon, alpha, dg_gaussian_type, dr )
  
  CALL dg_find_radix ( dr, cell_lengths, grid_b % npts )
  
  grid_b % bounds ( 1, : ) = - grid_b % npts / 2
  grid_b % bounds ( 2, : ) = + ( grid_b % npts - 1 ) / 2
  grid_s % npts ( : ) = grid_s % bounds ( 2, : ) - grid_s % bounds ( 1, : ) + 1
  
  CALL pw_find_cutoff ( grid_b % npts, box_b, cutoff )
  
  CALL pw_grid_setup ( box_b, grid_b, cutoff )
  
  CALL dg_find_basis ( grid_b % npts, box_b, unit_box )
  
  CALL dg_set_cell ( grid_s % npts, unit_box, box_s )
  
  CALL pw_grid_setup ( box_s, grid_s )
  
END SUBROUTINE dg_grid_setup

!******************************************************************************

SUBROUTINE dg_get_spacing ( npts, epsilon, alpha, dg_gaussian_type, dr )
  
  IMPLICIT NONE
  
! Arguments
  REAL ( dbl ), INTENT ( IN ) :: epsilon, alpha
  INTEGER, DIMENSION ( : ), INTENT ( IN ) :: npts
  CHARACTER ( LEN = * ), INTENT ( IN ) :: dg_gaussian_type
  
  REAL ( dbl ), DIMENSION ( : ), INTENT ( OUT ) :: dr
  
! Locals
  REAL ( dbl ) :: alphasq, norm
  
!------------------------------------------------------------------------------
  
  write(6,*) "dg_get_spacing, please check n/2 vs REAL(n)/2.0 !!!"
  
  SELECT CASE ( dg_gaussian_type )
     
  CASE ( "PME_GAUSS" )  ! use a Gaussian of form N exp(-2*alpha^2*r^2)
     alphasq = alpha ** 2
     norm = ( 2.0_dbl * alphasq / pi ) ** ( 1.5_dbl )
     dr ( : ) = SQRT ( -LOG ( epsilon / norm ) / alphasq / 2.0_dbl ) &
          / ( REAL ( npts ( : ), dbl ) / 2.0_dbl )
     
  CASE ( "S" )          ! use a Gaussian of form N exp(-2*alpha^2*r^2)
     alphasq = alpha ** 2
     norm = ( 2.0_dbl * alphasq / pi ) ** ( 1.5_dbl )
     dr ( : ) = SQRT ( -LOG ( epsilon / norm ) / alphasq / 2.0_dbl ) &
          / ( REAL ( npts ( : ), dbl ) / 2.0_dbl )
     
  CASE ( "P" )
     CALL stop_prg ( "dg_get_spacing", "'P; gaussian type not defined" )
     
  CASE DEFAULT
     CALL stop_prg ( "dg_get_spacing", "no suitable gaussian type specified" )
     
  END SELECT
  
END SUBROUTINE dg_get_spacing

!******************************************************************************

SUBROUTINE dg_find_radix ( dr, cell_lengths, npts )
  
  IMPLICIT NONE
  
! Arguments
  REAL ( dbl ), INTENT ( IN ) :: dr ( 3 )
  REAL ( dbl ), INTENT ( IN ) :: cell_lengths ( 3 )
  INTEGER, DIMENSION ( : ), INTENT ( OUT ) :: npts
  
!------------------------------------------------------------------------------
  
  npts ( : ) = NINT ( cell_lengths ( : ) / dr ( : ) )
  
END SUBROUTINE dg_find_radix

!******************************************************************************

SUBROUTINE dg_find_basis( npts, box, unit_box )
  IMPLICIT NONE

! Arguments
  INTEGER, DIMENSION ( : ), INTENT ( IN ) :: npts
  TYPE ( cell_type ), INTENT ( IN ) :: box
  TYPE ( cell_type ), INTENT ( OUT )  :: unit_box
 
! Locals
  INTEGER :: i
!------------------------------------------------------------------------------

  DO i = 1, 3
     unit_box % hmat ( :, i ) = box % hmat ( :, i ) / REAL ( npts ( : ), dbl )
  END DO

  CALL get_hinv ( unit_box )

END SUBROUTINE dg_find_basis

!******************************************************************************

!! Calculation of the basis on the mesh 'box'

SUBROUTINE dg_set_cell (npts, unit_box, box )
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, DIMENSION ( : ), INTENT ( IN ) :: npts
  TYPE ( cell_type ), INTENT ( OUT ) :: box
  TYPE ( cell_type ), INTENT ( IN ) :: unit_box
  
!------------------------------------------------------------------------------
! computing the unit vector along a, b, c and scaling it to length dr:

  
  box % hmat ( :, 1 ) = unit_box % hmat ( :, 1 ) * npts ( 1 )
  box % hmat ( :, 2 ) = unit_box % hmat ( :, 2 ) * npts ( 2 )
  box % hmat ( :, 3 ) = unit_box % hmat ( :, 3 ) * npts ( 3 )
  
  CALL get_hinv ( box )
  
END SUBROUTINE dg_set_cell

!******************************************************************************

SUBROUTINE dg_find_max_cutoff ( npts, grid_b )
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, DIMENSION ( : ), INTENT ( IN ) :: npts
  TYPE ( pw_grid_type ), DIMENSION ( : ), INTENT ( IN ) :: grid_b
  
!------------------------------------------------------------------------------
  
  CALL stop_prg ( "dg_find_max_cutoff", "not yet implemented" )
  
END SUBROUTINE dg_find_max_cutoff

!******************************************************************************

SUBROUTINE dg_get_strucfac ( box, r, npts_s, npts_b, centre, lb, ex, ey, ez )
  IMPLICIT NONE
  
! Arguments
  TYPE ( cell_type ), INTENT ( IN ) :: box
  REAL ( dbl ), DIMENSION ( : ), INTENT ( IN ) :: r
  INTEGER, DIMENSION ( : ), INTENT ( IN ) :: npts_s, npts_b
  INTEGER, INTENT ( IN ) :: lb ( 3 )
  
  INTEGER, INTENT ( OUT ) :: centre ( 3 )
  COMPLEX ( dbl ), DIMENSION ( lb ( 1 ): ), INTENT ( OUT ) :: ex
  COMPLEX ( dbl ), DIMENSION ( lb ( 2 ): ), INTENT ( OUT ) :: ey
  COMPLEX ( dbl ), DIMENSION ( lb ( 3 ): ), INTENT ( OUT ) :: ez
  
! Locals
  REAL ( dbl ) :: delta ( 3 )
  
!------------------------------------------------------------------------------
  
  CALL get_delta ( box, r, npts_s, npts_b, centre, delta)
  
  CALL structure_factor_evaluate ( delta, npts_s, lb, ex, ey, ez )
  
END SUBROUTINE dg_get_strucfac

!******************************************************************************

SUBROUTINE get_delta ( box, r, npts_s, npts_b, centre, delta )
  IMPLICIT NONE
  
! Arguments
  TYPE ( cell_type ), INTENT ( IN ) :: box
  REAL ( dbl ), DIMENSION ( : ), INTENT ( IN ) :: r
  INTEGER, DIMENSION ( : ), INTENT ( IN ) :: npts_s, npts_b
  
  INTEGER, DIMENSION ( : ), INTENT ( OUT ) :: centre
  REAL ( dbl ), DIMENSION ( : ), INTENT ( OUT ) :: delta
  
! Locals
  REAL ( dbl ), DIMENSION ( 3 ) :: s, grid_i
  
!------------------------------------------------------------------------------
  
! compute the scaled coordinate of atomi
  s = MATMUL ( box % h_inv, r )
  s = s - NINT ( s )
  
! find the continuous ``grid'' point (on big grid)
  grid_i ( 1:3 ) = REAL ( npts_b ( 1:3 ), dbl ) * s ( 1:3 )
  
! find the closest grid point (on big grid)
  centre ( : ) = NINT ( grid_i ( : ) )
  
! find the distance vector
  delta ( : ) = ( grid_i ( : ) - centre ( : ) ) / REAL ( npts_s ( : ), dbl )
  
END SUBROUTINE get_delta

!******************************************************************************

SUBROUTINE dg_get_1x_patch ( rho0, rhos_r, charge1, &
     lb, ex1, ey1, ez1 )
  
  IMPLICIT NONE
  
! Arguments
  REAL ( dbl ), INTENT ( IN ) :: charge1
  TYPE ( coeff_type ), INTENT ( IN ) :: rho0
  INTEGER, INTENT ( IN ) :: lb ( 3 )
  COMPLEX ( dbl ), DIMENSION ( lb(1): ), INTENT ( IN ) :: ex1
  COMPLEX ( dbl ), DIMENSION ( lb(2): ), INTENT ( IN ) :: ey1
  COMPLEX ( dbl ), DIMENSION ( lb(3): ), INTENT ( IN ) :: ez1
  
  TYPE ( coeff_type ), INTENT ( INOUT ), TARGET :: rhos_r
  
! Locals
  INTEGER :: i1, i2, i3
  COMPLEX ( dbl ) :: temp1a, temp2a, temp3a
  INTEGER, DIMENSION ( :, : ), POINTER :: bds
  INTEGER, DIMENSION ( : ), POINTER :: npts
  REAL ( dbl ), DIMENSION ( :, :, : ), POINTER :: rho
  COMPLEX ( dbl ), DIMENSION ( :, :, : ), POINTER :: rhos
  
!------------------------------------------------------------------------------
  
  rhos => rhos_r % pw % cc3d
  rho => rho0 % pw % cr3d
  bds => rhos_r % pw % pw_grid % bounds
  npts => rhos_r % pw % pw_grid % npts
  
  DO i3 = bds ( 1, 3 ), bds ( 2, 3 )
     temp3a = ez1 ( i3 ) * charge1
     
     DO i2 = bds ( 1, 2 ), bds ( 2, 2 )
        temp2a = temp3a * ey1 ( i2 )
        
        DO i1 = bds ( 1, 1 ), bds ( 2, 1 )
           temp1a = temp2a * ex1 ( i1 )
           
           rhos ( i1, i2, i3 ) = temp1a * rho ( i1, i2, i3 )
        END DO
     END DO
  END DO
  
! compute rho_s(r) for the current atoms
  CALL fft3d ( BWFFT, npts, rhos )
  
END SUBROUTINE dg_get_1x_patch

!******************************************************************************

!! The output comes in the real and imag parts of 'rhos_r'

SUBROUTINE dg_get_2x_patch ( rho0, rhos_r, charge1, charge2, &
     lb, ex1, ey1, ez1, ex2, ey2, ez2 )
  
  IMPLICIT NONE
  
! Arguments
  REAL ( dbl ), INTENT ( IN ) :: charge1, charge2
  TYPE ( coeff_type ), INTENT ( IN ), TARGET :: rho0
  INTEGER, INTENT ( IN ) :: lb ( 3 )
  COMPLEX ( dbl ), DIMENSION ( lb(1): ), INTENT ( IN ) :: ex1, ex2
  COMPLEX ( dbl ), DIMENSION ( lb(2): ), INTENT ( IN ) :: ey1, ey2
  COMPLEX ( dbl ), DIMENSION ( lb(3): ), INTENT ( IN ) :: ez1, ez2
  
  TYPE ( coeff_type ), INTENT ( INOUT ), TARGET :: rhos_r
  
! Locals
  INTEGER :: i1, i2, i3
  COMPLEX ( dbl ) :: temp1a, temp2a, temp3a, temp1b, temp2b, temp3b
  COMPLEX ( dbl ) :: ci
  INTEGER, DIMENSION ( :, : ), POINTER :: bds_s
  INTEGER, DIMENSION ( : ), POINTER :: npts_s
  REAL ( dbl ), DIMENSION ( :, :, : ), POINTER :: rho
  COMPLEX ( dbl ), DIMENSION ( :, :, : ), POINTER :: rhos
  
!------------------------------------------------------------------------------
  
  rhos => rhos_r % pw % cc3d
  rho => rho0 % pw % cr3d
  bds_s => rhos_r % pw % pw_grid % bounds
  npts_s => rhos_r % pw % pw_grid % npts
  
  ci = CMPLX ( 0.0_dbl, 1.0_dbl )
  
  DO i3 = bds_s ( 1, 3 ), bds_s ( 2, 3 )
     temp3a = ez1 ( i3 ) * charge1
     temp3b = ez2 ( i3 ) * charge2
     
     DO i2 = bds_s ( 1, 2 ), bds_s ( 2, 2 )
        temp2a = temp3a * ey1 ( i2 )
        temp2b = temp3b * ey2 ( i2 )
        
        DO i1 = bds_s ( 1, 1 ), bds_s ( 2, 1 )
           temp1a = temp2a * ex1 ( i1 )
           temp1b = temp2b * ex2 ( i1 )
           
           rhos ( i1, i2, i3 ) = ( temp1a + ci * temp1b ) * rho ( i1, i2, i3 )
        END DO
     END DO
  END DO
  
! compute rho_s(r) for the current atoms
  CALL fft3d ( BWFFT, npts_s, rhos )
  
END SUBROUTINE dg_get_2x_patch

!******************************************************************************

SUBROUTINE dg_sum_patch ( rhos, rhob, bounds_s, centre, &
     lb_pbc, pbc1_b, pbc2_b, pbc3_b, lb_rhob )
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, DIMENSION ( 2, 3 ), INTENT ( IN ) :: bounds_s
  REAL ( dbl ), DIMENSION ( bounds_s(1,1):bounds_s(2,1), &
       bounds_s(1,2):bounds_s(2,2), bounds_s(1,3):bounds_s(2,3) ), &
       INTENT ( IN ) :: rhos
  INTEGER, DIMENSION ( 3 ), INTENT ( IN ) :: centre
  INTEGER, DIMENSION ( 3 ), INTENT ( IN ) :: lb_pbc, lb_rhob
  INTEGER, DIMENSION ( lb_pbc ( 1 ) : ), INTENT ( IN ) :: pbc1_b
  INTEGER, DIMENSION ( lb_pbc ( 2 ) : ), INTENT ( IN ) :: pbc2_b
  INTEGER, DIMENSION ( lb_pbc ( 3 ) : ), INTENT ( IN ) :: pbc3_b
  
  REAL ( dbl ), DIMENSION ( lb_rhob(1):, lb_rhob(2):, lb_rhob(3): ), &
       INTENT ( INOUT ) :: rhob
  
! Locals
  INTEGER :: i1, i2, i3, ii1, ii2, ii3
  INTEGER :: lb1, lb2, lb3, ub1, ub2, ub3
  
!------------------------------------------------------------------------------
  
  lb1 = pbc1_b ( bounds_s ( 1, 1 ) + centre ( 1 ) ) - bounds_s ( 1, 1 )
  lb2 = pbc2_b ( bounds_s ( 1, 2 ) + centre ( 2 ) ) - bounds_s ( 1, 2 )
  lb3 = pbc3_b ( bounds_s ( 1, 3 ) + centre ( 3 ) ) - bounds_s ( 1, 3 )
  ub1 = pbc1_b ( bounds_s ( 2, 1 ) + centre ( 1 ) ) - bounds_s ( 1, 1 )
  ub2 = pbc2_b ( bounds_s ( 2, 2 ) + centre ( 2 ) ) - bounds_s ( 1, 2 )
  ub3 = pbc3_b ( bounds_s ( 2, 3 ) + centre ( 3 ) ) - bounds_s ( 1, 3 )
  
  IF ( lb3 < ub3 .AND. lb2 < ub2 .AND. lb1 < ub1 ) THEN
     DO i3 = bounds_s ( 1, 3 ), bounds_s ( 2, 3 )
        DO i2 = bounds_s ( 1, 2 ), bounds_s ( 2, 2 )
           DO i1 = bounds_s ( 1, 1 ), bounds_s ( 2, 1 )
              rhob ( lb1+i1, lb2+i2, lb3+i3 ) &
                   = rhob ( lb1+i1, lb2+i2, lb3+i3 ) &
                   + rhos ( i1, i2, i3 )
           END DO
        END DO
     END DO
  ELSE IF ( lb2 < ub2 .AND. lb1 < ub1 ) THEN
     DO i3 = bounds_s ( 1, 3 ), bounds_s ( 2, 3 )
        ii3 = pbc3_b ( i3 + centre ( 3 ) )
        DO i2 = bounds_s ( 1, 2 ), bounds_s ( 2, 2 )
           DO i1 = bounds_s ( 1, 1 ), bounds_s ( 2, 1 )
              rhob ( lb1+i1, lb2+i2, ii3 ) = rhob ( lb1+i1, lb2+i2, ii3 ) &
                   + rhos ( i1, i2, i3 )
           END DO
        END DO
     END DO
  ELSE IF ( lb3 < ub3 .AND. lb1 < ub1 ) THEN
     DO i3 = bounds_s ( 1, 3 ), bounds_s ( 2, 3 )
        DO i2 = bounds_s ( 1, 2 ), bounds_s ( 2, 2 )
           ii2 = pbc2_b ( i2 + centre ( 2 ) )
           DO i1 = bounds_s ( 1, 1 ), bounds_s ( 2, 1 )
              rhob ( lb1+i1, ii2, lb3+i3 ) = rhob ( lb1+i1, ii2, lb3+i3 ) &
                   + rhos ( i1, i2, i3 )
           END DO
        END DO
     END DO
  ELSE IF ( lb3 < ub3 .AND. lb2 < ub2 ) THEN
     DO i3 = bounds_s ( 1, 3 ), bounds_s ( 2, 3 )
        DO i2 = bounds_s ( 1, 2 ), bounds_s ( 2, 2 )
           DO i1 = bounds_s ( 1, 1 ), bounds_s ( 2, 1 )
              ii1 = pbc1_b ( i1 + centre ( 1 ) )
              rhob ( ii1, lb2+i2, lb3+i3 ) = rhob ( ii1, lb2+i2, lb3+i3 ) &
                   + rhos ( i1, i2, i3 )
           END DO
        END DO
     END DO
  ELSE IF ( lb1 < ub1 ) THEN
     DO i3 = bounds_s ( 1, 3 ), bounds_s ( 2, 3 )
        ii3 = pbc3_b ( i3 + centre ( 3 ) )
        DO i2 = bounds_s ( 1, 2 ), bounds_s ( 2, 2 )
           ii2 = pbc2_b ( i2 + centre ( 2 ) )
           DO i1 = bounds_s ( 1, 1 ), bounds_s ( 2, 1 )
              rhob ( lb1+i1, ii2, ii3 ) = rhob ( lb1+i1, ii2, ii3 ) &
                   + rhos ( i1, i2, i3 )
           END DO
        END DO
     END DO
  ELSE IF ( lb2 < ub2 ) THEN
     DO i3 = bounds_s ( 1, 3 ), bounds_s ( 2, 3 )
        ii3 = pbc3_b ( i3 + centre ( 3 ) )
        DO i2 = bounds_s ( 1, 2 ), bounds_s ( 2, 2 )
           DO i1 = bounds_s ( 1, 1 ), bounds_s ( 2, 1 )
              ii1 = pbc1_b ( i1 + centre ( 1 ) )
              rhob ( ii1, lb2+i2, ii3 ) = rhob ( ii1, lb2+i2, ii3 ) &
                   + rhos ( i1, i2, i3 )
           END DO
        END DO
     END DO
  ELSE IF ( lb1 < ub1 ) THEN
     DO i3 = bounds_s ( 1, 3 ), bounds_s ( 2, 3 )
        DO i2 = bounds_s ( 1, 2 ), bounds_s ( 2, 2 )
           ii2 = pbc2_b ( i2 + centre ( 2 ) )
           DO i1 = bounds_s ( 1, 1 ), bounds_s ( 2, 1 )
              ii1 = pbc1_b ( i1 + centre ( 1 ) )
              rhob ( ii1, ii2, lb3+i3 ) = rhob ( ii1, ii2, lb3+i3 ) &
                   + rhos ( i1, i2, i3 )
           END DO
        END DO
     END DO
  ELSE
     DO i3 = bounds_s ( 1, 3 ), bounds_s ( 2, 3 )
        ii3 = pbc3_b ( i3 + centre ( 3 ) )
        DO i2 = bounds_s ( 1, 2 ), bounds_s ( 2, 2 )
           ii2 = pbc2_b ( i2 + centre ( 2 ) )
           DO i1 = bounds_s ( 1, 1 ), bounds_s ( 2, 1 )
              ii1 = pbc1_b ( i1 + centre ( 1 ) )
              rhob ( ii1, ii2, ii3 ) = rhob ( ii1, ii2, ii3 ) &
                   + rhos ( i1, i2, i3 )
           END DO
        END DO
     END DO
  END IF
  
END SUBROUTINE dg_sum_patch

!******************************************************************************

SUBROUTINE dg_sum_patch_force ( rhos, rhob, f_xyz, bounds_s, centre, &
     lb_pbc, pbc1_b, pbc2_b, pbc3_b, lb_rhob )
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, DIMENSION ( 2, 3 ), INTENT ( IN ) :: bounds_s
  REAL ( dbl ), DIMENSION ( bounds_s(1,1):bounds_s(2,1), &
       bounds_s(1,2):bounds_s(2,2), bounds_s(1,3):bounds_s(2,3) ), &
       INTENT ( IN ) :: rhos
  INTEGER, DIMENSION ( 3 ), INTENT ( IN ) :: lb_rhob
  REAL ( dbl ), &
       DIMENSION ( 1:, lb_rhob ( 1 ) :, lb_rhob ( 2 ) :, lb_rhob ( 3 ) : ), &
       INTENT ( IN ) :: rhob
  INTEGER, DIMENSION ( 3 ), INTENT ( IN ) :: centre
  INTEGER, DIMENSION ( 3 ), INTENT ( IN ) :: lb_pbc
  INTEGER, DIMENSION ( lb_pbc ( 1 ) : ), INTENT ( IN ) :: pbc1_b
  INTEGER, DIMENSION ( lb_pbc ( 2 ) : ), INTENT ( IN ) :: pbc2_b
  INTEGER, DIMENSION ( lb_pbc ( 3 ) : ), INTENT ( IN ) :: pbc3_b
  
  REAL ( dbl ), DIMENSION ( : ), INTENT ( OUT ) :: f_xyz
  
! Locals
  INTEGER :: i1, i2, i3, ii1, ii2, ii3, ic
  INTEGER :: lb1, lb2, lb3, ub1, ub2, ub3
  
!------------------------------------------------------------------------------
  
  f_xyz ( : ) = 0.0_dbl

  lb1 = pbc1_b ( bounds_s ( 1, 1 ) + centre ( 1 ) ) - bounds_s ( 1, 1 )
  lb2 = pbc2_b ( bounds_s ( 1, 2 ) + centre ( 2 ) ) - bounds_s ( 1, 2 )
  lb3 = pbc3_b ( bounds_s ( 1, 3 ) + centre ( 3 ) ) - bounds_s ( 1, 3 )
  ub1 = pbc1_b ( bounds_s ( 2, 1 ) + centre ( 1 ) ) - bounds_s ( 1, 1 )
  ub2 = pbc2_b ( bounds_s ( 2, 2 ) + centre ( 2 ) ) - bounds_s ( 1, 2 )
  ub3 = pbc3_b ( bounds_s ( 2, 3 ) + centre ( 3 ) ) - bounds_s ( 1, 3 )
  
  IF ( lb3 < ub3 .AND. lb2 < ub2 .AND. lb1 < ub1 ) THEN
     DO i3 = bounds_s ( 1, 3 ), bounds_s ( 2, 3 )
        DO i2 = bounds_s ( 1, 2 ), bounds_s ( 2, 2 )
           DO i1 = bounds_s ( 1, 1 ), bounds_s ( 2, 1 )
              DO ic = 1, 3
                 f_xyz ( ic ) = f_xyz ( ic ) &
                      + rhob ( ic, lb1+i1, lb2+i2, lb3+i3 ) &
                      * rhos ( i1, i2, i3 )
              END DO
           END DO
        END DO
     END DO
  ELSE IF ( lb2 < ub2 .AND. lb1 < ub1 ) THEN
     DO i3 = bounds_s ( 1, 3 ), bounds_s ( 2, 3 )
        ii3 = pbc3_b ( i3 + centre ( 3 ) )
        DO i2 = bounds_s ( 1, 2 ), bounds_s ( 2, 2 )
           DO i1 = bounds_s ( 1, 1 ), bounds_s ( 2, 1 )
                 DO ic = 1, 3
              f_xyz ( ic ) = f_xyz ( ic ) &
                   + rhob ( ic, lb1+i1, lb2+i2, ii3 ) * rhos ( i1, i2, i3 )
           END DO
        END DO
     END DO
     END DO
  ELSE IF ( lb3 < ub3 .AND. lb1 < ub1 ) THEN
     DO i3 = bounds_s ( 1, 3 ), bounds_s ( 2, 3 )
        DO i2 = bounds_s ( 1, 2 ), bounds_s ( 2, 2 )
           ii2 = pbc2_b ( i2 + centre ( 2 ) )
           DO i1 = bounds_s ( 1, 1 ), bounds_s ( 2, 1 )
                 DO ic = 1, 3
              f_xyz ( ic ) = f_xyz ( ic ) &
                   + rhob ( ic, lb1+i1, ii2, lb3+i3 ) * rhos ( i1, i2, i3 )
           END DO
        END DO
     END DO
     END DO
  ELSE IF ( lb3 < ub3 .AND. lb2 < ub2 ) THEN
     DO i3 = bounds_s ( 1, 3 ), bounds_s ( 2, 3 )
        DO i2 = bounds_s ( 1, 2 ), bounds_s ( 2, 2 )
           DO i1 = bounds_s ( 1, 1 ), bounds_s ( 2, 1 )
              ii1 = pbc1_b ( i1 + centre ( 1 ) )
                 DO ic = 1, 3
              f_xyz ( ic ) = f_xyz ( ic ) &
                   + rhob ( ic, ii1, lb2+i2, lb3+i3 ) * rhos ( i1, i2, i3 )
           END DO
        END DO
     END DO
     END DO
  ELSE IF ( lb1 < ub1 ) THEN
     DO i3 = bounds_s ( 1, 3 ), bounds_s ( 2, 3 )
        ii3 = pbc3_b ( i3 + centre ( 3 ) )
        DO i2 = bounds_s ( 1, 2 ), bounds_s ( 2, 2 )
           ii2 = pbc2_b ( i2 + centre ( 2 ) )
           DO i1 = bounds_s ( 1, 1 ), bounds_s ( 2, 1 )
                 DO ic = 1, 3
              f_xyz ( ic ) = f_xyz ( ic ) &
                   + rhob ( ic, lb1+i1, ii2, ii3 ) * rhos ( i1, i2, i3 )
           END DO
        END DO
     END DO
     END DO
  ELSE IF ( lb2 < ub2 ) THEN
        DO i3 = bounds_s ( 1, 3 ), bounds_s ( 2, 3 )
           ii3 = pbc3_b ( i3 + centre ( 3 ) )
           DO i2 = bounds_s ( 1, 2 ), bounds_s ( 2, 2 )
              DO i1 = bounds_s ( 1, 1 ), bounds_s ( 2, 1 )
                 ii1 = pbc1_b ( i1 + centre ( 1 ) )
                 DO ic = 1, 3
                 f_xyz ( ic ) = f_xyz ( ic ) &
                      + rhob ( ic, ii1, lb2+i2, ii3 ) * rhos ( i1, i2, i3 )
              END DO
           END DO
        END DO
     END DO
  ELSE IF ( lb1 < ub1 ) THEN
     DO i3 = bounds_s ( 1, 3 ), bounds_s ( 2, 3 )
        DO i2 = bounds_s ( 1, 2 ), bounds_s ( 2, 2 )
           ii2 = pbc2_b ( i2 + centre ( 2 ) )
           DO i1 = bounds_s ( 1, 1 ), bounds_s ( 2, 1 )
              ii1 = pbc1_b ( i1 + centre ( 1 ) )
                 DO ic = 1, 3
              f_xyz ( ic ) = f_xyz ( ic ) &
                   + rhob ( ic, ii1, ii2, lb3+i3 ) * rhos ( i1, i2, i3 )
           END DO
        END DO
     END DO
     END DO
  ELSE
        DO i3 = bounds_s ( 1, 3 ), bounds_s ( 2, 3 )
           ii3 = pbc3_b ( i3 + centre ( 3 ) )
           DO i2 = bounds_s ( 1, 2 ), bounds_s ( 2, 2 )
              ii2 = pbc2_b ( i2 + centre ( 2 ) )
              DO i1 = bounds_s ( 1, 1 ), bounds_s ( 2, 1 )
                 ii1 = pbc1_b ( i1 + centre ( 1 ) )
                 DO ic = 1, 3
                 f_xyz ( ic ) = f_xyz ( ic ) &
                      + rhob ( ic, ii1, ii2, ii3 ) * rhos ( i1, i2, i3 )
              END DO
           END DO
        END DO
     END DO
  END IF
  
END SUBROUTINE dg_sum_patch_force

!******************************************************************************

SUBROUTINE dg_get_rho_tot ( rho_r, rhob_r )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( coeff_type ), INTENT ( INOUT ), TARGET :: rhob_r
  TYPE ( coeff_type ), INTENT ( INOUT ), TARGET :: rho_r
  
! Locals
  INTEGER :: lp1,mp1,np1,ln1,mn1,nn1
  INTEGER :: lp2,mp2,np2,ln2,mn2,nn2
  INTEGER :: gpt
  INTEGER, DIMENSION ( :, : ), POINTER :: bds1, bds2
  INTEGER, DIMENSION ( : ), POINTER :: ghat
  COMPLEX ( dbl ), DIMENSION ( :, :, : ), POINTER :: rho, rhob
  TYPE ( pw_grid_type ), POINTER :: grid1, grid2
  
!------------------------------------------------------------------------------
  
  grid1 => rhob_r % pw % pw_grid
  grid2 => rho_r % pw % pw_grid
  bds1 => grid1 % bounds
  bds2 => grid2 % bounds
  rhob => rhob_r % pw % cc3d
  rho => rho_r % pw % cc3d
  
  CALL fft_wrap ( rhob_r, dir = FWFFT, scale = grid1 % dvol )
  
  rho = 0.0_dbl
  
  DO gpt = 1, grid1 % ngpts_cut
    ghat => grid1 % g_hat ( :, gpt )
    
    lp1 = grid1 % mapl % pos ( ghat ( 1 ) )
    ln1 = grid1 % mapl % neg ( ghat ( 1 ) )
    mp1 = grid1 % mapm % pos ( ghat ( 2 ) )
    mn1 = grid1 % mapm % neg ( ghat ( 2 ) )
    np1 = grid1 % mapn % pos ( ghat ( 3 ) )
    nn1 = grid1 % mapn % neg ( ghat ( 3 ) )
    
    lp2 = grid2 % mapl % pos ( ghat ( 1 ) )
    ln2 = grid2 % mapl % neg ( ghat ( 1 ) )
    mp2 = grid2 % mapm % pos ( ghat ( 2 ) )
    mn2 = grid2 % mapm % neg ( ghat ( 2 ) )
    np2 = grid2 % mapn % pos ( ghat ( 3 ) )
    nn2 = grid2 % mapn % neg ( ghat ( 3 ) )
    
    !*apsi
    lp1 = lp1 + bds1 ( 1, 1 )
    mp1 = mp1 + bds1 ( 1, 2 )
    np1 = np1 + bds1 ( 1, 3 )
    ln1 = ln1 + bds1 ( 1, 1 )
    mn1 = mn1 + bds1 ( 1, 2 )
    nn1 = nn1 + bds1 ( 1, 3 )
    !*apsi
    lp2 = lp2 + bds1 ( 1, 1 )
    mp2 = mp2 + bds1 ( 1, 2 )
    np2 = np2 + bds1 ( 1, 3 )
    ln2 = ln2 + bds1 ( 1, 1 )
    mn2 = mn2 + bds1 ( 1, 2 )
    nn2 = nn2 + bds1 ( 1, 3 )
    
    rho ( lp2, mp2, np2 ) = rho ( lp2, mp2, np2 ) + rhob ( lp1, mp1, np1 )
    rho ( ln2, mn2, nn2 ) = rho ( ln2, mn2, nn2 ) + rhob ( ln1, mn1, nn1 )
    
  END DO
  
  CALL fft_wrap ( rho_r, dir = BWFFT )
  
END SUBROUTINE dg_get_rho_tot

!******************************************************************************

END MODULE dgs
