!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001 - 2002 CP2K developers group                           !
!-----------------------------------------------------------------------------!
!!****** cp2k/pw_xc [1.0] *
!!
!!   NAME
!!     pw_xc
!!
!!   FUNCTION
!!     Exchange and Correlation functional calculations for plane waves
!!
!!   AUTHOR
!!     JGH (24-Feb-2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE pw_xc
  
  USE kinds, ONLY: dbl
  USE message_passing, ONLY : mp_sum
  USE termination, ONLY : stop_program, stop_memory
  USE pw_grid_types, ONLY : pw_grid_type
  USE pw_grids, ONLY : pw_grid_compare
  USE pw_types, ONLY : pw_type, pw_allocate, pw_deallocate, pw_transfer, &
            pw_zero, pw_copy, pw_derive, pw_sumup_square, pw_sumup, &
            REALDATA1D, COMPLEXDATA1D, REALDATA3D, COMPLEXDATA3D, &
            NOSPACE, REALSPACE, RECIPROCALSPACE
  USE xc_functionals, ONLY : xc_lda
  USE dft_types, ONLY : dft_control_type
  
  IMPLICIT NONE
  
  INTERFACE pw_calculate_xc
     MODULE PROCEDURE pw_calculate_xcns, pw_calculate_xcspin
  END INTERFACE

!!*****
!******************************************************************************
  
CONTAINS

!******************************************************************************
!!****** xc/pw_calculate_xc [1.0] *
!!
!!   NAME
!!     pw_calculate_xc
!!
!!   FUNCTION
!!     Exchange and Correlation functional calculations for plane wave basis
!!
!!   AUTHOR
!!     JGH (24-Feb-2001)
!!
!!   MODIFICATION HISTORY
!!     JGH (7-Mar-2001) : Added spin polarised case and stress calculation
!!                        pw_grid is now a variable and can be different from
!!                        the grid used for density and potential
!!
!!*****
!******************************************************************************

  SUBROUTINE pw_calculate_xcns ( dft_control, density, xc_pot, &
                                 exc, vxc, xc_stress )

! Arguments
  TYPE ( dft_control_type ), INTENT ( IN ) :: dft_control
  TYPE ( pw_type ), INTENT ( INOUT ) :: density, xc_pot
  REAL ( dbl ), INTENT ( OUT ) :: exc, vxc
  REAL ( dbl ), DIMENSION ( 3, 3 ), INTENT ( OUT ) :: xc_stress

! Local
  TYPE ( pw_grid_type ), POINTER :: xc_grid
  TYPE ( pw_type ) :: rho, drho, v_xc, v_xcg, dv, dvg, vg
  TYPE ( pw_type ) :: rhog, drhog, drho_i ( 3 )
  INTEGER :: i, j, n ( 3 )

!------------------------------------------------------------------------------

! ! allocate local arrays
! CALL pw_allocate ( rho, xc_grid, REALDATA3D )
! rho % in_space = REALSPACE
! CALL pw_allocate ( v_xc, xc_grid, REALDATA3D )
! v_xc % in_space = REALSPACE

! CALL pw_transfer ( density, rho )

! ! gradient corrections
! IF ( dft_control % gradient_functional ) THEN
!   ! The density in G space
!   CALL pw_allocate ( rhog, xc_grid, COMPLEXDATA1D )
!   rhog % in_space = RECIPROCALSPACE
!   CALL pw_transfer ( rho, rhog )

!   ! drho will hold the absolut value of the derivative of the density
!   CALL pw_allocate ( drho, xc_grid, REALDATA3D )
!   drho % in_space = REALSPACE
!   CALL pw_zero ( drho )
!   ! derivative of the density
!   CALL pw_allocate ( drhog, xc_grid, COMPLEXDATA1D )
!   drhog % in_space = RECIPROCALSPACE
!   DO i = 1, 3
!     CALL pw_copy ( rhog, drhog )
!     n = 0
!     n ( i ) = 1
!     CALL pw_derive ( drhog, n )
!     CALL pw_allocate ( drho_i ( i ), xc_grid, REALDATA3D )
!     drho_i ( i ) % in_space = REALSPACE
!     CALL pw_transfer ( drhog, drho_i ( i ) )
!     CALL pw_sumup_square ( drho, drho_i ( i ) )
!   END DO
!   CALL pw_deallocate ( rhog )
!   drho % cr3d = SQRT ( drho % cr3d )
!   CALL pw_allocate ( v_xcg, xc_grid, REALDATA3D )
!   v_xcg % in_space = REALSPACE
! ELSE

!   IF ( dft_control%integration_grid%pw_default ) THEN
!   ELSE
!   END IF

! END IF

! CALL xc_lda ( dft_control%functional, dft_control%ftags, rho, drho, pot, 1 )

!!!! CALL becke88 ( rho % cr3d, drho % cr3d, exc, v_xc % cr3d, v_xcg % cr3d )
!
! IF ( xc_info % gradient_functional ) THEN
!   CALL pw_deallocate ( drho )
!   CALL pw_deallocate ( drhog )
! END IF

! ! contribution to stress from gradient corrections
! IF ( tcalc ( 3 ) ) THEN
!   xc_stress = 0._dbl
!   IF ( xc_info % gradient_functional ) THEN
!     DO i = 1, 3
!       DO j = i, 3
!         xc_stress ( i, j ) = SUM ( v_xcg % cr3d * drho_i ( i ) % cr3d * &
!            drho_i ( j ) % cr3d ) * xc_grid % dvol
!         xc_stress ( j, i ) = xc_stress ( i, j )
!       END DO
!     END DO
!   END IF
! END IF

! ! potential is needed
! vxc = 0._dbl
! IF ( tcalc ( 1 ) ) THEN
!    IF ( xc_info % gradient_functional ) THEN
!      CALL pw_allocate ( dv, xc_grid, REALDATA3D )
!      dv % in_space = REALSPACE
!      CALL pw_allocate ( dvg, xc_grid, COMPLEXDATA1D )
!      dvg % in_space = RECIPROCALSPACE
!      CALL pw_allocate ( vg, xc_grid, COMPLEXDATA1D )
!      vg % in_space = RECIPROCALSPACE
!      CALL pw_zero ( vg )
!      DO i = 1, 3
!        dv % cr3d = drho_i ( i ) % cr3d * v_xcg % cr3d
!        CALL pw_transfer ( dv, dvg )
!        n = 0
!        n ( i ) = 1
!        CALL pw_derive ( dvg, n )
!        CALL pw_sumup ( dvg, vg )
!      END DO
!      CALL pw_transfer ( vg, dv )
!      CALL pw_sumup ( dv, v_xc )
!      CALL pw_deallocate ( dv )
!      CALL pw_deallocate ( dvg )
!      CALL pw_deallocate ( vg )
!    END IF
!    CALL pw_transfer ( v_xc, xc_pot )
!    vxc = SUM ( v_xc % cr3d * rho % cr3d )
! END IF

! ! deallocate all local arrays
! CALL pw_deallocate ( rho )
! CALL pw_deallocate ( v_xc )
! IF ( xc_info % gradient_functional ) THEN
!   CALL pw_deallocate ( v_xcg )
!   DO i = 1, 3
!     CALL pw_deallocate ( drho_i ( i ) )
!   END DO
! END IF

! exc = exc * xc_grid % dvol
! vxc = vxc * xc_grid % dvol

! ! general contribution to stress (isotropic)
! IF ( tcalc ( 3 ) ) THEN
!   IF ( .NOT. tcalc ( 1 ) ) CALL stop_program ( "pw_calculate_xcns", &
!     "stress calculation needs full potential calculation" )
!   DO i = 1, 3
!     xc_stress ( i, i ) = xc_stress ( i, i ) + ( exc - vxc )
!   END DO
! END IF

! IF ( xc_grid % para % mode == 1 ) THEN
!   CALL mp_sum ( exc, xc_grid % para % group ) 
!   CALL mp_sum ( vxc, xc_grid % para % group ) 
!   IF ( tcalc ( 3 ) ) CALL mp_sum ( xc_stress, xc_grid % para % group ) 
! END IF

  END SUBROUTINE pw_calculate_xcns

!------------------------------------------------------------------------------

  SUBROUTINE pw_calculate_xcspin ( )

! Arguments
! TYPE ( xc_type ), INTENT ( IN ) :: xc_info
! LOGICAL, DIMENSION ( 3 ), INTENT ( IN ) :: tcalc
! TYPE ( pw_grid_type ), INTENT ( IN ) :: xc_grid
! TYPE ( pw_type ), INTENT ( INOUT ) :: density_a, density_b, &
!                                 xc_pot_a, xc_pot_b, dxc_pot_a, dxc_pot_b
! REAL ( dbl ), INTENT ( OUT ) :: exc, vxc
! REAL ( dbl ), DIMENSION ( 3, 3 ), INTENT ( OUT ) :: xc_stress

! Local
! TYPE ( pw_type ) :: rhoa, rhob, drhoa, drhob, v_xca, v_xcb, &
!                     v_xcga, v_xcgb, v_xcgab, dv, dvg, vg
! TYPE ( pw_type ) :: rhoga, rhogb, drhoga, drhogb, drho_ia ( 3 ), drho_ib ( 3 )
! INTEGER :: i, j, n ( 3 )

!------------------------------------------------------------------------------

! ! allocate local arrays
! CALL pw_allocate ( rhoa, xc_grid, REALDATA3D )
! rhoa % in_space = REALSPACE
! CALL pw_allocate ( rhob, xc_grid, REALDATA3D )
! rhoa % in_space = REALSPACE
! CALL pw_allocate ( v_xca, xc_grid, REALDATA3D )
! v_xca % in_space = REALSPACE
! CALL pw_allocate ( v_xcb, xc_grid, REALDATA3D )
! v_xcb % in_space = REALSPACE

! CALL pw_transfer ( density_a, rhoa )
! CALL pw_transfer ( density_b, rhob )

! ! gradient corrections
! IF ( xc_info % gradient_functional ) THEN
!   ! The density in G space
!   CALL pw_allocate ( rhoga, xc_grid, COMPLEXDATA1D )
!   rhoga % in_space = RECIPROCALSPACE
!   CALL pw_transfer ( rhoa, rhoga )
!   CALL pw_allocate ( rhogb, xc_grid, COMPLEXDATA1D )
!   rhogb % in_space = RECIPROCALSPACE
!   CALL pw_transfer ( rhob, rhogb )

!   ! drho will hold the absolut value of the derivative of the density
!   CALL pw_allocate ( drhoa, xc_grid, REALDATA3D )
!   drhoa % in_space = REALSPACE
!   CALL pw_allocate ( drhob, xc_grid, REALDATA3D )
!   drhob % in_space = REALSPACE
!   CALL pw_zero ( drhoa )
!   CALL pw_zero ( drhob )
!   ! derivative of the density
!   CALL pw_allocate ( drhoga, xc_grid, COMPLEXDATA1D )
!   drhoga % in_space = RECIPROCALSPACE
!   CALL pw_allocate ( drhogb, xc_grid, COMPLEXDATA1D )
!   drhogb % in_space = RECIPROCALSPACE

!   DO i = 1, 3
!     n = 0
!     n ( i ) = 1

!     ! alpha spin
!     CALL pw_copy ( rhoga, drhoga )
!     CALL pw_derive ( drhoga, n )
!     CALL pw_allocate ( drho_ia ( i ), xc_grid, REALDATA3D )
!     drho_ia ( i ) % in_space = REALSPACE
!     CALL pw_transfer ( drhoga, drho_ia ( i ) )
!     CALL pw_sumup_square ( drhoa, drho_ia ( i ) )

!     ! beta spin
!     CALL pw_copy ( rhogb, drhogb )
!     CALL pw_derive ( drhogb, n )
!     CALL pw_allocate ( drho_ib ( i ), xc_grid, REALDATA3D )
!     drho_ib ( i ) % in_space = REALSPACE
!     CALL pw_transfer ( drhogb, drho_ib ( i ) )
!     CALL pw_sumup_square ( drhob, drho_ib ( i ) )
!   END DO

!   CALL pw_deallocate ( rhoga )
!   CALL pw_deallocate ( rhogb )
!   drhoa % cr3d = SQRT ( drhoa % cr3d )
!   drhob % cr3d = SQRT ( drhob % cr3d )
!   CALL pw_allocate ( v_xcga, xc_grid, REALDATA3D )
!   v_xcga % in_space = REALSPACE
!   CALL pw_allocate ( v_xcgb, xc_grid, REALDATA3D )
!   v_xcgb % in_space = REALSPACE
!   CALL pw_allocate ( v_xcgab, xc_grid, REALDATA3D )
!   v_xcgab % in_space = REALSPACE
!   CALL pw_zero ( v_xcgab )
! END IF

! ! second derivatives
! IF ( tcalc ( 2 ) ) THEN
!   CALL stop_program ( "pw_calculate_xc", "not yet available" )
! END IF

! ! exchange functionals
! SELECT CASE ( xc_info % x_functional )
! CASE DEFAULT
!    CALL stop_program ( "calculate_xc", &
!         "no valid exchange functional" // xc_info % x_functional )
! CASE ( "NONE" )
! CASE ( "BECKE88" )
!   CALL becke88 ( rhoa % cr3d ,rhob % cr3d ,drhoa % cr3d ,drhob % cr3d ,&
!           exc, v_xca % cr3d, v_xcb % cr3d, v_xcga % cr3d, v_xcgb % cr3d)
! END SELECT

! ! correlation functionals
! SELECT CASE ( xc_info % c_functional )
! CASE DEFAULT
!    CALL stop_program ( "calculate_xc", &
!         "no valid correlation functional" // xc_info % c_functional )
! CASE ( "NONE" )
! END SELECT

! IF ( xc_info % gradient_functional ) THEN
!   CALL pw_deallocate ( drhoa )
!   CALL pw_deallocate ( drhob )
!   CALL pw_deallocate ( drhoga )
!   CALL pw_deallocate ( drhogb )
! END IF

! ! stress
! IF ( tcalc ( 3 ) ) THEN
!   xc_stress = 0._dbl
!   IF ( xc_info % gradient_functional ) THEN
!     ! contribution to stress from gradient corrections
!     DO i = 1, 3
!       DO j = i, 3
!         xc_stress ( i, j ) = SUM ( v_xcga % cr3d * drho_ia ( i ) % cr3d * &
!            drho_ia ( j ) % cr3d + v_xcgb % cr3d * drho_ib ( i ) % cr3d * &
!            drho_ib ( j ) % cr3d + v_xcgab % cr3d * drho_ib ( i ) % cr3d * &
!            drho_ia ( j ) % cr3d + v_xcgab % cr3d * drho_ia ( i ) % cr3d * &
!            drho_ib ( j ) % cr3d ) * xc_grid % dvol
!         xc_stress ( j, i ) = xc_stress ( i, j ) 
!       END DO
!     END DO
!   END IF
! END IF

! ! potential is needed
! IF ( tcalc ( 1 ) ) THEN
!    IF ( xc_info % gradient_functional ) THEN
!      CALL pw_allocate ( dv, xc_grid, REALDATA3D )
!      dv % in_space = REALSPACE
!      CALL pw_allocate ( dvg, xc_grid, COMPLEXDATA1D )
!      dvg % in_space = RECIPROCALSPACE
!      CALL pw_allocate ( vg, xc_grid, COMPLEXDATA1D )
!      vg % in_space = RECIPROCALSPACE

!      ! alpha spin
!      CALL pw_zero ( vg )
!      DO i = 1, 3
!        dv % cr3d = drho_ia ( i ) % cr3d * v_xcga % cr3d + &
!                    drho_ib ( i ) % cr3d * v_xcgab % cr3d
!        CALL pw_transfer ( dv, dvg )
!        n = 0
!        n ( i ) = 1
!        CALL pw_derive ( dvg, n )
!        CALL pw_sumup ( dvg, vg )
!      END DO
!      CALL pw_transfer ( vg, dv )
!      CALL pw_sumup ( dv, v_xca )

!      ! beta spin
!      CALL pw_zero ( vg )
!      DO i = 1, 3
!        dv % cr3d = drho_ib ( i ) % cr3d * v_xcgb % cr3d + &
!                    drho_ia ( i ) % cr3d * v_xcgab % cr3d
!        CALL pw_transfer ( dv, dvg )
!        n = 0
!        n ( i ) = 1
!        CALL pw_derive ( dvg, n )
!        CALL pw_sumup ( dvg, vg )
!      END DO
!      CALL pw_transfer ( vg, dv )
!      CALL pw_sumup ( dv, v_xcb )

!      CALL pw_deallocate ( dv )
!      CALL pw_deallocate ( dvg )
!      CALL pw_deallocate ( vg )
!    END IF
!    vxc = SUM ( v_xca % cr3d * rhoa % cr3d + &
!                v_xcb % cr3d * rhob % cr3d )
!    CALL pw_transfer ( v_xca, xc_pot_a )
!    CALL pw_transfer ( v_xcb, xc_pot_b )
! END IF

! ! deallocate all local arrays
! CALL pw_deallocate ( rhoa )
! CALL pw_deallocate ( rhob )
! CALL pw_deallocate ( v_xca )
! CALL pw_deallocate ( v_xcb )
! IF ( xc_info % gradient_functional ) THEN
!   CALL pw_deallocate ( v_xcga )
!   CALL pw_deallocate ( v_xcgb )
!   CALL pw_deallocate ( v_xcgab )
!   DO i = 1, 3
!     CALL pw_deallocate ( drho_ia ( i ) )
!     CALL pw_deallocate ( drho_ib ( i ) )
!   END DO
! END IF

! exc = exc * xc_grid % dvol
! vxc = vxc * xc_grid % dvol

! ! general contribution to stress (isotropic)
! IF ( tcalc ( 3 ) ) THEN
!   IF ( .NOT. tcalc ( 1 ) ) CALL stop_program ( "pw_calculate_xcspin", &
!     "stress calculation needs full potential calculation" )
!   DO i = 1, 3
!     xc_stress ( i, i ) = xc_stress ( i, i ) + ( exc - vxc )
!   END DO
! END IF

! IF ( xc_grid % para % mode == 1 ) THEN
!   CALL mp_sum ( exc, xc_grid % para % group ) 
!   CALL mp_sum ( vxc, xc_grid % para % group ) 
!   IF ( tcalc ( 3 ) ) CALL mp_sum ( xc_stress, xc_grid % para % group ) 
! END IF

  END SUBROUTINE pw_calculate_xcspin

!******************************************************************************

END MODULE pw_xc

!******************************************************************************
