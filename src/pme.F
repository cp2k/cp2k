!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/pme [2.0] *
!!
!!   NAME
!!     pme
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM and APSI
!!
!!   MODIFICATION HISTORY
!!     JGH (21-Mar-2001) : Complete rewrite
!!
!!   SOURCE
!******************************************************************************

MODULE pme

  USE coefficient_types, ONLY : coeff_type, coeff_allocate, &
       coeff_deallocate, coeff_zero, coeff_transform_space, &
       PW_COMPLEXDATA3D, PW_REALDATA3D, &
       PW_REALSPACE, PW_COMPLEXDATA1D, PW_RECIPROCALSPACE
  USE dgs, ONLY : dg_get_strucfac, dg_get_patch, dg_sum_patch, &
       dg_sum_patch_force_3d
  USE dg_types, ONLY : dg_type
  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE hartree, ONLY : calculate_hartree
  USE kinds, ONLY: dbl, sgl, dp_size, sp_size
  USE mathconstants, ONLY : fourpi, pi
  USE particle_types, ONLY : particle_type
  USE pw_grid_types, ONLY : pw_grid_type
  USE pws, ONLY : pw_poisson_solver, init_pw_poisson_solver, pw_integral_a2b, &
                  cleanup_pw_poisson_solver
  USE pw_types, ONLY : pw_transfer, pw_type, REALSPACE, RECIPROCALSPACE, NOSPACE
  USE realspace_grid_types, ONLY : realspace_grid_type, rs_grid_allocate, &
       rs_grid_deallocate, rs_grid_setup, rs_pw_transfer
  USE simulation_cell, ONLY : cell_type, get_hinv
  USE structure_factors, ONLY : structure_factor_allocate, &
       structure_factor_deallocate
  USE structure_factor_types, ONLY : structure_factor_type
  USE timings, ONLY: timeset, timestop

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: pme_evaluate

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** pme/pme_evaluate [2.0] *
!!
!!   NAME
!!     pme_evaluate
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM and APSI
!!
!!   MODIFICATION HISTORY
!!     JGH (15-Mar-2001) : New electrostatic calculation and pressure tensor
!!     JGH (21-Mar-2001) : Complete rewrite
!!     JGH (21-Mar-2001) : Introduced real space density type for future
!!                         parallelisation
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pme_evaluate ( dg, part, box, grid_s, grid_b, &
                    ewald_param, vg_coulomb, fg_coulomb, pv_g )

  IMPLICIT NONE

! Arguments
  TYPE ( dg_type ), INTENT ( IN ), TARGET :: dg
  TYPE ( particle_type ), DIMENSION ( : ), INTENT ( IN ) :: part
  TYPE ( cell_type ), INTENT ( IN ) :: box
  TYPE ( pw_grid_type ), INTENT ( IN ) :: grid_s
  TYPE ( pw_grid_type ), INTENT ( IN ) :: grid_b
  TYPE ( ewald_parameters_type ), INTENT ( IN ) :: ewald_param
  REAL ( dbl ), INTENT ( OUT ) :: vg_coulomb
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( OUT ), OPTIONAL :: fg_coulomb
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( OUT ), OPTIONAL :: pv_g

! Locals
  INTEGER :: handle, i, j, p1, p2, npart, ipart, jpart
  TYPE ( realspace_grid_type ) :: rden, drpot ( 3 )
  TYPE ( structure_factor_type ) :: exp_igr
  TYPE ( coeff_type ) :: rhob_r, phi_r, dphi_g ( 3 ), rhos1, rhos2
  REAL ( dbl ) :: inv_fourpi_eps0, ffa, dvols
  REAL ( dbl ), DIMENSION ( 3 ) :: fat
  REAL ( dbl ), DIMENSION ( 3, 3 ) :: h_stress, f_stress

!------------------------------------------------------------------------------

!-------------- INITIALISATION ---------------------

  CALL timeset ( 'PME', 'I', ' ', handle )

  inv_fourpi_eps0 = 1.0_dbl / ( fourpi * ewald_param % eps0 )
  npart = SIZE ( part )

  CALL structure_factor_allocate ( grid_s % bounds, npart, exp_igr, &
       allocate_centre = .TRUE. )

  CALL coeff_allocate ( rhos1, grid_s, use_data = PW_REALDATA3D )
  CALL coeff_allocate ( rhos2, grid_s, use_data = PW_REALDATA3D )

  CALL rs_grid_setup ( rden, grid_b, grid_s % npts )
  CALL rs_grid_allocate ( rden )
  rden % r = 0._dbl

  IF ( rden % parallel .AND. rden % group_dim ( 1 ) /= 1 ) THEN
    CALL get_center ( part, box, exp_igr % centre, grid_b % npts )
  END IF

!-------------- DENSITY CALCULATION ----------------

  ipart = 0
  jpart = 0
  DO 

     CALL set_list ( part, npart, exp_igr % centre, p1, p2, rden, &
          ipart, jpart )
     IF ( p1 == 0 .AND. p2 == 0 ) EXIT

     ! calculate function on small boxes (we use double packing in FFT)
     CALL get_patch ( dg, part, exp_igr, box, p1, p2, grid_b, grid_s, &
          rhos1, rhos2 )

     ! add boxes to real space grid (big box)
     CALL dg_sum_patch ( rden, rhos1, exp_igr % centre ( :, p1 ) )
     IF ( p2 /= 0 ) CALL dg_sum_patch ( rden, rhos2, &
                         exp_igr % centre ( :, p2 ) )

  END DO

  CALL coeff_allocate ( rhob_r, grid_b, use_data = PW_REALDATA3D )
  rhob_r % pw % in_space = REALSPACE

  CALL rs_pw_transfer ( rden, rhob_r % pw, "FORWARD" )

  CALL rs_grid_deallocate ( rden )

!-------------- ELECTROSTATIC CALCULATION -----------

  ! allocate intermediate arrays
  DO i = 1, 3
     CALL coeff_allocate ( dphi_g ( i ), grid_b, use_data = PW_COMPLEXDATA1D )
     dphi_g ( i ) % pw % in_space = RECIPROCALSPACE
  END DO
  CALL coeff_allocate ( phi_r, grid_b, use_data = PW_REALDATA3D )
  phi_r % pw % in_space = REALSPACE

  CALL init_pw_poisson_solver ( box % green )

  CALL calculate_hartree ( rhob_r, vg_coulomb, phi_r, dphi_g, h_stress )

  CALL coeff_deallocate ( phi_r )
  vg_coulomb = vg_coulomb * inv_fourpi_eps0

!---------- END OF ELECTROSTATIC CALCULATION --------

!------------- STRESS TENSOR CALCULATION ------------

  IF ( PRESENT ( pv_g ) ) THEN
    DO i = 1, 3
      DO j = i, 3
         f_stress ( i, j ) = pw_integral_a2b ( dphi_g ( i ) % pw, &
                                             dphi_g ( j ) % pw )
         f_stress ( j, i ) = f_stress ( i, j )
      END DO
    END DO
    ffa = ( 1.0_dbl / fourpi ) * ( 0.5_dbl / dg % dg_rho0 % zet ( 1 ) ) ** 2
    f_stress = -ffa * f_stress
    pv_g = inv_fourpi_eps0 * ( h_stress + f_stress )
  ENDIF

!--------END OF STRESS TENSOR CALCULATION -----------

  CALL rs_grid_setup ( drpot, grid_b , grid_s % npts)
  CALL rs_grid_allocate ( drpot ( 1 : 3 ) )

  DO i = 1, 3
     CALL coeff_transform_space ( dphi_g ( i ), rhob_r )
     CALL coeff_deallocate ( dphi_g ( i ) )
     CALL rs_pw_transfer ( drpot ( i ), rhob_r % pw, "BACKWARD" )
  END DO

  CALL coeff_deallocate ( rhob_r )

!----------------- FORCE CALCULATION ----------------

! initialize the forces
  IF ( PRESENT ( fg_coulomb ) ) THEN
    fg_coulomb = 0._dbl
    dvols = rhos1 % pw % pw_grid % dvol * inv_fourpi_eps0

    ipart = 0
    jpart = 0
    DO

       CALL set_list ( part, npart, exp_igr % centre, p1, p2, rden, ipart, jpart )
       IF ( p1 == 0 .AND. p2 == 0 ) EXIT

     ! calculate function on small boxes (we use double packing in FFT)
       CALL get_patch_again ( dg, part, exp_igr, p1, p2, rhos1, rhos2 )

     ! sum boxes on real space grids (big box)
       CALL dg_sum_patch_force_3d ( drpot, rhos1, exp_igr % centre ( :, p1 ), fat )
       fg_coulomb ( 1, p1 ) = fg_coulomb ( 1, p1 ) - fat ( 1 ) * dvols
       fg_coulomb ( 2, p1 ) = fg_coulomb ( 2, p1 ) - fat ( 2 ) * dvols
       fg_coulomb ( 3, p1 ) = fg_coulomb ( 3, p1 ) - fat ( 3 ) * dvols
       IF ( p2 /= 0 ) THEN
          CALL dg_sum_patch_force_3d ( drpot, rhos2, exp_igr % centre ( :, p2 ), fat )
          fg_coulomb ( 1, p2 ) = fg_coulomb ( 1, p2 ) - fat ( 1 ) * dvols
          fg_coulomb ( 2, p2 ) = fg_coulomb ( 2, p2 ) - fat ( 2 ) * dvols
          fg_coulomb ( 3, p2 ) = fg_coulomb ( 3, p2 ) - fat ( 3 ) * dvols
       END IF

    END DO
  ENDIF
!--------------END OF FORCE CALCULATION -------------

!------------------CLEANING UP ----------------------

  CALL rs_grid_deallocate ( drpot ( 1 : 3 ) )
  CALL coeff_deallocate ( rhos1 )
  CALL coeff_deallocate ( rhos2 )
  CALL structure_factor_deallocate ( exp_igr )
  CALL cleanup_pw_poisson_solver ( box % green )

  CALL timestop ( 0.0_dbl, handle )

END SUBROUTINE pme_evaluate

!!*****
!******************************************************************************
!!****** pme/get_patch [1.0] *
!!
!!   NAME
!!     get_patch
!!
!!   FUNCTION
!!     Calculates local density in a small box
!!
!!   AUTHOR
!!     JGH (21-Mar-2001)
!!
!!   MODIFICATION HISTORY
!!     JGH (23-Mar-2001) : Switch to integer from particle list pointers
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE get_patch ( dg, part, exp_igr, box, p1, p2, &
           grid_b, grid_s, rhos1, rhos2 )

  IMPLICIT NONE

! Arguments
  TYPE ( dg_type ), INTENT ( IN ), TARGET :: dg
  TYPE ( particle_type ), DIMENSION ( : ), INTENT ( IN ) :: part
  TYPE ( structure_factor_type ) :: exp_igr
  TYPE ( cell_type ), INTENT ( IN ) :: box
  INTEGER, INTENT ( IN ) :: p1, p2
  TYPE ( pw_grid_type ), INTENT ( IN ) :: grid_s
  TYPE ( pw_grid_type ), INTENT ( IN ) :: grid_b
  TYPE ( coeff_type ) :: rhos1
  TYPE ( coeff_type ) :: rhos2

! Local
  TYPE ( coeff_type ), POINTER :: rho0
  REAL ( dbl ), DIMENSION ( 3 ) :: r1, r2
  REAL ( dbl ) :: q1, q2
  COMPLEX ( dbl ), DIMENSION ( : ), POINTER :: ex1, ey1, ez1, ex2, ey2, ez2
  INTEGER, DIMENSION ( : ), POINTER :: center1, center2

!------------------------------------------------------------------------------

  rho0 => dg % dg_rho0 % density

  r1 = part ( p1 ) % r
  q1 = part ( p1 ) % prop % charge
  center1 => exp_igr % centre ( :, p1 )
  ex1 => exp_igr % ex ( :, p1 )
  ey1 => exp_igr % ey ( :, p1 )
  ez1 => exp_igr % ez ( :, p1 )
  CALL dg_get_strucfac ( box, r1, grid_s % npts, grid_b % npts, center1, &
       exp_igr % lb, ex1, ey1, ez1 )

  IF ( p2 /= 0 ) THEN
     r2 = part ( p2 ) % r
     q2 = part ( p2 ) % prop % charge
     center2 => exp_igr % centre ( :, p2 )
     ex2 => exp_igr % ex ( :, p2 )
     ey2 => exp_igr % ey ( :, p2 )
     ez2 => exp_igr % ez ( :, p2 )
     CALL dg_get_strucfac ( box, r2, grid_s % npts, grid_b % npts, center2, &
          exp_igr % lb, ex2, ey2, ez2 )
  END IF

  IF ( p2 == 0 ) THEN
     CALL dg_get_patch ( rho0, rhos1, q1, ex1, ey1, ez1 )
  ELSE
     CALL dg_get_patch ( rho0, rhos1, rhos2, q1, q2, &
          ex1, ey1, ez1, ex2, ey2, ez2 )
  END IF

END SUBROUTINE get_patch

!******************************************************************************

SUBROUTINE get_patch_again ( dg, part, exp_igr, p1, p2, rhos1, rhos2 )

  IMPLICIT NONE

! Arguments
  TYPE ( dg_type ), INTENT ( IN ), TARGET :: dg
  TYPE ( particle_type ), DIMENSION ( : ), INTENT ( IN ) :: part
  TYPE ( structure_factor_type ) :: exp_igr
  INTEGER, INTENT ( IN ) :: p1, p2
  TYPE ( coeff_type ) :: rhos1
  TYPE ( coeff_type ) :: rhos2

! Local
  TYPE ( coeff_type ), POINTER :: rho0
  REAL ( dbl ) :: q1, q2
  COMPLEX ( dbl ), DIMENSION ( : ), POINTER :: ex1, ey1, ez1, ex2, ey2, ez2

!------------------------------------------------------------------------------

  rho0 => dg % dg_rho0 % density

  q1 = part ( p1 ) % prop % charge
  ex1 => exp_igr % ex ( :, p1 )
  ey1 => exp_igr % ey ( :, p1 )
  ez1 => exp_igr % ez ( :, p1 )

  IF ( p2 /= 0 ) THEN
     q2 = part ( p2 ) % prop % charge
     ex2 => exp_igr % ex ( :, p2 )
     ey2 => exp_igr % ey ( :, p2 )
     ez2 => exp_igr % ez ( :, p2 )
  END IF

  IF ( p2 == 0 ) THEN
     CALL dg_get_patch ( rho0, rhos1, q1, ex1, ey1, ez1 )
  ELSE
     CALL dg_get_patch ( rho0, rhos1, rhos2, q1, q2, &
          ex1, ey1, ez1, ex2, ey2, ez2 )
  END IF

END SUBROUTINE get_patch_again

!!*****
!******************************************************************************
!!****** pme/set_list [1.0] *
!!
!!   NAME
!!     set_list
!!
!!   FUNCTION
!!     Radvances in the link list of particles to the next to be calculated
!!
!!   AUTHOR
!!     JGH (21-Mar-2001)
!!
!!   MODIFICATION HISTORY
!!     JGH (23-Mar-2001) : Switch to integer from particle list pointers
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE set_list ( part, npart, center, p1, p2, rs, ipart, jpart )

  IMPLICIT NONE

! Arguments
  TYPE ( particle_type ), DIMENSION ( : ), INTENT ( IN ) :: part
  INTEGER, INTENT ( IN ) :: npart
  INTEGER, DIMENSION ( :, : ), INTENT ( IN ) :: center
  INTEGER, INTENT ( OUT ) :: p1, p2
  TYPE ( realspace_grid_type ) :: rs
  INTEGER, INTENT ( INOUT ) :: ipart, jpart

! Local
  INTEGER :: dir, lb, ub, ndim, npos

!------------------------------------------------------------------------------

  p1 = 0
  p2 = 0
  IF ( rs % parallel .AND. rs % group_dim ( 1 ) /= 1 ) THEN
    dir = rs % direction
    lb = rs % lb_local ( dir ) + rs % border
    ub = rs % ub_local ( dir ) - rs % border
  END IF
  IF ( rs % parallel ) THEN
    ndim = rs % group_dim ( 2 )
    npos = rs % group_coor ( 2 )
  END IF

  DO
    ipart = ipart + 1
    IF ( ipart > npart ) EXIT
    IF ( part ( ipart ) % prop % charge == 0._dbl ) CYCLE
    IF ( rs % parallel ) THEN
      IF ( rs % group_dim ( 1 ) == 1 ) THEN
        ! All processor work on the same grid
        IF ( MOD ( ipart, ndim ) == npos ) THEN
          p1 = ipart
          EXIT
        END IF
      ELSE
        ! First check if this atom is on my grid
        IF ( in_slice ( center ( dir, ipart ), lb, ub ) ) THEN
          jpart = jpart + 1
          IF ( MOD ( jpart, ndim ) == npos ) THEN
            p1 = ipart
            EXIT
          END IF
        END IF
      END IF
    ELSE
      p1 = ipart
      EXIT
    END IF
  END DO
  DO
    ipart = ipart + 1
    IF ( ipart > npart ) EXIT
    IF ( part ( ipart ) % prop % charge == 0._dbl ) CYCLE
    IF ( rs % parallel ) THEN
      IF ( rs % group_dim ( 1 ) == 1 ) THEN
        ! All processor work on the same grid
        IF ( MOD ( ipart, ndim ) == npos ) THEN
          p2 = ipart
          EXIT
        END IF
      ELSE
        ! First check if this atom is on my grid
        IF ( in_slice ( center ( dir, ipart ), lb, ub ) ) THEN
          jpart = jpart + 1
          IF ( MOD ( jpart, ndim ) == npos ) THEN
            p2 = ipart
            EXIT
          END IF
        END IF
      END IF
    ELSE
      p2 = ipart
      EXIT
    END IF
  END DO

END SUBROUTINE set_list

!******************************************************************************

FUNCTION in_slice ( pos, lb, ub ) RESULT ( internal )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: pos, lb, ub
  LOGICAL :: internal

  IF ( pos >= lb .AND. pos <= ub ) THEN
    internal = .TRUE.
  ELSE
    internal = .FALSE.
  END IF

END FUNCTION in_slice

!!*****
!******************************************************************************

SUBROUTINE get_center ( part, box, centre, npts )

  IMPLICIT NONE

! Arguments
  TYPE ( particle_type ), DIMENSION ( : ), INTENT ( IN ) :: part
  TYPE ( cell_type ), INTENT ( IN ) :: box
  INTEGER, DIMENSION ( : ), INTENT ( IN ) :: npts
  INTEGER, DIMENSION ( :, : ), INTENT ( OUT ) :: centre

! Locals
  REAL ( dbl ), DIMENSION ( 3 ) :: s
  INTEGER :: ipart

!------------------------------------------------------------------------------

  DO ipart = 1, SIZE ( part )
    ! compute the scaled coordinate of atomi
    s = MATMUL ( box % h_inv, part ( ipart ) % r )
    s = s - NINT ( s )
    ! find the closest grid point (on big grid)
    centre ( :, ipart ) = NINT ( REAL ( npts ( 1:3 ), dbl ) * s ( 1:3 ) )
  END DO

END SUBROUTINE get_center

!******************************************************************************

END MODULE pme

!******************************************************************************
