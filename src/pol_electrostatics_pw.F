!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/pol_electrostatics_pw [1.0] *
!!
!!   NAME
!!     pol_electrostatics_pw
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM 
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE pol_electrostatics_pw

  USE coefficient_types, ONLY : coeff_type, coeff_allocate, &
       coeff_deallocate, coeff_zero, coeff_transform_space, &
       PW_COMPLEXDATA3D, PW_REALDATA3D, &
       PW_REALSPACE, PW_COMPLEXDATA1D, PW_RECIPROCALSPACE
  USE coefficient_lists, ONLY : coef_list_type     
  USE dgs, ONLY : dg_get_strucfac, dg_sum_patch, dg_sum_patch_force_3d, &
                  dg_sum_patch_force_1d, dg_get_patch
  USE dg_types, ONLY : dg_type
  USE md, ONLY : thermodynamic_type
  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE hartree, ONLY : calculate_hartree
  USE kinds, ONLY: dbl, sgl, dp_size, sp_size
  USE mathconstants, ONLY : fourpi, pi
  USE particle_lists, ONLY : particle_list_type
  USE particle_types, ONLY : particle_type
  USE pw_grid_types, ONLY : pw_grid_type, HALFSPACE, FULLSPACE
  USE pws, ONLY : pw_poisson_solver, init_pw_poisson_solver, pw_integral_ab
  USE pw_types, ONLY :  pw_type, pw_add, pw_allocate, pw_copy, pw_deallocate
  USE realspace_grid_types, ONLY : realspace_grid_type, rs_grid_allocate, &
       rs_grid_deallocate, rs_grid_setup, rs_pw_transfer
  USE simulation_cell, ONLY : cell_type
  USE structure_factors, ONLY : structure_factor_allocate, structure_factor_deallocate
  USE structure_factor_types, ONLY : structure_factor_type
  USE termination, ONLY : stop_memory
  USE timings, ONLY: timeset, timestop

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: electrostatics

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** pol_electrostatics_pw/electrostatics [1.0] *
!!
!!   NAME
!!     electrostatics
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM 
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE electrostatics ( dg_part, pw, part, box, grid_s,  &
                           grid_b, ewald_param, energy, f_coef, f_part )

  IMPLICIT NONE

! Arguments
  TYPE ( dg_type ), INTENT ( IN ), DIMENSION ( : ) :: dg_part
  TYPE ( pw_type ), INTENT ( INOUT ) :: pw
  TYPE ( particle_type ), DIMENSION ( : ), INTENT ( IN ) :: part
  TYPE ( cell_type ), INTENT ( IN ) :: box
  TYPE ( pw_grid_type ), INTENT ( IN ), DIMENSION ( : ) :: grid_s
  TYPE ( pw_grid_type ), INTENT ( IN ), DIMENSION ( : ) :: grid_b
  TYPE ( ewald_parameters_type ), INTENT ( IN ) :: ewald_param
  REAL ( dbl ), INTENT ( OUT ) :: energy
  COMPLEX ( dbl ), INTENT ( OUT ), DIMENSION ( : ), OPTIONAL :: f_coef
  REAL ( dbl ), INTENT ( OUT ), DIMENSION ( :, : ), OPTIONAL :: f_part

! Locals
  INTEGER :: handle
  INTEGER :: ig, allocstat, ngrids, ngt
  INTEGER :: npart, ncoef, gt, i, igtot
  TYPE ( structure_factor_type ), ALLOCATABLE, DIMENSION ( : ) :: exp_igr_part
  TYPE ( realspace_grid_type ) :: rho0_rs, drpot1 ( 3 )
  TYPE ( realspace_grid_type ) :: phi1_rs
  TYPE ( coeff_type ) :: rhos1, rhos2
  TYPE ( coeff_type ) :: rho0b_g, rho0b_r, phi1_r, dphi1_g ( 3 )
  TYPE ( coeff_type ) :: rho_tot, coeff_dummy 
  REAL ( dbl ) :: inv_fourpi_eps0
  REAL ( dbl ) :: e_rho_rho
  REAL ( dbl ) :: dvols
  LOGICAL :: first_time, match

!------------------------------------------------------------------------------

!-------------- INITIALISATION ---------------------

! Allocating the structure factor vector to the number of grid_types
  ngrids = SIZE ( dg_part )
  ALLOCATE ( exp_igr_part ( ngrids ), STAT=allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "pol_electrostatics_pw", &
   "exp_igr_part", ngrids )

  DO ig = 1, ngrids

    npart = dg_part ( ig ) % nparts
    
    gt = dg_part ( ig )% grid_index

    CALL structure_factor_allocate ( grid_s ( gt ) % bounds, npart,  &
          exp_igr_part ( ig ), allocate_centre = .TRUE. )
  ENDDO

! Allocate reference coefficient types ( e.g. with largest cut-off )

  CALL coeff_allocate ( rho_tot , pw % pw_grid, use_data = PW_COMPLEXDATA1D )
  rho_tot % pw % in_space = PW_RECIPROCALSPACE

  CALL coeff_allocate ( coeff_dummy , pw % pw_grid, use_data = PW_COMPLEXDATA1D )
  coeff_dummy % pw % in_space = PW_RECIPROCALSPACE
 
  CALL coeff_zero ( rho_tot )
  CALL coeff_zero ( coeff_dummy )

  IF ( PRESENT ( f_coef ) )  f_coef = 0._dbl
  IF ( PRESENT ( f_part ) )  f_part = 0._dbl

!-------------- DENSITY CALCULATION ---------------------

  ngt = SIZE ( grid_b )

  DO gt = 1, ngt
    first_time = .TRUE.
    igtot = SIZE ( dg_part )
    IF ( igtot == 0 ) EXIT
    match = .FALSE.
    DO ig = 1, igtot  ! RHO0
      IF ( gt == dg_part ( ig ) % grid_index ) THEN
        match = .TRUE.
        IF ( first_time ) THEN
! Allocate small coefficient types
          CALL coeff_allocate ( rhos1,  grid_s ( gt ), use_data = PW_REALDATA3D )
          CALL coeff_allocate ( rhos2,  grid_s ( gt ), use_data = PW_REALDATA3D )

! Allocate large real-space structure based on grid_b
          CALL rs_grid_setup ( rho0_rs, grid_b ( gt ), grid_s ( gt ) % npts )
          CALL rs_grid_allocate ( rho0_rs )
          rho0_rs % r = 0._dbl
          first_time = .FALSE.
        ENDIF
! Parrallelization
        IF (  rho0_rs % parallel .AND. rho0_rs % group_dim ( 1 ) /= 1 )  THEN
          CALL get_center ( part, box, exp_igr_part ( ig ) % centre,  &
          grid_b ( gt ) % npts )
        END IF
! Get density rho0
        CALL rho0_evaluate ( dg_part ( ig ), rho0_rs, part, box,  &
                  exp_igr_part ( ig ), grid_s ( gt ), grid_b ( gt ), rhos1, rhos2 ) 
      ELSE
        CYCLE
      ENDIF
    ENDDO
    IF ( match ) THEN
! Allocate large coefficient structure based on grid_b
      CALL coeff_allocate ( rho0b_r, grid_b ( gt ), use_data = PW_COMPLEXDATA3D )
      CALL coeff_allocate ( rho0b_g, grid_b ( gt ), use_data = PW_COMPLEXDATA1D )

      rho0b_r % pw % in_space = PW_REALSPACE

! Transfer from real-space -> coefficient_type
      CALL rs_pw_transfer ( rho0_rs, rho0b_r % pw, "FORWARD" )
 
! DO FFT from rho0b_r => rho0b_g 
      CALL coeff_transform_space ( rho0b_r, rho0b_g )

! Copy present rho_tot to coeff_dummy
      CALL pw_copy ( rho_tot % pw, coeff_dummy % pw )

! Sum to the reference grid
      CALL pw_add ( rho0b_g % pw, coeff_dummy % pw, rho_tot % pw )

! do a little cleaning
      CALL coeff_deallocate ( rho0b_g )
      CALL coeff_deallocate ( rho0b_r )
      CALL rs_grid_deallocate ( rho0_rs )
      CALL coeff_deallocate ( rhos2 )
      CALL coeff_deallocate ( rhos1 )
    ENDIF
  ENDDO

! Sum the contribution with coeff for the TOTAL density

! Copy present rho_tot to coeff_dummy
  CALL pw_copy ( rho_tot % pw, coeff_dummy % pw )
! Sum to the reference grid
  CALL pw_add ( pw, coeff_dummy % pw, rho_tot % pw )

  CALL coeff_deallocate ( coeff_dummy )

!-------------- ELECTROSTATIC CALCULATION -----------

  inv_fourpi_eps0 = 1.0_dbl / ( fourpi * ewald_param % eps0 )

! allocate intermediate coefficient arrays on the reference grid for
! the potential due to rho0 and drho separately

  CALL coeff_allocate ( phi1_r, pw % pw_grid, use_data = PW_COMPLEXDATA1D )
  phi1_r % pw % in_space = PW_RECIPROCALSPACE

  CALL init_pw_poisson_solver ( box % green )

  CALL calculate_hartree ( rho_tot, e_rho_rho, phi1_r )
  energy =  e_rho_rho * inv_fourpi_eps0

!-------------- ELECTROSTATIC CALCULATION DONE -----------


!-------------- FORCE CALCULATION ---------------------

  IF ( PRESENT ( f_coef ) )  THEN 

    f_coef = CMPLX ( 2._dbl * REAL ( phi1_r % pw % cc, dbl ), &
                     2._dbl * AIMAG ( phi1_r % pw % cc ), dbl )
  ENDIF 

  CALL coeff_deallocate ( phi1_r )
 
  IF ( PRESENT ( f_part ) ) THEN
    DO gt = 1, ngt
      first_time = .TRUE.
      igtot = SIZE ( dg_part )
      IF ( igtot == 0 ) EXIT
      match = .FALSE.
      DO ig = 1, igtot  ! RHO0
        IF ( gt == dg_part ( ig ) % grid_index ) THEN
          match = .TRUE.
          IF ( first_time ) THEN
! set up real-space structure for the derivative of the potential
! on the reference grid
            CALL rs_grid_setup ( drpot1, grid_b ( gt ) , grid_s ( gt ) % npts )

            CALL rs_grid_allocate ( drpot1 ( 1 : 3 ) )

            CALL coeff_allocate ( phi1_r, grid_b ( gt ), use_data = PW_COMPLEXDATA1D )
            phi1_r % pw % in_space = PW_RECIPROCALSPACE

! reallocate coeff with REALDATA3D to use as a dummy variable
! to transfer dphi -> drpot
            CALL coeff_allocate ( coeff_dummy , grid_b ( gt ), use_data = PW_REALDATA3D )
            coeff_dummy % pw % in_space = PW_REALSPACE

! allocate intermediate coefficient arrays on the reference grid for
! the derivative of the potential due to rho0 and drho separately
            DO i = 1, 3
               CALL coeff_allocate ( dphi1_g ( i ), grid_b ( gt ), use_data = PW_COMPLEXDATA1D )
               dphi1_g ( i ) % pw % in_space = PW_RECIPROCALSPACE
            END DO

            CALL calculate_hartree ( rho_tot, vhartree = phi1_r, dvhartree = dphi1_g )

            DO i = 1, 3
               CALL coeff_transform_space ( dphi1_g ( i ), coeff_dummy )
               CALL coeff_deallocate ( dphi1_g ( i ) )
               CALL rs_pw_transfer ( drpot1 ( i ), coeff_dummy % pw, "BACKWARD" )
            END DO

            CALL coeff_deallocate ( coeff_dummy )

            CALL coeff_allocate ( rhos1,  grid_s ( gt ), use_data = PW_REALDATA3D )
            CALL coeff_allocate ( rhos2,  grid_s ( gt ), use_data = PW_REALDATA3D )
            first_time = .FALSE.
          ENDIF
! force on particles
          CALL force_part_evaluate ( dg_part ( ig ), drpot1, part, &
             f_part, box, rhos1, rhos2, exp_igr_part ( ig ), grid_s ( gt ), &
             ewald_param )
          CALL structure_factor_deallocate ( exp_igr_part ( ig ) )
        ELSE
          CYCLE
        ENDIF
      ENDDO
      IF ( match ) THEN
        CALL coeff_deallocate ( rhos2 )
        CALL coeff_deallocate ( rhos1 )
        CALL coeff_deallocate ( phi1_r ) 
        CALL rs_grid_deallocate ( drpot1 ( 1 : 3 ) )
      ENDIF
    ENDDO
  ENDIF

! clean up
  CALL coeff_deallocate ( rho_tot )

  DEALLOCATE ( exp_igr_part, STAT=allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "pol_electrostatics_pw", &
   "exp_igr_part" )

END SUBROUTINE electrostatics

!******************************************************************************
!!****** pol_electrostatics_pw/rho0_evaluate [1.0] *
!!
!!   NAME
!!     rho0_evaluate
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM 
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE rho0_evaluate ( dg, rden, part, box, exp_igr, grid_s,  &
                             grid_b, rhos1, rhos2 ) 

  IMPLICIT NONE

! Arguments
  TYPE ( dg_type ), INTENT ( IN ), TARGET :: dg
  TYPE ( particle_type ), DIMENSION ( : ), INTENT ( IN ) :: part
  TYPE ( cell_type ), INTENT ( IN ) :: box
  TYPE ( pw_grid_type ), INTENT ( IN ) :: grid_s
  TYPE ( pw_grid_type ), INTENT ( IN ) :: grid_b
  TYPE ( realspace_grid_type ), INTENT ( INOUT ) :: rden
  TYPE ( structure_factor_type ), INTENT ( OUT ) :: exp_igr
  TYPE ( coeff_type ), INTENT ( INOUT ) :: rhos1, rhos2

! Locals
  INTEGER :: handle, p1, p2, ip1, ip2
  TYPE ( particle_list_type ), POINTER :: plist
  REAL ( dbl ) :: r1 ( 3 ), r2 ( 3 ), q1, q2

!------------------------------------------------------------------------------


  CALL timeset ( 'RHO0_EVALUATE', 'I', ' ', handle )



! Initialize the index of p1 and p2. This is to ensure the
! correct mapping between the particle and its array location
! in exp_igr
  ip1 = 0
  ip2 = 0

  plist => dg % plist
  DO 
     CALL set_list_rho0 ( plist, p1, p2, ip1, ip2 )
     IF ( p1 == 0 .AND. p2 == 0 )  EXIT
     r1 = part ( p1 ) % r
     q1 = 1._dbl 
     r2 = 0._dbl 
     q2 = 0._dbl 
     IF ( p2 /= 0 ) THEN 
       r2 = part ( p2 ) % r
       q2 = 1._dbl 
     ENDIF

     ! calculate function on small boxes (we use double packing in FFT)
     CALL get_patch ( dg, r1, r2, q1, q2, exp_igr, box, p1, p2, ip1, ip2, &
                     grid_b, grid_s, rhos1, rhos2 )
     ! add boxes to real space grid (big box)
     CALL dg_sum_patch ( rden, rhos1, exp_igr % centre ( :, ip1 ) )
     IF ( p2 /= 0 ) CALL dg_sum_patch ( rden, rhos2, exp_igr % centre ( :, ip2 ) )

     ip1 = ip2
  END DO

  CALL timestop ( 0.0_dbl, handle )

END SUBROUTINE rho0_evaluate

!******************************************************************************
!!****** pol_electrostatics_pw/force_part_evaluate [1.0] *
!!
!!   NAME
!!     force_evaluate_part
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM 
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE force_part_evaluate ( dg,  drpot1, part, force, box, &
                                rhos1, rhos2, exp_igr, grid_s, ewald_param )

! Arguments
  TYPE ( dg_type ), INTENT ( IN ) :: dg
  TYPE ( particle_type ), DIMENSION ( : ), INTENT ( IN ) :: part
  TYPE ( cell_type ), INTENT ( IN ) :: box
  TYPE ( pw_grid_type ), INTENT ( IN ) :: grid_s
  TYPE ( realspace_grid_type ), DIMENSION ( : ), INTENT ( INOUT ) :: drpot1
  TYPE ( structure_factor_type ), INTENT ( IN ) :: exp_igr
  TYPE ( coeff_type ), INTENT ( INOUT ) :: rhos1, rhos2
  TYPE ( ewald_parameters_type ), INTENT ( IN ) :: ewald_param
  REAL ( dbl ), INTENT (INOUT), DIMENSION (:,:) :: force

! Locals
  INTEGER :: handle, p1, p2, ip1, ip2
  REAL ( dbl ) :: dvols, inv_fourpi_eps0, q1, q2
  REAL ( dbl ), DIMENSION ( 3 ) :: fat
  TYPE ( particle_list_type ), POINTER :: plist

  CALL timeset ( 'FORCE_EVALUATE_PART', 'I', ' ', handle )

  inv_fourpi_eps0 = 1.0_dbl / ( fourpi * ewald_param % eps0 )

! initialize the forces
  dvols = rhos1 % pw % pw_grid % dvol * inv_fourpi_eps0

  plist => dg % plist
  ip1 = 0
  ip2 = 0
  DO

     CALL set_list_rho0 ( plist, p1, p2, ip1, ip2 )
     IF ( p1 == 0 .AND. p2 == 0 ) EXIT
     q1 = 1._dbl
     q2=0._dbl
     IF ( p2 /= 0 ) q2 = 1._dbl 
     ! calculate function on small boxes (we use double packing in FFT)
     CALL get_patch_again ( dg, q1, q2, exp_igr, p1, p2, ip1, ip2, &
                            rhos1, rhos2 )

     ! sum boxes on real space grids (big box) for forces due to drpot1:
     CALL dg_sum_patch_force_3d ( drpot1, rhos1, exp_igr % centre ( :, ip1 ), fat )
     force ( 1, p1 ) = force ( 1, p1 ) - fat ( 1 ) * dvols
     force ( 2, p1 ) = force ( 2, p1 ) - fat ( 2 ) * dvols
     force ( 3, p1 ) = force ( 3, p1 ) - fat ( 3 ) * dvols
     IF ( p2 /= 0 ) THEN
        CALL dg_sum_patch_force_3d ( drpot1, rhos2, exp_igr % centre ( :, ip2 ), fat )
        force ( 1, p2 ) = force ( 1, p2 ) - fat ( 1 ) * dvols
        force ( 2, p2 ) = force ( 2, p2 ) - fat ( 2 ) * dvols
        force ( 3, p2 ) = force ( 3, p2 ) - fat ( 3 ) * dvols
     END IF

     ip1 = ip2
  END DO

  CALL timestop ( 0.0_dbl, handle )

END SUBROUTINE force_part_evaluate

!******************************************************************************
!!****** pol_electrostatics_pw/set_list_rho0 [1.0] *
!!
!!   NAME
!!     set_list_rho0
!!
!!   FUNCTION
!!     Radvances in the link list of particles to the next to be calculated
!!
!!   AUTHOR
!!     CJM (18-April-2001)
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE set_list_rho0 ( plist, p1, p2, ip1, ip2 )

  IMPLICIT NONE

! Arguments
  TYPE ( particle_list_type ),  POINTER :: plist
  INTEGER, INTENT ( OUT ) :: p1, p2
  INTEGER, INTENT ( INOUT ) :: ip1, ip2

!------------------------------------------------------------------------------
  p1 = 0
  p2 = 0

  IF ( ASSOCIATED ( plist ) ) THEN
    ip1 = ip1 + 1
    ip2 = ip1
    p1 = plist % particle_index
    IF ( ASSOCIATED ( plist % next ) ) THEN
      ip2 = ip2 + 1
      plist => plist % next
      p2 = plist % particle_index
    ENDIF
    plist => plist % next
  END IF

END SUBROUTINE set_list_rho0

!******************************************************************************
!!****** pol_electrostatics_pw/get_patch [1.0] *
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

SUBROUTINE get_patch ( dg, r1, r2, q1, q2, exp_igr, box, p1, p2, &
           ip1, ip2, grid_b, grid_s, rhos1, rhos2 )

  IMPLICIT NONE

! Arguments
  TYPE ( dg_type ), INTENT ( IN ), TARGET :: dg
  TYPE ( structure_factor_type ) :: exp_igr
  TYPE ( cell_type ), INTENT ( IN ) :: box
  INTEGER, INTENT ( IN ) :: p1, p2, ip1, ip2
  TYPE ( pw_grid_type ), INTENT ( IN ) :: grid_s
  TYPE ( pw_grid_type ), INTENT ( IN ) :: grid_b
  TYPE ( coeff_type ) :: rhos1
  TYPE ( coeff_type ) :: rhos2
  REAL ( dbl ), DIMENSION ( 3 ), INTENT ( IN ) :: r1, r2
  REAL ( dbl ), INTENT ( IN ) :: q1, q2

! Local
  TYPE ( coeff_type ), POINTER :: rho0
  COMPLEX ( dbl ), DIMENSION ( : ), POINTER :: ex1, ey1, ez1, ex2, ey2, ez2
  INTEGER, DIMENSION ( : ), POINTER :: center1, center2

!------------------------------------------------------------------------------

  rho0 => dg % dg_rho0 % density

  center1 => exp_igr % centre ( :, ip1 )
  ex1 => exp_igr % ex ( :, ip1 )
  ey1 => exp_igr % ey ( :, ip1 )
  ez1 => exp_igr % ez ( :, ip1 )
  CALL dg_get_strucfac ( box, r1, grid_s % npts, grid_b % npts, center1, &
       exp_igr % lb, ex1, ey1, ez1 )

  IF ( p2 /= 0 ) THEN
     center2 => exp_igr % centre ( :, ip2 )
     ex2 => exp_igr % ex ( :, ip2 )
     ey2 => exp_igr % ey ( :, ip2 )
     ez2 => exp_igr % ez ( :, ip2 )
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

SUBROUTINE get_patch_again ( dg, q1, q2, exp_igr, p1, p2, ip1, ip2, rhos1, rhos2 )

  IMPLICIT NONE

! Arguments
  REAL ( dbl ), INTENT ( IN ) :: q1, q2
  TYPE ( dg_type ), INTENT ( IN ), TARGET :: dg
  TYPE ( structure_factor_type ) :: exp_igr
  INTEGER, INTENT ( IN ) :: p1, p2, ip1, ip2
  TYPE ( coeff_type ) :: rhos1
  TYPE ( coeff_type ) :: rhos2

! Local
  TYPE ( coeff_type ), POINTER :: rho0
  COMPLEX ( dbl ), DIMENSION ( : ), POINTER :: ex1, ey1, ez1, ex2, ey2, ez2

!------------------------------------------------------------------------------

  rho0 => dg % dg_rho0 % density

  ex1 => exp_igr % ex ( :, ip1 )
  ey1 => exp_igr % ey ( :, ip1 )
  ez1 => exp_igr % ez ( :, ip1 )

  IF ( p2 /= 0 ) THEN
     ex2 => exp_igr % ex ( :, ip2 )
     ey2 => exp_igr % ey ( :, ip2 )
     ez2 => exp_igr % ez ( :, ip2 )
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

END MODULE pol_electrostatics_pw

!******************************************************************************
