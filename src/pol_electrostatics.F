!-----------------------------------------------------------------------------
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/pol_electrostatics [1.0] *
!!
!!   NAME
!!     pol_electrostatics
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

MODULE pol_electrostatics

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
  USE simulation_cell, ONLY : cell_type, get_hinv
  USE structure_factors, ONLY : structure_factor_allocate, &
       structure_factor_deallocate
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
!!****** pol_electrostatics/electrostatics [1.0] *
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

SUBROUTINE electrostatics ( dg_part, dg_coef, part, coeff, box, grid_s,  &
                         grid_b, ewald_param, iref, energy, f_coef, f_part )

  IMPLICIT NONE

! Arguments
  TYPE ( dg_type ), INTENT ( IN ), DIMENSION ( : ) :: dg_part
  TYPE ( dg_type ), INTENT ( IN ), DIMENSION ( : ) :: dg_coef
  TYPE ( coeff_type ), INTENT ( IN ) :: coeff
  TYPE ( particle_type ), DIMENSION ( : ), INTENT ( IN ) :: part
  TYPE ( cell_type ), INTENT ( IN ) :: box
  TYPE ( pw_grid_type ), INTENT ( IN ), DIMENSION ( : ) :: grid_s
  TYPE ( pw_grid_type ), INTENT ( IN ), DIMENSION ( : ) :: grid_b
  TYPE ( ewald_parameters_type ), INTENT ( IN ) :: ewald_param
  INTEGER, INTENT ( IN ) :: iref
  REAL ( dbl ), INTENT ( OUT ) :: energy
  REAL ( dbl ), INTENT ( OUT ), DIMENSION ( : ), OPTIONAL :: f_coef
  REAL ( dbl ), INTENT ( OUT ), DIMENSION ( :, : ), OPTIONAL :: f_part

! Locals
  INTEGER :: handle
  INTEGER :: ig, allocstat, ngrids
  INTEGER :: npart, ncoef, gt, i
  TYPE ( realspace_grid_type ) :: rho0_rs, drpot1 ( 3 )
  TYPE ( realspace_grid_type ) :: drho_rs, drpot2 ( 3 )
  TYPE ( realspace_grid_type ) :: phi1_rs, phi2_rs
  TYPE ( structure_factor_type ), ALLOCATABLE, DIMENSION ( : ) :: exp_igr_part
  TYPE ( structure_factor_type ), ALLOCATABLE, DIMENSION ( : ) :: exp_igr_coef
  TYPE ( coeff_type ) :: rhos1, rhos2
  TYPE ( coeff_type ) :: rho0b_g, rho0b_r, phi1_r, dphi1_g ( 3 )
  TYPE ( coeff_type ) :: drhob_g, drhob_r, phi2_r, dphi2_g ( 3 )
  TYPE ( coeff_type ) :: rho0_ref, drho_ref, coeff_dummy 
  REAL ( dbl ) :: inv_fourpi_eps0
  REAL ( dbl ) :: e_rho0_rho0, e_drho_rho0, e_drho_drho

!------------------------------------------------------------------------------

!-------------- INITIALISATION ---------------------

! allocating the structure factor vector to the number of grid_types
  ngrids = SIZE ( dg_coef )
  ALLOCATE ( exp_igr_coef ( ngrids ), STAT=allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "pol_electrostatics", &
   "exp_igr_coef", ngrids )

  ngrids = SIZE ( dg_part )
  ALLOCATE ( exp_igr_part ( ngrids ), STAT=allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "pol_electrostatics", &
   "exp_igr_part", ngrids )

  DO ig = 1, SIZE ( dg_coef )

    ncoef = dg_coef ( ig ) % ncoefs

    gt = dg_coef (ig)% grid_index

    CALL structure_factor_allocate ( grid_s ( gt ) % bounds, ncoef,  &
          exp_igr_coef ( ig ), allocate_centre = .TRUE. )
  ENDDO

  DO ig = 1, SIZE ( dg_part )

    npart = dg_part ( ig ) % nparts
    
    gt = dg_part (ig)% grid_index

    CALL structure_factor_allocate ( grid_s ( gt ) % bounds, npart,  &
          exp_igr_part ( ig ), allocate_centre = .TRUE. )
  ENDDO

! Allocate reference coefficient types ( e.g. with largest cut-off )

  CALL coeff_allocate ( rho0_ref , grid_b ( iref ), use_data = PW_COMPLEXDATA1D )
  rho0_ref % pw % in_space = PW_RECIPROCALSPACE
  CALL coeff_allocate ( drho_ref , grid_b ( iref ), use_data = PW_COMPLEXDATA1D )
  drho_ref % pw % in_space = PW_RECIPROCALSPACE
  CALL coeff_allocate ( coeff_dummy , grid_b ( iref ), use_data = PW_COMPLEXDATA1D )
  coeff_dummy % pw % in_space = PW_RECIPROCALSPACE
 
  CALL coeff_zero ( rho0_ref )
  CALL coeff_zero ( drho_ref )
  CALL coeff_zero ( coeff_dummy )

  IF ( PRESENT ( f_coef ) )  f_coef = 0._dbl
  IF ( PRESENT ( f_part ) )  f_part = 0._dbl

!-------------- DENSITY CALCULATION ---------------------

  DO ig = 1, SIZE ( dg_coef )   ! DRHO
  
    gt = dg_coef ( ig ) % grid_index

! Allocate small coefficient types
    CALL coeff_allocate ( rhos1,  grid_s ( gt ), use_data = PW_REALDATA3D )
    CALL coeff_allocate ( rhos2,  grid_s ( gt ), use_data = PW_REALDATA3D )

! Allocate large real-space structure based on grid_b
    CALL rs_grid_setup ( drho_rs, grid_b ( gt ), grid_s ( gt ) % npts )
    CALL rs_grid_allocate ( drho_rs )
    drho_rs % r = 0._dbl

! Parrallelization
    IF (  drho_rs % parallel .AND. drho_rs % group_dim ( 1 ) /= 1 )  THEN
      CALL get_center ( part, box, exp_igr_coef ( ig ) % centre, &
                        grid_b ( gt ) % npts )
    END IF

! Get density drho
    CALL drho_evaluate ( dg_coef ( ig ), drho_rs, part, coeff, box,  &
              exp_igr_coef ( ig ), grid_s ( gt ), grid_b ( gt ), rhos1, rhos2 ) 

! Allocate large coefficient structure based on grid_b  
    CALL coeff_allocate ( drhob_r, grid_b ( gt ), use_data = PW_COMPLEXDATA3D )
    CALL coeff_allocate ( drhob_g, grid_b ( gt ), use_data = PW_COMPLEXDATA1D )

    drhob_r % pw % in_space = PW_REALSPACE

! Transfer from real-space -> coefficient_type
    CALL rs_pw_transfer ( drho_rs, drhob_r % pw, "FORWARD" )
 
! DO FFT from drhob_r => drhob_g 
    CALL coeff_transform_space ( drhob_r, drhob_g )

! Copy present drho_ref to coeff_dummy
    CALL pw_copy ( drho_ref % pw, coeff_dummy % pw )

! Sum to the reference grid
    CALL pw_add ( drhob_g % pw, coeff_dummy % pw, drho_ref % pw )

! do a little cleaning
    CALL coeff_deallocate ( drhob_g )
    CALL coeff_deallocate ( drhob_r )
    CALL rs_grid_deallocate ( drho_rs )
    CALL coeff_deallocate ( rhos2 )
    CALL coeff_deallocate ( rhos1 )

  ENDDO

  DO ig = 1, SIZE ( dg_part )  ! RHO0
  
    gt = dg_part ( ig ) % grid_index

! Allocate small coefficient types
    CALL coeff_allocate ( rhos1,  grid_s ( gt ), use_data = PW_REALDATA3D )
    CALL coeff_allocate ( rhos2,  grid_s ( gt ), use_data = PW_REALDATA3D )

! Allocate large real-space structure based on grid_b
    CALL rs_grid_setup ( rho0_rs, grid_b ( gt ), grid_s ( gt ) % npts )

    CALL rs_grid_allocate ( rho0_rs )

    rho0_rs % r = 0._dbl

! Parrallelization
    IF (  rho0_rs % parallel .AND. rho0_rs % group_dim ( 1 ) /= 1 )  THEN
      CALL get_center ( part, box, exp_igr_part ( ig ) % centre, grid_b ( gt ) % npts )
    END IF

! Get density rho0
    CALL rho0_evaluate ( dg_part ( ig ), rho0_rs, part, box, exp_igr_part ( ig ), &
                             grid_s ( gt ), grid_b ( gt ), rhos1, rhos2 ) 

! Allocate large coefficient structure based on grid_b
    CALL coeff_allocate ( rho0b_r, grid_b ( gt ), use_data = PW_COMPLEXDATA3D )
    CALL coeff_allocate ( rho0b_g, grid_b ( gt ), use_data = PW_COMPLEXDATA1D )

    rho0b_r % pw % in_space = PW_REALSPACE

! Transfer from real-space -> coefficient_type
    CALL rs_pw_transfer ( rho0_rs, rho0b_r % pw, "FORWARD" )
 
! DO FFT from rho0b_r => rho0b_g 
    CALL coeff_transform_space ( rho0b_r, rho0b_g )

! Copy present rho0_ref to coeff_dummy
    CALL pw_copy ( rho0_ref % pw, coeff_dummy % pw )

! Sum to the reference grid
    CALL pw_add ( rho0b_g % pw, coeff_dummy % pw, rho0_ref % pw )

! do a little cleaning
    CALL coeff_deallocate ( rho0b_g )
    CALL coeff_deallocate ( rho0b_r )
    CALL rs_grid_deallocate ( rho0_rs )
    CALL coeff_deallocate ( rhos2 )
    CALL coeff_deallocate ( rhos1 )

  ENDDO

  CALL coeff_deallocate ( coeff_dummy )

!-------------- ELECTROSTATIC CALCULATION -----------

  inv_fourpi_eps0 = 1.0_dbl / ( fourpi * ewald_param % eps0 )

! allocate intermediate coefficient arrays on the reference grid for
! the potential due to rho0 and drho separately

  CALL coeff_allocate ( phi1_r, grid_b ( iref ), use_data = PW_COMPLEXDATA1D )
  phi1_r % pw % in_space = PW_RECIPROCALSPACE

  CALL coeff_allocate ( phi2_r, grid_b ( iref ), use_data = PW_COMPLEXDATA1D )
  phi2_r % pw % in_space = PW_RECIPROCALSPACE

  CALL init_pw_poisson_solver ( box % green )

  CALL calculate_hartree ( rho0_ref, e_rho0_rho0, phi1_r )

  CALL calculate_hartree ( drho_ref, e_drho_drho, phi2_r )

  e_drho_rho0 = pw_integral_ab ( drho_ref % pw , phi1_r % pw )

  energy = ( e_rho0_rho0 + e_drho_drho + &
                                 e_drho_rho0 ) * inv_fourpi_eps0

  CALL coeff_deallocate ( phi2_r )
  CALL coeff_deallocate ( phi1_r )


!-------------- ELECTROSTATIC CALCULATION DONE -----------


!-------------- FORCE CALCULATION ---------------------

  IF ( PRESENT ( f_coef ) )  THEN 
    DO ig = 1, SIZE ( dg_coef )  ! DRHO
      gt = dg_coef ( ig ) % grid_index

! allocate intermediate coefficient arrays on the reference grid for
! the derivative of the potential due to rho0 and drho separately
      DO i = 1, 3
         CALL coeff_allocate ( dphi1_g ( i ), grid_b ( gt ), use_data = PW_COMPLEXDATA1D )
         dphi1_g ( i ) % pw % in_space = PW_RECIPROCALSPACE
 
         CALL coeff_allocate ( dphi2_g ( i ), grid_b ( gt ), use_data = PW_COMPLEXDATA1D )
         dphi2_g ( i ) % pw % in_space = PW_RECIPROCALSPACE
      END DO
      CALL coeff_allocate ( phi1_r, grid_b ( gt ), use_data = PW_COMPLEXDATA1D )
      phi1_r % pw % in_space = PW_RECIPROCALSPACE

      CALL calculate_hartree ( rho0_ref, vhartree = phi1_r, dvhartree = dphi1_g )

      CALL coeff_allocate ( phi2_r, grid_b ( gt ), use_data = PW_COMPLEXDATA1D )
      phi2_r % pw % in_space = PW_RECIPROCALSPACE

      CALL calculate_hartree ( drho_ref, vhartree = phi2_r, dvhartree = dphi2_g )

! Create real_space structure of phi1_r and phi2_r to use in the
! force calculation sum_patch_force only acts on real_space types

      CALL rs_grid_setup ( phi1_rs, grid_b ( gt ), grid_s ( gt ) % npts )
      CALL rs_grid_setup ( phi2_rs, grid_b ( gt ), grid_s ( gt ) % npts )

      CALL rs_grid_allocate ( phi1_rs )
      CALL rs_grid_allocate ( phi2_rs )

! reallocate coeff with REALDATA3D to use as a dummy variable
! to transfer dphi -> drpot and phi_r -> phi_rs
      CALL coeff_allocate ( coeff_dummy , grid_b ( gt ), use_data = PW_REALDATA3D )
      coeff_dummy % pw % in_space = PW_REALSPACE

      CALL coeff_transform_space ( phi1_r, coeff_dummy )
      CALL rs_pw_transfer ( phi1_rs, coeff_dummy % pw, "BACKWARD" )

      CALL coeff_transform_space ( phi2_r, coeff_dummy )
      CALL rs_pw_transfer ( phi2_rs, coeff_dummy % pw, "BACKWARD" )

! set up real-space structure for the derivative of the potential
! on the reference grid
      CALL rs_grid_setup ( drpot1, grid_b ( gt ) , grid_s ( gt ) % npts )
      CALL rs_grid_setup ( drpot2, grid_b ( gt ) , grid_s ( gt ) % npts )

      CALL rs_grid_allocate ( drpot1 ( 1 : 3 ) )
      CALL rs_grid_allocate ( drpot2 ( 1 : 3 ) )

      DO i = 1, 3
         CALL coeff_transform_space ( dphi1_g ( i ), coeff_dummy )
         CALL coeff_deallocate ( dphi1_g ( i ) )
         CALL rs_pw_transfer ( drpot1 ( i ), coeff_dummy % pw, "BACKWARD" )
      END DO

      DO i = 1, 3
         CALL coeff_transform_space ( dphi2_g ( i ), coeff_dummy )
         CALL coeff_deallocate ( dphi2_g ( i ) )
         CALL rs_pw_transfer ( drpot2 ( i ), coeff_dummy % pw, "BACKWARD" )
      END DO

      CALL coeff_deallocate ( coeff_dummy )

      CALL coeff_allocate ( rhos1,  grid_s ( gt ), use_data = PW_REALDATA3D )
      CALL coeff_allocate ( rhos2,  grid_s ( gt ), use_data = PW_REALDATA3D )

! force on coefs and ions
      CALL force_coef_evaluate ( dg_coef ( ig ), drpot1, drpot2, phi1_rs,  &
           phi2_rs, part, coeff, f_part, f_coef, box, rhos1, rhos2,  &
           exp_igr_coef ( ig ), grid_s ( gt ), ewald_param )

! clean up
      CALL coeff_deallocate ( rhos2 )
      CALL coeff_deallocate ( rhos1 )
      CALL rs_grid_deallocate ( drpot2 ( 1 : 3 ) )
      CALL rs_grid_deallocate ( drpot1 ( 1 : 3 ) )
      CALL rs_grid_deallocate ( phi2_rs )
      CALL rs_grid_deallocate ( phi1_rs )
      CALL structure_factor_deallocate ( exp_igr_coef ( ig ) )
    ENDDO
  ENDIF 
  IF ( PRESENT ( f_part ) ) THEN
    DO ig = 1, SIZE ( dg_part )  ! RHO0
      gt = dg_part ( ig ) % grid_index


! allocate intermediate coefficient arrays on the reference grid for
! the derivative of the potential due to rho0 and drho separately
      DO i = 1, 3
         CALL coeff_allocate ( dphi1_g ( i ), grid_b ( gt ), use_data = PW_COMPLEXDATA1D )
         dphi1_g ( i ) % pw % in_space = PW_RECIPROCALSPACE
 
         CALL coeff_allocate ( dphi2_g ( i ), grid_b ( gt ), use_data = PW_COMPLEXDATA1D )
         dphi2_g ( i ) % pw % in_space = PW_RECIPROCALSPACE
      END DO
      CALL coeff_allocate ( phi1_r, grid_b ( gt ), use_data = PW_COMPLEXDATA1D )
      phi1_r % pw % in_space = PW_RECIPROCALSPACE

      CALL calculate_hartree ( rho0_ref, vhartree = phi1_r, dvhartree = dphi1_g )

      CALL coeff_allocate ( phi2_r, grid_b ( gt ), use_data = PW_COMPLEXDATA1D )
      phi2_r % pw % in_space = PW_RECIPROCALSPACE

      CALL calculate_hartree ( drho_ref, vhartree = phi2_r, dvhartree = dphi2_g )

! set up real-space structure for the derivative of the potential
! on the reference grid
      CALL rs_grid_setup ( drpot1, grid_b ( gt ) , grid_s ( gt ) % npts )
      CALL rs_grid_setup ( drpot2, grid_b ( gt ) , grid_s ( gt ) % npts )

      CALL rs_grid_allocate ( drpot1 ( 1 : 3 ) )
      CALL rs_grid_allocate ( drpot2 ( 1 : 3 ) )

! reallocate coeff with REALDATA3D to use as a dummy variable
! to transfer dphi -> drpot
      CALL coeff_allocate ( coeff_dummy , grid_b ( gt ), use_data = PW_REALDATA3D )
      coeff_dummy % pw % in_space = PW_REALSPACE

      DO i = 1, 3
         CALL coeff_transform_space ( dphi1_g ( i ), coeff_dummy )
         CALL coeff_deallocate ( dphi1_g ( i ) )
         CALL rs_pw_transfer ( drpot1 ( i ), coeff_dummy % pw, "BACKWARD" )
      END DO

      DO i = 1, 3
         CALL coeff_transform_space ( dphi2_g ( i ), coeff_dummy )
         CALL coeff_deallocate ( dphi2_g ( i ) )
         CALL rs_pw_transfer ( drpot2 ( i ), coeff_dummy % pw, "BACKWARD" )
      END DO
      CALL coeff_deallocate ( coeff_dummy )

      CALL coeff_allocate ( rhos1,  grid_s ( gt ), use_data = PW_REALDATA3D )
      CALL coeff_allocate ( rhos2,  grid_s ( gt ), use_data = PW_REALDATA3D )

! force on coefs
      CALL force_part_evaluate ( dg_part ( ig ), drpot1, drpot2, part, &
         f_part, box, rhos1, rhos2, exp_igr_part ( ig ), grid_s ( gt ), &
         ewald_param )

      CALL coeff_deallocate ( rhos2 )
      CALL coeff_deallocate ( rhos1 )
      CALL rs_grid_deallocate ( drpot2 ( 1 : 3 ) )
      CALL rs_grid_deallocate ( drpot1 ( 1 : 3 ) )
      CALL structure_factor_deallocate ( exp_igr_part ( ig ) )
    ENDDO
  ENDIF

! clean up
  CALL coeff_deallocate ( drho_ref )
  CALL coeff_deallocate ( rho0_ref )

  DEALLOCATE ( exp_igr_part, STAT=allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "pol_electrostatics", &
   "exp_igr_part" )

  DEALLOCATE ( exp_igr_coef, STAT=allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "pol_electrostatics", &
   "exp_igr_coef" )

END SUBROUTINE electrostatics

!******************************************************************************

SUBROUTINE get_ref_density ( pw1, pw2 )
  IMPLICIT NONE 
  TYPE ( pw_type ), INTENT ( INOUT ), TARGET :: pw2
  TYPE ( pw_type ), INTENT ( IN ), TARGET :: pw1
! locals
  INTEGER :: l, m, n, ngpts, gpt, l0, m0, n0
  INTEGER, DIMENSION ( : ), POINTER :: mapl, mapm, mapn
  INTEGER, DIMENSION ( :, : ), POINTER :: ghat, bds
  INTEGER, PARAMETER :: IMPOSSIBLE = 10000

  mapl => pw2 % pw_grid % mapl % pos
  mapm => pw2 % pw_grid % mapm % pos
  mapn => pw2 % pw_grid % mapn % pos
  bds => pw2 % pw_grid % bounds

  ghat => pw1 % pw_grid % g_hat

  ngpts = SIZE ( pw1 % pw_grid % gsq ( : ) )

  IF ( -bds ( 1, 1 ) == bds ( 2, 1 ) ) THEN
    l0 = IMPOSSIBLE
  ELSE
    l0 = bds ( 1, 1 )
  END IF

  IF ( -bds ( 1, 2 ) == bds ( 2, 2 ) ) THEN
    m0 = IMPOSSIBLE
  ELSE
    m0 = bds ( 1, 2 )
  END IF

  IF ( -bds ( 1, 3 ) == bds ( 2, 3 ) ) THEN
    n0 = IMPOSSIBLE
  ELSE
    n0 = bds ( 1, 3 )
  END IF

  DO gpt = 1, ngpts

    l = mapl ( ghat ( 1, gpt ) ) + bds ( 1, 1 )
    m = mapm ( ghat ( 2, gpt ) ) + bds ( 1, 2 )
    n = mapn ( ghat ( 3, gpt ) ) + bds ( 1, 3 )
    pw2 % cc3d ( l, m, n ) = pw2 % cc3d ( l, m, n ) + pw1 % cc ( gpt ) 

    IF ( ghat ( 1, gpt ) == l0 .OR. ghat ( 2, gpt ) == m0 .OR. ghat ( 3, gpt ) == n0 ) THEN
      pw2 % cc3d ( l, m, n ) = 0.0_dbl
    END IF

  END DO


  IF ( pw1 % pw_grid % grid_span == HALFSPACE ) THEN

    mapl => pw2 % pw_grid % mapl % neg
    mapm => pw2 % pw_grid % mapm % neg
    mapn => pw2 % pw_grid % mapn % neg

    DO gpt = 1, ngpts

      l = mapl ( ghat ( 1, gpt ) ) + bds ( 1, 1 )
      m = mapm ( ghat ( 2, gpt ) ) + bds ( 1, 2 )
      n = mapn ( ghat ( 3, gpt ) ) + bds ( 1, 3 )
      pw2 % cc3d ( l, m, n ) = pw2 % cc3d ( l, m, n ) +  &
                               CONJG ( pw1 % cc ( gpt ) )

      IF ( ghat ( 1, gpt ) == l0 .OR. ghat ( 2, gpt ) == m0 .OR. ghat ( 3, gpt ) == n0 ) THEN
        pw2 % cc3d ( l, m, n ) = 0.0_dbl
      END IF

    END DO


  END IF

END SUBROUTINE get_ref_density 
                             
!******************************************************************************
!!****** pol_electrostatics/drho_evaluate [1.0] *
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

SUBROUTINE drho_evaluate ( dg, rden, part, coeff, box, exp_igr, grid_s,  &
                             grid_b, rhos1, rhos2 ) 

  IMPLICIT NONE

! Arguments
  TYPE ( dg_type ), INTENT ( IN ), TARGET :: dg
  TYPE ( particle_type ), DIMENSION ( : ), INTENT ( IN ) :: part
  TYPE ( coeff_type ), INTENT ( IN ) :: coeff
  TYPE ( cell_type ), INTENT ( IN ) :: box
  TYPE ( pw_grid_type ), INTENT ( IN ) :: grid_s
  TYPE ( pw_grid_type ), INTENT ( IN ) :: grid_b
  TYPE ( realspace_grid_type ), INTENT ( INOUT ) :: rden
  TYPE ( structure_factor_type ), INTENT ( OUT ) :: exp_igr
  TYPE ( coeff_type ), INTENT ( INOUT ) :: rhos1, rhos2

! Locals
  INTEGER :: handle, p1, p2, c1, c2, ip1, ip2
  TYPE ( coef_list_type ), POINTER :: clist
  REAL ( dbl ) :: r1 ( 3 ), r2 ( 3 ), q1, q2

!------------------------------------------------------------------------------


  CALL timeset ( 'DRHO_EVALUATE', 'I', ' ', handle )

! Initialize the index of p1 and p2. This is to ensure the
! correct mapping between the particle and its array location
! in exp_igr
  ip1 = 0
  ip2 = 0

  clist => dg % clist
  DO 
     CALL set_list_drho ( clist, p1, p2, ip1, ip2, c1, c2 )
     IF ( p1 == 0 .AND. p2 == 0 ) EXIT
     r1 = part ( p1 ) % r
     q1 = coeff % ao % cr ( c1 ) 
     r2 = 0._dbl
     q2 = 0._dbl
     IF ( p2 /= 0 ) THEN
       r2 = part ( p2 ) % r
       q2 = coeff % ao % cr ( c2 )
     ENDIF
     ! calculate function on small boxes (we use double packing in FFT)
     CALL get_patch ( dg, r1, r2, q1, q2, exp_igr, box, p1, p2, ip1, ip2, &
                      grid_b, grid_s, rhos1, rhos2 )

     ! add boxes to real space grid (big box)
     CALL dg_sum_patch ( rden, rhos1, exp_igr % centre ( :, ip1 ) )
     IF ( p2 /= 0 ) CALL dg_sum_patch ( rden, rhos2, &
                         exp_igr % centre ( :, ip2 ) )

     ip1 = ip2
  END DO

  CALL timestop ( 0.0_dbl, handle )

END SUBROUTINE drho_evaluate

!******************************************************************************
!!****** pol_electrostatics/rho0_evaluate [1.0] *
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
     IF ( p2 /= 0 ) CALL dg_sum_patch ( rden, rhos2, &
                         exp_igr % centre ( :, ip2 ) )

     ip1 = ip2
  END DO

  CALL timestop ( 0.0_dbl, handle )

END SUBROUTINE rho0_evaluate

!******************************************************************************
!!****** pol_electrostatics/force_part_evaluate [1.0] *
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
SUBROUTINE force_part_evaluate ( dg,  drpot1, drpot2, part, force, box, &
                                rhos1, rhos2, exp_igr, grid_s, ewald_param )

! Arguments
  TYPE ( dg_type ), INTENT ( IN ) :: dg
  TYPE ( particle_type ), DIMENSION ( : ), INTENT ( IN ) :: part
  TYPE ( cell_type ), INTENT ( IN ) :: box
  TYPE ( pw_grid_type ), INTENT ( IN ) :: grid_s
  TYPE ( realspace_grid_type ), DIMENSION ( : ), INTENT ( INOUT ) :: drpot1, drpot2
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

     ! sum boxes on real space grids (big box) for forces due to drpot2:
     CALL dg_sum_patch_force_3d ( drpot2, rhos1, exp_igr % centre ( :, ip1 ), fat )
     force ( 1, p1 ) = force ( 1, p1 ) - fat ( 1 ) * dvols
     force ( 2, p1 ) = force ( 2, p1 ) - fat ( 2 ) * dvols
     force ( 3, p1 ) = force ( 3, p1 ) - fat ( 3 ) * dvols
     IF ( p2 /= 0 ) THEN
        CALL dg_sum_patch_force_3d ( drpot2, rhos2, exp_igr % centre ( :, ip2 ), fat )
        force ( 1, p2 ) = force ( 1, p2 ) - fat ( 1 ) * dvols
        force ( 2, p2 ) = force ( 2, p2 ) - fat ( 2 ) * dvols
        force ( 3, p2 ) = force ( 3, p2 ) - fat ( 3 ) * dvols
     END IF

     ip1 = ip2
  END DO

  CALL timestop ( 0.0_dbl, handle )

END SUBROUTINE force_part_evaluate

!******************************************************************************
!!****** pol_electrostatics/force_coef_evaluate [1.0] *
!!
!!   NAME
!!     force_coef_evaluate
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
SUBROUTINE force_coef_evaluate ( dg, drpot1, drpot2, phi1_r, phi2_r, part, &
                coeff, fpart, fcoef, box, rhos1, rhos2, exp_igr, grid_s, ewald_param )

! Arguments
  TYPE ( dg_type ), INTENT ( IN ) :: dg
  TYPE ( coeff_type ), INTENT ( IN ) :: coeff
  TYPE ( particle_type ), DIMENSION ( : ), INTENT ( IN ) :: part
  REAL ( dbl ), INTENT ( INOUT ), DIMENSION ( :, : ) :: fpart
  REAL ( dbl ), INTENT ( INOUT ), DIMENSION ( : ) :: fcoef
  TYPE ( cell_type ), INTENT ( IN ) :: box
  TYPE ( pw_grid_type ), INTENT ( IN ) :: grid_s
  TYPE ( realspace_grid_type ), INTENT ( INOUT ) :: phi1_r, phi2_r
  TYPE ( realspace_grid_type ), DIMENSION ( : ), INTENT ( INOUT ) :: drpot1, drpot2
  TYPE ( structure_factor_type ), INTENT ( IN ) :: exp_igr
  TYPE ( coeff_type ), INTENT ( IN ) :: rhos1, rhos2
  TYPE ( ewald_parameters_type ), INTENT ( IN ) :: ewald_param

! Locals
  INTEGER :: handle, p1, p2, ip1, ip2, c1, c2
  TYPE ( coef_list_type ), POINTER :: clist
  REAL ( dbl ) :: fc, dvols, inv_fourpi_eps0, q1, q2
  REAL ( dbl ), DIMENSION (3) :: fat
  CALL timeset ( 'FORCE_EVALUATE_COEF', 'I', ' ', handle )

  inv_fourpi_eps0 = 1.0_dbl / ( fourpi * ewald_param % eps0 )

! initialize the forces
  dvols = rhos1 % pw % pw_grid % dvol * inv_fourpi_eps0

  clist => dg % clist
  ip1 = 0
  ip2 = 0
  DO

     CALL set_list_drho ( clist, p1, p2, ip1, ip2, c1, c2 )
     IF ( p1 == 0 .AND. p2 == 0 ) EXIT
     q1 =  coeff % ao % cr ( c1 )
     q2 = 0._dbl
     IF ( p2 /= 0 ) q2 =  coeff % ao % cr ( c2 )

     ! calculate function on small boxes (we use double packing in FFT)
     CALL get_patch_again ( dg, 1._dbl, 1._dbl, exp_igr, p1, p2, ip1, ip2, &
                            rhos1, rhos2 )

! forces on the ions

     ! sum boxes on real space grids (big box) due to drpot1
     CALL dg_sum_patch_force_3d ( drpot1, rhos1, exp_igr % centre ( :, ip1 ), fat )
     fpart ( 1, p1 ) = fpart ( 1, p1 ) - q1 * fat ( 1 ) * dvols
     fpart ( 2, p1 ) = fpart ( 2, p1 ) - q1 * fat ( 2 ) * dvols
     fpart ( 3, p1 ) = fpart ( 3, p1 ) - q1 * fat ( 3 ) * dvols
     IF ( p2 /= 0 ) THEN
        CALL dg_sum_patch_force_3d ( drpot1, rhos2, exp_igr % centre ( :, ip2 ), fat )
        fpart ( 1, p2 ) = fpart ( 1, p2 ) - q2 * fat ( 1 ) * dvols
        fpart ( 2, p2 ) = fpart ( 2, p2 ) - q2 * fat ( 2 ) * dvols
        fpart ( 3, p2 ) = fpart ( 3, p2 ) - q2 * fat ( 3 ) * dvols
     END IF

     ! sum boxes on real space grids (big box) due to drpot2
     CALL dg_sum_patch_force_3d ( drpot2, rhos1, exp_igr % centre ( :, ip1 ), fat )
     fpart ( 1, p1 ) = fpart ( 1, p1 ) - q1 * fat ( 1 ) * dvols
     fpart ( 2, p1 ) = fpart ( 2, p1 ) - q1 * fat ( 2 ) * dvols
     fpart ( 3, p1 ) = fpart ( 3, p1 ) - q1 * fat ( 3 ) * dvols
     IF ( p2 /= 0 ) THEN
        CALL dg_sum_patch_force_3d ( drpot2, rhos2, exp_igr % centre ( :, ip2 ), fat )
        fpart ( 1, p2 ) = fpart ( 1, p2 ) - q2 * fat ( 1 ) * dvols
        fpart ( 2, p2 ) = fpart ( 2, p2 ) - q2 * fat ( 2 ) * dvols
        fpart ( 3, p2 ) = fpart ( 3, p2 ) - q2 * fat ( 3 ) * dvols
     END IF

! forces on the coeffs

     ! sum boxes on real space grids (big box) due to phi1_r
     
     CALL dg_sum_patch_force_1d ( phi1_r, rhos1, exp_igr % centre ( :, ip1 ), fc )
     fcoef ( c1 ) = fcoef ( c1 ) - fc * dvols
     IF ( p2 /= 0 ) THEN
        CALL dg_sum_patch_force_1d ( phi1_r, rhos2, exp_igr % centre ( :, ip2 ), fc )
        fcoef ( c2 ) = fcoef ( c2 ) - fc * dvols
     END IF

     ! sum boxes on real space grids (big box) due to phi2_r
     
     CALL dg_sum_patch_force_1d ( phi2_r, rhos1, exp_igr % centre ( :, ip1 ), fc )
     fcoef ( c1 ) = fcoef ( c1 ) - fc * dvols
     IF ( p2 /= 0 ) THEN
        CALL dg_sum_patch_force_1d ( phi2_r, rhos2, exp_igr % centre ( :, ip2 ), fc )
        fcoef ( c2 ) = fcoef ( c2 ) - fc * dvols
     END IF

     ip1 = ip2 
  END DO 
  CALL timestop ( 0.0_dbl, handle )

END SUBROUTINE force_coef_evaluate

!******************************************************************************
!!****** pol_electrostatics/set_list_drho [1.0] *
!!
!!   NAME
!!     set_list_drho
!!
!!   FUNCTION
!!     Radvances in the link list of particles to the next to be calculated
!!
!!   AUTHOR
!!     CJM (18-April-2001)
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE set_list_drho ( clist, p1, p2, ip1, ip2, c1, c2 )

  IMPLICIT NONE

! Arguments
  TYPE ( coef_list_type ),  POINTER :: clist
  INTEGER, INTENT ( OUT ) :: p1, p2, c1, c2
  INTEGER, INTENT ( INOUT ) :: ip1, ip2

!------------------------------------------------------------------------------
  p1 = 0
  p2 = 0
  c1 = 0
  c2 = 0

  IF ( ASSOCIATED ( clist ) ) THEN
    ip1 = ip1 + 1
    ip2 = ip1
    p1 = clist % particle_index
    c1 = clist % coef_index
    IF ( ASSOCIATED ( clist % next ) ) THEN
      ip2 = ip2 + 1
      clist => clist % next 
      p2 = clist % particle_index
      c2 = clist % coef_index
    END IF
    clist => clist % next 
  END IF  
 

END SUBROUTINE set_list_drho

!******************************************************************************
!!****** pol_electrostatics/set_list_rho0 [1.0] *
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

END MODULE pol_electrostatics

!******************************************************************************
