!-----------------------------------------------------------------------------
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/pol_electrostatics_ao [1.0] *
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

  USE coefficient_lists,               ONLY: coef_list_type
  USE coefficient_types,               ONLY: PW_COMPLEXDATA1D,&
                                             PW_COMPLEXDATA3D,&
                                             PW_REALDATA3D,&
                                             PW_REALSPACE,&
                                             PW_RECIPROCALSPACE,&
                                             coeff_allocate,&
                                             coeff_deallocate,&
                                             coeff_transform_space,&
                                             coeff_type,&
                                             coeff_zero
  USE dg_types,                        ONLY: dg_type
  USE dgs,                             ONLY: dg_get_patch,&
                                             dg_get_strucfac,&
                                             dg_sum_patch,&
                                             dg_sum_patch_force_1d,&
                                             dg_sum_patch_force_3d
  USE ewald_parameters_types,          ONLY: ewald_parameters_type
  USE hartree,                         ONLY: calculate_hartree
  USE kinds,                           ONLY: dbl
  USE mathconstants,                   ONLY: fourpi
  USE particle_lists,                  ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_types,                        ONLY: pw_add,&
                                             pw_copy
  USE pws,                             ONLY: init_pw_poisson_solver,&
                                             pw_integral_ab
  USE realspace_grid_types,            ONLY: realspace_grid_type,&
                                             rs_grid_allocate,&
                                             rs_grid_deallocate,&
                                             rs_grid_setup,&
                                             rs_pw_transfer
  USE simulation_cell,                 ONLY: cell_type
  USE structure_factor_types,          ONLY: structure_factor_type
  USE structure_factors,               ONLY: structure_factor_allocate,&
                                             structure_factor_deallocate
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: electrostatics

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** pol_electrostatics_ao/electrostatics [1.0] *
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


    TYPE(dg_type), DIMENSION(:), INTENT(IN)  :: dg_part, dg_coef
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: part
    TYPE(coeff_type), INTENT(IN)             :: coeff
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(pw_grid_type), DIMENSION(:), &
      INTENT(IN)                             :: grid_s, grid_b
    TYPE(ewald_parameters_type), INTENT(IN)  :: ewald_param
    INTEGER, INTENT(IN)                      :: iref
    REAL(dbl), INTENT(OUT)                   :: energy
    REAL(dbl), DIMENSION(:), INTENT(OUT), &
      OPTIONAL                               :: f_coef
    REAL(dbl), DIMENSION(:, :), &
      INTENT(OUT), OPTIONAL                  :: f_part

    INTEGER                                  :: allocstat, gt, i, ig, igtot, &
                                                ncoef, ngrids, ngt, npart
    LOGICAL                                  :: first_time, match
    REAL(dbl)                                :: e_drho_drho, e_drho_rho0, &
                                                e_rho0_rho0, inv_fourpi_eps0
    TYPE(coeff_type) :: coeff_dummy, dphi1_g( 3 ), dphi2_g( 3 ), drho_ref, &
      drhob_g, drhob_r, phi1_r, phi2_r, rho0_ref, rho0b_g, rho0b_r, rhos1, &
      rhos2
    TYPE(realspace_grid_type)                :: drho_rs, drpot1( 3 ), &
                                                drpot2( 3 ), phi1_rs, &
                                                phi2_rs, rho0_rs
    TYPE(structure_factor_type), &
      ALLOCATABLE, DIMENSION(:)              :: exp_igr_coef, exp_igr_part

!------------------------------------------------------------------------------
!-------------- INITIALISATION ---------------------
! allocating the structure factor vector to the number of grid_types

  ngrids = SIZE ( dg_coef )
  ALLOCATE ( exp_igr_coef ( ngrids ), STAT=allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "pol_electrostatics_ao", &
   "exp_igr_coef", ngrids )

  ngrids = SIZE ( dg_part )
  ALLOCATE ( exp_igr_part ( ngrids ), STAT=allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "pol_electrostatics_ao", &
   "exp_igr_part", ngrids )
  DO ig = 1, SIZE ( dg_coef )

    ncoef = dg_coef ( ig ) % ncoefs

    gt = dg_coef ( ig )% grid_index

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

  ngt = SIZE ( grid_b )

! coefficients
  DO gt = 1, ngt
    first_time = .TRUE.
    igtot = SIZE ( dg_coef )
    IF ( igtot == 0 ) EXIT
    match = .FALSE.
    DO ig = 1, igtot
      IF ( gt == dg_coef ( ig ) % grid_index ) THEN
        match = .TRUE.
        IF ( first_time ) THEN
! Allocate small coefficient types
          CALL coeff_allocate ( rhos1,  grid_s ( gt ), use_data = PW_REALDATA3D )
          CALL coeff_allocate ( rhos2,  grid_s ( gt ), use_data = PW_REALDATA3D )

! Allocate large real-space structure based on grid_b
          CALL rs_grid_setup ( drho_rs, grid_b ( gt ), -1 )
          CALL rs_grid_allocate ( drho_rs )
          drho_rs % r = 0._dbl
          first_time = .FALSE.
        ENDIF
        IF (  drho_rs % parallel .AND. drho_rs % group_dim ( 1 ) /= 1 )  THEN
          CALL get_center ( part, box, exp_igr_coef ( ig ) % centre, &
                        grid_b ( gt ) % npts )
        END IF
! Get density drho
        CALL drho_evaluate ( dg_coef ( ig ), drho_rs, part, coeff, box,  &
                  exp_igr_coef ( ig ), grid_s ( gt ), grid_b ( gt ), rhos1, rhos2 ) 
      ELSE
        CYCLE
      ENDIF
    ENDDO
! Allocate large coefficient structure based on grid_b  
    IF ( match ) THEN
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
    ENDIF
  ENDDO

! particles
  DO gt = 1, ngt
    first_time = .TRUE.
    igtot = SIZE ( dg_part )
    IF ( igtot == 0 ) EXIT
    match = .FALSE.
    DO ig = 1, igtot
      IF ( gt == dg_part ( ig ) % grid_index ) THEN
        match = .TRUE.
        IF ( first_time ) THEN
! Allocate small coefficient types
          CALL coeff_allocate ( rhos1,  grid_s ( gt ), use_data = PW_REALDATA3D )
          CALL coeff_allocate ( rhos2,  grid_s ( gt ), use_data = PW_REALDATA3D )

! Allocate large real-space structure based on grid_b
          CALL rs_grid_setup ( rho0_rs, grid_b ( gt ), -1 )
          CALL rs_grid_allocate ( rho0_rs )
          rho0_rs % r = 0._dbl
          first_time = .FALSE.
        ENDIF
        IF (  rho0_rs % parallel .AND. rho0_rs % group_dim ( 1 ) /= 1 )  THEN
          CALL get_center ( part, box, exp_igr_part ( ig ) % centre, &
                        grid_b ( gt ) % npts )
        END IF
! Get density drho
        CALL rho0_evaluate ( dg_part ( ig ), rho0_rs, part, box, exp_igr_part ( ig ), &
                             grid_s ( gt ), grid_b ( gt ), rhos1, rhos2 ) 
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

! Copy present rho_ref to coeff_dummy
      CALL pw_copy ( rho0_ref % pw, coeff_dummy % pw )

! Sum to the reference grid
      CALL pw_add ( rho0b_g % pw, coeff_dummy % pw, rho0_ref % pw )
! do a little cleaning
      CALL coeff_deallocate ( rho0b_g )
      CALL coeff_deallocate ( rho0b_r )
      CALL rs_grid_deallocate ( rho0_rs )
      CALL coeff_deallocate ( rhos2 )
      CALL coeff_deallocate ( rhos1 )
    ENDIF
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

  energy = ( e_rho0_rho0 + e_drho_drho + e_drho_rho0 ) * inv_fourpi_eps0

  CALL coeff_deallocate ( phi2_r )
  CALL coeff_deallocate ( phi1_r )


!-------------- ELECTROSTATIC CALCULATION DONE -----------


!-------------- FORCE CALCULATION ---------------------

  IF ( PRESENT ( f_coef ) )  THEN 
    DO gt = 1, ngt
      first_time = .TRUE.
      igtot = SIZE ( dg_coef )
      IF ( igtot == 0 ) EXIT
      match = .FALSE.
      DO ig = 1, igtot
        IF ( gt == dg_coef ( ig ) % grid_index ) THEN
          match = .TRUE. 
          IF ( first_time ) THEN
! set up real-space structure for the derivative of the potential
! on the reference grid
            CALL rs_grid_setup ( drpot1, grid_b ( gt ) , -1 )
            CALL rs_grid_setup ( drpot2, grid_b ( gt ) , -1 )

            CALL rs_grid_allocate ( drpot1 ( 1 : 3 ) )
            CALL rs_grid_allocate ( drpot2 ( 1 : 3 ) )

! Create real_space structure of phi1_r and phi2_r to use in the
! force calculation sum_patch_force only acts on real_space types

            CALL rs_grid_setup ( phi1_rs, grid_b ( gt ), -1 )
            CALL rs_grid_setup ( phi2_rs, grid_b ( gt ), -1 )

            CALL rs_grid_allocate ( phi1_rs )
            CALL rs_grid_allocate ( phi2_rs )

            CALL coeff_allocate ( phi1_r, grid_b ( gt ), use_data = PW_COMPLEXDATA1D )
            phi1_r % pw % in_space = PW_RECIPROCALSPACE

            CALL coeff_allocate ( phi2_r, grid_b ( gt ), use_data = PW_COMPLEXDATA1D )
            phi2_r % pw % in_space = PW_RECIPROCALSPACE

! reallocate coeff with REALDATA3D to use as a dummy variable
! to transfer dphi -> drpot and phi_r -> phi_rs
            CALL coeff_allocate ( coeff_dummy , grid_b ( gt ), use_data = PW_REALDATA3D )
            coeff_dummy % pw % in_space = PW_REALSPACE

! allocate intermediate coefficient arrays on the reference grid for
! the derivative of the potential due to rho0 and drho separately
            DO i = 1, 3
               CALL coeff_allocate ( dphi1_g ( i ), grid_b ( gt ), use_data = PW_COMPLEXDATA1D )
               dphi1_g ( i ) % pw % in_space = PW_RECIPROCALSPACE
 
               CALL coeff_allocate ( dphi2_g ( i ), grid_b ( gt ), use_data = PW_COMPLEXDATA1D )
               dphi2_g ( i ) % pw % in_space = PW_RECIPROCALSPACE
            END DO

            CALL calculate_hartree ( rho0_ref, vhartree = phi1_r, dvhartree = dphi1_g )

            CALL calculate_hartree ( drho_ref, vhartree = phi2_r, dvhartree = dphi2_g )

            CALL coeff_transform_space ( phi1_r, coeff_dummy )
            CALL rs_pw_transfer ( phi1_rs, coeff_dummy % pw, "BACKWARD" )

            CALL coeff_transform_space ( phi2_r, coeff_dummy )
            CALL rs_pw_transfer ( phi2_rs, coeff_dummy % pw, "BACKWARD" )

            DO i = 1, 3
               CALL coeff_transform_space ( dphi2_g ( i ), coeff_dummy )
               CALL coeff_deallocate ( dphi2_g ( i ) )
               CALL rs_pw_transfer ( drpot2 ( i ), coeff_dummy % pw, "BACKWARD" )
            END DO

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
! force on coefs and ions
          CALL force_coef_evaluate ( dg_coef ( ig ), drpot1, drpot2, phi1_rs,  &
          phi2_rs, part, coeff, f_part, f_coef, box, rhos1, rhos2,  &
          exp_igr_coef ( ig ), grid_s ( gt ), ewald_param )
          CALL structure_factor_deallocate ( exp_igr_coef ( ig ) )
        ELSE
         CYCLE
        ENDIF
      ENDDO
! clean up
      IF ( match ) THEN
        CALL coeff_deallocate ( rhos2 )
        CALL coeff_deallocate ( rhos1 )
        CALL coeff_deallocate ( phi2_r )
        CALL coeff_deallocate ( phi1_r )
        CALL rs_grid_deallocate ( phi2_rs )
        CALL rs_grid_deallocate ( phi1_rs )
        CALL rs_grid_deallocate ( drpot2 ( 1 : 3 ) )
        CALL rs_grid_deallocate ( drpot1 ( 1 : 3 ) )
      ENDIF
    ENDDO
  ENDIF 
 
  IF ( PRESENT ( f_part ) ) THEN
    DO gt = 1, ngt
      first_time = .TRUE.
      igtot = SIZE ( dg_part )
      IF ( igtot == 0 ) EXIT
      match = .FALSE. 
      DO ig = 1, igtot 
        IF ( gt == dg_part ( ig ) % grid_index ) THEN
          match = .TRUE.
          IF ( first_time ) THEN

! set up real-space structure for the derivative of the potential
! on the reference grid
             CALL rs_grid_setup ( drpot1, grid_b ( gt ) , -1 )
             CALL rs_grid_setup ( drpot2, grid_b ( gt ) , -1 )

             CALL rs_grid_allocate ( drpot1 ( 1 : 3 ) )
             CALL rs_grid_allocate ( drpot2 ( 1 : 3 ) )

             CALL coeff_allocate ( phi1_r, grid_b ( gt ), use_data = PW_COMPLEXDATA1D )
             phi1_r % pw % in_space = PW_RECIPROCALSPACE

             CALL coeff_allocate ( phi2_r, grid_b ( gt ), use_data = PW_COMPLEXDATA1D )
             phi2_r % pw % in_space = PW_RECIPROCALSPACE

! reallocate coeff with REALDATA3D to use as a dummy variable
! to transfer dphi -> drpot
             CALL coeff_allocate ( coeff_dummy , grid_b ( gt ), use_data = PW_REALDATA3D )
             coeff_dummy % pw % in_space = PW_REALSPACE

! allocate intermediate coefficient arrays on the reference grid for
! the derivative of the potential due to rho0 and drho separately
             DO i = 1, 3
                CALL coeff_allocate ( dphi1_g ( i ), grid_b ( gt ), use_data = PW_COMPLEXDATA1D )
                dphi1_g ( i ) % pw % in_space = PW_RECIPROCALSPACE
 
                CALL coeff_allocate ( dphi2_g ( i ), grid_b ( gt ), use_data = PW_COMPLEXDATA1D )
                dphi2_g ( i ) % pw % in_space = PW_RECIPROCALSPACE
             END DO

             CALL calculate_hartree ( rho0_ref, vhartree = phi1_r, dvhartree = dphi1_g )

             CALL calculate_hartree ( drho_ref, vhartree = phi2_r, dvhartree = dphi2_g )

             DO i = 1, 3
                CALL coeff_transform_space ( dphi2_g ( i ), coeff_dummy )
                CALL coeff_deallocate ( dphi2_g ( i ) )
                CALL rs_pw_transfer ( drpot2 ( i ), coeff_dummy % pw, "BACKWARD" )
             END DO

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

! force on coefs
          CALL force_part_evaluate ( dg_part ( ig ), drpot1, drpot2, part, &
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
        CALL coeff_deallocate ( phi2_r ) 
        CALL coeff_deallocate ( phi1_r ) 
        CALL rs_grid_deallocate ( drpot2 ( 1 : 3 ) )
        CALL rs_grid_deallocate ( drpot1 ( 1 : 3 ) )
      ENDIF
    ENDDO
  ENDIF

! clean up
  CALL coeff_deallocate ( drho_ref )
  CALL coeff_deallocate ( rho0_ref )

  DEALLOCATE ( exp_igr_part, STAT=allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "pol_electrostatics_ao", &
   "exp_igr_part" )

  DEALLOCATE ( exp_igr_coef, STAT=allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "pol_electrostatics_ao", &
   "exp_igr_coef" )

END SUBROUTINE electrostatics

!******************************************************************************
!!****** pol_electrostatics_ao/drho_evaluate [1.0] *
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


    TYPE(dg_type), INTENT(IN), TARGET        :: dg
    TYPE(realspace_grid_type), INTENT(INOUT) :: rden
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: part
    TYPE(coeff_type), INTENT(IN)             :: coeff
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(structure_factor_type), INTENT(OUT) :: exp_igr
    TYPE(pw_grid_type), INTENT(IN)           :: grid_s, grid_b
    TYPE(coeff_type), INTENT(INOUT)          :: rhos1, rhos2

    INTEGER                                  :: c1, c2, handle, ip1, ip2, p1, &
                                                p2
    REAL(dbl)                                :: q1, q2, r1( 3 ), r2( 3 )
    TYPE(coef_list_type), POINTER            :: clist

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
!!****** pol_electrostatics_ao/rho0_evaluate [1.0] *
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


    TYPE(dg_type), INTENT(IN), TARGET        :: dg
    TYPE(realspace_grid_type), INTENT(INOUT) :: rden
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: part
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(structure_factor_type), INTENT(OUT) :: exp_igr
    TYPE(pw_grid_type), INTENT(IN)           :: grid_s, grid_b
    TYPE(coeff_type), INTENT(INOUT)          :: rhos1, rhos2

    INTEGER                                  :: handle, ip1, ip2, p1, p2
    REAL(dbl)                                :: q1, q2, r1( 3 ), r2( 3 )
    TYPE(particle_list_type), POINTER        :: plist

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
!!****** pol_electrostatics_ao/force_part_evaluate [1.0] *
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

    TYPE(dg_type), INTENT(IN)                :: dg
    TYPE(realspace_grid_type), &
      DIMENSION(:), INTENT(INOUT)            :: drpot1, drpot2
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: part
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: force
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(coeff_type), INTENT(INOUT)          :: rhos1, rhos2
    TYPE(structure_factor_type), INTENT(IN)  :: exp_igr
    TYPE(pw_grid_type), INTENT(IN)           :: grid_s
    TYPE(ewald_parameters_type), INTENT(IN)  :: ewald_param

    INTEGER                                  :: handle, ip1, ip2, p1, p2
    REAL(dbl)                                :: dvols, inv_fourpi_eps0, q1, q2
    REAL(dbl), DIMENSION(3)                  :: fat
    TYPE(particle_list_type), POINTER        :: plist

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
!!****** pol_electrostatics_ao/force_coef_evaluate [1.0] *
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

    TYPE(dg_type), INTENT(IN)                :: dg
    TYPE(realspace_grid_type), &
      DIMENSION(:), INTENT(INOUT)            :: drpot1, drpot2
    TYPE(realspace_grid_type), INTENT(INOUT) :: phi1_r, phi2_r
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: part
    TYPE(coeff_type), INTENT(IN)             :: coeff
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: fpart
    REAL(dbl), DIMENSION(:), INTENT(INOUT)   :: fcoef
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(coeff_type), INTENT(IN)             :: rhos1, rhos2
    TYPE(structure_factor_type), INTENT(IN)  :: exp_igr
    TYPE(pw_grid_type), INTENT(IN)           :: grid_s
    TYPE(ewald_parameters_type), INTENT(IN)  :: ewald_param

    INTEGER                                  :: c1, c2, handle, ip1, ip2, p1, &
                                                p2
    REAL(dbl)                                :: dvols, fc, inv_fourpi_eps0, &
                                                q1, q2
    REAL(dbl), DIMENSION(3)                  :: fat
    TYPE(coef_list_type), POINTER            :: clist

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
!!****** pol_electrostatics_ao/set_list_drho [1.0] *
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


    TYPE(coef_list_type), POINTER            :: clist
    INTEGER, INTENT(OUT)                     :: p1, p2
    INTEGER, INTENT(INOUT)                   :: ip1, ip2
    INTEGER, INTENT(OUT)                     :: c1, c2

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
!!****** pol_electrostatics_ao/set_list_rho0 [1.0] *
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


    TYPE(particle_list_type), POINTER        :: plist
    INTEGER, INTENT(OUT)                     :: p1, p2
    INTEGER, INTENT(INOUT)                   :: ip1, ip2

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


    TYPE(dg_type), INTENT(IN), TARGET        :: dg
    REAL(dbl), DIMENSION(3), INTENT(IN)      :: r1, r2
    REAL(dbl), INTENT(IN)                    :: q1, q2
    TYPE(structure_factor_type)              :: exp_igr
    TYPE(cell_type), INTENT(IN)              :: box
    INTEGER, INTENT(IN)                      :: p1, p2, ip1, ip2
    TYPE(pw_grid_type), INTENT(IN)           :: grid_b, grid_s
    TYPE(coeff_type)                         :: rhos1, rhos2

    COMPLEX(dbl), DIMENSION(:), POINTER      :: ex1, ex2, ey1, ey2, ez1, ez2
    INTEGER, DIMENSION(:), POINTER           :: center1, center2
    TYPE(coeff_type), POINTER                :: rho0

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


    TYPE(dg_type), INTENT(IN), TARGET        :: dg
    REAL(dbl), INTENT(IN)                    :: q1, q2
    TYPE(structure_factor_type)              :: exp_igr
    INTEGER, INTENT(IN)                      :: p1, p2, ip1, ip2
    TYPE(coeff_type)                         :: rhos1, rhos2

    COMPLEX(dbl), DIMENSION(:), POINTER      :: ex1, ex2, ey1, ey2, ez1, ez2
    TYPE(coeff_type), POINTER                :: rho0

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


    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: part
    TYPE(cell_type), INTENT(IN)              :: box
    INTEGER, DIMENSION(:, :), INTENT(OUT)    :: centre
    INTEGER, DIMENSION(:), INTENT(IN)        :: npts

    INTEGER                                  :: ipart
    REAL(dbl), DIMENSION(3)                  :: s

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
