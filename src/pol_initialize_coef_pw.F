!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/pol_initialize_coef_pw [1.0] *
!!
!!   NAME
!!     pol_initialize_coef_pw
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     gt
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE pol_initialize_coef_pw

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
                                             dg_sum_patch
  USE ewald_parameters_types,          ONLY: ewald_parameters_type
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_types,                        ONLY: pw_add,&
                                             pw_copy,&
                                             pw_type
  USE realspace_grid_types,            ONLY: realspace_grid_type,&
                                             rs_grid_allocate,&
                                             rs_grid_deallocate,&
                                             rs_grid_setup,&
                                             rs_pw_transfer
  USE simulation_cell,                 ONLY: cell_type
  USE structure_factor_types,          ONLY: structure_factor_type
  USE structure_factors,               ONLY: structure_factor_allocate
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: initialize_coef_pw

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** pol_initialize_coef_pw/initialize_coef_pw [1.0] *
!!
!!   NAME
!!     initialize_coef_pw
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

SUBROUTINE initialize_coef_pw ( dg_proj, pw, part, box, grid_s,  &
                           grid_b, ewald_param )


    TYPE(dg_type), DIMENSION(:), INTENT(IN)  :: dg_proj
    TYPE(pw_type), INTENT(INOUT)             :: pw
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: part
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(pw_grid_type), DIMENSION(:), &
      INTENT(IN)                             :: grid_s, grid_b
    TYPE(ewald_parameters_type), INTENT(IN)  :: ewald_param

    INTEGER                                  :: allocstat, gt, ig, igtot, &
                                                ngrids, ngt, nproj
    LOGICAL                                  :: first_time, match
    TYPE(coeff_type)                         :: coeff_dummy, drhob_g, &
                                                drhob_r, rhos1, rhos2
    TYPE(realspace_grid_type)                :: drho_rs
    TYPE(structure_factor_type), &
      ALLOCATABLE, DIMENSION(:)              :: exp_igr_proj

!------------------------------------------------------------------------------

 ngrids = SIZE ( dg_proj )

 ALLOCATE ( exp_igr_proj ( ngrids ), STAT=allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "pol_electrostatics_pw", &
   "exp_igr_proj", ngrids )

 DO ig = 1, ngrids

    nproj = dg_proj ( ig ) % ncoefs

    gt = dg_proj ( ig )% grid_index

    CALL structure_factor_allocate ( grid_s ( gt ) % bounds, nproj,  &
          exp_igr_proj ( ig ), allocate_centre = .TRUE. )
  ENDDO

  CALL coeff_allocate ( coeff_dummy , pw % pw_grid, use_data = PW_COMPLEXDATA1D )
  coeff_dummy % pw % in_space = PW_RECIPROCALSPACE

  CALL coeff_zero ( coeff_dummy )

!-------------- DENSITY CALCULATION ---------------------

  ngt = SIZE ( grid_b )

!initialize the pws
  DO gt = 1, ngt
    first_time = .TRUE.
    igtot = SIZE ( dg_proj )
    IF ( igtot == 0 ) EXIT
    match = .FALSE.
    DO ig = 1, igtot   ! DRHO in AO basis

      IF ( gt == dg_proj ( ig ) % grid_index ) THEN
        match = .TRUE.
       IF ( first_time ) THEN
! Allocate small coefficient types
          CALL coeff_allocate ( rhos1,  grid_s ( gt ), use_data = PW_REALDATA3D )
          CALL coeff_allocate ( rhos2,  grid_s ( gt ), use_data = PW_REALDATA3D )
! Allocate large real-space structure based on grid_b
          CALL rs_grid_setup ( drho_rs, grid_b ( gt ), -1 )
          CALL rs_grid_allocate ( drho_rs )
          drho_rs % r = 0.0_dp
          first_time = .FALSE.
        ENDIF
        IF (  drho_rs % parallel .AND. drho_rs % group_dim ( 1 ) /= 1 )  THEN
          CALL get_center ( part, box, exp_igr_proj ( ig ) % centre, &
                        grid_b ( gt ) % npts )
        END IF
! Get density drho
        CALL drho_evaluate ( dg_proj ( ig ), drho_rs, part,  box,  &
                  exp_igr_proj ( ig ), grid_s ( gt ), grid_b ( gt ), rhos1, rhos2 )
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
      CALL pw_copy ( pw, coeff_dummy % pw )

! Sum to the reference grid
      CALL pw_add ( drhob_g % pw, coeff_dummy % pw,  pw )
! do a little cleaning
      CALL coeff_deallocate ( drhob_g )
      CALL coeff_deallocate ( drhob_r )
      CALL rs_grid_deallocate ( drho_rs )
      CALL coeff_deallocate ( rhos2 )
      CALL coeff_deallocate ( rhos1 )
    ENDIF
  ENDDO

 DEALLOCATE ( exp_igr_proj, STAT=allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "pol_initialize_coef_pw", &
   "exp_igr_proj" )

END SUBROUTINE initialize_coef_pw
!******************************************************************************

SUBROUTINE drho_evaluate ( dg, rden, part, box, exp_igr, grid_s,  &
                             grid_b, rhos1, rhos2 )


    TYPE(dg_type), INTENT(IN), TARGET        :: dg
    TYPE(realspace_grid_type), INTENT(INOUT) :: rden
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: part
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(structure_factor_type), INTENT(OUT) :: exp_igr
    TYPE(pw_grid_type), INTENT(IN)           :: grid_s, grid_b
    TYPE(coeff_type), INTENT(INOUT)          :: rhos1, rhos2

    INTEGER                                  :: c1, c2, handle, ip1, ip2, p1, &
                                                p2
    REAL(KIND=dp)                                :: q1, q2, r1( 3 ), r2( 3 )
    TYPE(coef_list_type), POINTER            :: clist

!-----------------------------------------------------------------

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
     q1 = 1.0_dp
     r2 = 0.0_dp
     q2 = 0.0_dp
     IF ( p2 /= 0 ) THEN
       r2 = part ( p2 ) % r
       q2 = 1.0_dp
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

  CALL timestop ( 0.0_dp, handle )

END SUBROUTINE drho_evaluate
!******************************************************************************

SUBROUTINE get_patch ( dg, r1, r2, q1, q2, exp_igr, box, p1, p2, &
           ip1, ip2, grid_b, grid_s, rhos1, rhos2 )


    TYPE(dg_type), INTENT(IN), TARGET        :: dg
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)      :: r1, r2
    REAL(KIND=dp), INTENT(IN)                    :: q1, q2
    TYPE(structure_factor_type)              :: exp_igr
    TYPE(cell_type), INTENT(IN)              :: box
    INTEGER, INTENT(IN)                      :: p1, p2, ip1, ip2
    TYPE(pw_grid_type), INTENT(IN)           :: grid_b, grid_s
    TYPE(coeff_type)                         :: rhos1, rhos2

    COMPLEX(KIND=dp), DIMENSION(:), POINTER      :: ex1, ex2, ey1, ey2, ez1, ez2
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
!*****
!******************************************************************************

SUBROUTINE get_center ( part, box, centre, npts )


    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: part
    TYPE(cell_type), INTENT(IN)              :: box
    INTEGER, DIMENSION(:, :), INTENT(OUT)    :: centre
    INTEGER, DIMENSION(:), INTENT(IN)        :: npts

    INTEGER                                  :: ipart
    REAL(KIND=dp), DIMENSION(3)                  :: s

!------------------------------------------------------------------------------

  DO ipart = 1, SIZE ( part )
    ! compute the scaled coordinate of atomi
    s = MATMUL ( box % h_inv, part ( ipart ) % r )
    s = s - NINT ( s )
    ! find the closest grid point (on big grid)
    centre ( :, ipart ) = NINT ( REAL ( npts ( 1:3 ),KIND=dp) * s ( 1:3 ) )
  END DO

END SUBROUTINE get_center
!******************************************************************************

END MODULE pol_initialize_coef_pw

!******************************************************************************
