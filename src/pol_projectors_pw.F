!-----------------------------------------------------------------------------
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/pol_projectors_pw [1.0] *
!!
!!   NAME
!!     pol_projectors_pw
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

MODULE pol_projectors_pw

  USE ao_types,                        ONLY: ao_type
  USE coefficient_lists,               ONLY: coef_list_type
  USE coefficient_types,               ONLY: PW_COMPLEXDATA1D,&
                                             PW_COMPLEXDATA3D,&
                                             PW_REALDATA3D,&
                                             PW_REALSPACE,&
                                             PW_RECIPROCALSPACE,&
                                             coeff_allocate,&
                                             coeff_deallocate,&
                                             coeff_transform_space,&
                                             coeff_type
  USE dg_types,                        ONLY: dg_type
  USE dgs,                             ONLY: dg_get_patch,&
                                             dg_get_strucfac,&
                                             dg_sum_patch_force_1d,&
                                             dg_sum_patch_force_3d
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_types,                        ONLY: pw_copy,&
                                             pw_derive,&
                                             pw_type
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
  
  TYPE  pol_projector_type
    INTEGER :: ipart
    REAL (KIND=dp), DIMENSION ( 3 ) :: force
  END TYPE pol_projector_type

  PRIVATE
  PUBLIC :: projectors, pol_projector_type

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** pol_projector_pw/projectors [1.0] *
!!
!!   NAME
!!     projectors
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

SUBROUTINE projectors ( dg_proj, ao, pw, part, box, grid_s,  &
                        grid_b, fproj_part )


    TYPE(dg_type), DIMENSION(:), INTENT(IN)  :: dg_proj
    TYPE(ao_type), INTENT(INOUT)             :: ao
    TYPE(pw_type), INTENT(IN)                :: pw
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: part
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(pw_grid_type), DIMENSION(:), &
      INTENT(IN)                             :: grid_s, grid_b
    TYPE(pol_projector_type), OPTIONAL       :: fproj_part( : )

    INTEGER                                  :: allocstat, gt, i, ig, igtot, &
                                                n( 3 ), ngrids, ngt, nproj
    LOGICAL                                  :: first_time, match
    TYPE(coeff_type)                         :: coeff_dummy, dphi1_g( 3 ), &
                                                drhob_r, rhos1, rhos2
    TYPE(realspace_grid_type)                :: drho_rs, drpot1( 3 )
    TYPE(structure_factor_type), &
      ALLOCATABLE, DIMENSION(:)              :: exp_igr_proj

!------------------------------------------------------------------------------
!-------------- INITIALISATION ---------------------
! allocating the structure factor vector to the number of grid_types

  ngrids = SIZE ( dg_proj )
  ALLOCATE ( exp_igr_proj ( ngrids ), STAT=allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "pol_projector_pw", &
   "exp_igr_proj", ngrids )

  DO ig = 1, SIZE ( dg_proj )

    nproj = dg_proj ( ig ) % ncoefs

    gt = dg_proj ( ig )% grid_index

    CALL structure_factor_allocate ( grid_s ( gt ) % bounds, nproj,  &
          exp_igr_proj ( ig ), allocate_centre = .TRUE. )
  ENDDO

!-------------- PROJECTOR AND FORCE CALCULATION ---------------------

  ngt = SIZE ( grid_b )
 

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
          CALL coeff_allocate ( drhob_r, grid_b ( gt ), use_data = PW_COMPLEXDATA3D )
          drhob_r % pw % in_space = PW_REALSPACE

! Allocate 1D COMPLEX coefficient type 
          CALL coeff_allocate ( coeff_dummy, grid_b ( gt ), use_data = PW_COMPLEXDATA1D )
          coeff_dummy % pw % in_space = PW_RECIPROCALSPACE
! Copy the PW-coefficients to coeff_dummy
          CALL pw_copy ( pw, coeff_dummy % pw )
! DO FFT from coeff_dummy => rhob_r 
          CALL coeff_transform_space ( coeff_dummy, drhob_r )
! Allocate large real-space structure based on grid_b
          CALL rs_grid_setup ( drho_rs, grid_b ( gt ), -1 )
          CALL rs_grid_allocate ( drho_rs )
          drho_rs % r = 0.0_dp
! Transfer from  real-space <- coefficient_type
          CALL rs_pw_transfer ( drho_rs, drhob_r % pw, "BACKWARD" )

          IF ( .NOT. PRESENT ( fproj_part ) ) CALL coeff_deallocate ( coeff_dummy )
          IF ( PRESENT ( fproj_part ) ) THEN 
            CALL rs_grid_setup ( drpot1, grid_b ( gt ) , -1 )
            CALL rs_grid_allocate ( drpot1 ( 1 : 3 ) )
! allocate intermediate coefficient arrays on the reference grid for
! the derivative of the potential of drho 
            DO i = 1, 3
               CALL coeff_allocate ( dphi1_g ( i ), grid_b ( gt ), use_data = PW_COMPLEXDATA1D )
               dphi1_g ( i ) % pw % in_space = PW_RECIPROCALSPACE
            END DO

            DO i = 1, 3
              n = 0
              n ( i ) = 1
              CALL pw_derive ( coeff_dummy % pw, dphi1_g ( i ) % pw, n )
            END DO

            CALL coeff_deallocate ( coeff_dummy )

! Reallocate coeff with REALDATA3D to use as a dummy variable
! to transfer dphi -> drpot and phi_r -> phi_rs

            CALL coeff_allocate ( coeff_dummy , grid_b ( gt ), use_data = PW_REALDATA3D )
            coeff_dummy % pw % in_space = PW_REALSPACE

            DO i = 1, 3
               CALL coeff_transform_space ( dphi1_g ( i ), coeff_dummy )
               CALL coeff_deallocate ( dphi1_g ( i ) )
               CALL rs_pw_transfer ( drpot1 ( i ), coeff_dummy % pw, "BACKWARD" )
            END DO
            CALL coeff_deallocate ( coeff_dummy )
          ENDIF
          first_time = .FALSE. 
        ENDIF
! Parrallelization
        IF (  drho_rs % parallel .AND. drho_rs % group_dim ( 1 ) /= 1 )  THEN
          CALL get_center ( part, box, exp_igr_proj ( ig ) % centre, grid_b ( gt ) % npts )
        END IF
! Get get force on particle due to projector

! Get projector
        CALL proj_evaluate ( dg_proj ( ig ), drho_rs, part, ao, box,  &
               exp_igr_proj ( ig ), grid_s ( gt ), grid_b ( gt ), rhos1, &
               rhos2, fproj_part, drpot1 ) 
        CALL structure_factor_deallocate ( exp_igr_proj ( ig ) )
      ELSE
        CYCLE
      ENDIF
    ENDDO
! do a little cleaning
    IF ( match ) THEN
      IF ( PRESENT ( fproj_part ) ) CALL rs_grid_deallocate ( drpot1 ( 1 : 3 ) )
      CALL coeff_deallocate ( drhob_r )
      CALL rs_grid_deallocate ( drho_rs )
      CALL coeff_deallocate ( rhos2 )
      CALL coeff_deallocate ( rhos1 )
    ENDIF
  ENDDO
 
  DEALLOCATE ( exp_igr_proj, STAT=allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "pol_projectors_pw", &
   "exp_igr_proj", ngrids )

END SUBROUTINE projectors

!******************************************************************************
!!****** pol_projectors_pw/proj_evaluate [1.0] *
!!
!!   NAME
!!     proj_evaluate
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
SUBROUTINE proj_evaluate ( dg, rden, part, ao, box, exp_igr, grid_s,  &
                           grid_b, rhos1, rhos2, fproj_part, drpot1 ) 


    TYPE(dg_type), INTENT(IN), TARGET        :: dg
    TYPE(realspace_grid_type), INTENT(INOUT) :: rden
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: part
    TYPE(ao_type), INTENT(INOUT)             :: ao
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(structure_factor_type), INTENT(OUT) :: exp_igr
    TYPE(pw_grid_type), INTENT(IN)           :: grid_s, grid_b
    TYPE(coeff_type), INTENT(INOUT)          :: rhos1, rhos2
    TYPE(pol_projector_type), DIMENSION(:), &
      INTENT(INOUT), OPTIONAL                :: fproj_part
    TYPE(realspace_grid_type), &
      DIMENSION(:), INTENT(INOUT), OPTIONAL  :: drpot1

    INTEGER                                  :: c1, c2, handle, ip1, ip2, p1, &
                                                p2
    REAL(KIND=dp)                                :: dvols, ivol, q1, q2, r1( 3 ), &
                                                r2( 3 ), sum, vol
    REAL(KIND=dp), DIMENSION(3)                  :: fat
    TYPE(coef_list_type), POINTER            :: clist

!------------------------------------------------------------------------------

  CALL timeset ( 'PROJECTOR_EVALUATE', 'I', ' ', handle )

  dvols = rhos1 % pw % pw_grid % dvol
  vol = box%deth 
  ivol = 1.0_dp/vol

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

! Projector calculation
     CALL dg_sum_patch_force_1d ( rden, rhos1, exp_igr % centre ( :, ip1 ), sum )
     ao % cr ( c1 ) = sum * dvols * ivol * ao % norm ( c1 )
     IF ( p2 /= 0 ) THEN  
       CALL dg_sum_patch_force_1d ( rden, rhos2, exp_igr % centre ( :, ip2 ), sum )
       ao % cr ( c2 ) = sum * dvols * ivol * ao % norm ( c2 )
     ENDIF

     IF ( PRESENT ( fproj_part ) ) THEN
! forces on the ions
     ! sum boxes on real space grids (big box) due to drpot1
       CALL dg_sum_patch_force_3d ( drpot1, rhos1, exp_igr % centre ( :, ip1 ), fat )
       fproj_part ( c1 ) % ipart = p1
       fproj_part ( c1 ) % force ( 1 ) = -fat ( 1 ) * dvols * ivol * ao % norm ( c1 )
       fproj_part ( c1 ) % force ( 2 ) = -fat ( 2 ) * dvols * ivol * ao % norm ( c1 )
       fproj_part ( c1 ) % force ( 3 ) = -fat ( 3 ) * dvols * ivol * ao % norm ( c1 )
       IF ( p2 /= 0 ) THEN
          CALL dg_sum_patch_force_3d ( drpot1, rhos2, exp_igr % centre ( :, ip2 ), fat )
          fproj_part ( c2 ) % ipart = p2
          fproj_part ( c2 ) % force ( 1 ) =  -fat ( 1 ) * dvols * ivol * ao % norm ( c2 )
          fproj_part ( c2 ) % force ( 2 ) =  -fat ( 2 ) * dvols * ivol * ao % norm ( c2 )
          fproj_part ( c2 ) % force ( 3 ) =  -fat ( 3 ) * dvols * ivol * ao % norm ( c2 )
       END IF
     ENDIF
     ip1 = ip2
  END DO

  CALL timestop ( 0.0_dp, handle )

END SUBROUTINE proj_evaluate

!******************************************************************************
!!****** pol_projectors_pw/set_list_drho [1.0] *
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
!!****** pol_projectors_pw/get_patch [1.0] *
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

END MODULE pol_projectors_pw

!******************************************************************************
