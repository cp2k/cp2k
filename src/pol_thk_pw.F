!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/pol_thk_pw [1.0] *
!!
!!   NAME
!!     pol_thk_pw
!!
!!   FUNCTION
!!     Calculation of the first order contribution of the Hohenbeg-Kohn
!!     functional t
!!
!!     E (hk) =   c_i* thk_i
!!
!!   AUTHOR
!!     gloria (06 Jul 2001)     
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE pol_thk_pw

! *****************************************************************************

  USE ao_types, ONLY : ao_type
  USE coefficient_lists, ONLY : coef_list_type
  USE coefficient_types, ONLY : coeff_type, coeff_allocate, &
       coeff_deallocate, coeff_transform_space, coeff_zero, &
       PW_COMPLEXDATA3D, PW_REALDATA3D, &
       PW_REALSPACE, PW_COMPLEXDATA1D, PW_RECIPROCALSPACE
  USE dgs, ONLY : dg_get_strucfac, dg_sum_patch, dg_get_patch
  USE dg_types, ONLY : dg_type
  USE empirical_parameters,     ONLY : empirical_parameter_type
  USE kinds, ONLY: dbl        
  USE particle_types, ONLY : particle_type
  USE pol_projectors_pw, ONLY : pol_projector_type
  USE pw_grid_types, ONLY : pw_grid_type
  USE pws, ONLY : pw_poisson_solver, init_pw_poisson_solver, pw_integral_ab
  USE pw_types, ONLY :  pw_type, pw_add, pw_allocate, pw_copy, pw_deallocate
  USE realspace_grid_types, ONLY : realspace_grid_type, rs_grid_allocate, &
       rs_grid_deallocate, rs_grid_setup, rs_pw_transfer
  USE simulation_cell, ONLY : cell_type
  USE structure_factors, ONLY : structure_factor_allocate, structure_factor_deallocate
  USE structure_factor_types, ONLY : structure_factor_type
  USE termination, ONLY : stop_memory, stop_program
  USE timings, ONLY : timeset, timestop


  IMPLICIT NONE

  PRIVATE
 
  PUBLIC :: force_thk

!*****************************************************************************

CONTAINS

!*****************************************************************************
  SUBROUTINE force_thk ( dg, ao, pw_grid, part, box, empparm, grid_s, grid_b,  &
                        energy, f_coef, f_part, fproj_part )

  TYPE ( dg_type ), INTENT ( IN ), DIMENSION ( : ) :: dg
  TYPE ( ao_type ), INTENT(INOUT) :: ao 
  TYPE ( empirical_parameter_type ), DIMENSION ( : ), INTENT ( IN ) :: empparm
  TYPE ( pw_grid_type ), INTENT ( IN ) :: pw_grid
  REAL ( dbl ), INTENT(OUT) :: energy
  TYPE ( particle_type ), DIMENSION ( : ), INTENT ( IN ) :: part
  TYPE ( cell_type ), INTENT ( IN ) :: box
  TYPE ( pw_grid_type ), INTENT ( IN ), DIMENSION ( : ) :: grid_s
  TYPE ( pw_grid_type ), INTENT ( IN ), DIMENSION ( : ) :: grid_b
  COMPLEX ( dbl ), DIMENSION ( : ), INTENT( OUT ), OPTIONAL :: f_coef
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( OUT ), OPTIONAL :: f_part
  TYPE ( pol_projector_type ), DIMENSION ( : ), INTENT ( IN ), OPTIONAL :: fproj_part

! locals

  INTEGER :: iproj, nproj, ikind, iset, allocstat, ig, gt, ngt, igtot, ngrids
  REAL ( dbl ), ALLOCATABLE, DIMENSION ( : ) :: thk
  TYPE ( realspace_grid_type ) :: rho_rs
  TYPE ( structure_factor_type ), ALLOCATABLE, DIMENSION ( : ) :: exp_igr
  TYPE ( coeff_type ) :: rhos1, rhos2
  TYPE ( coeff_type ) :: rhob_r, rhob_g
  TYPE ( coeff_type ) :: coeff_dummy, rho_tot
  LOGICAL :: match, first_time
  REAL (dbl) :: vol

  ngt = SIZE ( grid_b )
  vol = box % deth

! allocating the structure factor vector to the number of grid_types
  ngrids = SIZE ( dg )
  ALLOCATE ( exp_igr ( ngrids ), STAT=allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "pol_thk_pw", &
   "exp_igr", ngrids )

  DO ig = 1, SIZE ( dg )

    nproj = dg ( ig ) % ncoefs

    gt = dg ( ig ) % grid_index

    CALL structure_factor_allocate ( grid_s ( gt ) % bounds, nproj,  &
          exp_igr ( ig ), allocate_centre = .TRUE. )
  ENDDO

  nproj = SIZE ( ao % cr )
  ALLOCATE ( thk ( nproj ), STAT=allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "pol_thk_pw", &
   "thk", nproj )

  IF ( PRESENT ( f_coef ) ) f_coef ( : ) = 0.0_dbl
  IF ( PRESENT ( f_part ) ) f_part ( :, : ) = 0.0_dbl


  CALL coeff_allocate ( rho_tot, pw_grid, use_data = PW_COMPLEXDATA1D )
  rho_tot % pw % in_space = PW_RECIPROCALSPACE
  CALL coeff_zero ( rho_tot )

  CALL coeff_allocate ( coeff_dummy, pw_grid, use_data = PW_COMPLEXDATA1D )
  coeff_dummy % pw % in_space = PW_RECIPROCALSPACE
  CALL coeff_zero ( coeff_dummy )

  energy = 0._dbl
  DO iproj = 1, nproj
     ikind = ao % coef_to_basis ( iproj )
     iset = ao % coef_to_set ( iproj )
     thk ( iproj ) = empparm ( ikind ) % hk_param ( iset ) * ao % norm ( iproj )
     energy = energy + ao % cr ( iproj ) * thk ( iproj )
  END DO

  IF ( PRESENT ( f_coef ) .OR. PRESENT ( f_part ) ) THEN
    DO gt = 1, ngt
      first_time = .TRUE.
      igtot = SIZE ( dg )
      IF ( igtot == 0 ) EXIT
      match = .FALSE.
      DO ig = 1, igtot  ! DRHO in AO basis
        IF ( gt == dg ( ig ) % grid_index ) THEN
          match = .TRUE.
          IF ( first_time ) THEN

! Allocate small coefficient types
            CALL coeff_allocate ( rhos1,  grid_s ( gt ), use_data = PW_REALDATA3D )
            CALL coeff_allocate ( rhos2,  grid_s ( gt ), use_data = PW_REALDATA3D )
! Allocate large real-space structure based on grid_b
            CALL rs_grid_setup ( rho_rs, grid_b ( gt ), -1 )
            CALL rs_grid_allocate ( rho_rs )
            rho_rs % r = 0._dbl
            first_time = .FALSE.
          ENDIF
          IF (  rho_rs % parallel .AND. rho_rs % group_dim ( 1 ) /= 1 )  THEN
            CALL get_center ( part, box, exp_igr ( ig ) % centre, &
                          grid_b ( gt ) % npts )
          END IF
! denisty and particle_forces
          CALL rho_evaluate ( dg ( ig ), rho_rs, part, box,  &
              exp_igr ( ig ), grid_s ( gt ), grid_b ( gt ), rhos1, rhos2, &
              thk, ao % norm, fproj_part, f_part ) 
          CALL structure_factor_deallocate ( exp_igr ( ig ) )
        ELSE
          CYCLE
        ENDIF
      ENDDO
! Allocate large coefficient structure based on grid_b  
      IF ( match ) THEN
        CALL coeff_allocate ( rhob_r, grid_b ( gt ), use_data = PW_COMPLEXDATA3D )
        CALL coeff_allocate ( rhob_g, grid_b ( gt ), use_data = PW_COMPLEXDATA1D )
        rhob_r % pw % in_space = PW_REALSPACE

! Transfer from real-space -> coefficient_type
        CALL rs_pw_transfer ( rho_rs, rhob_r % pw, "FORWARD" )
 
! DO FFT from rhob_r => rhob_g 
        CALL coeff_transform_space ( rhob_r, rhob_g )

! Copy present drho_ref to coeff_dummy
        CALL pw_copy ( rho_tot % pw, coeff_dummy % pw )

! Sum to the reference grid
        CALL pw_add ( rhob_g % pw, coeff_dummy % pw, rho_tot % pw )

! do a little cleaning
        CALL coeff_deallocate ( rhob_g )
        CALL coeff_deallocate ( rhob_r )
        CALL rs_grid_deallocate ( rho_rs )
        CALL coeff_deallocate ( rhos2 )
        CALL coeff_deallocate ( rhos1 )
      ENDIF
    ENDDO
  ENDIF

  CALL coeff_deallocate ( coeff_dummy )
  DEALLOCATE ( exp_igr, STAT=allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "pol_projectors_pw", &
   "exp_igr", ngrids )
  DEALLOCATE ( thk, STAT=allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "pol_projectors_pw", &
   "thk", nproj )
  
  IF ( PRESENT ( f_coef ) ) THEN
   
    rho_tot % pw % cc = CONJG ( rho_tot % pw % cc )
    f_coef = CMPLX (-2._dbl * REAL ( rho_tot % pw % cc, dbl ), &
                    2._dbl * AIMAG ( rho_tot % pw % cc ), dbl )
! renormalize to 0 the force on g=0
    f_coef (1) =  CMPLX ( 0.0_dbl, 0.0_dbl, dbl )
  ENDIF
  CALL coeff_deallocate ( rho_tot )
 END SUBROUTINE force_thk

!******************************************************************************
!!****** pol_force_thk_pw/rho_evaluate [1.0] *
!!
!!   NAME
!!     rho_evaluate
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

SUBROUTINE rho_evaluate ( dg, rden, part, box, exp_igr, grid_s, &
    grid_b, rhos1, rhos2, thk, norm, fproj_part, f_part )

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
  REAL ( dbl ), DIMENSION ( : ), INTENT ( IN ) :: thk
  REAL ( dbl ), DIMENSION ( : ), INTENT ( IN ) :: norm
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( INOUT ), OPTIONAL :: f_part
  TYPE ( pol_projector_type ), DIMENSION ( : ), INTENT ( IN ), OPTIONAL :: fproj_part

! Locals
  INTEGER :: handle, p1, p2, c1, c2, ip1, ip2
  TYPE ( coef_list_type ), POINTER :: clist
  REAL ( dbl ) :: r1 ( 3 ), r2 ( 3 ), q1, q2

!------------------------------------------------------------------------------


  CALL timeset ( 'FORCE_EVALUATE', 'I', ' ', handle )

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
     q1 = thk ( c1 ) * norm ( c1 )
     r2 = 0._dbl
     q2 = 0._dbl
     IF ( p2 /= 0 ) THEN
       r2 = part ( p2 ) % r
       q2 = thk ( c2 ) * norm ( c2 )
     ENDIF
! calculate function on small boxes (we use double packing in FFT)
     CALL get_patch ( dg, r1, r2, q1, q2, exp_igr, box, p1, p2, ip1, ip2, &
                      grid_b, grid_s, rhos1, rhos2 )
! add boxes to real space grid (big box)
     CALL dg_sum_patch ( rden, rhos1, exp_igr % centre ( :, ip1 ) )
     IF ( p2 /= 0 ) CALL dg_sum_patch ( rden, rhos2, exp_igr % centre ( :, ip2 ) )

     IF ( PRESENT ( f_part ) ) THEN
! forces on the ions
       f_part ( 1, p1 ) = f_part ( 1, p1 ) + thk ( c1 ) * fproj_part ( c1 ) % force ( 1 )
       f_part ( 2, p1 ) = f_part ( 2, p1 ) + thk ( c1 ) * fproj_part ( c1 ) % force ( 2 ) 
       f_part ( 3, p1 ) = f_part ( 3, p1 ) + thk ( c1 ) * fproj_part ( c1 ) % force ( 3 )
       IF ( p2 /= 0 ) THEN
         f_part ( 1, p2 ) = f_part ( 1, p2 ) + thk ( c2 ) * fproj_part ( c2 ) % force ( 1 )
         f_part ( 2, p2 ) = f_part ( 2, p2 ) + thk ( c2 ) * fproj_part ( c2 ) % force ( 2 ) 
         f_part ( 3, p2 ) = f_part ( 3, p2 ) + thk ( c2 ) * fproj_part ( c2 ) % force ( 3 ) 
       END IF
     ENDIF
     ip1 = ip2
  END DO

  CALL timestop ( 0.0_dbl, handle )

END SUBROUTINE rho_evaluate
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

! *****************************************************************************

END MODULE pol_thk_pw

! *****************************************************************************
