!-----------------------------------------------------------------------------
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/pol_projectors_pw [1.0] *
!!
!!   NAME
!!     pol_fo_kernel_pw
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

MODULE pol_fo_kernel_pw

  USE ao_types, ONLY : ao_type
  USE coefficient_types, ONLY : coeff_type, coeff_allocate, &
       coeff_deallocate, coeff_zero, coeff_transform_space, &
       PW_COMPLEXDATA3D, PW_REALDATA3D, &
       PW_REALSPACE, PW_COMPLEXDATA1D, PW_RECIPROCALSPACE
  USE coefficient_lists, ONLY : coef_list_type     
  USE dgs, ONLY : dg_get_strucfac, dg_sum_patch_force_3d, &
                  dg_sum_patch_force_1d, dg_get_patch
  USE dg_types, ONLY : dg_type
  USE kinds, ONLY: dbl, sgl
  USE particle_lists, ONLY : particle_list_type
  USE particle_types, ONLY : particle_type
  USE pw_grid_types, ONLY : pw_grid_type, HALFSPACE, FULLSPACE
  USE pws, ONLY : pw_integral_ab
  USE pw_types, ONLY :  pw_type, pw_add, pw_allocate,  &
                        pw_copy, pw_deallocate, pw_derive
  USE realspace_grid_types, ONLY : realspace_grid_type, rs_grid_allocate, &
       rs_grid_deallocate, rs_grid_setup, rs_pw_transfer
  USE simulation_cell, ONLY : cell_type
  USE structure_factors, ONLY : structure_factor_allocate, &
       structure_factor_deallocate
  USE structure_factor_types, ONLY : structure_factor_type
  USE termination, ONLY : stop_memory, stop_program
  USE timings, ONLY: timeset, timestop

  IMPLICIT NONE 
  
  PRIVATE
  PUBLIC :: force_fo_kernel

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** pol_fo_kernel_pw/fo_kernel [1.0] *
!!
!!   NAME
!!     force_fo_kernel_pw
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

SUBROUTINE force_fo_kernel ( dg_part, pw, rho0, part, box, grid_s,  &
                        grid_b, energy, f_coef, f_part )

  IMPLICIT NONE

! Arguments
  TYPE ( dg_type ), INTENT ( IN ), DIMENSION ( : ) :: dg_part
  TYPE ( pw_type ), INTENT ( IN ) :: pw
  TYPE ( coeff_type ), INTENT ( IN ) :: rho0
  TYPE ( particle_type ), DIMENSION ( : ), INTENT ( IN ) :: part
  TYPE ( cell_type ), INTENT ( IN ) :: box
  TYPE ( pw_grid_type ), INTENT ( IN ), DIMENSION ( : ) :: grid_s
  TYPE ( pw_grid_type ), INTENT ( IN ), DIMENSION ( : ) :: grid_b
  REAL ( dbl ), INTENT(OUT) :: energy
  COMPLEX ( dbl ), DIMENSION ( : ), INTENT( OUT ), OPTIONAL :: f_coef
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( OUT ), OPTIONAL :: f_part

! Locals
  INTEGER :: handle
  INTEGER :: ig, allocstat, ngrids, igtot, ngt
  INTEGER :: npart,  gt, i, n ( 3 )
  REAL ( dbl ) :: ivol, b, a
  TYPE ( realspace_grid_type ) ::  dpotk1_rs ( 3 )
  TYPE ( structure_factor_type ), ALLOCATABLE, DIMENSION ( : ) :: exp_igr_part
  TYPE ( coeff_type ) :: rhos1, rhos2
  TYPE ( coeff_type ) :: k1_g, dk1_g, drhodk1_g, dpotk1_g ( 3 )
  TYPE ( coeff_type ) :: k1_r, dk1_r, rho0_r     
  TYPE ( coeff_type ) :: coeff_dummy, coeff_dummy2 
  LOGICAL :: first_time, match

!------------------------------------------------------------------------------

!-------------- INITIALISATION ---------------------
  ivol = 1._dbl / box % deth
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

  CALL coeff_allocate (k1_g, pw % pw_grid, use_data = PW_COMPLEXDATA1D )
  k1_g % pw % in_space = PW_RECIPROCALSPACE
  CALL coeff_zero ( k1_g )

  CALL coeff_allocate ( dk1_g, pw % pw_grid, use_data = PW_COMPLEXDATA1D )
  dk1_g % pw % in_space = PW_RECIPROCALSPACE
  CALL coeff_zero ( dk1_g )

  CALL coeff_allocate ( coeff_dummy, pw % pw_grid, use_data = PW_COMPLEXDATA1D )
  coeff_dummy % pw % in_space = PW_RECIPROCALSPACE
  CALL coeff_zero ( coeff_dummy )

  CALL coeff_allocate ( rho0_r , pw % pw_grid, use_data = PW_REALDATA3D)
  rho0_r % pw % in_space = PW_REALSPACE
  CALL coeff_zero ( rho0_r )

  CALL coeff_allocate ( k1_r , pw % pw_grid, use_data = PW_REALDATA3D)
  k1_r % pw % in_space = PW_REALSPACE
  CALL coeff_zero ( k1_r )

  CALL coeff_allocate ( dk1_r , pw % pw_grid, use_data = PW_REALDATA3D)
  dk1_r % pw % in_space = PW_REALSPACE
  CALL coeff_zero ( dk1_r )

  CALL coeff_transform_space (rho0 % pw, rho0_r  )

! calculate thomas-fermi + dirac exchange
  CALL calc_tfx_pot ( rho0_r, k1_r, dk1_r )
! transform the potentials to gspace
  CALL coeff_transform_space ( k1_r, k1_g )
  CALL coeff_transform_space ( dk1_r, dk1_g )
  coeff_dummy % pw % cc = CONJG ( pw % cc )

  energy = ivol * pw_integral_ab ( coeff_dummy % pw, k1_g % pw )
  IF ( PRESENT ( f_coef ) )  THEN
    k1_g % pw % cc = CONJG (k1_g % pw % cc )
    f_coef = - 1._dbl *  &
               CMPLX ( 2._dbl * REAL ( k1_g % pw % cc, dbl ), &
                       2._dbl * AIMAG ( k1_g % pw % cc ), dbl )
  ENDIF
  CALL coeff_deallocate ( k1_g )
  CALL coeff_deallocate ( k1_r )
  CALL coeff_deallocate ( dk1_r )
  CALL coeff_deallocate ( rho0_r )
  CALL coeff_deallocate ( coeff_dummy )

!--------------  FORCE CALCULATION ---------------------

  IF ( PRESENT ( f_part ) ) THEN 
     ngt = SIZE ( grid_b )
     DO gt = 1, ngt
       first_time = .TRUE.
       igtot = SIZE ( dg_part )
       IF ( igtot == 0 ) EXIT
       match = .FALSE.
       DO ig = 1, igtot   ! DRHO in AO basis
  
         IF ( gt == dg_part ( ig ) % grid_index ) THEN
           match = .TRUE. 
           IF ( first_time ) THEN
! Allocate small coefficient types
             CALL coeff_allocate ( rhos1,  grid_s ( gt ), use_data = PW_REALDATA3D )
             CALL coeff_allocate ( rhos2,  grid_s ( gt ), use_data = PW_REALDATA3D )
! Allocate 1D COMPLEX coefficient type 
             CALL coeff_allocate ( coeff_dummy, grid_b ( gt ), use_data = PW_COMPLEXDATA1D )
             coeff_dummy % pw % in_space = PW_RECIPROCALSPACE
! Copy the PW-coefficients to coeff_dummy
             CALL pw_copy ( pw, coeff_dummy % pw )
             coeff_dummy % pw % cc = CONJG ( coeff_dummy %  pw % cc )
! Allocate 1D COMPLEX coefficient type 
             CALL coeff_allocate ( coeff_dummy2, grid_b ( gt ), use_data = PW_COMPLEXDATA1D )
             coeff_dummy2 % pw % in_space = PW_RECIPROCALSPACE
             CALL pw_copy ( dk1_g % pw, coeff_dummy2 % pw )
             CALL coeff_allocate ( drhodk1_g, grid_b ( gt ), use_data = PW_COMPLEXDATA1D )
             drhodk1_g % pw % in_space = PW_RECIPROCALSPACE
             DO i = 1, SIZE ( drhodk1_g % pw % cc )   
                a = REAL ( coeff_dummy % pw % cc (i)) * REAL ( coeff_dummy2 % pw % cc (i) ) &
                   - AIMAG ( coeff_dummy % pw % cc (i)) * AIMAG ( coeff_dummy2 % pw % cc (i)) 
                b = REAL ( coeff_dummy % pw % cc (i)) * AIMAG ( coeff_dummy2 % pw % cc (i) ) &
                   + AIMAG ( coeff_dummy % pw % cc (i)) * REAL ( coeff_dummy2 % pw % cc (i)) 
               drhodk1_g % pw % cc (i) = CMPLX ( a, b, dbl )
             END DO
             CALL coeff_deallocate ( coeff_dummy )
             CALL coeff_deallocate ( coeff_dummy2 )

! Allocate large real-space structure based on grid_b
             CALL rs_grid_setup ( dpotk1_rs, grid_b ( gt ), -1 )
             CALL rs_grid_allocate ( dpotk1_rs (1 : 3) )

! allocate intermediate coefficient arrays on the reference grid for
! the derivative of the potential of drho 
             DO i = 1, 3
                CALL coeff_allocate ( dpotk1_g ( i ), grid_b ( gt ), use_data = PW_COMPLEXDATA1D )
                dpotk1_g ( i ) % pw % in_space = PW_RECIPROCALSPACE
             END DO

             DO i = 1, 3
               n = 0
               n ( i ) = 1
               CALL pw_copy ( drhodk1_g % pw, dpotk1_g ( i ) % pw )
               CALL pw_derive ( dpotk1_g ( i ) % pw, n )
             END DO

             CALL coeff_deallocate ( drhodk1_g )

! Reallocate coeff with REALDATA3D to use as a dummy variable
! to transfer dphi -> drpot and phi_r -> phi_rs

             CALL coeff_allocate ( coeff_dummy , grid_b ( gt ), use_data = PW_REALDATA3D )
             coeff_dummy % pw % in_space = PW_REALSPACE

             DO i = 1, 3
                CALL coeff_transform_space ( dpotk1_g ( i ), coeff_dummy )
                CALL coeff_deallocate ( dpotk1_g ( i ) )
                CALL rs_pw_transfer ( dpotk1_rs ( i ), coeff_dummy % pw, "BACKWARD" )
             END DO
             CALL coeff_deallocate ( coeff_dummy )
             first_time = .FALSE. 
         ENDIF
! Parrallelization
         IF (  dpotk1_rs (1) % parallel .AND. dpotk1_rs (1) % group_dim ( 1 ) /= 1 )  THEN
           CALL get_center ( part, box, exp_igr_part ( ig ) % centre, grid_b ( gt ) % npts )
         END IF
! Get get force on particle 
         CALL force_part_evaluate ( dg_part ( ig ), part, f_part, box,  &
                       rhos1, rhos2, exp_igr_part ( ig ), grid_b ( gt ),&
                       grid_s ( gt ), dpotk1_rs )
         CALL structure_factor_deallocate ( exp_igr_part ( ig ) )
       ELSE
         CYCLE
       ENDIF
     ENDDO
! do a little cleaning
     IF ( match ) THEN
       CALL rs_grid_deallocate ( dpotk1_rs ( 1 : 3 ) )
       CALL coeff_deallocate ( rhos2 )
       CALL coeff_deallocate ( rhos1 )
     ENDIF
  ENDDO
 END IF
 
 CALL coeff_deallocate ( dk1_g )
 DEALLOCATE ( exp_igr_part, STAT=allocstat )
 IF ( allocstat /= 0 ) CALL stop_memory ( "pol_fo_kernel", &
  "exp_igr_part", ngrids )

END SUBROUTINE force_fo_kernel

!******************************************************************************
!!****** pol_fo_kernel/force_part_evaluate [1.0] *
!!
!!   NAME
!!     force_part_evaluate
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
SUBROUTINE force_part_evaluate ( dg, part, force, box, rhos1, rhos2,  &
                                 exp_igr, grid_b, grid_s, drpot1 )

! Arguments
  TYPE ( dg_type ), INTENT ( IN ) :: dg
  TYPE ( particle_type ), DIMENSION ( : ), INTENT ( IN ) :: part
  TYPE ( cell_type ), INTENT ( IN ) :: box
  TYPE ( pw_grid_type ), INTENT ( IN ) :: grid_b
  TYPE ( pw_grid_type ), INTENT ( IN ) :: grid_s
  TYPE ( structure_factor_type ), INTENT ( IN ) :: exp_igr
  TYPE ( coeff_type ), INTENT ( INOUT ) :: rhos1, rhos2
  REAL ( dbl ), INTENT (INOUT), DIMENSION (:,:) :: force
  TYPE ( realspace_grid_type ), DIMENSION ( : ), INTENT ( INOUT ) :: drpot1

! Locals
  INTEGER :: handle, p1, p2, ip1, ip2
  REAL ( dbl ) :: ivol, dvols, q1, q2
  REAL ( dbl ) :: r1 ( 3 ), r2 ( 3 )
  REAL ( dbl ), DIMENSION ( 3 ) :: fat
  TYPE ( particle_list_type ), POINTER :: plist

  CALL timeset ( 'FORCE_EVALUATE_PART', 'I', ' ', handle )

  ivol = 1._dbl / box % deth

! initialize the forces
  dvols = rhos1 % pw % pw_grid % dvol 

  plist => dg % plist
  ip1 = 0
  ip2 = 0
  DO

     CALL set_list_rho0 ( plist, p1, p2, ip1, ip2 )
     IF ( p1 == 0 .AND. p2 == 0 ) EXIT
     r1 = part ( p1 ) % r
     q1 = 1._dbl
     r2=0._dbl
     q2=0._dbl
     IF ( p2 /= 0 ) THEN
       r2 = part ( p2 ) % r
       q2 = 1._dbl                  
     ENDIF
     ! calculate function on small boxes (we use double packing in FFT)
! calculate function on small boxes (we use double packing in FFT)
     CALL get_patch ( dg, r1, r2, q1, q2, exp_igr, box, p1, p2, ip1, ip2, &
                      grid_b, grid_s, rhos1, rhos2 )

     ! sum boxes on real space grids (big box) for forces due to drpot1:
     CALL dg_sum_patch_force_3d ( drpot1, rhos1, exp_igr % centre ( :, ip1 ), fat )
     force ( 1, p1 ) = force ( 1, p1 ) - fat ( 1 ) * dvols * ivol
     force ( 2, p1 ) = force ( 2, p1 ) - fat ( 2 ) * dvols * ivol 
     force ( 3, p1 ) = force ( 3, p1 ) - fat ( 3 ) * dvols * ivol
     IF ( p2 /= 0 ) THEN
        CALL dg_sum_patch_force_3d ( drpot1, rhos2, exp_igr % centre ( :, ip2 ), fat )
        force ( 1, p2 ) = force ( 1, p2 ) - fat ( 1 ) * dvols * ivol
        force ( 2, p2 ) = force ( 2, p2 ) - fat ( 2 ) * dvols * ivol
        force ( 3, p2 ) = force ( 3, p2 ) - fat ( 3 ) * dvols * ivol
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


!---------------------------------------------------------------------------
SUBROUTINE calc_tfx_pot ( rho, pot, dpot )
!---------------------------------------------------------------------------

  IMPLICIT NONE

! Arguments
  TYPE ( coeff_type ) , INTENT(IN), TARGET :: rho
  TYPE ( coeff_type ) , INTENT(INOUT)      :: pot
  TYPE ( coeff_type ) , INTENT(INOUT)      :: dpot

! *** Local variables ***
  INTEGER, DIMENSION ( :, : ), POINTER :: bs
  INTEGER :: nd ( 3 ), ub ( 3 ), lb ( 3 )
  INTEGER :: n, ierr, i, j, k
  REAL(dbl), PARAMETER :: pi = 3.14159265358979323846264338_dbl
  REAL(dbl), PARAMETER :: rsfac = 0.6203504908994000166680065_dbl ! (4*pi/3)^(-1/3)
  REAL(dbl), PARAMETER :: eps_rho = 1.E-12_dbl              
  REAL(dbl), PARAMETER :: f13 = 1.0_dbl/3.0_dbl, &
                          f23 = 2.0_dbl*f13, &
                          f43 = 4.0_dbl*f13, &
                          f53 = 5.0_dbl*f13
  REAL(dbl) :: ctf, cxa, xparam, f, f2, fxa, fxa2
  REAL(dbl), DIMENSION(:,:,:), ALLOCATABLE :: r13
!   ---------------------------------------------------------------------------

   ctf = 0.3_dbl*(3._dbl*pi*pi)**f23
   xparam = 2._dbl / 3._dbl
   cxa = -9._dbl/8._dbl * xparam * (3._dbl/pi)**f13

   bs => rho % pw % pw_grid % bounds
   lb = bs (1,:)
   ub = bs (2,:)
   nd = rho % pw % pw_grid % npts
   n = PRODUCT ( nd )
   ALLOCATE ( r13 (lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)), STAT = ierr )
   IF ( ierr /= 0 ) call stop_memory ( "calc_tfx_pot","r13",n)
! initialize to 0
   r13=0._dbl

   DO i = bs ( 1, 1 ), bs ( 2, 1 )
      DO j = bs ( 1, 2 ), bs ( 2, 2 )
         DO k = bs ( 1, 3 ), bs ( 2, 3 )
            IF (  ABS (rho % pw % cr3d (i, j, k)) < eps_rho ) THEN
               r13(i,j,k) = 0.0_dbl
             ELSE
               IF ( rho % pw % cr3d (i, j, k) < - eps_rho ) THEN
                  r13(i,j,k) =  - ABS (rho % pw % cr3d (i, j, k))**f13 
               ELSE IF ( rho % pw % cr3d (i, j, k) >   eps_rho ) THEN
                  r13(i,j,k) =    ABS (rho % pw % cr3d (i, j, k))**f13 
               END IF
             END IF
         END DO
      END DO
   END DO

!V_thomas_fermi + V_xc (first derivatives of sr energy)

  f = f53 * ctf
  fxa = f43 * cxa

   DO i = bs ( 1, 1 ), bs ( 2, 1 )
      DO j = bs ( 1, 2 ), bs ( 2, 2 )
         DO k = bs ( 1, 3 ), bs ( 2, 3 )
            IF (  ABS ( rho % pw % cr3d (i, j, k) ) > eps_rho ) THEN
              pot % pw % cr3d (i,j,k) = f * r13(i,j,k)*r13(i,j,k) +  fxa * r13(i,j,k)
            END IF
         END DO
      END DO
   END DO

! derivative of potential

  f2 = f23 * f53 * ctf
  fxa2 = f13 * f43 * cxa

   DO i = bs ( 1, 1 ), bs ( 2, 1 )
      DO j = bs ( 1, 2 ), bs ( 2, 2 )
         DO k = bs ( 1, 3 ), bs ( 2, 3 )
            IF (  ABS ( rho % pw % cr3d (i, j, k) ) > eps_rho ) THEN
               dpot % pw % cr3d (i,j,k) =  f2 / r13(i,j,k)  + fxa2 * r13(i,j,k) / rho% pw% cr3d (i,j,k)
            END IF
         END DO
      END DO
   END DO
   DEALLOCATE (r13, STAT=ierr )
END SUBROUTINE calc_tfx_pot

!******************************************************************************

 END MODULE pol_fo_kernel_pw     

!******************************************************************************
