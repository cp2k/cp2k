!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000, 2006  CP2K developers group                           !
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

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type
  USE coefficient_types,               ONLY: coeff_type
  USE dg_rho0_types,                   ONLY: dg_rho0_type
  USE dg_types,                        ONLY: dg_get,&
                                             dg_type
  USE dgs,                             ONLY: dg_get_patch,&
                                             dg_get_strucfac,&
                                             dg_sum_patch,&
                                             dg_sum_patch_force_3d
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE ewald_pw_types,                  ONLY: ewald_pw_get,&
                                             ewald_pw_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fourpi
  USE particle_types,                  ONLY: particle_type
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_poisson_methods,              ONLY: pw_poisson_solve
  USE pw_poisson_types,                ONLY: pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_coeff,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_init_coeff,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_integral_a2b,&
                                             pw_p_type,&
                                             pw_transfer,&
                                             pw_type
  USE realspace_grid_types,            ONLY: &
       pw2rs, realspace_grid_p_type, realspace_grid_type, rs2pw, &
       rs_grid_create, rs_grid_p_create, rs_grid_p_release, rs_grid_release, &
       rs_pw_transfer
  USE structure_factor_types,          ONLY: structure_factor_type
  USE structure_factors,               ONLY: structure_factor_allocate,&
                                             structure_factor_deallocate
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: pme_evaluate
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pme'

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

  SUBROUTINE pme_evaluate ( ewald_env, ewald_pw, particle_set, vg_coulomb, fg_coulomb,&
       pv_g, force_env_section, error )
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), INTENT(OUT)               :: vg_coulomb
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT), OPTIONAL                  :: fg_coulomb, pv_g
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pme_evaluate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: group, handle, i, ipart, j, &
                                                jpart, npart, p1, p2
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: alpha, dvols, eps0, ffa, &
                                                inv_fourpi_eps0
    REAL(KIND=dp), DIMENSION(3)              :: fat
    REAL(KIND=dp), DIMENSION(3, 3)           :: f_stress, h_stress
    TYPE(cell_type), POINTER                 :: box
    TYPE(coeff_type)                         :: rhos1, rhos2
    TYPE(dg_rho0_type), POINTER              :: dg_rho0
    TYPE(dg_type), POINTER                   :: dg
    TYPE(pw_grid_type), POINTER              :: grid_b, grid_s
    TYPE(pw_p_type), DIMENSION(3)            :: dphi_g
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_type), POINTER              :: pw_big_pool, pw_small_pool
    TYPE(pw_type), POINTER                   :: phi_r, rhob_r
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: drpot
    TYPE(realspace_grid_type), POINTER       :: rden
    TYPE(structure_factor_type)              :: exp_igr

    NULLIFY( box, poisson_env, rden, drpot)
    failure = .FALSE.
    CALL timeset ( 'PME', 'I', ' ', handle )
    CALL ewald_env_get (ewald_env, alpha=alpha, eps0 = eps0, group = group )
    CALL ewald_pw_get (ewald_pw, pw_big_pool=pw_big_pool, & 
                       pw_small_pool=pw_small_pool, cell=box, &
                       poisson_env=poisson_env, dg = dg )

    grid_b => pw_big_pool % pw_grid
    grid_s => pw_small_pool % pw_grid

    CALL dg_get ( dg, dg_rho0 = dg_rho0 )

    inv_fourpi_eps0 = 1.0_dp / ( fourpi * eps0 )
    npart = SIZE ( particle_set )

    CALL structure_factor_allocate ( grid_s % bounds, npart, exp_igr, &
         allocate_centre = .TRUE. )

    CALL pw_pool_init_coeff ( pw_small_pool, rhos1, &
                              use_data = REALDATA3D )
    CALL pw_pool_init_coeff ( pw_small_pool, rhos2, &
                              use_data = REALDATA3D )

    CALL rs_grid_create( rden, grid_b, -MAXVAL ( grid_s % npts ),&
         error=error)
    rden % r = 0.0_dp

    IF ( rden % parallel .AND. rden % group_dim ( 1 ) /= 1 ) THEN
       CALL get_center ( particle_set, box, exp_igr % centre, grid_b % npts )
    END IF

    !-------------- DENSITY CALCULATION ----------------

    ipart = 0
    jpart = 0
    DO 

       CALL set_list ( particle_set, npart, exp_igr % centre, p1, p2, rden, &
            ipart, jpart )
       IF ( p1 == 0 .AND. p2 == 0 ) EXIT
       
       ! calculate function on small boxes (we use double packing in FFT)
       CALL get_patch ( dg, particle_set, exp_igr, box, p1, p2, grid_b, grid_s, &
            rhos1, rhos2 )
       
       ! add boxes to real space grid (big box)
       CALL dg_sum_patch ( rden, rhos1, exp_igr % centre ( :, p1 ) )
       IF ( p2 /= 0 ) CALL dg_sum_patch ( rden, rhos2, exp_igr % centre ( :, p2 ) )
       
    END DO
    
    CALL pw_pool_create_pw ( pw_big_pool, rhob_r, use_data = REALDATA3D, in_space = REALSPACE )
    CALL rs_pw_transfer ( rden, rhob_r, rs2pw )

    !-------------- ELECTROSTATIC CALCULATION -----------

    ! allocate intermediate arrays
    DO i = 1, 3
       NULLIFY(dphi_g ( i )%pw)
       CALL pw_pool_create_pw ( pw_big_pool, dphi_g ( i )%pw, &
                                use_data = COMPLEXDATA1D, &
                                in_space = RECIPROCALSPACE )
    END DO
    CALL pw_pool_create_pw ( pw_big_pool, phi_r, &
                             use_data = REALDATA3D, &
                             in_space = REALSPACE )

    CALL pw_poisson_solve ( poisson_env, rhob_r, vg_coulomb, phi_r, dphi_g, h_stress )
    
    CALL pw_pool_give_back_pw ( pw_big_pool, phi_r )
    
    vg_coulomb = vg_coulomb * inv_fourpi_eps0

    !---------- END OF ELECTROSTATIC CALCULATION --------

    !------------- STRESS TENSOR CALCULATION ------------

    IF ( PRESENT ( pv_g ) ) THEN
       DO i = 1, 3
          DO j = i, 3
             f_stress ( i, j ) = pw_integral_a2b ( dphi_g ( i ) % pw, dphi_g ( j ) % pw )
             f_stress ( j, i ) = f_stress ( i, j )
          END DO
       END DO
       ffa = ( 1.0_dp / fourpi ) * ( 0.5_dp / dg_rho0 % zet ( 1 ) ) ** 2
       f_stress = -ffa * f_stress
       pv_g = inv_fourpi_eps0 * ( h_stress + f_stress )
    ENDIF
    
    !--------END OF STRESS TENSOR CALCULATION -----------

    CALL rs_grid_p_create ( drpot, 3, grid_b , -MAXVAL ( grid_s % npts ),&
          error=error )

    DO i = 1, 3
       CALL pw_transfer ( dphi_g ( i )%pw, rhob_r )
       CALL pw_pool_give_back_pw ( pw_big_pool, dphi_g ( i )%pw )
       CALL rs_pw_transfer ( drpot ( i )%rs_grid, rhob_r, pw2rs )
    END DO
    
    CALL pw_pool_give_back_pw ( pw_big_pool, rhob_r )
    
    !----------------- FORCE CALCULATION ----------------

    ! initialize the forces
    IF ( PRESENT ( fg_coulomb ) ) THEN
       fg_coulomb = 0.0_dp
       dvols = rhos1 % pw % pw_grid % dvol * inv_fourpi_eps0
       
       ipart = 0
       jpart = 0
       DO
          
          CALL set_list ( particle_set, npart, exp_igr % centre, p1, p2, rden, ipart, jpart )
          IF ( p1 == 0 .AND. p2 == 0 ) EXIT
          
          ! calculate function on small boxes (we use double packing in FFT)
          CALL get_patch_again ( dg, particle_set, exp_igr, p1, p2, rhos1, rhos2 )
          
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

    CALL rs_grid_release ( rden, error)
    CALL rs_grid_p_release ( drpot, error)
    CALL pw_pool_give_back_coeff ( pw_small_pool, rhos1 )
    CALL pw_pool_give_back_coeff ( pw_small_pool, rhos2 )
    CALL structure_factor_deallocate ( exp_igr )

    CALL timestop ( 0.0_dp, handle )
    
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

SUBROUTINE get_patch ( dg, particle_set, exp_igr, box, p1, p2, &
           grid_b, grid_s, rhos1, rhos2 )


    TYPE(dg_type), POINTER                   :: dg
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(structure_factor_type)              :: exp_igr
    TYPE(cell_type), INTENT(IN)              :: box
    INTEGER, INTENT(IN)                      :: p1, p2
    TYPE(pw_grid_type), INTENT(IN)           :: grid_b, grid_s
    TYPE(coeff_type)                         :: rhos1, rhos2

    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: ex1, ex2, ey1, ey2, ez1, ez2
    INTEGER, DIMENSION(:), POINTER           :: center1, center2
    REAL(KIND=dp)                            :: q1, q2
    REAL(KIND=dp), DIMENSION(3)              :: r1, r2
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(coeff_type), POINTER                :: rho0
    TYPE(dg_rho0_type), POINTER              :: dg_rho0

!------------------------------------------------------------------------------

  CALL dg_get ( dg, dg_rho0 = dg_rho0 )
  rho0 => dg_rho0 % density

  r1 = particle_set ( p1 ) % r
  
  atomic_kind => particle_set (p1) % atomic_kind
  CALL get_atomic_kind (atomic_kind=atomic_kind,qeff=q1) 

  center1 => exp_igr % centre ( :, p1 )
  ex1 => exp_igr % ex ( :, p1 )
  ey1 => exp_igr % ey ( :, p1 )
  ez1 => exp_igr % ez ( :, p1 )
  CALL dg_get_strucfac ( box, r1, grid_s % npts, grid_b % npts, center1, &
       exp_igr % lb, ex1, ey1, ez1 )

  IF ( p2 /= 0 ) THEN
     r2 = particle_set ( p2 ) % r
     atomic_kind => particle_set (p2) % atomic_kind
     CALL get_atomic_kind (atomic_kind=atomic_kind,qeff=q2) 
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

SUBROUTINE get_patch_again ( dg, particle_set, exp_igr, p1, p2, rhos1, rhos2 )


    TYPE(dg_type), POINTER                   :: dg
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(structure_factor_type)              :: exp_igr
    INTEGER, INTENT(IN)                      :: p1, p2
    TYPE(coeff_type)                         :: rhos1, rhos2

    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: ex1, ex2, ey1, ey2, ez1, ez2
    REAL(KIND=dp)                            :: q1, q2
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(coeff_type), POINTER                :: rho0
    TYPE(dg_rho0_type), POINTER              :: dg_rho0

!------------------------------------------------------------------------------

  CALL dg_get ( dg, dg_rho0 = dg_rho0 )
  rho0 => dg_rho0 % density

  atomic_kind => particle_set (p1) % atomic_kind
  CALL get_atomic_kind (atomic_kind=atomic_kind,qeff=q1) 
  ex1 => exp_igr % ex ( :, p1 )
  ey1 => exp_igr % ey ( :, p1 )
  ez1 => exp_igr % ez ( :, p1 )

  IF ( p2 /= 0 ) THEN
     atomic_kind => particle_set (p2) % atomic_kind
     CALL get_atomic_kind (atomic_kind=atomic_kind,qeff=q2) 
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

SUBROUTINE set_list ( particle_set, npart, center, p1, p2, rs, ipart, jpart )


    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: npart
    INTEGER, DIMENSION(:, :), INTENT(IN)     :: center
    INTEGER, INTENT(OUT)                     :: p1, p2
    TYPE(realspace_grid_type)                :: rs
    INTEGER, INTENT(INOUT)                   :: ipart, jpart

    INTEGER                                  :: dir, lb, ndim, npos, ub
    REAL(KIND=dp)                            :: charge
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

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
    atomic_kind => particle_set ( ipart ) % atomic_kind
    CALL get_atomic_kind ( atomic_kind=atomic_kind, qeff=charge )
    IF ( charge == 0.0_dp ) CYCLE
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
    atomic_kind => particle_set ( ipart ) % atomic_kind
    CALL get_atomic_kind ( atomic_kind=atomic_kind, qeff=charge )
    IF ( charge == 0.0_dp ) CYCLE
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


    INTEGER, INTENT(IN)                      :: pos, lb, ub
    LOGICAL                                  :: internal

  IF ( pos >= lb .AND. pos <= ub ) THEN
    internal = .TRUE.
  ELSE
    internal = .FALSE.
  END IF

END FUNCTION in_slice

!!*****
!******************************************************************************

SUBROUTINE get_center ( particle_set, box, centre, npts )


    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), INTENT(IN)              :: box
    INTEGER, DIMENSION(:, :), INTENT(OUT)    :: centre
    INTEGER, DIMENSION(:), INTENT(IN)        :: npts

    INTEGER                                  :: ipart
    REAL(KIND=dp), DIMENSION(3)              :: s

!------------------------------------------------------------------------------

  DO ipart = 1, SIZE ( particle_set )
    ! compute the scaled coordinate of atomi
    s = MATMUL ( box % h_inv, particle_set ( ipart ) % r )
    s = s - NINT ( s )
    ! find the closest grid point (on big grid)
    centre ( :, ipart ) = NINT ( REAL ( npts ( 1:3 ),KIND=dp) * s ( 1:3 ) )
  END DO

END SUBROUTINE get_center

!******************************************************************************

END MODULE pme

!******************************************************************************
