!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/spme [1.0] *
!!
!!   NAME
!!     spme
!!
!!   FUNCTION
!!     Calculate the electrostatic energy by the Smooth Particle Ewald method
!!
!!   AUTHOR
!!     JGH (21-Mar-2001)
!!
!!   MODIFICATION HISTORY
!!     JGH (03-May-2001) : first correctly working version
!!
!!   SOURCE
!******************************************************************************

MODULE spme

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE coefficient_types,               ONLY: coeff_transform_space,&
                                             coeff_type
  USE dgs,                             ONLY: dg_sum_patch,&
                                             dg_sum_patch_force_3d
  USE ewald_environment_types,         ONLY: ewald_environment_type, &
                                             ewald_env_get
  USE ewald_pw_types,                  ONLY: ewald_pw_type, &
                                             ewald_pw_get
  USE greens_fn,                       ONLY: greens_function_type
  USE hartree,                         ONLY: calculate_hartree
  USE kinds,                           ONLY: dbl
  USE mathconstants,                   ONLY: fourpi
  USE particle_types,                  ONLY: particle_type
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_pool_types,                   ONLY: pw_pool_type, &
                                             pw_pool_init_coeff, &
                                             pw_pool_give_back_coeff
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             COMPLEXDATA3D,&
                                             REALDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE, &
                                             pw_type
  USE pws,                             ONLY: cleanup_pw_poisson_solver,&
                                             init_pw_poisson_solver,&
                                             pw_integral_a2b
  USE realspace_grid_types,            ONLY: realspace_grid_type,&
                                             rs_grid_allocate,&
                                             rs_grid_deallocate,&
                                             rs_grid_setup,&
                                             rs_pw_transfer
  USE simulation_cell,                 ONLY: cell_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: spme_evaluate

!!***
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** spme/spme_evaluate [1.0] *
!!
!!   NAME
!!     spme_evaluate
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH (21-Mar-2001)
!!
!!   MODIFICATION HISTORY
!!     JGH (03-May-2001) : SPME with charge definition
!!
!!*** *************************************************************************

SUBROUTINE spme_evaluate ( ewald_env, ewald_pw, particle_set, &
                           fg_coulomb, vg_coulomb, pv_g )


    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: fg_coulomb
    REAL(dbl), INTENT(OUT)                   :: vg_coulomb
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: pv_g
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE ( ewald_pw_type ), POINTER          :: ewald_pw

    INTEGER                                  :: handle, i, ierr, ipart, j, &
                                                jpart, n, ndim( 3 ), npart, p1
    INTEGER                                  :: o_spline, group
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: centre
    REAL(dbl)                                :: dvols, ffa, inv_fourpi_eps0
    REAL(dbl), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: rhos
    REAL(dbl), DIMENSION(3)                  :: fat
    REAL(dbl), DIMENSION(3, 3)               :: f_stress, h_stress
    REAL ( dbl )                             :: alpha, eps0
    TYPE(coeff_type)                         :: dphi_g( 3 ), phi_g, rhob_g, &
                                                rhob_r
    TYPE(realspace_grid_type)                :: drpot( 3 ), rden
    TYPE(pw_grid_type), POINTER              :: grid_spme
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(cell_type), POINTER           :: box

!------------------------------------------------------------------------------
!-------------- INITIALISATION ---------------------

  CALL timeset ( 'SPME', 'I', ' ', handle )
  CALL ewald_env_get ( ewald_env, alpha=alpha, eps0 = eps0, &
                       o_spline = o_spline, group = group )
  CALL ewald_pw_get ( ewald_pw, pw_big_pool=pw_pool, cell = box )
  grid_spme => pw_pool % pw_grid

  inv_fourpi_eps0 = 1.0_dbl / ( fourpi * eps0 )
  npart = SIZE ( particle_set )
 
  n = o_spline
  ALLOCATE ( rhos ( n, n, n ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "spme_evaluate", "rhos", n**3 )
  ndim = n
  CALL rs_grid_setup ( rden, grid_spme, -MAXVAL ( ndim ) )
  CALL rs_grid_allocate ( rden )
  rden % r = 0._dbl

  ALLOCATE ( centre ( 3, npart ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "spme_evaluate", "centre", 3*npart )
  CALL get_center ( particle_set, box, centre, grid_spme % npts, n )

!-------------- DENSITY CALCULATION ----------------

  ipart = 0
  jpart = 0
  DO

     CALL set_list ( particle_set, npart, centre, p1, rden, ipart, jpart )
     IF ( p1 == 0 ) EXIT

     ! calculate function on small boxes
     CALL get_patch ( particle_set, box, grid_spme % npts, p1, rhos )

     ! add boxes to real space grid (big box)
     CALL dg_sum_patch ( rden, rhos, centre(:,p1) )

  END DO

  CALL pw_pool_init_coeff ( pw_pool, rhob_r,  &
                           use_data = REALDATA3D, &
                           in_space = REALSPACE )

  CALL rs_pw_transfer ( rden, rhob_r % pw, "FORWARD" )

  CALL rs_grid_deallocate ( rden )

  ! transform density to G space and add charge function
  CALL pw_pool_init_coeff ( pw_pool, rhob_g,  &
                           use_data = COMPLEXDATA1D, &
                           in_space = RECIPROCALSPACE )
  CALL coeff_transform_space ( rhob_r, rhob_g )
  ! update charge function
  CALL p3m_update ( box % green )
  rhob_g % pw % cc = rhob_g % pw % cc * box % green % p3m_charge % cr

!-------------- ELECTROSTATIC CALCULATION -----------

  ! allocate intermediate arrays
  DO i = 1, 3
     CALL pw_pool_init_coeff ( pw_pool, dphi_g ( i ),  &
                              use_data = COMPLEXDATA1D,   &
                              in_space = RECIPROCALSPACE )
  END DO
  CALL pw_pool_init_coeff ( pw_pool, phi_g,  &
                           use_data = COMPLEXDATA1D,   &
                           in_space = RECIPROCALSPACE )

  CALL init_pw_poisson_solver ( box % green )

  CALL calculate_hartree ( rhob_g, vg_coulomb, phi_g, dphi_g, h_stress )

  CALL pw_pool_give_back_coeff ( pw_pool, phi_g )
  CALL pw_pool_give_back_coeff ( pw_pool, rhob_g )
  vg_coulomb = vg_coulomb * inv_fourpi_eps0

!---------- END OF ELECTROSTATIC CALCULATION --------

!------------- STRESS TENSOR CALCULATION ------------

  DO i = 1, 3
    DO j = i, 3
       f_stress ( i, j ) = pw_integral_a2b ( dphi_g ( i ) % pw, &
                                             dphi_g ( j ) % pw )
       f_stress ( j, i ) = f_stress ( i, j )
    END DO
  END DO
  ffa = ( 1.0_dbl / fourpi ) * &
        ( 0.5_dbl / alpha ) ** 2
  f_stress = -ffa * f_stress
  pv_g = inv_fourpi_eps0 * ( h_stress + f_stress )

!--------END OF STRESS TENSOR CALCULATION -----------

  ! move derivative of potential to real space grid and
  ! multiply by charge function in g-space
  CALL rs_grid_setup ( drpot, grid_spme, -MAXVAL ( ndim ) )
  CALL rs_grid_allocate ( drpot ( 1 : 3 ) )
  DO i = 1, 3
     dphi_g ( i ) % pw % cc = dphi_g ( i ) % pw % cc * &
                              box % green % p3m_charge % cr
     CALL coeff_transform_space ( dphi_g ( i ), rhob_r )
     CALL pw_pool_give_back_coeff ( pw_pool, dphi_g ( i ) )
     CALL rs_pw_transfer ( drpot ( i ), rhob_r % pw, "BACKWARD" )
  END DO

  CALL pw_pool_give_back_coeff ( pw_pool, rhob_r )

!----------------- FORCE CALCULATION ----------------

! initialize the forces
  fg_coulomb = 0._dbl
  dvols = grid_spme % dvol * inv_fourpi_eps0

  ipart = 0
  jpart = 0
  DO

     CALL set_list ( particle_set, npart, centre, p1, rden, ipart, jpart )
     IF ( p1 == 0 ) EXIT

     ! calculate function on small boxes
     CALL get_patch ( particle_set, box, grid_spme % npts, p1, rhos )

     ! add boxes to real space grid (big box)
     CALL dg_sum_patch_force_3d ( drpot, rhos, centre(:,p1), fat )
     fg_coulomb ( 1, p1 ) = fg_coulomb ( 1, p1 ) - fat ( 1 ) * dvols
     fg_coulomb ( 2, p1 ) = fg_coulomb ( 2, p1 ) - fat ( 2 ) * dvols
     fg_coulomb ( 3, p1 ) = fg_coulomb ( 3, p1 ) - fat ( 3 ) * dvols

  END DO

!--------------END OF FORCE CALCULATION -------------

!------------------CLEANING UP ----------------------

  CALL rs_grid_deallocate ( drpot ( 1 : 3 ) )
  DEALLOCATE ( rhos, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "spme_evaluate", "rhos" )
  DEALLOCATE ( centre, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "spme_evaluate", "centre" )
  CALL cleanup_pw_poisson_solver ( box % green )

  CALL timestop ( 0.0_dbl, handle )

END SUBROUTINE spme_evaluate

!******************************************************************************

SUBROUTINE p3m_update ( gf )


    TYPE(greens_function_type), &
      INTENT(INOUT), TARGET                  :: gf

    INTEGER                                  :: ig, l, m, n
    REAL(dbl)                                :: arg, novol
    REAL(dbl), DIMENSION(:, :), POINTER      :: bm2
    TYPE(pw_grid_type), POINTER              :: grid
    TYPE(pw_type), POINTER                   :: pc

!------------------------------------------------------------------------------

  grid => gf % influence_function % pw_grid

  ! check if charge function is consistent with current box volume
  IF ( gf % volume /= grid % vol ) THEN

    pc => gf % p3m_charge
    bm2 => gf % p3m_bm2
    arg = 1._dbl / ( 8._dbl * gf % p3m_alpha ** 2 )
    novol = REAL ( grid % ngpts, dbl ) / grid % vol
    DO ig = 1, grid % ngpts_cut_local
       l = grid % g_hat ( 1, ig )
       m = grid % g_hat ( 2, ig )
       n = grid % g_hat ( 3, ig )
       pc % cr ( ig ) = novol * exp ( -arg * grid % gsq ( ig ) ) * &
                        bm2 ( 1, l ) * bm2 ( 2, m ) * bm2 ( 3, n )
    END DO

  END IF

END SUBROUTINE p3m_update

!******************************************************************************
!!****** spme/get_patch [1.0] *
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
!!     none
!!
!!*** *************************************************************************

SUBROUTINE get_patch ( part, box, npts, p, rhos )


    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: part
    TYPE(cell_type), INTENT(IN)              :: box
    INTEGER, DIMENSION(3), INTENT(IN)        :: npts
    INTEGER, INTENT(IN)                      :: p
    REAL(dbl), DIMENSION(:, :, :), &
      INTENT(OUT)                            :: rhos

    INTEGER                                  :: nbox
    REAL(dbl)                                :: q
    REAL(dbl), DIMENSION(3)                  :: delta, r

!------------------------------------------------------------------------------

  nbox = SIZE ( rhos, 1 )
  r = part(p)% r
   
  CALL get_atomic_kind ( atomic_kind=part(p)%atomic_kind, qeff=q)

  CALL get_delta ( box, r, npts, delta, nbox )

  CALL spme_get_patch ( rhos, nbox, delta, q, box % green % p3m_coeff )

END SUBROUTINE get_patch 

!******************************************************************************
!!****** spme/spme_get_patch [1.0] *
!!
!!   NAME
!!     spme_get_patch
!!
!!   FUNCTION
!!     Calculates SPME charge assignment
!!
!!   AUTHOR
!!     JGH (22-Mar-2001)
!!
!!   MODIFICATION HISTORY
!!     DG (29-Mar-2001) : code implemented
!!
!!*** *************************************************************************

SUBROUTINE spme_get_patch ( rhos, n, delta, q, coeff )

    REAL(dbl), DIMENSION(:, :, :), &
      INTENT(OUT)                            :: rhos
    INTEGER, INTENT(IN)                      :: n
    REAL(dbl), DIMENSION(3), INTENT(IN)      :: delta
    REAL(dbl), INTENT(IN)                    :: q
    REAL(dbl), &
      DIMENSION(-(n-1):n-1, 0:n-1), &
      INTENT(IN)                             :: coeff

    INTEGER, PARAMETER                       :: nmax = 12

    INTEGER                                  :: i, i1, i2, i3, j, l
    REAL(dbl)                                :: r2, r3
    REAL(dbl), DIMENSION(3, -nmax:nmax)      :: w_assign
    REAL(dbl), DIMENSION(3, 0:nmax-1)        :: deltal
    REAL(dbl), DIMENSION(3, 1:nmax)          :: f_assign

!------------------------------------------------------------------------------

  IF ( n > nmax ) THEN
    CALL stop_program ( "spme_get_patch","nmax value too small" )
  END IF
! calculate the assignment function values and
! the charges on the grid (small box)

  deltal ( 1, 0 ) = 1._dbl
  deltal ( 2, 0 ) = 1._dbl
  deltal ( 3, 0 ) = 1._dbl
  DO l = 1, n-1
     deltal ( 1, l ) = deltal ( 1, l-1 ) * delta ( 1 )
     deltal ( 2, l ) = deltal ( 2, l-1 ) * delta ( 2 )
     deltal ( 3, l ) = deltal ( 3, l-1 ) * delta ( 3 )
  END DO

  w_assign = 0._dbl
  DO j = -(n-1), n-1, 2
    DO l = 0, n-1
      w_assign ( 1, j ) =  w_assign ( 1, j ) + coeff ( j, l ) * deltal ( 1, l )
      w_assign ( 2, j ) =  w_assign ( 2, j ) + coeff ( j, l ) * deltal ( 2, l )
      w_assign ( 3, j ) =  w_assign ( 3, j ) + coeff ( j, l ) * deltal ( 3, l )
    END DO
  END DO
  DO i = 1, n
    j = n + 1 - 2 * i
    f_assign (1, i ) = w_assign ( 1, j )
    f_assign (2, i ) = w_assign ( 2, j )
    f_assign (3, i ) = w_assign ( 3, j )
  END DO

  DO i3 = 1, n
    r3 = q * f_assign ( 3, i3 )
    DO i2 = 1, n
      r2 = r3 * f_assign ( 2, i2 )
      DO i1 = 1, n
        rhos ( i1, i2, i3 ) = r2 * f_assign ( 1, i1 )
      END DO
    END DO
  END DO

END SUBROUTINE spme_get_patch

!******************************************************************************
!!****** spme/set_list [1.0] *
!!
!!   NAME
!!     set_list
!!
!!   FUNCTION
!!     Advances in the link list of particles to the next to be calculated
!!
!!   AUTHOR
!!     JGH (21-Mar-2001)
!!
!!   MODIFICATION HISTORY
!!     JGH (30-Mar-2001) : spme version adapted from pme
!!
!!*** *************************************************************************

SUBROUTINE set_list ( part, npart, center, p1, rs, ipart, jpart )


    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: part
    INTEGER, INTENT(IN)                      :: npart
    INTEGER, DIMENSION(:, :), INTENT(IN)     :: center
    INTEGER, INTENT(OUT)                     :: p1
    TYPE(realspace_grid_type)                :: rs
    INTEGER, INTENT(INOUT)                   :: ipart, jpart

    INTEGER                                  :: dir, lb, ndim, npos, ub
    REAL(dbl)                                :: charge
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

!------------------------------------------------------------------------------

  p1 = 0
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
    atomic_kind => part(ipart)%atomic_kind
    CALL get_atomic_kind (atomic_kind=atomic_kind,qeff=charge)
    IF ( charge == 0._dbl ) CYCLE
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

!******************************************************************************

SUBROUTINE get_center ( part, box, centre, npts, n )


    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: part
    TYPE(cell_type), INTENT(IN)              :: box
    INTEGER, DIMENSION(:, :), INTENT(OUT)    :: centre
    INTEGER, DIMENSION(:), INTENT(IN)        :: npts
    INTEGER, INTENT(IN)                      :: n

    INTEGER                                  :: ipart, mp
    REAL(dbl)                                :: rmp
    REAL(dbl), DIMENSION(3)                  :: gp, s

!------------------------------------------------------------------------------

  mp = MAXVAL ( npts )
  rmp = REAL ( mp, dbl )
  DO ipart = 1, SIZE ( part )
    ! compute the scaled coordinate of atomi
    s = MATMUL ( box % h_inv, part ( ipart ) % r )
    s = s - NINT ( s )
    gp = REAL ( npts, dbl ) * s
    ! find the closest grid point (on big grid)
    IF ( MOD ( n, 2 ) == 0 ) THEN
      centre ( :, ipart ) = INT ( gp + rmp ) - mp
    ELSE
      centre ( :, ipart ) = NINT ( gp )
    END IF
  END DO

END SUBROUTINE get_center

!******************************************************************************

SUBROUTINE get_delta ( box, r, npts, delta, n )

    TYPE(cell_type), INTENT(IN)              :: box
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: r
    INTEGER, DIMENSION(:), INTENT(IN)        :: npts
    REAL(dbl), DIMENSION(:), INTENT(OUT)     :: delta
    INTEGER, INTENT(IN)                      :: n

    INTEGER                                  :: mp
    INTEGER, DIMENSION(3)                    :: centre
    REAL(dbl)                                :: rmp
    REAL(dbl), DIMENSION(3)                  :: ca, grid_i, s

!------------------------------------------------------------------------------

  mp = MAXVAL ( npts )
  rmp = REAL ( mp, dbl )
! compute the scaled coordinate of atomi
  s = MATMUL ( box % h_inv, r )
  s = s - REAL ( NINT ( s ), dbl )

! find the continuous ``grid'' point
  grid_i ( 1:3 ) = REAL ( npts ( 1:3 ), dbl ) * s ( 1:3 )

! find the closest grid point

  IF ( MOD ( n, 2 ) == 0 ) THEN
    centre ( : ) = INT ( grid_i ( : ) + rmp ) - mp
    ca ( : ) = REAL ( centre ( : ) ) + 0.5_dbl
  ELSE
    centre ( : ) = NINT ( grid_i ( : ) )
    ca ( : ) = REAL ( centre ( : ) )
  END IF

! find the distance vector
  delta ( : ) = grid_i ( : ) - ca ( : )

END SUBROUTINE get_delta


END MODULE spme

!******************************************************************************
