!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/greens_fn [1.0] *
!!
!!   NAME
!!     greens_fn
!!
!!   FUNCTION
!!     Routines to calculate the Green's function for the Poisson equation '''
!!
!!   AUTHOR
!!     JGH (19-Dec-2000) 
!! 
!!   MODIFICATION HISTORY
!!     JGH (9-Mar-2001) : include influence_function into greens_function_type
!!                        add cell volume as indicator for updates
!!     JGH (30-Mar-2001) : Added B-spline routines
!!
!!   SOURCE
!******************************************************************************

MODULE greens_fn
  USE bessel_lib,                      ONLY: bessj0,&
                                             bessj1,&
                                             bessk0,&
                                             bessk1
  USE kinds,                           ONLY: dbl
  USE mathconstants,                   ONLY: fourpi,&
                                             pi,&
                                             twopi
  USE pw_grid_types,                   ONLY: HALFSPACE,&
                                             pw_grid_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             COMPLEXDATA3D,&
                                             NOSPACE,&
                                             REALDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_allocate,&
                                             pw_deallocate,&
                                             pw_transfer,&
                                             pw_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: pw_green_fn, pw_green_fn_init, greens_function_type
  PUBLIC :: PERIODIC3D, ANALYTIC2D, ANALYTIC1D, ANALYTIC0D, &
            MT2D, MT1D, MT0D, HOCKNEY2D, HOCKNEY1D, HOCKNEY0D
  
  TYPE greens_function_type
    INTEGER :: method
    INTEGER :: special_dimension
    INTEGER :: id_nr
    REAL ( dbl ) :: radius
    REAL ( dbl ) :: slab_size
    REAL ( dbl ) :: alpha
    REAL ( dbl ) :: volume
    LOGICAL :: p3m
    INTEGER :: p3m_order
    REAL ( dbl ) :: p3m_alpha
    REAL ( dbl ), DIMENSION ( :, : ), POINTER :: p3m_coeff
    REAL ( dbl ), DIMENSION ( :, : ), POINTER :: p3m_bm2
    TYPE ( pw_type ), POINTER :: influence_function
    TYPE ( pw_type ), POINTER :: p3m_charge
  END TYPE greens_function_type

  INTEGER, PARAMETER :: PERIODIC3D = 1000, ANALYTIC2D = 1001,  &
                        ANALYTIC1D = 1002, ANALYTIC0D = 1003,  &
                        MT2D = 1101, MT1D = 1102, MT0D = 1103, &
                        HOCKNEY2D = 1201, HOCKNEY1D = 1202, HOCKNEY0D = 1203
  INTEGER, SAVE, PRIVATE :: last_greens_fn_id_nr=0

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** greens_fn/pw_green_fn_init [1.0] *
!!
!!   NAME
!!     pw_green_fn_init
!!
!!   SYNOPSIS
!!     Subroutine pw_green_fn_init(gftype, grid, p3m, alpha)
!!       Implicit None
!!       Type(greens_function_type), Intent (INOUT):: gftype
!!       Type(pw_grid_type), Intent (IN):: grid
!!       Integer, Optional, Intent (IN):: p3m
!!       Real(Kind=dbl), Optional, Intent (IN):: alpha
!!     End Subroutine pw_green_fn_init
!!
!!   FUNCTION
!!     Initialisation
!!
!!   AUTHOR
!!     JGH (12-Mar-2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_green_fn_init ( gftype, grid, p3m, alpha )

  IMPLICIT NONE

! Arguments
  TYPE ( greens_function_type ), INTENT ( INOUT ) :: gftype
  TYPE ( pw_grid_type ), INTENT ( IN )  :: grid
  INTEGER, OPTIONAL, INTENT ( IN ) :: p3m
  REAL ( dbl ), OPTIONAL, INTENT ( IN ) :: alpha

! Locals
  INTEGER :: use_data, ierr, n

!------------------------------------------------------------------------------

  SELECT CASE ( gftype % method )
  CASE ( PERIODIC3D, ANALYTIC2D, ANALYTIC1D, ANALYTIC0D, MT2D, MT1D, MT0D )
     ALLOCATE ( gftype % influence_function, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "pw_green_fn_init", &
       "gftype % influence_function", 0 )
     use_data = REALDATA1D
     CALL pw_allocate ( gftype % influence_function, grid, use_data )
     gftype % influence_function % in_space = RECIPROCALSPACE
     IF ( PRESENT ( p3m ) ) THEN
        gftype % p3m = .TRUE.
        gftype % p3m_order = p3m
        IF ( PRESENT ( alpha ) ) THEN
          gftype % p3m_alpha = alpha
        ELSE
          CALL stop_program ( "pw_green_fn_init", "P3M needs alpha PARAMETER")
        END IF
        ALLOCATE ( gftype % p3m_charge, STAT = ierr )
        IF ( ierr /= 0 ) CALL stop_memory ( "pw_green_fn_init", &
           "gftype % p3m_charge", 0 )
        CALL pw_allocate ( gftype % p3m_charge, grid, REALDATA1D )
        gftype % p3m_charge % in_space = RECIPROCALSPACE
        CALL influence_factor ( gftype )
     ELSE
        gftype % p3m = .FALSE.
     END IF
  CASE ( HOCKNEY2D, HOCKNEY1D, HOCKNEY0D )
     CALL stop_program ( "pw_green_fn_init", " Method not implemented ")
  CASE DEFAULT
     CALL stop_program ( "pw_green_fn_init", " Method not implemented ")
  END SELECT
  last_greens_fn_id_nr = last_greens_fn_id_nr+1
  gftype % id_nr = last_greens_fn_id_nr

END SUBROUTINE pw_green_fn_init

!!*****
!******************************************************************************
!!****** greens_fn/pw_green_fn [1.0] *
!!
!!   NAME
!!     pw_green_fn
!!
!!   SYNOPSIS
!!     Subroutine pw_green_fn(gftype)
!!       Implicit None
!!       Type(greens_function_type), Intent (INOUT):: gftype
!!     End Subroutine pw_green_fn
!!
!!   FUNCTION
!!     Calculates the Green's function in reciprocal space  '''
!!
!!   AUTHOR
!!     JGH (19-Dec-2000)
!!
!!   MODIFICATION HISTORY
!!     JGH (9-Mar-2001) : update check through volume
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_green_fn ( gftype )

  IMPLICIT NONE

! Arguments
  TYPE ( greens_function_type ), INTENT ( INOUT ) :: gftype

! Locals
  TYPE ( pw_type ), POINTER :: gf
  TYPE ( pw_grid_type ), POINTER :: grid
  INTEGER :: method, ig, iz, nz, l, m, n
  REAL ( dbl ) :: rlength, zlength, g3d, gxy, g2, j0g, j1g, k0g, k1g, gg

!------------------------------------------------------------------------------

  method = gftype % method
  gf => gftype % influence_function
  grid => gftype % influence_function % pw_grid

  ! some tests
  SELECT CASE ( method )
  CASE ( PERIODIC3D, ANALYTIC2D, ANALYTIC1D, ANALYTIC0D, MT2D, MT1D, MT0D )
    IF ( gf % in_space /= RECIPROCALSPACE .OR. &
         gf % in_use /= REALDATA1D ) &
      CALL stop_program ( "greens_fn", "gf wrong space or DATA TYPE" )
  CASE ( HOCKNEY2D, HOCKNEY1D, HOCKNEY0D )
     CALL stop_program ( "pw_green_fn", " Method not implemented ")
  CASE DEFAULT
     CALL stop_program ( "pw_green_fn", " Method not implemented ")
  END SELECT

  ! check if influence function is consistent with current box volume
  IF ( gftype % volume /= grid % vol ) THEN

    gftype % volume = grid % vol

    SELECT CASE ( method )

    CASE ( PERIODIC3D )

      DO ig = grid % first_gne0, grid % ngpts_cut_local
        g2 = grid % gsq ( ig )
        gf % cr ( ig ) = fourpi / g2
      END DO
      IF ( grid % have_g0 ) gf % cr ( 1 ) = 0._dbl

    CASE ( ANALYTIC2D )

      iz = gftype % special_dimension ! iz is the direction with NO PBC
      zlength = gftype % slab_size    ! zlength is the thickness of the cell
      DO ig = grid % first_gne0, grid % ngpts_cut_local
        nz = grid % g_hat ( iz, ig )
        g2 = grid % gsq ( ig )
        g3d = fourpi / g2
        gg = 0.5_dbl * SQRT ( g2 )
        gf % cr ( ig ) = g3d * ( 1._dbl - (-1._dbl)**nz * EXP ( - gg * zlength ) )
      END DO
      IF ( grid % have_g0 ) gf % cr ( 1 ) = 0._dbl

    CASE ( ANALYTIC1D )

      ! iz is the direction of the PBC ( can be 1,2,3 -> x,y,z )
      iz = gftype % special_dimension
      ! rlength is the radius of the tube
      rlength = gftype % radius
      DO ig = grid % first_gne0, grid % ngpts_cut_local
        g2 = grid % gsq ( ig )
        g3d = fourpi / g2
        gxy = SQRT ( g2 - grid % g(iz,ig) * grid % g(iz,ig) )
        j0g = bessj0 ( rlength * gxy )
        j1g = bessj1 ( rlength * gxy )
        k0g = bessk0 ( rlength * grid % g(iz,ig) )
        k1g = bessk1 ( rlength * grid % g(iz,ig) )
        gf % cr ( ig ) = g3d * ( 1._dbl - rlength * &
          ( gxy * j1g * k0g - grid % g(iz,ig) * j0g * k1g ) )
      END DO
      IF ( grid % have_g0 ) gf % cr ( 1 ) = 0._dbl

    CASE ( ANALYTIC0D )

      rlength = gftype % radius   ! rlength is the radius of the sphere
      DO ig = grid % first_gne0, grid % ngpts_cut_local
        g2 = grid % gsq ( ig )
        gg = 0.5_dbl * SQRT ( g2 )
        g3d = fourpi / g2
        gf % cr ( ig ) = g3d * ( 1._dbl - COS ( rlength * gg ) )
      END DO
      IF ( grid % have_g0 ) &
        gf % cr ( 1 ) = 0.5_dbl * fourpi * rlength * rlength

    CASE ( MT2D )
       CALL stop_program ( "pw_green_fn", " Method not implemented ")
    CASE ( MT1D )
       CALL stop_program ( "pw_green_fn", " Method not implemented ")
    CASE ( MT0D )
       CALL stop_program ( "pw_green_fn", " Method not implemented ")
    CASE ( HOCKNEY2D )
       CALL stop_program ( "pw_green_fn", " Method not implemented ")
    CASE ( HOCKNEY1D )
       CALL stop_program ( "pw_green_fn", " Method not implemented ")
    CASE ( HOCKNEY0D )
       CALL stop_program ( "pw_green_fn", " Method not implemented ")
    CASE DEFAULT
       CALL stop_program ( "pw_green_fn", " Method not implemented ")
    END SELECT

  END IF

END SUBROUTINE pw_green_fn

!!*****
!******************************************************************************
!!****** greens_fn/influence_factor [1.0] *
!!
!!   NAME
!!     influence_factor
!!
!!   FUNCTION
!!     Calculates the influence_factor for the
!!     SPME Green's function in reciprocal space'''
!!
!!   AUTHOR
!!     DH (29-Mar-2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE influence_factor ( gftype )

  IMPLICIT NONE

! Arguments
  TYPE ( greens_function_type ), TARGET, INTENT ( INOUT ) :: gftype

! Locals
  INTEGER :: dim, ierr, j, k, l, n, pt
  INTEGER, DIMENSION ( : ), POINTER :: lb, ub
  INTEGER, DIMENSION ( 3 ) :: npts
  REAL ( dbl ) :: l_arg, prod_arg, val
  REAL ( dbl ), DIMENSION ( : ), ALLOCATABLE :: m_assign
  COMPLEX ( dbl ) :: exp_m, sum_m, b_m

!------------------------------------------------------------------------------

  n = gftype % p3m_order

! calculate the assignment function values

  lb => gftype % influence_function % pw_grid % bounds (1, : )
  ub => gftype % influence_function % pw_grid % bounds (2, : )
  ALLOCATE ( gftype % p3m_bm2 ( 3, MINVAL(lb):MAXVAL(ub) ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( 'influence_factor', &
     'gftype % p3m_bm2', 3*(MAXVAL(ub)-MINVAL(lb)+1) )

  ALLOCATE ( m_assign ( 0:n-2 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( 'influence_factor', 'm_assign', n-1 )
  m_assign = 0._dbl
  DO k = 0, n-2
     j = -(n-1) + 2 * k
     DO l = 0, n-1
        l_arg = 0.5_dbl ** l
        prod_arg = gftype % p3m_coeff ( j, l ) * l_arg
        m_assign ( k ) =  m_assign ( k ) + prod_arg
     END DO
  END DO

! calculate the absolute b values

  npts ( : ) = ub ( : ) - lb ( : ) + 1
  DO dim = 1, 3
     DO pt = lb (dim), ub (dim)
        val = twopi * ( REAL ( pt, dbl ) / REAL ( npts ( dim ), dbl ) )
        exp_m = CMPLX ( COS ( val ), -SIN ( val ), dbl )
        sum_m = CMPLX ( 0._dbl, 0._dbl, dbl )
        DO k = 0, n-2
           sum_m  =  sum_m + m_assign ( k ) * exp_m ** k
        END DO
        b_m = exp_m ** ( n - 1 ) / sum_m
        gftype % p3m_bm2 ( dim, pt ) = SQRT ( REAL ( b_m * CONJG ( b_m ), dbl ) )
     END DO
  END DO

  DEALLOCATE ( m_assign, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( 'influence_factor', 'm_assign', 0 )

END SUBROUTINE influence_factor

!!*****
!******************************************************************************

END MODULE greens_fn

!******************************************************************************
