!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/radial_grids [1.0] *
!!
!!   NAME
!!     radial_grids
!!
!!   FUNCTION
!!     Radial grid transformation from
!!       D. Andrae and J. Hinze, IJQC Vol. 63 pp 65 (1997)
!!
!!   AUTHOR
!!     JGH 17-NOV-2000
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE radial_grids

  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: uppercase
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: radial_grid_type, init_radial_grid

  TYPE radial_grid_type
    CHARACTER (LEN=30) :: TYPE
    CHARACTER (LEN=60) :: form
    INTEGER :: n
    REAL(KIND=dp) :: h
    LOGICAL :: closed
    REAL(KIND=dp) :: z,b(5)
    REAL(KIND=dp) :: rmax
    REAL(KIND=dp) :: q
    REAL(KIND=dp), DIMENSION ( : ), POINTER :: r, w, wdd
  END TYPE radial_grid_type

!!*****
!-----------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------!
!!****** radial_grids/init_radial_grid [1.0] *
!!
!!   NAME
!!     init_radial_grid
!!
!!   FUNCTION
!!     Initialization routine for radial grids
!!
!!   AUTHOR
!!     JGH 17-NOV-2000
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!-----------------------------------------------------------------------------!

  SUBROUTINE init_radial_grid ( rg, TYPE, n, z, b, rmax )

    TYPE(radial_grid_type), INTENT(OUT)      :: rg
    CHARACTER(LEN=*), INTENT(IN)             :: TYPE
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), INTENT(IN)                :: z
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(IN), OPTIONAL                   :: b
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: rmax

    INTEGER                                  :: i, ierr, m
    REAL(KIND=dp)                            :: ab, abt, p, s, t, t0, tm, &
                                                ts2, tt

!------------------------------------------------------------------------------

  ALLOCATE ( rg % r ( 0:n+1 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "init_radial_grid", "rg % r", n+2 )
  ALLOCATE ( rg % w ( 0:n+1 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "init_radial_grid", "rg % w", n+2 )
  ALLOCATE ( rg % wdd ( 0:n+1 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "init_radial_grid", "rg % wdd", n+2 )

  rg % type = TYPE
  CALL uppercase ( rg % type )
  rg % n = n
  rg % z = z
  rg % h = 1.0_dp / REAL ( n + 1,KIND=dp)
  IF ( PRESENT ( rmax ) ) THEN
    rg % rmax = rmax
  ELSE
    rg % rmax = 100.0_dp
  END IF
  IF ( PRESENT ( b ) ) THEN
    rg % b = 0.0_dp
    m = SIZE ( b )
    rg % b ( 1:m ) = b ( 1:m )
  ELSE
    rg % b = 1.0_dp
  END IF

  SELECT CASE ( rg % type )

  CASE DEFAULT
    CALL stop_program ( "init_radial_grid", "Variable transformation unknown" )

  CASE ( "RATIONAL" )
    rg % form = "(b/Z)[Ts/(1-Ts)]"
    t0 = 0.0_dp
    tm = rg % z * rg % rmax / ( rg % z * rg % rmax + rg % b(1) )
    t = tm - t0
    ab = rg % z / rg % b(1)
    abt = SQRT ( ab / t )
    DO i = 0, n+1
      s = REAL ( i,KIND=dp) / REAL ( n + 1,KIND=dp)
      rg % r ( i ) = ( t * s ) / ( 1.0_dp - t * s ) / ab
      rg % w ( i ) = abt * ( 1.0_dp - t * s )
    END DO
    rg % wdd ( 0:n+1 ) = 0.0_dp
    rg % q = 1.0_dp
    rg % closed = .TRUE.

  CASE ( "ALGEBRAIC SQUARE ROOT" )
    rg % form = "(b/Z)[Ts/SQRT(1-(Ts)^2)]"
    ab = rg % z / rg % b(1)
    t0 = 0.0_dp
    tm = SIN ( ATAN ( ab * rg % rmax ) )
    t = tm - t0
    abt = SQRT ( ab / t )
    DO i = 0, n+1
      s = REAL ( i,KIND=dp) / REAL ( n + 1,KIND=dp)
      ts2 = t * t * s * s
      rg % r ( i ) = ( t * s ) / SQRT ( 1.0_dp - ts2 ) / ab
      rg % w ( i ) = abt * ( 1.0_dp - ts2 )**(0.75_dp)
      rg % wdd ( i ) = -0.75_dp * t * t * abt * &
       ( 2.0_dp - ts2 ) / ( 1.0_dp - ts2 )**(1.25_dp)
    END DO
    rg % q = 1.0_dp
    rg % closed = .TRUE.

  CASE ( "LINEAR" )
    rg % form = "Rmax*s"
    DO i = 0, n+1
      s = REAL ( i,KIND=dp) / REAL ( n + 1,KIND=dp)
      rg % r ( i ) = s * rmax
    END DO
    rg % w ( 0:n+1 ) = SQRT ( 1.0_dp / rg % rmax )
    rg % wdd ( 0:n+1 ) = 0.0_dp
    rg % q = 1.0_dp
    rg % closed = .FALSE.

  CASE ( "POWER" )
    rg % form = "(b/Z)[(Ts/b^p + 1)^(1/p) - 1]"
    p = rg % b(2)
    ab = rg % z / rg % b(1)
    t0 = rg % b(1) ** p
    tm = ( z * rg % rmax + rg % b(1) ) ** p
    t = tm - t0
    abt = SQRT ( p * rg % b(1)**p * ab / t )
!bug : check for the correct value of s!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CALL stop_program("set the value of s","radial_grids.F")
    tt = ( ( t * s ) / rg % b(1) ** p + 1.0_dp )
    DO i = 0, n+1
      s = REAL ( i,KIND=dp) / REAL ( n + 1,KIND=dp)
      rg % r ( i ) = ( tt ** (1.0_dp/p) - 1.0_dp ) / ab
      rg % w ( i ) = abt * tt ** ((p-1.0_dp)/(2.0_dp*p))
      rg % wdd ( i ) = -0.25_dp * abt * ( p - 1.0_dp ) * ( p + 1.0_dp ) * &
              tt ** ((p-1.0_dp)/(2.0_dp*p)-2.0_dp)
    END DO
    rg % q = 1.0_dp
    rg % closed = .FALSE.

  CASE ( "LOG" )
    rg % form = "(b/Z)[EXP(Ts) - 1]"
    ab = rg % z / rg % b(1)
    t0 = LOG ( rg % b(1) )
    tm = LOG ( z * rg % rmax + rg % b(1) )
    t = tm - t0
    abt = SQRT ( ab / t )
    DO i = 0, n+1
      s = REAL ( i,KIND=dp) / REAL ( n + 1,KIND=dp)
      rg % r ( i ) = ( EXP ( t * s ) - 1.0_dp ) / ab
      rg % w ( i ) = abt * EXP ( -0.5_dp * t * s )
      rg % wdd ( i ) = 0.25_dp * t * t * abt * EXP ( -0.5_dp * t * s )
    END DO
    rg % q = 1.0_dp
    rg % closed = .FALSE.

  END SELECT

  END SUBROUTINE init_radial_grid

!!*****
!******************************************************************************

END MODULE radial_grids

!******************************************************************************
