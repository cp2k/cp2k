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

  USE kinds,                           ONLY: dbl
  USE string_utilities,                ONLY: uppercase
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: radial_grid_type, init_radial_grid

  TYPE radial_grid_type
    CHARACTER (LEN=30) :: type
    CHARACTER (LEN=60) :: form
    INTEGER :: n
    REAL(dbl) :: h
    LOGICAL :: closed
    REAL(dbl) :: z,b(5)
    REAL(dbl) :: rmax
    REAL(dbl) :: q
    REAL(dbl), DIMENSION ( : ), POINTER :: r, w, wdd
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

  SUBROUTINE init_radial_grid ( rg, type, n, z, b, rmax )

    TYPE(radial_grid_type), INTENT(OUT)      :: rg
    CHARACTER(LEN=*), INTENT(IN)             :: type
    INTEGER, INTENT(IN)                      :: n
    REAL(dbl), INTENT(IN)                    :: z
    REAL(dbl), DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: b
    REAL(dbl), INTENT(IN), OPTIONAL          :: rmax

    INTEGER                                  :: i, ierr, m
    REAL(dbl)                                :: ab, abt, p, s, t, t0, tm, &
                                                ts2, tt

!------------------------------------------------------------------------------

  ALLOCATE ( rg % r ( 0:n+1 ), STAT = ierr )
  IF ( ierr /= 0 ) call stop_memory ( "init_radial_grid", "rg % r", n+2 )
  ALLOCATE ( rg % w ( 0:n+1 ), STAT = ierr )
  IF ( ierr /= 0 ) call stop_memory ( "init_radial_grid", "rg % w", n+2 )
  ALLOCATE ( rg % wdd ( 0:n+1 ), STAT = ierr )
  IF ( ierr /= 0 ) call stop_memory ( "init_radial_grid", "rg % wdd", n+2 )

  rg % type = type
  CALL uppercase ( rg % type )
  rg % n = n
  rg % z = z
  rg % h = 1._dbl / REAL ( n + 1, dbl )
  IF ( PRESENT ( rmax ) ) THEN
    rg % rmax = rmax
  ELSE
    rg % rmax = 100._dbl
  END IF
  IF ( PRESENT ( b ) ) THEN
    rg % b = 0._dbl
    m = SIZE ( b )
    rg % b ( 1:m ) = b ( 1:m )
  ELSE
    rg % b = 1._dbl
  END IF

  SELECT CASE ( rg % type )

  CASE DEFAULT
    call stop_program ( "init_radial_grid", "Variable transformation unknown" )

  CASE ( "RATIONAL" )
    rg % form = "(b/Z)[Ts/(1-Ts)]"
    t0 = 0._dbl
    tm = rg % z * rg % rmax / ( rg % z * rg % rmax + rg % b(1) )
    t = tm - t0
    ab = rg % z / rg % b(1)
    abt = SQRT ( ab / t )
    DO i = 0, n+1
      s = REAL ( i, dbl ) / REAL ( n + 1, dbl )
      rg % r ( i ) = ( t * s ) / ( 1._dbl - t * s ) / ab
      rg % w ( i ) = abt * ( 1._dbl - t * s )
    END DO
    rg % wdd ( 0:n+1 ) = 0._dbl
    rg % q = 1._dbl
    rg % closed = .TRUE.

  CASE ( "ALGEBRAIC SQUARE ROOT" )
    rg % form = "(b/Z)[Ts/SQRT(1-(Ts)^2)]"
    ab = rg % z / rg % b(1)
    t0 = 0._dbl
    tm = SIN ( ATAN ( ab * rg % rmax ) )
    t = tm - t0
    abt = SQRT ( ab / t )
    DO i = 0, n+1
      s = REAL ( i, dbl ) / REAL ( n + 1, dbl )
      ts2 = t * t * s * s
      rg % r ( i ) = ( t * s ) / SQRT ( 1._dbl - ts2 ) / ab
      rg % w ( i ) = abt * ( 1._dbl - ts2 )**(0.75_dbl)
      rg % wdd ( i ) = -0.75_dbl * t * t * abt * &
       ( 2._dbl - ts2 ) / ( 1._dbl - ts2 )**(1.25_dbl)
    END DO
    rg % q = 1._dbl
    rg % closed = .TRUE.

  CASE ( "LINEAR" )
    rg % form = "Rmax*s"
    DO i = 0, n+1
      s = REAL ( i, dbl ) / REAL ( n + 1, dbl )
      rg % r ( i ) = s * rmax
    END DO
    rg % w ( 0:n+1 ) = SQRT ( 1._dbl / rg % rmax )
    rg % wdd ( 0:n+1 ) = 0._dbl
    rg % q = 1._dbl
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
    tt = ( ( t * s ) / rg % b(1) ** p + 1._dbl )
    DO i = 0, n+1
      s = REAL ( i, dbl ) / REAL ( n + 1, dbl )
      rg % r ( i ) = ( tt ** (1._dbl/p) - 1._dbl ) / ab
      rg % w ( i ) = abt * tt ** ((p-1._dbl)/(2._dbl*p))
      rg % wdd ( i ) = -0.25_dbl * abt * ( p - 1._dbl ) * ( p + 1._dbl ) * &
              tt ** ((p-1._dbl)/(2._dbl*p)-2._dbl)
    END DO
    rg % q = 1._dbl
    rg % closed = .FALSE.

  CASE ( "LOG" )
    rg % form = "(b/Z)[EXP(Ts) - 1]"
    ab = rg % z / rg % b(1)
    t0 = log ( rg % b(1) )
    tm = log ( z * rg % rmax + rg % b(1) )
    t = tm - t0
    abt = SQRT ( ab / t )
    DO i = 0, n+1
      s = REAL ( i, dbl ) / REAL ( n + 1, dbl )
      rg % r ( i ) = ( EXP ( t * s ) - 1._dbl ) / ab
      rg % w ( i ) = abt * EXP ( -0.5_dbl * t * s )
      rg % wdd ( i ) = 0.25_dbl * t * t * abt * EXP ( -0.5_dbl * t * s )
    END DO
    rg % q = 1._dbl
    rg % closed = .FALSE.

  END SELECT

  END SUBROUTINE init_radial_grid

!!*****
!******************************************************************************

END MODULE radial_grids

!******************************************************************************
