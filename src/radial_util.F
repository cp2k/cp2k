!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/radial_util [1.0] *
!!
!!   NAME
!!     radial_util
!!
!!   FUNCTION
!!     Utility functions for radial grids
!!
!!   AUTHOR
!!     JGH 17-NOV-2000
!!
!!   MODIFICATION HISTORY
!!     TCH 20-JUN-2002, added a function radial_normalize
!!
!!   SOURCE
!******************************************************************************

MODULE radial_util

  USE kinds, ONLY : dbl
  USE radial_grids, ONLY : radial_grid_type
  USE termination, ONLY : stop_memory
  USE mathconstants, ONLY : pi

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: radial_diff, radial_integrate, radial_normalize

  INTERFACE radial_diff
     MODULE PROCEDURE radial_diff_r, radial_diff_s
  END INTERFACE

  INTERFACE radial_integrate
     MODULE PROCEDURE radial_int_r, radial_int_s
  END INTERFACE

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************

SUBROUTINE radial_int_s ( n, f, h, result, ierr )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: n
  REAL ( dbl ), DIMENSION (1:n), INTENT ( IN ) :: f
  REAL ( dbl ), INTENT ( IN ) :: h
  REAL ( dbl ), INTENT ( OUT ) :: result
  INTEGER, INTENT ( OUT ) :: ierr

! Locals
  REAL ( dbl ), PARAMETER :: p1 = 0.375_dbl, &
                             p2 = 1.16666666666666666667_dbl, &
                             p3 = 0.95833333333333333333_dbl

! Simpson's integration rule '''
!------------------------------------------------------------------------------

  IF ( n < 6 ) THEN
    result = HUGE ( 0._dbl )
    ierr = -1
  ELSE
    result = p1 * ( f ( 1 ) + f ( n ) ) + p2 * ( f ( 2 ) + f ( n-1 ) ) + &
             p3 * ( f ( 3 ) + f ( n-2 ) )
    result = result + SUM ( f ( 4:n-3 ) )
    result = result * h
    ierr = 0
  ENDIF

END SUBROUTINE radial_int_s

!******************************************************************************

SUBROUTINE radial_int_r ( n, f, rg, result, ierr )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: n
  REAL ( dbl ), DIMENSION (1:n), INTENT ( IN ) :: f
  TYPE ( radial_grid_type ), INTENT ( IN ) :: rg
  REAL ( dbl ), INTENT ( OUT ) :: result
  INTEGER, INTENT ( OUT ) :: ierr

! Locals
  REAL ( dbl ), DIMENSION ( : ), ALLOCATABLE :: fs
  INTEGER :: stat

!------------------------------------------------------------------------------

  ALLOCATE ( fs ( 1:n ), STAT = stat )
  IF ( stat /= 0 ) call stop_memory ( "radial_int_r", "fs", n )

  fs(1:n) = f(1:n) / rg%w(1:n) ** 2
  call radial_int_s ( n, fs, rg % h, result, ierr )

  DEALLOCATE ( fs, STAT = stat )
  IF ( stat /= 0 ) call stop_memory ( "radial_int_r", "fs" )

END SUBROUTINE radial_int_r

!******************************************************************************

SUBROUTINE radial_diff_s ( n, f, df, order, stat )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: n
  REAL ( dbl ), DIMENSION (1:n), INTENT ( IN ) :: f
  REAL ( dbl ), DIMENSION (1:n), INTENT ( OUT ) :: df
  INTEGER, INTENT ( IN ) :: order
  INTEGER, INTENT ( OUT ) :: stat

! Locals
  REAL ( dbl ), PARAMETER :: t = 0.33333333333333333333_dbl, &
                             z = 0.08333333333333333333_dbl
  INTEGER :: i, o

! Differentiation using Lagrange's 5 point formula    '''
!------------------------------------------------------------------------------

    o = order

  IF ( n < 5 ) THEN
    df = 0._dbl
    stat = -1
  ELSE
! differentiation wrt to s ( the current grid )
    SELECT CASE ( o )
    CASE DEFAULT
      df = 0._dbl
      stat = 1
    CASE ( 0 )
      df = f
      stat = 0
    CASE ( 1 )
      df ( 1 ) = -25._dbl * z * f ( 1 ) + 4._dbl * f ( 2 ) &
                 - 3._dbl * f ( 3 ) + 4._dbl * t * f ( 4 ) &
                 - 0.25_dbl * f ( 5 )
      df ( 2 ) = -0.25_dbl * f ( 1 ) - 10._dbl * z * f ( 2 ) &
                 + 1.5_dbl * f ( 3 ) - 0.5_dbl * f ( 4 ) &
                 + z * f ( 5 )
      DO i = 3, n-2
        df ( i ) = z * f ( i-2 ) - 2._dbl * t * f ( i-1 ) + &
                   2._dbl * t * f ( i+1 ) - z * f ( i+2 )
      END DO
      df ( n-1 ) = -z * f ( n-4 ) + 0.5_dbl * f ( n-3 ) &
                   - 1.5_dbl * f ( n-2 ) + 10._dbl * z * f ( n-1 ) &
                   + 0.25_dbl * f ( n )
      df ( n ) = 0.25_dbl * f ( n-4 ) - 4._dbl * t * f ( n-3 ) &
                 + 3._dbl * f ( n-2 ) - 4._dbl * f ( n-1 ) &
                 + 25._dbl * z * f ( n )
      stat = 0
    CASE ( 2 )
      df ( 1 ) = 35._dbl * z * f ( 1 ) - 26._dbl * t * f ( 2 ) &
                 + 9.5_dbl * f ( 3 ) - 14._dbl * t * f ( 4 ) &
                 + 11._dbl * z * f ( 5 )
      df ( 2 ) = 11._dbl * z * f ( 1 ) - 5._dbl * t * f ( 2 ) &
                 + 0.5_dbl * f ( 3 ) + t * f ( 4 ) - z * f ( 5 )
      DO i = 3, n-2
        df ( i ) = -z * f ( i-2 ) + 4._dbl * t * f ( i-1 ) &
                   - 2.5_dbl * f ( i ) + 4._dbl * t * f ( i+1 ) &
                   - z * f ( i+2 )
      END DO
      df ( n-1 ) = -z * f ( n-4 ) + t * f ( n-3 ) + 0.5_dbl * f ( n-2 ) &
                   - 5._dbl * t * f ( n-1 ) + 11._dbl * t * f ( n )
      df ( n ) = 11._dbl * z * f ( n-4 ) - 14._dbl * t * f ( n-3 ) &
                 + 9.5_dbl * f ( n-2 ) - 26._dbl * t * f ( n-1 ) &
                 + 35._dbl * z * f ( n )
      stat = 0
    END SELECT
  ENDIF

END SUBROUTINE radial_diff_s

!******************************************************************************

SUBROUTINE radial_diff_r ( n, f, df, order, rg, stat )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: n
  REAL ( dbl ), DIMENSION (1:n), INTENT ( IN ) :: f
  REAL ( dbl ), DIMENSION (1:n), INTENT ( OUT ) :: df
  INTEGER, INTENT ( IN ) :: order
  TYPE ( radial_grid_type ), INTENT ( IN ) :: rg
  INTEGER, INTENT ( OUT ) :: stat

! Locals
  INTEGER :: i, o, ierr
  REAL ( dbl ) :: oh
  REAL ( dbl ), DIMENSION ( : ), ALLOCATABLE :: fs

!------------------------------------------------------------------------------

  o = order

  IF ( n < 5 ) THEN
    df = 0._dbl
    stat = -1
  ELSE
! differentiation wrt to r
    oh = 1._dbl / rg % h
    SELECT CASE ( o )
    CASE DEFAULT
      df = 0._dbl
      stat = 1
    CASE ( 0 )
      df = f
      stat = 0
    CASE ( 1 )
      CALL radial_diff_s ( n, f, df, 1, STAT=stat )
      IF ( stat == 0 ) THEN
        df = oh * df * rg % w ** 2
      ENDIF
    CASE ( 2 )
      ALLOCATE ( fs ( 1:n ), STAT = ierr )
      IF ( ierr /= 0 ) call stop_memory ( "radial_diff_r", "fs", n )
      fs = f * rg % w
      CALL radial_diff_s ( n, fs, df, 2, STAT=stat )
      IF ( stat == 0 ) THEN
        df = oh * oh * rg % w ** 3 * ( df - rg % wdd * f )
      ENDIF
      DEALLOCATE ( fs, STAT = ierr )
      IF ( ierr /= 0 ) call stop_memory ( "radial_diff_r", "fs" )
    END SELECT
  ENDIF

END SUBROUTINE radial_diff_r

! **************************************************************************

!!****f* cp2k/radial_util/radial_normalize *
!!
!!   NAME
!!     radial_normalize
!!
!!   FUNCTION
!!     Normalizes a function that dependa only on r (the magnitude
!!     of a vector in 3-dim space) to unity
!!
!!   ARGUMENTS
!!     f - REAL(DBL), DIMENSION(:) - the function (vector) to be normalized
!!     rg - TYPE(RADIAL_GRID_TYPE) - the radial grid on which the function
!!                                   is defined
!!
!!   MODIFICATION HISTORY
!!     20-JUN-2002, TCH, working version
!!
!!   SEE ALSO
!!     radial_grids
!!
!!****
  SUBROUTINE radial_normalize(f, rg)

    IMPLICIT NONE

    REAL ( DBL ), DIMENSION ( 1: ), INTENT ( INOUT ) :: f
    TYPE ( RADIAL_GRID_TYPE ), INTENT ( IN ) :: rg

    REAL ( DBL ) :: res
    INTEGER      :: errno, np

    np = rg%n

    CALL radial_integrate ( np, rg%r(1:np)**2 * f(1:np)**2, rg, res, errno )
    ! should do some errno checking
    CALL DSCAL(np, 1._dbl/SQRT(4._dbl*pi*res), f(1:np), 1)
    
  END SUBROUTINE radial_normalize

!******************************************************************************

END MODULE radial_util

!******************************************************************************
