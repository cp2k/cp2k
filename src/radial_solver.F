!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/radial_solver [1.0] *
!!
!!   NAME
!!     radial_solver
!!
!!   FUNCTION
!!     Solvers for radial differential equations (second order)
!!
!!   AUTHOR
!!     JGH 17-NOV-2000
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE radial_solver

  USE kinds, ONLY : dbl
  USE termination, ONLY : stop_memory, stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: numerov, diff5p

!!*****
!-----------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------!
!!****** radial_solver/numerov [1.0] *
!!
!!   NAME
!!     numerov
!!
!!   FUNCTION
!!     Numerov solver
!!
!!   AUTHOR
!!     JGH 17-NOV-2000
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!-----------------------------------------------------------------------------!

  SUBROUTINE numerov ( f, g, w, corr )

!Arguments
  REAL ( dbl ), DIMENSION ( 0: ), INTENT ( IN ) :: f
  REAL ( dbl ), DIMENSION ( 0: ), INTENT ( INOUT ) :: g
  REAL ( dbl ), DIMENSION ( 0: ), INTENT ( IN ) :: w
  REAL ( dbl ), DIMENSION ( : ), INTENT ( IN ) :: corr

!Locals
  INTEGER :: n, ierr, info, i
  REAL ( dbl ) :: e, ww
  INTEGER, DIMENSION ( : ), ALLOCATABLE :: ipiv
  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: ab
  REAL ( dbl ), DIMENSION ( : ), ALLOCATABLE :: s

!------------------------------------------------------------------------------

  n = SIZE ( f ) - 2
  ALLOCATE ( ab ( 1:4, 1:n ), STAT = ierr )
  IF ( ierr /= 0 ) call stop_memory ( "numerov", "ab", 4*n )
  ALLOCATE ( s ( 1:n ), STAT = ierr )
  IF ( ierr /= 0 ) call stop_memory ( "numerov", "s", n )
  ALLOCATE ( ipiv ( 1:n ), STAT = ierr )
  IF ( ierr /= 0 ) call stop_memory ( "numerov", "ipiv", n )

  ab = 0._dbl
  DO i = 1, n
    s ( i ) = g ( i - 1 ) + 10._dbl * g ( i ) + g ( i + 1 )
    e = 1._dbl + f ( i )
    ab ( 2, i ) = e
    ab ( 3, i ) = 10._dbl * e - 12._dbl
    ab ( 4, i ) = e
  END DO
  ab ( 3, 1 ) = ab ( 3, 1 ) + corr ( 1 )
  ab ( 2, 2 ) = ab ( 2, 2 ) + corr ( 2 )
  s ( n ) = s ( n ) - corr ( 3 ) * ( 1._dbl + f ( n+1 ) )

  DO i = 1, n
    ww = w ( i )
    ab ( 2, i ) = ab ( 2, i ) * ww
    ab ( 3, i ) = ab ( 3, i ) * ww
    ab ( 4, i ) = ab ( 4, i ) * ww
  END DO

  CALL dgbsv ( n, 1, 1, 1, ab, 4, ipiv, s, n, info )
  IF ( info /= 0 ) call stop_program ( "numerov", "DGBSV: info /= 0" )
  g ( 1:n ) = s ( 1:n )

  DEALLOCATE ( ab, STAT = ierr )
  IF ( ierr /= 0 ) call stop_memory ( "numerov", "ab" )
  DEALLOCATE ( s, STAT = ierr )
  IF ( ierr /= 0 ) call stop_memory ( "numerov", "s" )
  DEALLOCATE ( ipiv, STAT = ierr )
  IF ( ierr /= 0 ) call stop_memory ( "numerov", "ipiv" )

  END SUBROUTINE numerov

!!*****
!******************************************************************************
!!****** radial_solver/diff5p [1.0] *
!!
!!   NAME
!!     diff5p
!!
!!   FUNCTION
!!     Lagrange differentiation solver
!!
!!   AUTHOR
!!     JGH 17-NOV-2000
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!-----------------------------------------------------------------------------!

  SUBROUTINE diff5p ( f, g, w )

!Arguments
  REAL ( dbl ), DIMENSION ( 0: ), INTENT ( IN ) :: f
  REAL ( dbl ), DIMENSION ( 0: ), INTENT ( INOUT ) :: g
  REAL ( dbl ), DIMENSION ( 0: ), INTENT ( IN ) :: w

!Locals
  INTEGER :: n, ierr, info, i
  REAL ( dbl ) :: ww
  INTEGER, DIMENSION ( : ), ALLOCATABLE :: ipiv
  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: ab

!------------------------------------------------------------------------------

  n = SIZE ( f ) - 2

  END SUBROUTINE diff5p

!!*****
!******************************************************************************

END MODULE radial_solver

!******************************************************************************
