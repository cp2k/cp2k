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
  USE stop_program, ONLY : stop_prg, stop_memory
  
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

  CALL DGBSV ( n, 1, 1, 1, ab, 4, ipiv, s, n, info )
  IF ( info /= 0 ) call stop_prg ( "numerov", "DGBSV: info /= 0" )
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
!deb
if(1 == 0) then
!deb
  ALLOCATE ( ab ( 1:7, 1:n ), STAT = ierr )
  IF ( ierr /= 0 ) call stop_memory ( "diff5p", "ab", 7*n )
  ALLOCATE ( ipiv ( 1:n ), STAT = ierr )
  IF ( ierr /= 0 ) call stop_memory ( "diff5p", "ipiv", n )

  ab = 0._dbl
  DO i = 1, n
    ab ( 3, i ) = -5._dbl
    ab ( 4, i ) = 80._dbl
    ab ( 5, i ) = -150._dbl + f ( i )
    ab ( 6, i ) = 80._dbl
    ab ( 7, i ) = -5._dbl
  END DO

  ab ( 5, 1 ) = f ( 1 ) + 720._dbl
  ab ( 4, 2 ) = 60._dbl + 10._dbl * ab ( 5, 2 )
  ab ( 3, 3 ) = 870._dbl + ab ( 5, 3)
  ab ( 5, n ) = f ( n ) + 720._dbl
  ab ( 6, n - 1 ) = 60._dbl + 10._dbl * ab ( 5, n - 1 )
  ab ( 7, n - 2 ) = 870._dbl + ab ( 5, n - 2 )

  g ( 1 ) = g ( 1 ) + 10._dbl * g ( 2 ) + g ( 3 )
  g ( n-1) = g ( n-1) - 5._dbl * w(n-1)
  g ( n ) = g ( n ) + 10._dbl * g ( n - 1 ) + g ( n - 2 )

  DO i = 1, n
    ww = w ( i )
    ab ( 3, i ) = ab ( 3, i ) * ww
    ab ( 4, i ) = ab ( 4, i ) * ww
    ab ( 5, i ) = ab ( 5, i ) * ww
    ab ( 6, i ) = ab ( 6, i ) * ww
    ab ( 7, i ) = ab ( 7, i ) * ww
  END DO

  CALL DGBSV ( n, 2, 2, 1, ab, 7, ipiv, g(1:n), n, info )
  IF ( info /= 0 ) call stop_prg ( "diff5p", "DGBSV: info /= 0" )

!deb
else
!deb
  ALLOCATE ( ab ( 1:13, 1:n ), STAT = ierr )
  IF ( ierr /= 0 ) call stop_memory ( "diff5p", "ab", 13*n )
  ALLOCATE ( ipiv ( 1:n ), STAT = ierr )
  IF ( ierr /= 0 ) call stop_memory ( "diff5p", "ipiv", n )

  ab = 0._dbl
  DO i = 1, n
    ab ( 7, i ) = -5._dbl
    ab ( 8, i ) = 80._dbl
    ab ( 9, i ) = -150._dbl + f ( i )
    ab (10, i ) = 80._dbl
    ab (11, i ) = -5._dbl
  END DO

  ab ( 9, 1 ) = 175._dbl + f ( 1 )
  ab ( 8, 2 ) = -520._dbl
  ab ( 7, 3 ) = 570._dbl
  ab ( 6, 4 ) = -280._dbl
  ab ( 5, 5 ) = 55._dbl
  ab (10, 1 ) = 55._dbl
  ab ( 9, 2 ) = -100._dbl + f ( 2 )
  ab ( 8, 3 ) = 30._dbl
  ab ( 7, 4 ) = 20._dbl
  ab ( 6, 5 ) = -5._dbl

  ab ( 8, n )   = 55._dbl
  ab ( 9, n-1 ) = -100._dbl + f ( n-1 )
  ab (10, n-2 ) = 30._dbl
  ab (11, n-3 ) = 20._dbl
  ab (12, n-4 ) = -5._dbl
  ab ( 9, n )   = 175._dbl + f ( n )
  ab (10, n-1 ) = -520._dbl
  ab (11, n-2 ) = 570._dbl
  ab (12, n-3 ) = -280._dbl
  ab (13, n-4 ) = 55._dbl

  DO i = 1, n
    ww = w ( i )
    ab ( 5, i ) = ab ( 5, i ) * ww
    ab ( 6, i ) = ab ( 6, i ) * ww
    ab ( 7, i ) = ab ( 7, i ) * ww
    ab ( 8, i ) = ab ( 8, i ) * ww
    ab ( 9, i ) = ab ( 9, i ) * ww
    ab (10, i ) = ab (10, i ) * ww
    ab (11, i ) = ab (11, i ) * ww
    ab (12, i ) = ab (12, i ) * ww
    ab (13, i ) = ab (13, i ) * ww
  END DO

!deb
  write(6,*) ab (13,n-4),ab(12,n-3),ab(11,n-2),ab(10,n-1),ab(9,n)
!deb

  CALL DGBSV ( n, 4, 4, 1, ab,13, ipiv, g(1:n), n, info )
  IF ( info /= 0 ) call stop_prg ( "diff5p", "DGBSV: info /= 0" )
!deb
ENDIF
!deb


  DEALLOCATE ( ab, STAT = ierr )
  IF ( ierr /= 0 ) call stop_memory ( "diff5p", "ab" )
  DEALLOCATE ( ipiv, STAT = ierr )
  IF ( ierr /= 0 ) call stop_memory ( "diff5p", "ipiv" )

  END SUBROUTINE diff5p

!!*****
!******************************************************************************

END MODULE radial_solver

!******************************************************************************
