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

  USE kinds,                           ONLY: dbl
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: numerov, diff6p

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
!!     TCH 17-OCT-2001: - changed the way, how f and g are passed 
!!                        to the solver. they must no longer be 
!!                        multiplied by the factor h^2/12
!!
!!   SOURCE
!-----------------------------------------------------------------------------!

  SUBROUTINE numerov ( f, g, w, corr )


    REAL(dbl), DIMENSION(0:), INTENT(IN)     :: f
    REAL(dbl), DIMENSION(0:), INTENT(INOUT)  :: g
    REAL(dbl), DIMENSION(0:), INTENT(IN)     :: w
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: corr

    INTEGER                                  :: i, ierr, info, n
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ipiv
    REAL(dbl)                                :: e, h212, ww
    REAL(dbl), ALLOCATABLE, DIMENSION(:)     :: s
    REAL(dbl), ALLOCATABLE, DIMENSION(:, :)  :: ab

! the potential (in s-form)
! first the inhomegenity (in s-form) 
! later the solution
! the weights W
! corrections for the boundaries
! square of stepsize divided by 12
! the coefficient matrix
! solution
!------------------------------------------------------------------------------

    n = SIZE ( f ) - 2

    ALLOCATE ( ab ( 1:4, 1:n ), STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "numerov", "ab", 4*n )
    ALLOCATE ( s ( 1:n ), STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "numerov", "s", n )
    ALLOCATE ( ipiv ( 1:n ), STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "numerov", "ipiv", n )

    ! the stepsize is 1/(n+2-1)
    h212 = ((1._dbl/DBLE(n+1))**2)/12._dbl

    ab = 0._dbl
    DO i = 1, n
       
       ! g(0) and g(n+1) must be 0 for this!
       s ( i ) = (g(i-1) + 10._dbl*g(i) + g(i+1))*h212
       e = 1._dbl + f(i)*h212

       ! (i,j) -> (3+i-j,j)
       ab ( 2, i ) = e
       ab ( 3, i ) = 10._dbl * e - 12._dbl
       ab ( 4, i ) = e

    END DO

    ! elements that do not exist
    ab(2,1) = 0._dbl
    ab(4,n) = 0._dbl

    ! elements that need corrections
    ab ( 3, 1 ) = ab ( 3, 1 ) + corr ( 1 )
    ab ( 2, 2 ) = ab ( 2, 2 ) + corr ( 2 )
    s(n) = s(n) - corr(3)*(1._dbl + h212*f(n+1))

    ! finally multiply each column of a with the corresponding weight
    DO i = 1, n
       ww = w ( i )
       ab ( 2, i ) = ab ( 2, i ) * ww
       ab ( 3, i ) = ab ( 3, i ) * ww
       ab ( 4, i ) = ab ( 4, i ) * ww
    END DO

    CALL dgbsv ( n, 1, 1, 1, ab, 4, ipiv, s, n, info )
    IF ( info /= 0 ) CALL stop_program ( "numerov", "DGBSV: info /= 0" )
    g ( 1:n ) = s ( 1:n )

    DEALLOCATE ( ab, STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "numerov", "ab" )
    DEALLOCATE ( s, STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "numerov", "s" )
    DEALLOCATE ( ipiv, STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "numerov", "ipiv" )

  END SUBROUTINE numerov

!!*****
!******************************************************************************
!!****** radial_solver/diff6p [1.0] *
!!
!!   NAME
!!     diff6p
!!
!!   FUNCTION
!!     6 point Lagrange differentiation solver
!!     solves (in s-space):
!!        [d^2/ds^2 + f]wP = g for P.
!!     The solution must be 0 at the origin! If the solution
!!     is known to decay only slowly, then a correction is required
!!     to compensate for a finite value of the solution at the 
!!     upper boundary.
!!
!!   AUTHOR
!!     JGH 17-NOV-2000
!!
!!   MODIFICATION HISTORY
!!     TCH 09-OCT-2001: - implemented the algo 
!!                      - changed name from diff5p to diff6p
!!     TCH 16-OCT-2001: - implemented the upper boundary correction
!!                      - added some on function description
!!     TCH 17-OCT-2001: - added the correction term (for long range
!!                        potentials) and made it an optional argument
!!
!!   SOURCE
!-----------------------------------------------------------------------------!

  SUBROUTINE diff6p ( f, g, w, corr )


    REAL(DBL), DIMENSION(0:), INTENT(in)     :: f
    REAL(DBL), DIMENSION(0:), INTENT(INOUT)  :: g
    REAL(DBL), DIMENSION(0:), INTENT(in)     :: w
    REAL(DBL), INTENT(IN), OPTIONAL          :: corr

    INTEGER                                  :: i, ierr, info, j, n
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ipiv
    REAL(DBL)                                :: h2_60, ww
    REAL(DBL), ALLOCATABLE, DIMENSION(:, :)  :: ab

! as in the equation above (in s-form)
! on input an inhomogenity on output
! the solution 
! weights from the variable transformation
! correction for the case that
! y(n+1) /= 0
! value should be w(n+1)*y(n+1)
! loop variables (i: row, j: column)
! factor needed in calc
! vector of pivot elements for dgbsv
! matrix to solve
!------------------------------------------------------------------------------

    n = SIZE(f) - 2

    ALLOCATE(ab(1:7,1:n), stat=ierr)
    IF (ierr /= 0) CALL stop_memory("radial_solver/diff6p", "ab", 7*n)

    ALLOCATE(ipiv(1:n), stat=ierr)
    IF (ierr /= 0) CALL stop_memory("radial_solver/diff6p", "ipiv", n)
    
    ab = 0._dbl
    ipiv = 0

    ! in s-space the interval is from 0 to 1
    h2_60 = ((1._dbl/DBLE(n+1))**2)*60._dbl

    ! ab(11+i-j,j) <- a(i,j)
    DO i=3, n-2
       ab(7,i-2) = -5._dbl 
       ab(6,i-1) = 80._dbl
       ab(5,i)   = -150._dbl
       ab(4,i+1) = 80._dbl
       ab(3,i+2) = -5._dbl
    END DO
    
    ab(6,1)=80._dbl; ab(5,2)=-150._dbl; ab(4,3)=80._dbl; ab(3,4)=-5._dbl 
    ab(4,n)=80._dbl; ab(5,n-1)=-150._dbl; ab(6,n-2)=80._dbl; ab(7,n-3)=-5._dbl

    ! special treatment of first and last row to make the matrix banded 
    ! with 2 lower and upper diagonals
    ab(5,1) = 720._dbl
    ab(4,2) = -1440._dbl + 10._dbl*h2_60*f(2)
    ab(3,3) = 720._dbl + h2_60*f(3)

    ab(5,n) = 720._dbl
    ab(6,n-1) = -1440._dbl + 10._dbl*h2_60*f(n-1)
    ab(7,n-2) = 720._dbl + h2_60*f(n-2)

    g(1) = g(1) + 10._dbl*g(2) + g(3)
    g(n) = g(n) + 10._dbl*g(n-1) + g(n-2)

    ! add the function f to the diagonal element (i=j) and multiply by correct factor
    DO j=1, n
       ab(5,j) = ab(5,j) + h2_60*f(j)
       g(j) = h2_60*g(j)
    END DO

    ! apply the correction
    IF (PRESENT(corr)) g(n-1) = g(n-1) + 5._dbl*corr

    ! now multiply each column of the matrix a with the corresponding weight
    DO j=1, n
       ww = w(j)
       DO i=j-2, j+2
          ab(5+i-j,j) = ab(5+i-j,j) * ww
       END DO
    END DO

    ! call LAPACK routine here
    CALL dgbsv(n,2,2,1,ab,7,ipiv,g(1:n),n,info)
    IF (info /= 0) CALL stop_program("radial_solver/diff6p", "DGBSV: info /= 0")

    ! g(1:n) now contains the solution to the problem

    DEALLOCATE(ab, stat=ierr)
    IF (ierr /= 0) CALL stop_memory("radial_solver/diff6p", "ab")

    DEALLOCATE(ipiv, stat=ierr)
    IF (ierr /= 0) CALL stop_memory("radial_solver/diff6p", "ipiv")

  END SUBROUTINE diff6p

!!*****
!******************************************************************************

END MODULE radial_solver

!******************************************************************************
