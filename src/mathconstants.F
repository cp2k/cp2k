!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/mathconstants [1.0] *
!!
!!   NAME
!!     mathconstants
!!
!!   FUNCTION
!!     Definition of mathematical constants and functions.
!!
!!   AUTHOR
!!     Matthias Krack
!!
!!   MODIFICATION HISTORY
!!     Adapted for use in CP2K (JGH)
!!
!!   SOURCE
!******************************************************************************

MODULE mathconstants
  
  USE kinds, ONLY : dbl
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: init_mathcon
  PUBLIC :: maxfac, pi, pio2, twopi, fourpi, root2, rootpi, oorootpi, &
       zero, one, half, degree, radians, gaussi, fac, dfac
  
  INTEGER, PARAMETER :: maxfac = 30
  
  REAL ( dbl ), PARAMETER :: pi = 3.14159265358979323846264338_dbl
  REAL ( dbl ), PARAMETER :: pio2 = 1.57079632679489661923132169_dbl
  REAL ( dbl ), PARAMETER :: twopi = 6.28318530717958647692528677_dbl
  REAL ( dbl ), PARAMETER :: fourpi = 12.56637061435917295385057353_dbl
  REAL ( dbl ), PARAMETER :: root2 = 1.41421356237309504880168872_dbl
  REAL ( dbl ), PARAMETER :: rootpi = 1.77245387556702677717522568_dbl
  REAL ( dbl ), PARAMETER :: oorootpi = 1._dbl/rootpi
  
  REAL ( dbl ), PARAMETER :: zero = 0._dbl, one = 1._dbl, half = 0.5_dbl
  REAL ( dbl ), PARAMETER :: degree = 180._dbl/pi, radians = one/degree
  COMPLEX ( dbl ) :: gaussi
  
  REAL ( dbl ), DIMENSION (0:maxfac) :: fac
  REAL ( dbl ), DIMENSION (-1:2*maxfac+1) :: dfac

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** cp2k/init_mathcon [1.0] *
!!
!!   NAME
!!     init_mathcon
!!
!!   FUNCTION
!!     Definition of mathematical constants and functions.
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
  
SUBROUTINE init_mathcon
  
  IMPLICIT NONE
  
! Locals
  INTEGER :: i
  
!------------------------------------------------------------------------------
  
  gaussi = CMPLX ( zero, one )
  
! 
!    *** Factorial function ***
! 
  fac(0) = one
  DO i = 1, maxfac
     fac(i) = real(i,dbl)*fac(i-1)
  END DO
! 
!    *** Double factorial function ***
! 
  dfac ( -1 ) = one
  dfac (  0 ) = one
  dfac (  1 ) = one
  DO i = 3, 2 * maxfac + 1
     dfac(i) = REAL ( i, dbl ) * dfac ( i - 2 )
  END DO
  
END SUBROUTINE init_mathcon

!!*****
!******************************************************************************

END MODULE mathconstants

!******************************************************************************
