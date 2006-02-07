!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2004 CP2K developers group                           !
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
!!     JGH (16-06-2002) : Added Gamma functions
!!     JGH (10-08-2004) : Added Euler constant (gamma)
!!
!!   SOURCE
!******************************************************************************

MODULE mathconstants

  USE kinds,                           ONLY: dp

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: init_mathcon
  PUBLIC :: maxfac, pi, pio2, twopi, fourpi, root2, rootpi, oorootpi, &
       zero, one, half, degree, radians, gaussi, fac, dfac, gamma0, gamma1, &
       euler

#if defined(__SGL)
  INTEGER, PARAMETER :: maxfac = 24
#else
  INTEGER, PARAMETER :: maxfac = 30
#endif

  REAL (KIND=dp), PARAMETER :: pi = 3.14159265358979323846264338_dp
  REAL (KIND=dp), PARAMETER :: pio2 = 1.57079632679489661923132169_dp
  REAL (KIND=dp), PARAMETER :: twopi = 6.28318530717958647692528677_dp
  REAL (KIND=dp), PARAMETER :: fourpi = 12.56637061435917295385057353_dp
  REAL (KIND=dp), PARAMETER :: root2 = 1.41421356237309504880168872_dp
  REAL (KIND=dp), PARAMETER :: rootpi = 1.77245387556702677717522568_dp
  REAL (KIND=dp), PARAMETER :: oorootpi = 1.0_dp/rootpi

  REAL (KIND=dp), PARAMETER :: euler=0.577215664901532860606512_dp

  REAL (KIND=dp), PARAMETER :: zero = 0.0_dp, one = 1.0_dp, half = 0.5_dp
  REAL (KIND=dp), PARAMETER :: degree = 180.0_dp/pi, radians = one/degree
  COMPLEX (KIND=dp) :: gaussi

  REAL (KIND=dp), DIMENSION (0:maxfac) :: fac
  REAL (KIND=dp), DIMENSION (-1:2*maxfac+1) :: dfac
  REAL (KIND=dp), DIMENSION (0:maxfac) :: gamma0, gamma1

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
!!*****
!******************************************************************************

SUBROUTINE init_mathcon

    INTEGER                                  :: i

!------------------------------------------------------------------------------

  gaussi = CMPLX ( zero, one,KIND=dp)

!
!    *** Factorial function ***
!
  fac(0) = one
  DO i = 1, maxfac
     fac(i) = REAL(i,dp)*fac(i-1)
  END DO
!
!    *** Double factorial function ***
!
  dfac ( -1 ) = one
  dfac (  0 ) = one
  dfac (  1 ) = one
  DO i = 2, 2 * maxfac + 1
     dfac(i) = REAL ( i,KIND=dp) * dfac ( i - 2 )
  END DO
!
! Gamma functions 
! gamma ( n ) = gamma0 ( n ) = (n-1)!
! gamma ( n + 1/2 ) = gamma1 ( n ) = SQRT(pi)/2^n (2n-1)!!
!
  gamma0 ( 0 ) = 0
  gamma1 ( 0 ) = rootpi
  DO i = 1, maxfac
     gamma0 ( i ) = fac ( i-1 )
     gamma1 ( i ) = rootpi / 2**i * dfac ( 2*i-1 )
  END DO

END SUBROUTINE init_mathcon

!******************************************************************************

END MODULE mathconstants

!******************************************************************************
