!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/functionals_utilities [1.0] *
!!
!!   NAME
!!     functionals_utilities
!!
!!   FUNCTION
!!     Utility routines for the functional calculations
!!
!!   AUTHOR
!!     JGH (15.02.2002)
!!
!!   MODIFICATION HISTORY
!!     - none
!!
!!   SOURCE
!******************************************************************************

MODULE functionals_utilities

! *****************************************************************************

  USE kinds, ONLY: dbl
  USE termination, ONLY: stop_program

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  REAL(dbl), PARAMETER :: pi = 3.14159265358979323846264338_dbl
  REAL(dbl), PARAMETER :: f13 = 1.0_dbl/3.0_dbl, &
                          f23 = 2.0_dbl*f13, &
                          f43 = 4.0_dbl*f13, &
                          f53 = 5.0_dbl*f13

  PUBLIC :: set_util, calc_rs, calc_rho13, calc_fx

  REAL(dbl) :: eps_rho

! *****************************************************************************

CONTAINS

! *****************************************************************************

    SUBROUTINE set_util(cutoff)

    REAL(dbl) :: cutoff

    eps_rho = cutoff

    END SUBROUTINE set_util

! *****************************************************************************

    SUBROUTINE calc_rs ( rho, rs )

!   ---------------------------------------------------------------------------
!   rs parameter : f*rho**(-1/3)
!   ---------------------------------------------------------------------------

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho
    REAL(dbl), DIMENSION(:), INTENT(OUT) :: rs

!   *** Local variables ***

    INTEGER :: n, ip
    REAL(dbl) :: rsfac

!   ---------------------------------------------------------------------------

    rsfac = (f43*pi)**(-f13)
    n = SIZE ( rho )
    IF ( n > SIZE(rs) ) &
      CALL stop_program ( "calc_rs", "Incompatible array sizes" )

    DO ip = 1, n
      IF ( rho(ip) < eps_rho ) THEN
         rs(ip) = 0.0_dbl
      ELSE
         rs(ip) = rsfac * rho(ip)**(-f13)
      END IF
    END DO

  END SUBROUTINE calc_rs

! *****************************************************************************

    SUBROUTINE calc_rho13 ( rho, r13 )

!   ---------------------------------------------------------------------------

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho
    REAL(dbl), DIMENSION(:), INTENT(OUT) :: r13

!   *** Local variables ***

    INTEGER :: n, ip

!   ---------------------------------------------------------------------------

    n = SIZE ( rho )
    IF ( n > SIZE(r13) ) &
      CALL stop_program ( "calc_rho13", "Incompatible array sizes" )

    DO ip = 1, n
      IF ( rho(ip) < eps_rho ) THEN
         r13(ip) = 0.0_dbl
      ELSE
         r13(ip) = rho(ip)**f13
      END IF
    END DO

  END SUBROUTINE calc_rho13

! *****************************************************************************

  SUBROUTINE calc_fx ( rhoa, rhob, fx, m )

!   ---------------------------------------------------------------------------
!   spin interpolation function and derivatives
!   
!   f(x) = ( (1+x)^(4/3) + (1-x)^(4/3) - 2 ) / (2^(4/3)-2)
!   df(x) = (4/3)( (1+x)^(1/3) - (1-x)^(1/3) ) / (2^(4/3)-2)
!   d2f(x) = (4/9)( (1+x)^(-2/3) + (1-x)^(-2/3) ) / (2^(4/3)-2)
!   d3f(x) = (-8/27)( (1+x)^(-5/3) - (1-x)^(-5/3) ) / (2^(4/3)-2)
!
!   ---------------------------------------------------------------------------

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rhoa, rhob
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: fx
    INTEGER, INTENT(IN) :: m

!   *** Local variables ***

    INTEGER :: n, na, nb, ip
    REAL(dbl) :: fac, x, rhoab

!   ---------------------------------------------------------------------------

    fac = 1._dbl/(2._dbl**f43-2._dbl)
    na = SIZE ( rhoa )
    nb = SIZE ( rhob )
    n = MAX ( na, nb)
    IF ( n > SIZE(fx,1) .OR. na /= nb ) &
      CALL stop_program ( "calc_fx", "Incompatible array sizes" )
    IF ( m > SIZE(fx,2) ) &
      CALL stop_program ( "calc_fx", "Incompatible array sizes" )

    DO ip = 1, n
      rhoab = rhoa(ip) + rhob(ip)
      IF ( rhoab < eps_rho ) THEN
         fx(ip,1:m) = 0.0_dbl
      ELSE
         x = (rhoa(ip) - rhob(ip)) / rhoab
         IF ( x < -1._dbl ) THEN
           IF ( m >= 0 ) fx(ip,1) = 1._dbl
           IF ( m >= 1 ) fx(ip,1) = -f43*fac*2._dbl**f13
           IF ( m >= 2 ) fx(ip,3) = f13*f43*fac/2._dbl**f23
           IF ( m >= 3 ) fx(ip,4) = f23*f13*f43*fac/2._dbl**f53
         ELSE IF ( x > 1._dbl ) THEN
           IF ( m >= 0 ) fx(ip,1) = 1._dbl
           IF ( m >= 1 ) fx(ip,2) = f43*fac*2._dbl**f13
           IF ( m >= 2 ) fx(ip,3) = f13*f43*fac/2._dbl**f23
           IF ( m >= 3 ) fx(ip,4) = -f23*f13*f43*fac/2._dbl**f53
         ELSE
           IF ( m >= 0 ) &
              fx(ip,1) = ( (1._dbl+x)**f43 + (1._dbl-x)**f43 - 2._dbl ) * fac
           IF ( m >= 1 ) &
              fx(ip,2) = ( (1._dbl+x)**f13 - (1._dbl-x)**f13 ) * fac * f43
           IF ( m >= 2 ) &
              fx(ip,3) = ( (1._dbl+x)**(-f23) + (1._dbl-x)**(-f23) ) * &
                         fac * f43 * f13
           IF ( m >= 3 ) &
              fx(ip,4) = ( (1._dbl+x)**(-f53) - (1._dbl-x)**(-f53) ) * &
                         fac * f43 * f13 * (-f23)
         END IF
      END IF
    END DO

  END SUBROUTINE calc_fx

! *****************************************************************************

END MODULE functionals_utilities

! *****************************************************************************
