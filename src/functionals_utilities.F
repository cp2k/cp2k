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
!!     JGH (20.02.2001) : Added setup routine
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

  PUBLIC :: set_util, setup_calculation, calc_rs, calc_rho13, calc_fx, &
            calc_wave_vector, calc_srs, calc_rs_d, calc_z

  REAL(dbl) :: eps_rho

! *****************************************************************************

CONTAINS

! *****************************************************************************

    SUBROUTINE set_util(cutoff)

    REAL(dbl) :: cutoff

    eps_rho = cutoff

    END SUBROUTINE set_util

! *****************************************************************************

    SUBROUTINE setup_calculation(order,m,calc,tag)

    INTEGER, INTENT(IN) :: order
    INTEGER, INTENT(OUT) :: m(0:,:)
    LOGICAL, INTENT(OUT) :: calc(0:)
    INTEGER, INTENT(IN) :: tag

    IF ( ABS(order) > 3 ) &
       CALL stop_program("setup_calculation","Order of derivative too high")

    SELECT CASE (tag)
    CASE ( 100 )
          m(0,1) = 1
          m(0,2) = 1
          m(1,1) = 2
          m(1,2) = 2
          m(2,1) = 3
          m(2,2) = 3
          m(3,1) = 4
          m(3,2) = 4
    CASE ( 110 )
          m(0,1) = 1
          m(0,2) = 1
          m(1,1) = 2
          m(1,2) = 3
          m(2,1) = 4
          m(2,2) = 6
          m(3,1) = 7
          m(3,2) = 10
    CASE ( 200 )
          m(0,1) = 1
          m(0,2) = 1
          m(1,1) = 2
          m(1,2) = 3
          m(2,1) = 4
          m(2,2) = 5
          m(3,1) = 6
          m(3,2) = 7
    CASE ( 201 )
          m(0,1) = 1
          m(0,2) = 1
          m(1,1) = 2
          m(1,2) = 3
          m(2,1) = 4
          m(2,2) = 6
          m(3,1) = 7
          m(3,2) = 10
    CASE ( 210 )
          m(0,1) = 1
          m(0,2) = 1
          m(1,1) = 2
          m(1,2) = 5
          m(2,1) = 6
          m(2,2) = 11
          m(3,1) = 12
          m(3,2) = 19
    CASE ( 211 )
          m(0,1) = 1
          m(0,2) = 1
          m(1,1) = 2
          m(1,2) = 6
          m(2,1) = 7
          m(2,2) = 21
          m(3,1) = 22
          m(3,2) = 57
    CASE DEFAULT
       CALL stop_program("setup_calculation","Invalid tag")
    END SELECT

    calc = .FALSE.
    IF ( order >= 0 ) THEN
      calc(0:order) = .TRUE.
    ELSE
      calc(-order) = .TRUE.
      m(0:3,2) = m(0:3,2) - m(0:3,1) + 1
      m(0:3,1) = 1
    END IF

    END SUBROUTINE setup_calculation

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

    SUBROUTINE calc_srs ( rho, x )

!   ---------------------------------------------------------------------------
!   rs parameter : f*rho**(-1/3)
!   x = sqrt(rs)
!   ---------------------------------------------------------------------------

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho
    REAL(dbl), DIMENSION(:), INTENT(OUT) :: x

!   *** Local variables ***

    INTEGER :: n, ip

!   ---------------------------------------------------------------------------

    CALL calc_rs ( rho, x )
    n = SIZE ( rho )

    DO ip = 1, n
      x(ip) = sqrt(x(ip))
    END DO

  END SUBROUTINE calc_srs

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

  SUBROUTINE calc_wave_vector ( tag, rho, grho, s )

!   ---------------------------------------------------------------------------
!   wave vector s = |nabla rho| / (2(3pi^2)^1/3 * rho^4/3)
!   ---------------------------------------------------------------------------

    CHARACTER (len=*), INTENT(IN)        :: tag
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho
    REAL(dbl), DIMENSION(:), INTENT(IN)  :: grho
    REAL(dbl), DIMENSION(:), INTENT(OUT) :: s

!   *** Local variables ***

    INTEGER :: n, ip
    REAL(dbl) :: fac

!   ---------------------------------------------------------------------------

    fac = 1._dbl / (2._dbl*(3._dbl*pi*pi)**f13)
    IF ( tag(1:1)=="u" .OR. tag(1:1)=="U" ) fac = fac*(2._dbl)**f13

    n = SIZE ( rho )
    IF ( n > SIZE(s) ) &
      CALL stop_program ( "calc_wave_vector", "Incompatible array sizes" )
    IF ( n > SIZE(grho) ) &
      CALL stop_program ( "calc_wave_vector", "Incompatible array sizes" )

    DO ip = 1, n
      IF ( rho(ip) < eps_rho ) THEN
         s(ip) = 0.0_dbl
      ELSE
         s(ip) = fac*grho(ip)*rho(ip)**(-f43)
      END IF
    END DO

  END SUBROUTINE calc_wave_vector

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
    IF ( m > 3 ) CALL stop_program ( "calc_fx", "Order too high" )

    DO ip = 1, n
      rhoab = rhoa(ip) + rhob(ip)
      IF ( rhoab < eps_rho ) THEN
         fx(ip,1:m) = 0.0_dbl
      ELSE
         x = (rhoa(ip) - rhob(ip)) / rhoab
         IF ( x < -1._dbl ) THEN
           IF ( m >= 0 ) fx(ip,1) = 1._dbl
           IF ( m >= 1 ) fx(ip,2) = -f43*fac*2._dbl**f13
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

  SUBROUTINE calc_rs_d(rho, rs, order)

    IMPLICIT NONE

    REAL(DBL), DIMENSION(1:), INTENT(IN) :: rho
    REAL(DBL), DIMENSION(1:,0:), INTENT(OUT) :: rs
    INTEGER, INTENT(IN) :: order

    CALL calc_rs(rho, rs(:,0))
    IF (order >= 1) rs(:,1) = (-1._dbl/3._dbl) * rs(:,0) / rho(:)
    IF (order >= 2) rs(:,2) = (-4._dbl/3._dbl) * rs(:,1) / rho(:)
    IF (order >= 3) rs(:,3) = (-7._dbl/3._dbl) * rs(:,2) / rho(:) 

  END SUBROUTINE calc_rs_d

! *****************************************************************************

  SUBROUTINE calc_z ( a, b, z, order )

    IMPLICIT NONE

    REAL(DBL), INTENT(IN) :: a, b
    REAL(DBL), DIMENSION(0:,0:), INTENT(OUT) :: z
    INTEGER, INTENT(IN) :: order

    REAL(DBL) :: c,d

    c = a+b

    z(0,0) = (a-b)/c
    IF (order >= 1) THEN
       d = c**2
       z(1,0) = 2._dbl*b/d
       z(0,1) = -2._dbl*a/d
    END IF
    IF (order >= 2) THEN
       d = c**3
       z(2,0) = -4._dbl*b/d
       z(1,1) = 2._dbl*(a-b)/d
       z(0,2) = 4._dbl*a/d
    END IF
    IF (order >= 3) THEN
       d = c**4
       z(3,0) = 12._dbl*b/d
       z(2,1) = -4._dbl*(a-2._dbl*b)/d
       z(1,2) = -4._dbl*(2._dbl*a-b)/d
       z(0,3) = -12._dbl*a/d
    END IF

  END SUBROUTINE calc_z
    
! *****************************************************************************

END MODULE functionals_utilities

! *****************************************************************************
