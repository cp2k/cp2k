!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/xc_functionals_utilities [1.0] *
!!
!!   NAME
!!     xc_functionals_utilities
!!
!!   FUNCTION
!!     Utility routines for the functional calculations
!!
!!   AUTHOR
!!     JGH (15.02.2002)
!!
!!   MODIFICATION HISTORY
!!     JGH (20.02.2001) : Added setup routine
!!     JGH (26.02.2003) : OpenMP enabled
!!
!!   SOURCE
!******************************************************************************

MODULE xc_functionals_utilities

! *****************************************************************************

  USE kinds, ONLY: dbl
  USE termination, ONLY: stop_program

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  REAL(dbl), PARAMETER :: pi = 3.14159265358979323846264338_dbl
  REAL(dbl), PARAMETER :: rsfac = 0.6203504908994000166680065_dbl ! (4*pi/3)^(-1/3)
  REAL(dbl), PARAMETER :: f13 = 1.0_dbl/3.0_dbl, &
                          f23 = 2.0_dbl*f13, &
                          f43 = 4.0_dbl*f13, &
                          f53 = 5.0_dbl*f13
  REAL(dbl), PARAMETER :: fxfac = 1.923661050931536319759455_dbl ! 1/(2^(4/3) - 2)

  CHARACTER(len=*), PARAMETER :: moduleN="xc_functionals_utilities"

  PUBLIC :: set_util, setup_calculation, calc_rs, calc_rho13, calc_fx, &
            calc_wave_vector, calc_srs, calc_z, calc_rs_pw

  INTERFACE calc_fx
     MODULE PROCEDURE calc_fx_array, calc_fx_single
  END INTERFACE

  INTERFACE calc_rs
     MODULE PROCEDURE calc_rs_array, calc_rs_single
  END INTERFACE

  REAL(dbl) :: eps_rho

!!***
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
    CASE ( 100 )          ! LDA
          m(0,1) = 1
          m(0,2) = 1
          m(1,1) = 2
          m(1,2) = 2
          m(2,1) = 3
          m(2,2) = 3
          m(3,1) = 4
          m(3,2) = 4
    CASE ( 110 )          ! GGA non-spin polarized
          m(0,1) = 1
          m(0,2) = 1
          m(1,1) = 2
          m(1,2) = 3
          m(2,1) = 4
          m(2,2) = 6
          m(3,1) = 7
          m(3,2) = 10
    CASE ( 200 )          ! LSD no crossterms
          m(0,1) = 1
          m(0,2) = 1
          m(1,1) = 2
          m(1,2) = 3
          m(2,1) = 4
          m(2,2) = 5
          m(3,1) = 6
          m(3,2) = 7
    CASE ( 201 )          ! LSD with crossterms
          m(0,1) = 1
          m(0,2) = 1
          m(1,1) = 2
          m(1,2) = 3
          m(2,1) = 4
          m(2,2) = 6
          m(3,1) = 7
          m(3,2) = 10
    CASE ( 210 )          ! GGA spin-polarized, no crossterms
          m(0,1) = 1
          m(0,2) = 1
          m(1,1) = 2
          m(1,2) = 5
          m(2,1) = 6
          m(2,2) = 11
          m(3,1) = 12
          m(3,2) = 19
    CASE ( 211 )          ! GGA spin-polarized, with crossterms
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

    SUBROUTINE calc_rs_single ( rho, rs )

!   ---------------------------------------------------------------------------
!   rs parameter : f*rho**(-1/3)
!   ---------------------------------------------------------------------------

    REAL(dbl), INTENT(IN)  :: rho
    REAL(dbl), INTENT(OUT) :: rs

!   ---------------------------------------------------------------------------

    IF ( rho < eps_rho ) THEN
       rs = 0.0_dbl
    ELSE
       rs = rsfac * rho**(-f13)
    END IF

  END SUBROUTINE calc_rs_single

! *****************************************************************************

  SUBROUTINE calc_rs_array ( rho, rs )
    
!   ---------------------------------------------------------------------------
!   rs parameter : f*rho**(-1/3)
!   ---------------------------------------------------------------------------

    REAL(dbl), DIMENSION(:), INTENT(IN)  :: rho
    REAL(dbl), INTENT(OUT), DIMENSION(:) :: rs

    INTEGER :: k

!   ---------------------------------------------------------------------------

    IF (SIZE(rs) < SIZE(rho)) THEN
       CALL stop_program("functionals_utilities/calc_rs_array","size of array rs too small")
    END IF

!$omp parallel do private(k)
    DO k=1, SIZE(rs)
       IF ( rho(k) < eps_rho ) THEN
          rs(k) = 0.0_dbl
       ELSE
          rs(k) = rsfac * rho(k)**(-f13)
       END IF
    END DO

  END SUBROUTINE calc_rs_array

! *****************************************************************************

  SUBROUTINE calc_rs_pw ( rho, rs, n )
    
!   ---------------------------------------------------------------------------
!   rs parameter : f*rho**(-1/3)
!   ---------------------------------------------------------------------------

    REAL(KIND=dbl), DIMENSION(*), INTENT(IN)  :: rho
    REAL(KIND=dbl), DIMENSION(*), INTENT(OUT) :: rs
    INTEGER, INTENT(IN) :: n

    INTEGER :: k

!   ---------------------------------------------------------------------------

!$omp parallel do private(k)
    DO k=1, n
       IF ( rho(k) < eps_rho ) THEN
          rs(k) = 0.0_dbl
       ELSE
          rs(k) = rsfac * rho(k)**(-f13)
       END IF
    END DO

  END SUBROUTINE calc_rs_pw

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

!$omp parallel do private(ip)
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

!$omp parallel do private(ip)
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

!   TAGS: U: total density, spin wave vector
!         R: spin density, total density wave vector
!   ---------------------------------------------------------------------------

    fac = 1._dbl / (2._dbl*(3._dbl*pi*pi)**f13)
    IF ( tag(1:1)=="u" .OR. tag(1:1)=="U" ) fac = fac*(2._dbl)**f13
    IF ( tag(1:1)=="r" .OR. tag(1:1)=="R" ) fac = fac*(2._dbl)**f13

    n = SIZE ( rho )
    IF ( n > SIZE(s) ) &
      CALL stop_program ( "calc_wave_vector", "Incompatible array sizes" )
    IF ( n > SIZE(grho) ) &
      CALL stop_program ( "calc_wave_vector", "Incompatible array sizes" )

!$omp parallel do private(ip)
    DO ip = 1, n
      IF ( rho(ip) < eps_rho ) THEN
         s(ip) = 0.0_dbl
      ELSE
         s(ip) = fac*grho(ip)*rho(ip)**(-f43)
      END IF
    END DO

  END SUBROUTINE calc_wave_vector

! *****************************************************************************

  SUBROUTINE calc_fx_array(n, rhoa, rhob, fx, m )

!   ---------------------------------------------------------------------------
!   spin interpolation function and derivatives
!   
!   f(x) = ( (1+x)^(4/3) + (1-x)^(4/3) - 2 ) / (2^(4/3)-2)
!   df(x) = (4/3)( (1+x)^(1/3) - (1-x)^(1/3) ) / (2^(4/3)-2)
!   d2f(x) = (4/9)( (1+x)^(-2/3) + (1-x)^(-2/3) ) / (2^(4/3)-2)
!   d3f(x) = (-8/27)( (1+x)^(-5/3) - (1-x)^(-5/3) ) / (2^(4/3)-2)
!
!   ---------------------------------------------------------------------------

    INTEGER, INTENT(IN)                 :: n ! number of points
    REAL(dbl), DIMENSION(*), INTENT(IN) :: rhoa, rhob
    REAL(dbl), DIMENSION(:,:), INTENT(OUT) :: fx
    INTEGER, INTENT(IN)                 :: m ! order of derivatives

!   *** Local variables ***

    INTEGER :: ierr, ip
    REAL(dbl) :: x, rhoab

!   *** Parameters ***
    CHARACTER(len=*), PARAMETER :: routineN="calc_fx_array", &
                                   routineP=moduleN//"/"//routineN

!   ---------------------------------------------------------------------------

    IF ( m > 3 ) CALL stop_program (routineP, "Order too high" )
!!    IF (.NOT.ASSOCIATED(fx)) THEN
!!       ALLOCATE(fx(n,m+1), STAT=ierr)
!!       IF (ierr /= 0) CALL stop_memory(routineP, "fx", n*m)
!!    ELSE
       IF (SIZE(fx,1)<n) CALL stop_program(routineP, "SIZE(fx,1) too small")
       IF (SIZE(fx,2)<m) CALL stop_program(routineP, "SIZE(fx,2) too small")
!!    END IF
       
!$omp parallel do private(ip,x,rhoab)
    DO ip = 1, n
      rhoab = rhoa(ip) + rhob(ip)
      IF ( rhoab < eps_rho ) THEN
         fx(ip,1:m) = 0.0_dbl
      ELSE
         x = (rhoa(ip) - rhob(ip)) / rhoab
         IF ( x < -1._dbl ) THEN
           IF ( m >= 0 ) fx(ip,1) = 1._dbl
           IF ( m >= 1 ) fx(ip,2) = -f43*fxfac*2._dbl**f13
           IF ( m >= 2 ) fx(ip,3) = f13*f43*fxfac/2._dbl**f23
           IF ( m >= 3 ) fx(ip,4) = f23*f13*f43*fxfac/2._dbl**f53
         ELSE IF ( x > 1._dbl ) THEN
           IF ( m >= 0 ) fx(ip,1) = 1._dbl
           IF ( m >= 1 ) fx(ip,2) = f43*fxfac*2._dbl**f13
           IF ( m >= 2 ) fx(ip,3) = f13*f43*fxfac/2._dbl**f23
           IF ( m >= 3 ) fx(ip,4) = -f23*f13*f43*fxfac/2._dbl**f53
         ELSE
           IF ( m >= 0 ) &
              fx(ip,1) = ( (1._dbl+x)**f43 + (1._dbl-x)**f43 - 2._dbl ) * fxfac
           IF ( m >= 1 ) &
              fx(ip,2) = ( (1._dbl+x)**f13 - (1._dbl-x)**f13 ) * fxfac * f43
           IF ( m >= 2 ) &
              fx(ip,3) = ( (1._dbl+x)**(-f23) + (1._dbl-x)**(-f23) ) * &
                         fxfac * f43 * f13
           IF ( m >= 3 ) &
              fx(ip,4) = ( (1._dbl+x)**(-f53) - (1._dbl-x)**(-f53) ) * &
                         fxfac * f43 * f13 * (-f23)
         END IF
      END IF
    END DO

  END SUBROUTINE calc_fx_array

! *****************************************************************************

  SUBROUTINE calc_fx_single ( rhoa, rhob, fx, m )

!   ---------------------------------------------------------------------------
!   spin interpolation function and derivatives
!   
!   f(x) = ( (1+x)^(4/3) + (1-x)^(4/3) - 2 ) / (2^(4/3)-2)
!   df(x) = (4/3)( (1+x)^(1/3) - (1-x)^(1/3) ) / (2^(4/3)-2)
!   d2f(x) = (4/9)( (1+x)^(-2/3) + (1-x)^(-2/3) ) / (2^(4/3)-2)
!   d3f(x) = (-8/27)( (1+x)^(-5/3) - (1-x)^(-5/3) ) / (2^(4/3)-2)
!
!   ---------------------------------------------------------------------------

    REAL(dbl), INTENT(IN)  :: rhoa, rhob
    REAL(dbl), DIMENSION(:), INTENT(OUT) :: fx
    INTEGER, INTENT(IN) :: m

!   *** Local variables ***

    INTEGER :: n, na, nb, ip
    REAL(dbl) :: x, rhoab

!   ---------------------------------------------------------------------------

    rhoab = rhoa + rhob
    IF ( rhoab < eps_rho ) THEN
       fx(1:m) = 0.0_dbl
    ELSE
       x = (rhoa - rhob) / rhoab
       IF ( x < -1._dbl ) THEN
          IF ( m >= 0 ) fx(1) = 1._dbl
          IF ( m >= 1 ) fx(2) = -f43*fxfac*2._dbl**f13
          IF ( m >= 2 ) fx(3) = f13*f43*fxfac/2._dbl**f23
          IF ( m >= 3 ) fx(4) = f23*f13*f43*fxfac/2._dbl**f53
       ELSE IF ( x > 1._dbl ) THEN
          IF ( m >= 0 ) fx(1) = 1._dbl
          IF ( m >= 1 ) fx(2) = f43*fxfac*2._dbl**f13
          IF ( m >= 2 ) fx(3) = f13*f43*fxfac/2._dbl**f23
          IF ( m >= 3 ) fx(4) = -f23*f13*f43*fxfac/2._dbl**f53
       ELSE
          IF ( m >= 0 ) &
               fx(1) = ( (1._dbl+x)**f43 + (1._dbl-x)**f43 - 2._dbl ) * fxfac
          IF ( m >= 1 ) &
               fx(2) = ( (1._dbl+x)**f13 - (1._dbl-x)**f13 ) * fxfac * f43
          IF ( m >= 2 ) &
               fx(3) = ( (1._dbl+x)**(-f23) + (1._dbl-x)**(-f23) ) * &
               fxfac * f43 * f13
          IF ( m >= 3 ) &
               fx(4) = ( (1._dbl+x)**(-f53) - (1._dbl-x)**(-f53) ) * &
               fxfac * f43 * f13 * (-f23)
       END IF
    END IF
    
  END SUBROUTINE calc_fx_single

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
       d = c*c
       z(1,0) = 2._dbl*b/d
       z(0,1) = -2._dbl*a/d
    END IF
    IF (order >= 2) THEN
       d = d*c
       z(2,0) = -4._dbl*b/d
       z(1,1) = 2._dbl*(a-b)/d
       z(0,2) = 4._dbl*a/d
    END IF
    IF (order >= 3) THEN
       d = d*c
       z(3,0) = 12._dbl*b/d
       z(2,1) = -4._dbl*(a-2._dbl*b)/d
       z(1,2) = -4._dbl*(2._dbl*a-b)/d
       z(0,3) = -12._dbl*a/d
    END IF

  END SUBROUTINE calc_z
    
! *****************************************************************************

END MODULE xc_functionals_utilities

! *****************************************************************************
