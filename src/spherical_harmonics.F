!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/spherical_harmonics [1.0] *
!!
!!   NAME
!!     spherical_harmonics
!!
!!   FUNCTION
!!     Calculate spherical harmonics
!!
!!   AUTHOR
!!     JGH 6-Oct-2000
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!     Numerical Stability up to L=15
!!     Accuracy > 1.E-12 up to L=15 tested
!!
!!   SOURCE
!******************************************************************************

MODULE spherical_harmonics

! *****************************************************************************

  USE kinds,         ONLY: dbl

  USE mathconstants, ONLY: pi, fac
  USE termination,   ONLY: stop_memory, stop_program

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: y_lm

  INTERFACE y_lm
     MODULE PROCEDURE rvy_lm, rry_lm, cvy_lm, ccy_lm
  END INTERFACE


!!*****
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE rvy_lm ( r, y, l, m )
!
! Real Spherical Harmonics
!                   _                   _
!                  |  [(2l+1)(l-|m|)!]   |1/2 m         cos(m p)   m>=0
!  Y_lm ( t, p ) = |---------------------|   P_l(cos(t))
!                  |[2Pi(1+d_m0)(l+|m|)!]|              sin(|m| p) m<0
!
! Input: r == (x,y,z) : normalised    x^2 + y^2 + z^2 = 1
!
  REAL(dbl), DIMENSION (:,:), INTENT(IN) :: r
  REAL(dbl), DIMENSION (:), INTENT(OUT) :: y
  INTEGER, INTENT(IN) :: l, m

  INTEGER :: i
  REAL(dbl) :: pf, lpm, lmm, plm, z, rxy, t, cp, sp

  SELECT CASE ( l )
  CASE ( : -1 )
    call stop_program("y_lm","Negative l value")
  CASE ( 0 )
    pf = SQRT ( 1._dbl / ( 4._dbl * pi ) )
    IF ( m /= 0 )  call stop_program("y_lm","l = 0 and m value out of bounds")
    y(:) = pf
  CASE ( 1 )
    SELECT CASE ( m )
    CASE DEFAULT
      call stop_program("y_lm","l = 1 and m value out of bounds")
    CASE ( 1 )
      pf = - SQRT ( 3._dbl / ( 4._dbl * pi ) )
      y ( : ) = pf * r ( 1, : )
    CASE ( 0 )
      pf = SQRT ( 3._dbl / ( 4._dbl * pi ) )
      y ( : ) = pf * r ( 3, : )
    CASE ( -1 )
      pf = - SQRT ( 3._dbl / ( 4._dbl * pi ) )
      y ( : ) = pf * r ( 2, : )
    END SELECT
  CASE ( 2 )
    SELECT CASE ( m )
    CASE DEFAULT
      call stop_program("y_lm","l = 2 and m value out of bounds")
    CASE ( 2 )
      pf = SQRT ( 15._dbl / ( 16._dbl * pi ) )
      y ( : ) = pf * ( r ( 1, : ) * r ( 1, : ) - r ( 2, : ) * r ( 2, : ) )
    CASE ( 1 )
      pf = - SQRT ( 15._dbl / ( 4._dbl * pi ) )
      y ( : ) = pf * r ( 3, : ) * r ( 1, : )
    CASE ( 0 )
      pf = SQRT ( 5._dbl / ( 16._dbl * pi ) )
      y ( : ) = pf * ( 3._dbl * r ( 3, : ) * r ( 3, : ) - 1._dbl )
    CASE ( -1 )
      pf = - SQRT ( 15._dbl / ( 4._dbl * pi ) )
      y ( : ) = pf * r ( 3, : ) * r ( 2, : )
    CASE ( -2 )
      pf = SQRT ( 15._dbl / ( 16._dbl * pi ) )
      y ( : ) = pf * 2._dbl * r ( 1, : ) * r ( 2, : )
    END SELECT
  CASE ( 3 )
    SELECT CASE ( m )
    CASE DEFAULT
      call stop_program("y_lm","l = 3 and m value out of bounds")
    CASE ( 3 )
      pf = - SQRT ( 35._dbl / ( 32._dbl * pi ) )
      y ( : ) = pf * r ( 1, : ) * ( r ( 1, : )**2 - 3._dbl * r ( 2, : )**2 )
    CASE ( 2 )
      pf = SQRT ( 105._dbl / ( 16._dbl * pi ) )
      y ( : ) = pf * r ( 3, : ) * ( r ( 1, : )**2 - r ( 2, : )**2 )
    CASE ( 1 )
      pf = - SQRT ( 21._dbl / ( 32._dbl * pi ) )
      y ( : ) = pf * r ( 1, : ) * ( 5._dbl * r ( 3, : ) * r ( 3, : ) - 1._dbl )
    CASE ( 0 )
      pf = SQRT ( 7._dbl / ( 16._dbl * pi ) )
      y ( : ) = pf * r ( 3, : ) * ( 5._dbl * r ( 3, : ) * r ( 3, : ) - 3._dbl )
    CASE ( -1 )
      pf = - SQRT ( 21._dbl / ( 32._dbl * pi ) )
      y ( : ) = pf * r ( 2, : ) * ( 5._dbl * r ( 3, : ) * r ( 3, : ) - 1._dbl )
    CASE ( -2 )
      pf = SQRT ( 105._dbl / ( 16._dbl * pi ) )
      y ( : ) = pf * 2._dbl * r ( 1, : ) * r ( 2, : ) * r ( 3, : )
    CASE ( -3 )
      pf = - SQRT ( 35._dbl / ( 32._dbl * pi ) )
      y ( : ) = pf * r ( 2, : ) * ( 3._dbl * r ( 1, : )**2 - r ( 2, : )**2 )
    END SELECT
  CASE DEFAULT
    IF ( m < -l .OR. m > l ) call stop_program("y_lm","m value out of bounds")
    lpm = fac ( l + ABS ( m ) )
    lmm = fac ( l - ABS ( m ) )
    IF ( m == 0 ) THEN
      t = 4._dbl * pi
    ELSE
      t = 2._dbl * pi
    END IF
    IF ( ABS ( lpm ) < EPSILON ( 1._dbl ) ) THEN
      pf = REAL ( 2*l+1, dbl )  / t
    ELSE
      pf = ( REAL ( 2*l+1, dbl ) * lmm ) / ( t * lpm )
    ENDIF
    pf = sqrt ( pf )
    DO i = 1, SIZE ( r ( 1, : ) )
      z = r ( 3, i )
      plm = legendre ( z, l, m )
      IF ( m == 0 ) THEN
        y ( i ) = pf * plm
      ELSE
        rxy = SQRT ( r ( 1, i ) ** 2  +  r ( 2, i ) ** 2 )
        IF ( rxy < EPSILON ( 1._dbl ) ) THEN
          y ( i ) = 0._dbl
        ELSE
          cp = r ( 1, i ) / rxy
          sp = r ( 2, i ) / rxy
          IF ( m > 0 ) THEN
            y ( i ) = pf * plm * cosn ( m, cp, sp )
          ELSE
            y ( i ) = pf * plm * sinn ( -m , cp, sp )
          END IF
        END IF
      END IF
    END DO
  END SELECT

  END SUBROUTINE rvy_lm

! *****************************************************************************

  SUBROUTINE rry_lm ( r, y, l, m )
!
! Real Spherical Harmonics
!                   _                   _
!                  |  [(2l+1)(l-|m|)!]   |1/2 m         cos(m p)   m>=0
!  Y_lm ( t, p ) = |---------------------|   P_l(cos(t))
!                  |[2Pi(1+d_m0)(l+|m|)!]|              sin(|m| p) m<0
!
! Input: r == (x,y,z) : normalised    x^2 + y^2 + z^2 = 1
!
  REAL(dbl), DIMENSION (:), INTENT(IN) :: r
  REAL(dbl), INTENT(OUT) :: y
  INTEGER, INTENT(IN) :: l, m

  INTEGER :: i
  REAL(dbl) :: pf, lpm, lmm, plm, z, rxy, t, cp, sp

  SELECT CASE ( l )
  CASE ( : -1 )
    call stop_program("y_lm","Negative l value")
  CASE ( 0 )
    pf = SQRT ( 1._dbl / ( 4._dbl * pi ) )
    IF ( m /= 0 )  call stop_program("y_lm","l = 0 and m value out of bounds")
    y = pf
  CASE ( 1 )
    SELECT CASE ( m )
    CASE DEFAULT
      call stop_program("y_lm","l = 1 and m value out of bounds")
    CASE ( 1 )
      pf = - SQRT ( 3._dbl / ( 4._dbl * pi ) )
      y = pf * r ( 1 )
    CASE ( 0 )
      pf = SQRT ( 3._dbl / ( 4._dbl * pi ) )
      y = pf * r ( 3 )
    CASE ( -1 )
      pf = - SQRT ( 3._dbl / ( 4._dbl * pi ) )
      y = pf * r ( 2 )
    END SELECT
  CASE ( 2 )
    SELECT CASE ( m )
    CASE DEFAULT
      call stop_program("y_lm","l = 2 and m value out of bounds")
    CASE ( 2 )
      pf = SQRT ( 15._dbl / ( 16._dbl * pi ) )
      y = pf * ( r ( 1 ) * r ( 1 ) - r ( 2 ) * r ( 2 ) )
    CASE ( 1 )
      pf = - SQRT ( 15._dbl / ( 4._dbl * pi ) )
      y = pf * r ( 3 ) * r ( 1 )
    CASE ( 0 )
      pf = SQRT ( 5._dbl / ( 16._dbl * pi ) )
      y = pf * ( 3._dbl * r ( 3 ) * r ( 3 ) - 1._dbl )
    CASE ( -1 )
      pf = - SQRT ( 15._dbl / ( 4._dbl * pi ) )
      y = pf * r ( 3 ) * r ( 2 )
    CASE ( -2 )
      pf = SQRT ( 15._dbl / ( 16._dbl * pi ) )
      y = pf * 2._dbl * r ( 1 ) * r ( 2 )
    END SELECT
  CASE ( 3 )
    SELECT CASE ( m )
    CASE DEFAULT
      call stop_program("y_lm","l = 3 and m value out of bounds")
    CASE ( 3 )
      pf = - SQRT ( 35._dbl / ( 32._dbl * pi ) )
      y = pf * r ( 1 ) * ( r ( 1 )**2 - 3._dbl * r ( 2 )**2 )
    CASE ( 2 )
      pf = SQRT ( 105._dbl / ( 16._dbl * pi ) )
      y = pf * r ( 3 ) * ( r ( 1 )**2 - r ( 2 )**2 )
    CASE ( 1 )
      pf = - SQRT ( 21._dbl / ( 32._dbl * pi ) )
      y = pf * r ( 1 ) * ( 5._dbl * r ( 3 ) * r ( 3 ) - 1._dbl )
    CASE ( 0 )
      pf = SQRT ( 7._dbl / ( 16._dbl * pi ) )
      y = pf * r ( 3 ) * ( 5._dbl * r ( 3 ) * r ( 3 ) - 3._dbl )
    CASE ( -1 )
      pf = - SQRT ( 21._dbl / ( 32._dbl * pi ) )
      y = pf * r ( 2 ) * ( 5._dbl * r ( 3 ) * r ( 3 ) - 1._dbl )
    CASE ( -2 )
      pf = SQRT ( 105._dbl / ( 16._dbl * pi ) )
      y = pf * 2._dbl * r ( 1 ) * r ( 2 ) * r ( 3 )
    CASE ( -3 )
      pf = - SQRT ( 35._dbl / ( 32._dbl * pi ) )
      y = pf * r ( 2 ) * ( 3._dbl * r ( 1 )**2 - r ( 2 )**2 )
    END SELECT
  CASE DEFAULT
    IF ( m < -l .OR. m > l ) call stop_program("y_lm","m value out of bounds")
    lpm = fac ( l + ABS ( m ) )
    lmm = fac ( l - ABS ( m ) )
    IF ( m == 0 ) THEN
      t = 4._dbl * pi
    ELSE
      t = 2._dbl * pi
    END IF
    IF ( ABS ( lpm ) < EPSILON ( 1._dbl ) ) THEN
      pf = REAL ( 2*l+1, dbl )  / t
    ELSE
      pf = ( REAL ( 2*l+1, dbl ) * lmm ) / ( t * lpm )
    ENDIF
    pf = sqrt ( pf )
    z = r ( 3 )
    plm = legendre ( z, l, m )
    IF ( m == 0 ) THEN
      y = pf * plm
    ELSE
      rxy = SQRT ( r ( 1 ) ** 2  +  r ( 2 ) ** 2 )
      IF ( rxy < EPSILON ( 1._dbl ) ) THEN
        y = 0._dbl
      ELSE
        cp = r ( 1 ) / rxy
        sp = r ( 2 ) / rxy
        IF ( m > 0 ) THEN
          y = pf * plm * cosn ( m, cp, sp )
        ELSE
          y = pf * plm * sinn ( -m , cp, sp )
        END IF
      END IF
    END IF
  END SELECT

  END SUBROUTINE rry_lm

! *****************************************************************************

  SUBROUTINE cvy_lm ( r, y, l, m )
!
! Real Spherical Harmonics
!                   _                _
!                  | [(2l+1)(l-|m|)!] |1/2 m
!  Y_lm ( t, p ) = |------------------|   P_l(cos(t)) Exp[ i m p ]
!                  |  [4Pi(l+|m|)!]|  |
!
! Input: r == (x,y,z) : normalised    x^2 + y^2 + z^2 = 1
!
  REAL(dbl), DIMENSION (:,:), INTENT(IN) :: r
  COMPLEX(dbl), DIMENSION (:), INTENT(OUT) :: y
  INTEGER, INTENT(IN) :: l, m

  INTEGER :: i, n
  REAL(dbl) :: pf, lpm, lmm, plm, z, rxy, t, cp, sp, yp, ym

  n = SIZE ( r ( 1, : ) )
  SELECT CASE ( l )
  CASE ( : -1 )
    call stop_program("y_lm","Negative l value")
  CASE ( 0 )
    pf = SQRT ( 1._dbl / ( 4._dbl * pi ) )
    IF ( m /= 0 )  call stop_program("y_lm","l = 0 and m value out of bounds")
    y(:) = pf
  CASE ( 1 )
    SELECT CASE ( m )
    CASE DEFAULT
      call stop_program("y_lm","l = 1 and m value out of bounds")
    CASE ( 1 )
      pf = - SQRT ( 3._dbl / ( 8._dbl * pi ) )
      DO i = 1, n
        yp = r ( 1, i )
        ym = r ( 2, i )
        y ( i ) = pf * CMPLX ( yp, ym, dbl )
      END DO
    CASE ( 0 )
      pf = SQRT ( 3._dbl / ( 4._dbl * pi ) )
      y ( : ) = pf * r ( 3, : )
    CASE ( -1 )
      pf = - SQRT ( 3._dbl / ( 8._dbl * pi ) )
      DO i = 1, n
        yp = r ( 1, i )
        ym = r ( 2, i )
        y ( i ) = pf * CMPLX ( yp, -ym, dbl )
      END DO
    END SELECT
  CASE ( 2 )
    SELECT CASE ( m )
    CASE DEFAULT
      call stop_program("y_lm","l = 2 and m value out of bounds")
    CASE ( 2 )
      pf = SQRT ( 15._dbl / ( 32._dbl * pi ) )
      DO i = 1, n
        yp = ( r ( 1, i ) * r ( 1, i ) - r ( 2, i ) * r ( 2, i ) )
        ym = 2._dbl * r ( 1, i ) * r ( 2, i )
        y ( i ) = pf * CMPLX ( yp, ym, dbl )
      END DO
    CASE ( 1 )
      pf = - SQRT ( 15._dbl / ( 8._dbl * pi ) )
      DO i = 1, n
        yp = r ( 3, i ) * r ( 1, i )
        ym = r ( 3, i ) * r ( 2, i )
        y ( i ) = pf * CMPLX ( yp, ym, dbl )
      END DO
    CASE ( 0 )
      pf = SQRT ( 5._dbl / ( 16._dbl * pi ) )
      y ( : ) = pf * ( 3._dbl * r ( 3, : ) * r ( 3, : ) - 1._dbl )
    CASE ( -1 )
      pf = - SQRT ( 15._dbl / ( 8._dbl * pi ) )
      DO i = 1, n
        yp = r ( 3, i ) * r ( 1, i )
        ym = r ( 3, i ) * r ( 2, i )
        y ( i ) = pf * CMPLX ( yp, -ym, dbl )
      END DO
    CASE ( -2 )
      pf = SQRT ( 15._dbl / ( 32._dbl * pi ) )
      DO i = 1, n
        yp = ( r ( 1, i ) * r ( 1, i ) - r ( 2, i ) * r ( 2, i ) )
        ym = 2._dbl * r ( 1, i ) * r ( 2, i )
        y ( i ) = pf * CMPLX ( yp, -ym, dbl )
      END DO
    END SELECT
  CASE ( 3 )
    SELECT CASE ( m )
    CASE DEFAULT
      call stop_program("y_lm","l = 3 and m value out of bounds")
    CASE ( 3 )
      pf = - SQRT ( 35._dbl / ( 64._dbl * pi ) )
      DO i = 1, n
        yp = r ( 1, i ) * ( r ( 1, i )**2 - 3._dbl * r ( 2, i )**2 )
        ym = r ( 2, i ) * ( 3._dbl * r ( 1, i )**2 - r ( 2, i )**2 )
        y ( i ) = pf * CMPLX ( yp, ym, dbl )
      END DO
    CASE ( 2 )
      pf = SQRT ( 105._dbl / ( 32._dbl * pi ) )
      DO i = 1, n
        yp = r ( 3, i ) * ( r ( 1, i )**2 - r ( 2, i )**2 )
        ym = 2._dbl * r ( 1, i ) * r ( 2, i ) * r ( 3, i )
        y ( i ) = pf * CMPLX ( yp, ym, dbl )
      END DO
    CASE ( 1 )
      pf = - SQRT ( 21._dbl / ( 64._dbl * pi ) )
      DO i = 1, n
        yp = r ( 1, i ) * ( 5._dbl * r ( 3, i ) * r ( 3, i ) - 1._dbl )
        ym = r ( 2, i ) * ( 5._dbl * r ( 3, i ) * r ( 3, i ) - 1._dbl )
        y ( i ) = pf * CMPLX ( yp, ym, dbl )
      END DO
    CASE ( 0 )
      pf = SQRT ( 7._dbl / ( 16._dbl * pi ) )
      y ( : ) = pf * r ( 3, : ) * ( 5._dbl * r ( 3, : ) * r ( 3, : ) - 3._dbl )
    CASE ( -1 )
      pf = - SQRT ( 21._dbl / ( 64._dbl * pi ) )
      DO i = 1, n
        yp = r ( 1, i ) * ( 5._dbl * r ( 3, i ) * r ( 3, i ) - 1._dbl )
        ym = r ( 2, i ) * ( 5._dbl * r ( 3, i ) * r ( 3, i ) - 1._dbl )
        y ( i ) = pf * CMPLX ( yp, -ym, dbl )
      END DO
    CASE ( -2 )
      pf = SQRT ( 105._dbl / ( 32._dbl * pi ) )
      DO i = 1, n
        yp = r ( 3, i ) * ( r ( 1, i )**2 - r ( 2, i )**2 )
        ym = 2._dbl * r ( 1, i ) * r ( 2, i ) * r ( 3, i )
        y ( i ) = pf * CMPLX ( yp, -ym, dbl )
      END DO
    CASE ( -3 )
      pf = - SQRT ( 35._dbl / ( 64._dbl * pi ) )
      DO i = 1, n
        yp = r ( 1, i ) * ( r ( 1, i )**2 - 3._dbl * r ( 2, i )**2 )
        ym = r ( 2, i ) * ( 3._dbl * r ( 1, i )**2 - r ( 2, i )**2 )
        y ( i ) = pf * CMPLX ( yp, -ym, dbl )
      END DO
    END SELECT
  CASE DEFAULT
    IF ( m < -l .OR. m > l ) call stop_program("y_lm","m value out of bounds")
    lpm = fac ( l + ABS ( m ) )
    lmm = fac ( l - ABS ( m ) )
    t = 4._dbl * pi
    IF ( ABS ( lpm ) < EPSILON ( 1._dbl ) ) THEN
      pf = REAL ( 2*l+1, dbl )  / t
    ELSE
      pf = ( REAL ( 2*l+1, dbl ) * lmm ) / ( t * lpm )
    ENDIF
    pf = sqrt ( pf )
    DO i = 1, n
      z = r ( 3, i )
      plm = legendre ( z, l, m )
      IF ( m == 0 ) THEN
        y ( i ) = pf * plm
      ELSE
        rxy = SQRT ( r ( 1, i ) ** 2  +  r ( 2, i ) ** 2 )
        IF ( rxy < EPSILON ( 1._dbl ) ) THEN
          y ( i ) = 0._dbl
        ELSE
          cp = r ( 1, i ) / rxy
          sp = r ( 2, i ) / rxy
          IF ( m > 0 ) THEN
            yp = cosn ( m, cp, sp )
            ym = sinn ( m, cp, sp )
            y ( i ) = pf * plm * CMPLX( yp, ym, dbl )
          ELSE
            yp = cosn ( -m, cp, sp )
            ym = sinn ( -m, cp, sp )
            y ( i ) = pf * plm * CMPLX( yp, -ym, dbl )
          END IF
        END IF
      END IF
    END DO
  END SELECT

  END SUBROUTINE cvy_lm

! *****************************************************************************

  SUBROUTINE ccy_lm ( r, y, l, m )
!
! Real Spherical Harmonics
!                   _                _
!                  | [(2l+1)(l-|m|)!] |1/2 m
!  Y_lm ( t, p ) = |------------------|   P_l(cos(t)) Exp[ i m p ]
!                  |  [4Pi(l+|m|)!]|  |
!
! Input: r == (x,y,z) : normalised    x^2 + y^2 + z^2 = 1
!
  REAL(dbl), DIMENSION (:), INTENT(IN) :: r
  COMPLEX(dbl), INTENT(OUT) :: y
  INTEGER, INTENT(IN) :: l, m

  INTEGER :: i
  REAL(dbl) :: pf, lpm, lmm, plm, z, rxy, t, cp, sp, yp, ym

  SELECT CASE ( l )
  CASE ( : -1 )
    call stop_program("y_lm","Negative l value")
  CASE ( 0 )
    pf = SQRT ( 1._dbl / ( 4._dbl * pi ) )
    IF ( m /= 0 )  call stop_program("y_lm","l = 0 and m value out of bounds")
    y = pf
  CASE ( 1 )
    SELECT CASE ( m )
    CASE DEFAULT
      call stop_program("y_lm","l = 1 and m value out of bounds")
    CASE ( 1 )
      pf = - SQRT ( 3._dbl / ( 8._dbl * pi ) )
      yp = r ( 1 )
      ym = r ( 2 )
      y = pf * CMPLX ( yp, ym, dbl )
    CASE ( 0 )
      pf = SQRT ( 3._dbl / ( 4._dbl * pi ) )
      y = pf * r ( 3 )
    CASE ( -1 )
      pf = - SQRT ( 3._dbl / ( 8._dbl * pi ) )
      yp = r ( 1 )
      ym = r ( 2 )
      y = pf * CMPLX ( yp, -ym, dbl )
    END SELECT
  CASE ( 2 )
    SELECT CASE ( m )
    CASE DEFAULT
      call stop_program("y_lm","l = 2 and m value out of bounds")
    CASE ( 2 )
      pf = SQRT ( 15._dbl / ( 32._dbl * pi ) )
      yp = ( r ( 1 ) * r ( 1 ) - r ( 2 ) * r ( 2 ) )
      ym = 2._dbl * r ( 1 ) * r ( 2 )
      y = pf * CMPLX ( yp, ym, dbl )
    CASE ( 1 )
      pf = - SQRT ( 15._dbl / ( 8._dbl * pi ) )
      yp = r ( 3 ) * r ( 1 )
      ym = r ( 3 ) * r ( 2 )
      y = pf * CMPLX ( yp, ym, dbl )
    CASE ( 0 )
      pf = SQRT ( 5._dbl / ( 16._dbl * pi ) )
      y = pf * ( 3._dbl * r ( 3 ) * r ( 3 ) - 1._dbl )
    CASE ( -1 )
      pf = - SQRT ( 15._dbl / ( 8._dbl * pi ) )
      yp = r ( 3 ) * r ( 1 )
      ym = r ( 3 ) * r ( 2 )
      y = pf * CMPLX ( yp, -ym, dbl )
    CASE ( -2 )
      pf = SQRT ( 15._dbl / ( 32._dbl * pi ) )
      yp = ( r ( 1 ) * r ( 1 ) - r ( 2 ) * r ( 2 ) )
      ym = 2._dbl * r ( 1 ) * r ( 2 )
      y = pf * CMPLX ( yp, -ym, dbl )
    END SELECT
  CASE ( 3 )
    SELECT CASE ( m )
    CASE DEFAULT
      call stop_program("y_lm","l = 3 and m value out of bounds")
    CASE ( 3 )
      pf = - SQRT ( 35._dbl / ( 64._dbl * pi ) )
      yp = r ( 1 ) * ( r ( 1 )**2 - 3._dbl * r ( 2 )**2 )
      ym = r ( 2 ) * ( 3._dbl * r ( 1 )**2 - r ( 2 )**2 )
      y = pf * CMPLX ( yp, ym, dbl )
    CASE ( 2 )
      pf = SQRT ( 105._dbl / ( 32._dbl * pi ) )
      yp = r ( 3 ) * ( r ( 1 )**2 - r ( 2 )**2 )
      ym = 2._dbl * r ( 1 ) * r ( 2 ) * r ( 3 )
      y = pf * CMPLX ( yp, ym, dbl )
    CASE ( 1 )
      pf = - SQRT ( 21._dbl / ( 64._dbl * pi ) )
      yp = r ( 1 ) * ( 5._dbl * r ( 3 ) * r ( 3 ) - 1._dbl )
      ym = r ( 2 ) * ( 5._dbl * r ( 3 ) * r ( 3 ) - 1._dbl )
      y = pf * CMPLX ( yp, ym, dbl )
    CASE ( 0 )
      pf = SQRT ( 7._dbl / ( 16._dbl * pi ) )
      y = pf * r ( 3 ) * ( 5._dbl * r ( 3 ) * r ( 3 ) - 3._dbl )
    CASE ( -1 )
      pf = - SQRT ( 21._dbl / ( 64._dbl * pi ) )
      yp = r ( 1 ) * ( 5._dbl * r ( 3 ) * r ( 3 ) - 1._dbl )
      ym = r ( 2 ) * ( 5._dbl * r ( 3 ) * r ( 3 ) - 1._dbl )
      y = pf * CMPLX ( yp, -ym, dbl )
    CASE ( -2 )
      pf = SQRT ( 105._dbl / ( 32._dbl * pi ) )
      yp = r ( 3 ) * ( r ( 1 )**2 - r ( 2 )**2 )
      ym = 2._dbl * r ( 1 ) * r ( 2 ) * r ( 3 )
      y = pf * CMPLX ( yp, -ym, dbl )
    CASE ( -3 )
      pf = - SQRT ( 35._dbl / ( 64._dbl * pi ) )
      yp = r ( 1 ) * ( r ( 1 )**2 - 3._dbl * r ( 2 )**2 )
      ym = r ( 2 ) * ( 3._dbl * r ( 1 )**2 - r ( 2 )**2 )
      y = pf * CMPLX ( yp, -ym, dbl )
    END SELECT
  CASE DEFAULT
    IF ( m < -l .OR. m > l ) call stop_program("y_lm","m value out of bounds")
    lpm = fac ( l + ABS ( m ) )
    lmm = fac ( l - ABS ( m ) )
    t = 4._dbl * pi
    IF ( ABS ( lpm ) < EPSILON ( 1._dbl ) ) THEN
      pf = REAL ( 2*l+1, dbl )  / t
    ELSE
      pf = ( REAL ( 2*l+1, dbl ) * lmm ) / ( t * lpm )
    ENDIF
    pf = sqrt ( pf )
    z = r ( 3 )
    plm = legendre ( z, l, m )
    IF ( m == 0 ) THEN
      y = pf * plm
    ELSE
      rxy = SQRT ( r ( 1 ) ** 2  +  r ( 2 ) ** 2 )
      IF ( rxy < EPSILON ( 1._dbl ) ) THEN
        y = 0._dbl
      ELSE
        cp = r ( 1 ) / rxy
        sp = r ( 2 ) / rxy
        IF ( m > 0 ) THEN
          yp = cosn ( m, cp, sp )
          ym = sinn ( m , cp, sp )
          y = pf * plm * CMPLX( yp, ym, dbl )
        ELSE
          yp = cosn ( -m, cp, sp )
          ym = sinn ( -m , cp, sp )
          y = pf * plm * CMPLX( yp, -ym, dbl )
        END IF
      END IF
    END IF
  END SELECT

  END SUBROUTINE ccy_lm

! *****************************************************************************

  FUNCTION legendre ( x, l, m ) RESULT ( plm )

  REAL(dbl), INTENT(IN) :: x
  INTEGER, INTENT(IN) :: l, m
  REAL(dbl) :: plm
  INTEGER :: mm, err, im, il
  REAL(dbl) :: pmm, somx2, fact, pmmp1, pll

  IF ( ABS ( x ) > 1._dbl ) call stop_program("legendre","x value > 1")
  SELECT CASE ( l )
  CASE ( : -1 )
    call stop_program("legendre","Negative l value")
  CASE ( 0 )
    plm = 1._dbl
  CASE ( 1 )
    SELECT CASE ( ABS ( m ) )
    CASE DEFAULT
      call stop_program("legendre","l = 1 and m value out of bounds")
    CASE ( 1 )
      plm = - SQRT ( 1._dbl - x * x )
    CASE ( 0 )
      plm = x
    END SELECT
  CASE ( 2 )
    SELECT CASE ( ABS ( m ) )
    CASE DEFAULT
      call stop_program("legendre","l = 2 and m value out of bounds")
    CASE ( 2 )
      plm = 3._dbl * ( 1._dbl - x * x )
    CASE ( 1 )
      plm = - 3._dbl * x * SQRT ( 1._dbl - x * x )
    CASE ( 0 )
      plm = 1.5_dbl * x * x - 0.5_dbl
    END SELECT
  CASE ( 3 )
    SELECT CASE ( ABS ( m ) )
    CASE DEFAULT
      call stop_program("legendre","l = 3 and m value out of bounds")
    CASE ( 3 )
      plm = - 15._dbl * ( 1._dbl - x * x ) ** 1.5_dbl
    CASE ( 2 )
      plm = 15._dbl * x * ( 1._dbl - x * x )
    CASE ( 1 )
      plm = ( -7.5_dbl * x * x + 1.5_dbl ) * SQRT ( 1._dbl - x * x )
    CASE ( 0 )
      plm = 2.5_dbl * x**3 - 1.5_dbl * x
    END SELECT
  CASE ( 4 )
    SELECT CASE ( ABS ( m ) )
    CASE DEFAULT
      call stop_program("legendre","l = 4 and m value out of bounds")
    CASE ( 4 )
      plm = 105._dbl * ( 1._dbl - x * x ) ** 2
    CASE ( 3 )
      plm = - 105._dbl * x * ( 1._dbl - x * x ) ** 1.5_dbl
    CASE ( 2 )
      plm = ( 52.5_dbl * x * x - 7.5_dbl ) * ( 1._dbl - x * x )
    CASE ( 1 )
      plm = ( -17.5_dbl * x ** 3 + 7.5_dbl * x ) * SQRT ( 1._dbl - x * x )
    CASE ( 0 )
      plm = 4.375_dbl * x ** 4 - 3.75_dbl * x ** 2 + 0.375_dbl
    END SELECT
  CASE ( 5 )
    SELECT CASE ( ABS ( m ) )
    CASE DEFAULT
      call stop_program("legendre","l = 5 and m value out of bounds")
    CASE ( 5 )
      plm = - 945._dbl * ( 1._dbl - x * x ) ** 2.5_dbl
    CASE ( 4 )
      plm = 945._dbl * x * ( 1._dbl - x * x ) ** 2
    CASE ( 3 )
      plm = ( -472.5_dbl * x * x + 52.5_dbl ) * ( 1._dbl - x * x ) ** 1.5_dbl
    CASE ( 2 )
      plm = ( 157.5_dbl * x ** 3 - 52.5_dbl * x ) * ( 1._dbl - x * x )
    CASE ( 1 )
      plm = ( -39.375_dbl * x ** 4 + 26.25_dbl * x * x - &
              1.875_dbl ) * SQRT ( 1._dbl - x * x )
    CASE ( 0 )
      plm = 7.875_dbl * x ** 5 - 8.75_dbl * x ** 3 + 1.875_dbl * x
    END SELECT
  CASE ( 6 )
    SELECT CASE ( ABS ( m ) )
    CASE DEFAULT
      call stop_program("legendre","l = 6 and m value out of bounds")
    CASE ( 6 )
      plm = 10395._dbl * ( 1._dbl - x * x ) ** 3
    CASE ( 5 )
      plm = - 10395._dbl * x * ( 1._dbl - x * x ) ** 2.5_dbl
    CASE ( 4 )
      plm = ( 5197.5_dbl * x * x - 472.5_dbl ) * ( 1._dbl - x * x ) ** 2
    CASE ( 3 )
      plm = ( - 1732.5_dbl * x ** 3 + 472.5_dbl * x ) * &
            ( 1._dbl - x * x ) ** 1.5_dbl
    CASE ( 2 )
      plm = ( 433.125_dbl * x ** 4 - 236.25_dbl * x ** 2 + &
              13.125_dbl ) * ( 1._dbl - x * x )
    CASE ( 1 )
      plm = ( -86.625_dbl * x ** 5 + 78.75_dbl * x ** 3 - &
              13.125_dbl * x ) * SQRT ( 1._dbl - x * x )
    CASE ( 0 )
      plm = 14.4375_dbl * x ** 6 - 19.6875_dbl * x ** 4 + &
            6.5625_dbl * x ** 2 - 0.3125_dbl
    END SELECT
  CASE DEFAULT
    mm = ABS ( m )
    IF ( mm > l ) call stop_program ( " legendre ","m out of bounds" )
! use recurence from numerical recipies
    pmm = 1._dbl
    IF ( mm > 0 ) THEN
      somx2 = sqrt((1._dbl-x)*(1._dbl+x))
      fact = 1._dbl
      DO im = 1, mm
        pmm = pmm * fact * somx2
        fact = fact + 2._dbl
      END DO
    END IF
    IF ( l == mm ) THEN
      plm = pmm
    ELSE
      pmmp1 = x * REAL ( 2 * mm + 1, dbl ) * pmm
      IF ( l == mm + 1 ) THEN
        plm = pmmp1
      ELSE
        DO il = mm + 2, l
          pll = ( x * REAL ( 2 * il - 1, dbl ) * pmmp1 - &
                      REAL ( il + mm - 1, dbl ) * pmm ) / REAL ( il - mm, dbl )
          pmm = pmmp1
          pmmp1 = pll
        END DO
        plm = pll
      END IF
    END IF
  END SELECT

  END FUNCTION legendre

! *****************************************************************************

  FUNCTION choose ( n, k )

  INTEGER, INTENT ( IN ) :: n, k
  REAL ( dbl ) :: choose

  IF ( n >= k ) THEN
    choose = REAL ( NINT ( fac ( n ) / ( fac ( k ) * fac ( n - k ) ) ), dbl )
  ELSE
    choose = 0._dbl
  ENDIF

  END FUNCTION choose

! *****************************************************************************

  FUNCTION cosn ( n, c, s )

  INTEGER, INTENT ( IN ) :: n
  REAL(dbl), INTENT ( IN ) :: c, s
  REAL(dbl) :: cosn
  INTEGER :: i, j

  cosn = 0._dbl

  DO i = n, 0, -2
    j = n - i
    IF ( MOD(j/2,2) == 0 ) THEN
      cosn = cosn + choose ( n, j ) * c**i * s**j
    ELSE
      cosn = cosn - choose ( n, j ) * c**i * s**j
    END IF
  END DO

  END FUNCTION cosn

! *****************************************************************************

  FUNCTION sinn ( n, c, s )

  INTEGER, INTENT ( IN ) :: n
  REAL(dbl), INTENT ( IN ) :: c, s
  REAL(dbl) :: sinn
  INTEGER :: i, j

  sinn = 0._dbl

  DO i = n-1, 0, -2
    j = n - i
    IF ( MOD(j/2,2) == 0 ) THEN
      sinn = sinn + choose ( n, j ) * c**i * s**j
    ELSE
      sinn = sinn - choose ( n, j ) * c**i * s**j
    END IF
  END DO

  END FUNCTION sinn

! *****************************************************************************

END MODULE spherical_harmonics

