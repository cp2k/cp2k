!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculate spherical harmonics
!> \note
!>      Spherical Harmonics
!>          Numerical Stability up to L=15
!>          Accuracy > 1.E-12 up to L=15 tested
!>          Definition is consistent with orbital_transformation_matrices
!> \todo
!>      1) Check if this definition is consistent with the
!>         Slater-Koster module
!> \par History
!>      JGH 28-Feb-2002 : Change of sign convention (-1^m)
!>      JGH  1-Mar-2002 : Clebsch-Gordon Coefficients
!>      - Clebsch-Gordon coefficients and Wigner 3-j symbols added as generic routines using the
!>        standard normalization (19.09.2022, MK)
!>      - Clebsch-Gordon coefficients (which were really Gaunt coefficients of the
!>        real spherical harmonics) moved to gaunt_coefficients, and the unused
!>        generic Clebsch-Gordon and Wigner 3-j routines dropped, in favour of the
!>        exact evaluation in libwignernj (09.2026, Susi Lehtola)
!> \author JGH 6-Oct-2000, MK
! **************************************************************************************************
MODULE spherical_harmonics

   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fac,&
                                              pi
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: y_lm, dy_lm
   PUBLIC :: legendre, dlegendre

   INTERFACE y_lm
      MODULE PROCEDURE rvy_lm, rry_lm, cvy_lm, ccy_lm
   END INTERFACE

   INTERFACE dy_lm
      MODULE PROCEDURE dry_lm, dcy_lm
   END INTERFACE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'spherical_harmonics'

CONTAINS

! Calculation of Spherical Harmonics

! **************************************************************************************************
!> \brief ...
!> \param r ...
!> \param y ...
!> \param l ...
!> \param m ...
! **************************************************************************************************
   SUBROUTINE rvy_lm(r, y, l, m)
!
! Real Spherical Harmonics
!                   _                   _
!                  |  [(2l+1)(l-|m|)!]   |1/2 m         cos(m p)   m>=0
!  Y_lm ( t, p ) = |---------------------|   P_l(cos(t))
!                  |[2Pi(1+d_m0)(l+|m|)!]|              sin(|m| p) m<0
!
! Input: r == (x,y,z) : normalised    x^2 + y^2 + z^2 = 1
!
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: r
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: y
      INTEGER, INTENT(IN)                                :: l, m

      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: cp, lmm, lpm, pf, plm, rxy, sp, t, z

      SELECT CASE (l)
      CASE (:-1)
         CPABORT("Negative l value")
      CASE (0)
         pf = SQRT(1.0_dp/(4.0_dp*pi))
         IF (m /= 0) CPABORT("l = 0 and m value out of bounds")
         y(:) = pf
      CASE (1)
         SELECT CASE (m)
         CASE DEFAULT
            CPABORT("l = 1 and m value out of bounds")
         CASE (1)
            pf = SQRT(3.0_dp/(4.0_dp*pi))
            y(:) = pf*r(1, :)
         CASE (0)
            pf = SQRT(3.0_dp/(4.0_dp*pi))
            y(:) = pf*r(3, :)
         CASE (-1)
            pf = SQRT(3.0_dp/(4.0_dp*pi))
            y(:) = pf*r(2, :)
         END SELECT
      CASE (2)
         SELECT CASE (m)
         CASE DEFAULT
            CPABORT("l = 2 and m value out of bounds")
         CASE (2)
            pf = SQRT(15.0_dp/(16.0_dp*pi))
            y(:) = pf*(r(1, :)*r(1, :) - r(2, :)*r(2, :))
         CASE (1)
            pf = SQRT(15.0_dp/(4.0_dp*pi))
            y(:) = pf*r(3, :)*r(1, :)
         CASE (0)
            pf = SQRT(5.0_dp/(16.0_dp*pi))
            y(:) = pf*(3.0_dp*r(3, :)*r(3, :) - 1.0_dp)
         CASE (-1)
            pf = SQRT(15.0_dp/(4.0_dp*pi))
            y(:) = pf*r(3, :)*r(2, :)
         CASE (-2)
            pf = SQRT(15.0_dp/(16.0_dp*pi))
            y(:) = pf*2.0_dp*r(1, :)*r(2, :)
         END SELECT
      CASE (3)
         SELECT CASE (m)
         CASE DEFAULT
            CPABORT("l = 3 and m value out of bounds")
         CASE (3)
            pf = SQRT(35.0_dp/(32.0_dp*pi))
            y(:) = pf*r(1, :)*(r(1, :)**2 - 3.0_dp*r(2, :)**2)
         CASE (2)
            pf = SQRT(105.0_dp/(16.0_dp*pi))
            y(:) = pf*r(3, :)*(r(1, :)**2 - r(2, :)**2)
         CASE (1)
            pf = SQRT(21.0_dp/(32.0_dp*pi))
            y(:) = pf*r(1, :)*(5.0_dp*r(3, :)*r(3, :) - 1.0_dp)
         CASE (0)
            pf = SQRT(7.0_dp/(16.0_dp*pi))
            y(:) = pf*r(3, :)*(5.0_dp*r(3, :)*r(3, :) - 3.0_dp)
         CASE (-1)
            pf = SQRT(21.0_dp/(32.0_dp*pi))
            y(:) = pf*r(2, :)*(5.0_dp*r(3, :)*r(3, :) - 1.0_dp)
         CASE (-2)
            pf = SQRT(105.0_dp/(16.0_dp*pi))
            y(:) = pf*2.0_dp*r(1, :)*r(2, :)*r(3, :)
         CASE (-3)
            pf = SQRT(35.0_dp/(32.0_dp*pi))
            y(:) = pf*r(2, :)*(3.0_dp*r(1, :)**2 - r(2, :)**2)
         END SELECT
      CASE DEFAULT
         IF (m < -l .OR. m > l) CPABORT("m value out of bounds")
         lpm = fac(l + ABS(m))
         lmm = fac(l - ABS(m))
         IF (m == 0) THEN
            t = 4.0_dp*pi
         ELSE
            t = 2.0_dp*pi
         END IF
         IF (ABS(lpm) < EPSILON(1.0_dp)) THEN
            pf = REAL(2*l + 1, KIND=dp)/t
         ELSE
            pf = (REAL(2*l + 1, KIND=dp)*lmm)/(t*lpm)
         END IF
         pf = SQRT(pf)
         DO i = 1, SIZE(r, 2)
            z = r(3, i)
!      plm = legendre ( z, l, m )
            plm = legendre(z, l, ABS(m))
            IF (m == 0) THEN
               y(i) = pf*plm
            ELSE
               rxy = SQRT(r(1, i)**2 + r(2, i)**2)
               IF (rxy < EPSILON(1.0_dp)) THEN
                  y(i) = 0.0_dp
               ELSE
                  cp = r(1, i)/rxy
                  sp = r(2, i)/rxy
                  IF (m > 0) THEN
                     y(i) = pf*plm*cosn(m, cp, sp)
                  ELSE
                     y(i) = pf*plm*sinn(-m, cp, sp)
                  END IF
               END IF
            END IF
         END DO
      END SELECT

   END SUBROUTINE rvy_lm

! **************************************************************************************************
!> \brief ...
!> \param r ...
!> \param y ...
!> \param l ...
!> \param m ...
! **************************************************************************************************
   SUBROUTINE rry_lm(r, y, l, m)
!
! Real Spherical Harmonics
!                   _                   _
!                  |  [(2l+1)(l-|m|)!]   |1/2 m         cos(m p)   m>=0
!  Y_lm ( t, p ) = |---------------------|   P_l(cos(t))
!                  |[2Pi(1+d_m0)(l+|m|)!]|              sin(|m| p) m<0
!
! Input: r == (x,y,z) : normalised    x^2 + y^2 + z^2 = 1
!
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: r
      REAL(KIND=dp), INTENT(OUT)                         :: y
      INTEGER, INTENT(IN)                                :: l, m

      REAL(KIND=dp)                                      :: cp, lmm, lpm, pf, plm, rxy, sp, t, z

      SELECT CASE (l)
      CASE (:-1)
         CPABORT("Negative l value")
      CASE (0)
         pf = SQRT(1.0_dp/(4.0_dp*pi))
         IF (m /= 0) CPABORT("l = 0 and m value out of bounds")
         y = pf
      CASE (1)
         SELECT CASE (m)
         CASE DEFAULT
            CPABORT("l = 1 and m value out of bounds")
         CASE (1)
            pf = SQRT(3.0_dp/(4.0_dp*pi))
            y = pf*r(1)
         CASE (0)
            pf = SQRT(3.0_dp/(4.0_dp*pi))
            y = pf*r(3)
         CASE (-1)
            pf = SQRT(3.0_dp/(4.0_dp*pi))
            y = pf*r(2)
         END SELECT
      CASE (2)
         SELECT CASE (m)
         CASE DEFAULT
            CPABORT("l = 2 and m value out of bounds")
         CASE (2)
            pf = SQRT(15.0_dp/(16.0_dp*pi))
            y = pf*(r(1)*r(1) - r(2)*r(2))
         CASE (1)
            pf = SQRT(15.0_dp/(4.0_dp*pi))
            y = pf*r(3)*r(1)
         CASE (0)
            pf = SQRT(5.0_dp/(16.0_dp*pi))
            y = pf*(3.0_dp*r(3)*r(3) - 1.0_dp)
         CASE (-1)
            pf = SQRT(15.0_dp/(4.0_dp*pi))
            y = pf*r(3)*r(2)
         CASE (-2)
            pf = SQRT(15.0_dp/(16.0_dp*pi))
            y = pf*2.0_dp*r(1)*r(2)
         END SELECT
      CASE (3)
         SELECT CASE (m)
         CASE DEFAULT
            CPABORT("l = 3 and m value out of bounds")
         CASE (3)
            pf = SQRT(35.0_dp/(32.0_dp*pi))
            y = pf*r(1)*(r(1)**2 - 3.0_dp*r(2)**2)
         CASE (2)
            pf = SQRT(105.0_dp/(16.0_dp*pi))
            y = pf*r(3)*(r(1)**2 - r(2)**2)
         CASE (1)
            pf = SQRT(21.0_dp/(32.0_dp*pi))
            y = pf*r(1)*(5.0_dp*r(3)*r(3) - 1.0_dp)
         CASE (0)
            pf = SQRT(7.0_dp/(16.0_dp*pi))
            y = pf*r(3)*(5.0_dp*r(3)*r(3) - 3.0_dp)
         CASE (-1)
            pf = SQRT(21.0_dp/(32.0_dp*pi))
            y = pf*r(2)*(5.0_dp*r(3)*r(3) - 1.0_dp)
         CASE (-2)
            pf = SQRT(105.0_dp/(16.0_dp*pi))
            y = pf*2.0_dp*r(1)*r(2)*r(3)
         CASE (-3)
            pf = SQRT(35.0_dp/(32.0_dp*pi))
            y = pf*r(2)*(3.0_dp*r(1)**2 - r(2)**2)
         END SELECT
      CASE DEFAULT
         IF (m < -l .OR. m > l) CPABORT("m value out of bounds")
         lpm = fac(l + ABS(m))
         lmm = fac(l - ABS(m))
         IF (m == 0) THEN
            t = 4.0_dp*pi
         ELSE
            t = 2.0_dp*pi
         END IF
         IF (ABS(lpm) < EPSILON(1.0_dp)) THEN
            pf = REAL(2*l + 1, KIND=dp)/t
         ELSE
            pf = (REAL(2*l + 1, KIND=dp)*lmm)/(t*lpm)
         END IF
         pf = SQRT(pf)
         z = r(3)
         plm = legendre(z, l, m)
         IF (m == 0) THEN
            y = pf*plm
         ELSE
            rxy = SQRT(r(1)**2 + r(2)**2)
            IF (rxy < EPSILON(1.0_dp)) THEN
               y = 0.0_dp
            ELSE
               cp = r(1)/rxy
               sp = r(2)/rxy
               IF (m > 0) THEN
                  y = pf*plm*cosn(m, cp, sp)
               ELSE
                  y = pf*plm*sinn(-m, cp, sp)
               END IF
            END IF
         END IF
      END SELECT

   END SUBROUTINE rry_lm

! **************************************************************************************************
!> \brief ...
!> \param r ...
!> \param y ...
!> \param l ...
!> \param m ...
! **************************************************************************************************
   SUBROUTINE cvy_lm(r, y, l, m)
!
! Complex Spherical Harmonics
!                   _                _
!                  | [(2l+1)(l-|m|)!] |1/2 m
!  Y_lm ( t, p ) = |------------------|   P_l(cos(t)) Exp[ i m p ]
!                  |  [4Pi(l+|m|)!]|  |
!
! Input: r == (x,y,z) : normalised    x^2 + y^2 + z^2 = 1
!
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: r
      COMPLEX(KIND=dp), DIMENSION(:), INTENT(OUT)        :: y
      INTEGER, INTENT(IN)                                :: l, m

      INTEGER                                            :: i, n
      REAL(KIND=dp)                                      :: cp, lmm, lpm, pf, plm, rxy, sp, t, ym, &
                                                            yp, z

      n = SIZE(r, 2)
      SELECT CASE (l)
      CASE (:-1)
         CPABORT("Negative l value")
      CASE (0)
         pf = SQRT(1.0_dp/(4.0_dp*pi))
         IF (m /= 0) CPABORT("l = 0 and m value out of bounds")
         y(:) = pf
      CASE (1)
         SELECT CASE (m)
         CASE DEFAULT
            CPABORT("l = 1 and m value out of bounds")
         CASE (1)
            pf = SQRT(3.0_dp/(8.0_dp*pi))
            DO i = 1, n
               yp = r(1, i)
               ym = r(2, i)
               y(i) = pf*CMPLX(yp, ym, KIND=dp)
            END DO
         CASE (0)
            pf = SQRT(3.0_dp/(4.0_dp*pi))
            y(:) = pf*r(3, :)
         CASE (-1)
            pf = SQRT(3.0_dp/(8.0_dp*pi))
            DO i = 1, n
               yp = r(1, i)
               ym = r(2, i)
               y(i) = pf*CMPLX(yp, -ym, KIND=dp)
            END DO
         END SELECT
      CASE (2)
         SELECT CASE (m)
         CASE DEFAULT
            CPABORT("l = 2 and m value out of bounds")
         CASE (2)
            pf = SQRT(15.0_dp/(32.0_dp*pi))
            DO i = 1, n
               yp = (r(1, i)*r(1, i) - r(2, i)*r(2, i))
               ym = 2.0_dp*r(1, i)*r(2, i)
               y(i) = pf*CMPLX(yp, ym, KIND=dp)
            END DO
         CASE (1)
            pf = SQRT(15.0_dp/(8.0_dp*pi))
            DO i = 1, n
               yp = r(3, i)*r(1, i)
               ym = r(3, i)*r(2, i)
               y(i) = pf*CMPLX(yp, ym, KIND=dp)
            END DO
         CASE (0)
            pf = SQRT(5.0_dp/(16.0_dp*pi))
            y(:) = pf*(3.0_dp*r(3, :)*r(3, :) - 1.0_dp)
         CASE (-1)
            pf = SQRT(15.0_dp/(8.0_dp*pi))
            DO i = 1, n
               yp = r(3, i)*r(1, i)
               ym = r(3, i)*r(2, i)
               y(i) = pf*CMPLX(yp, -ym, KIND=dp)
            END DO
         CASE (-2)
            pf = SQRT(15.0_dp/(32.0_dp*pi))
            DO i = 1, n
               yp = (r(1, i)*r(1, i) - r(2, i)*r(2, i))
               ym = 2.0_dp*r(1, i)*r(2, i)
               y(i) = pf*CMPLX(yp, -ym, KIND=dp)
            END DO
         END SELECT
      CASE (3)
         SELECT CASE (m)
         CASE DEFAULT
            CPABORT("l = 3 and m value out of bounds")
         CASE (3)
            pf = SQRT(35.0_dp/(64.0_dp*pi))
            DO i = 1, n
               yp = r(1, i)*(r(1, i)**2 - 3.0_dp*r(2, i)**2)
               ym = r(2, i)*(3.0_dp*r(1, i)**2 - r(2, i)**2)
               y(i) = pf*CMPLX(yp, ym, KIND=dp)
            END DO
         CASE (2)
            pf = SQRT(105.0_dp/(32.0_dp*pi))
            DO i = 1, n
               yp = r(3, i)*(r(1, i)**2 - r(2, i)**2)
               ym = 2.0_dp*r(1, i)*r(2, i)*r(3, i)
               y(i) = pf*CMPLX(yp, ym, KIND=dp)
            END DO
         CASE (1)
            pf = SQRT(21.0_dp/(64.0_dp*pi))
            DO i = 1, n
               yp = r(1, i)*(5.0_dp*r(3, i)*r(3, i) - 1.0_dp)
               ym = r(2, i)*(5.0_dp*r(3, i)*r(3, i) - 1.0_dp)
               y(i) = pf*CMPLX(yp, ym, KIND=dp)
            END DO
         CASE (0)
            pf = SQRT(7.0_dp/(16.0_dp*pi))
            y(:) = pf*r(3, :)*(5.0_dp*r(3, :)*r(3, :) - 3.0_dp)
         CASE (-1)
            pf = SQRT(21.0_dp/(64.0_dp*pi))
            DO i = 1, n
               yp = r(1, i)*(5.0_dp*r(3, i)*r(3, i) - 1.0_dp)
               ym = r(2, i)*(5.0_dp*r(3, i)*r(3, i) - 1.0_dp)
               y(i) = pf*CMPLX(yp, -ym, KIND=dp)
            END DO
         CASE (-2)
            pf = SQRT(105.0_dp/(32.0_dp*pi))
            DO i = 1, n
               yp = r(3, i)*(r(1, i)**2 - r(2, i)**2)
               ym = 2.0_dp*r(1, i)*r(2, i)*r(3, i)
               y(i) = pf*CMPLX(yp, -ym, KIND=dp)
            END DO
         CASE (-3)
            pf = SQRT(35.0_dp/(64.0_dp*pi))
            DO i = 1, n
               yp = r(1, i)*(r(1, i)**2 - 3.0_dp*r(2, i)**2)
               ym = r(2, i)*(3.0_dp*r(1, i)**2 - r(2, i)**2)
               y(i) = pf*CMPLX(yp, -ym, KIND=dp)
            END DO
         END SELECT
      CASE DEFAULT
         IF (m < -l .OR. m > l) CPABORT("m value out of bounds")
         lpm = fac(l + ABS(m))
         lmm = fac(l - ABS(m))
         t = 4.0_dp*pi
         IF (ABS(lpm) < EPSILON(1.0_dp)) THEN
            pf = REAL(2*l + 1, KIND=dp)/t
         ELSE
            pf = (REAL(2*l + 1, KIND=dp)*lmm)/(t*lpm)
         END IF
         pf = SQRT(pf)
         DO i = 1, n
            z = r(3, i)
            plm = legendre(z, l, m)
            IF (m == 0) THEN
               y(i) = pf*plm
            ELSE
               rxy = SQRT(r(1, i)**2 + r(2, i)**2)
               IF (rxy < EPSILON(1.0_dp)) THEN
                  y(i) = 0.0_dp
               ELSE
                  cp = r(1, i)/rxy
                  sp = r(2, i)/rxy
                  IF (m > 0) THEN
                     yp = cosn(m, cp, sp)
                     ym = sinn(m, cp, sp)
                     y(i) = pf*plm*CMPLX(yp, ym, KIND=dp)
                  ELSE
                     yp = cosn(-m, cp, sp)
                     ym = sinn(-m, cp, sp)
                     y(i) = pf*plm*CMPLX(yp, -ym, KIND=dp)
                  END IF
               END IF
            END IF
         END DO
      END SELECT

   END SUBROUTINE cvy_lm

! **************************************************************************************************
!> \brief ...
!> \param r ...
!> \param y ...
!> \param l ...
!> \param m ...
! **************************************************************************************************
   SUBROUTINE ccy_lm(r, y, l, m)
!
! Complex Spherical Harmonics
!                   _                _
!                  | [(2l+1)(l-|m|)!] |1/2 m
!  Y_lm ( t, p ) = |------------------|   P_l(cos(t)) Exp[ i m p ]
!                  |  [4Pi(l+|m|)!]|  |
!
! Input: r == (x,y,z) : normalised    x^2 + y^2 + z^2 = 1
!
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: r
      COMPLEX(KIND=dp), INTENT(OUT)                      :: y
      INTEGER, INTENT(IN)                                :: l, m

      REAL(KIND=dp)                                      :: cp, lmm, lpm, pf, plm, rxy, sp, t, ym, &
                                                            yp, z

      SELECT CASE (l)
      CASE (:-1)
         CPABORT("Negative l value")
      CASE (0)
         pf = SQRT(1.0_dp/(4.0_dp*pi))
         IF (m /= 0) CPABORT("l = 0 and m value out of bounds")
         y = pf
      CASE (1)
         SELECT CASE (m)
         CASE DEFAULT
            CPABORT("l = 1 and m value out of bounds")
         CASE (1)
            pf = SQRT(3.0_dp/(8.0_dp*pi))
            yp = r(1)
            ym = r(2)
            y = pf*CMPLX(yp, ym, KIND=dp)
         CASE (0)
            pf = SQRT(3.0_dp/(4.0_dp*pi))
            y = pf*r(3)
         CASE (-1)
            pf = SQRT(3.0_dp/(8.0_dp*pi))
            yp = r(1)
            ym = r(2)
            y = pf*CMPLX(yp, -ym, KIND=dp)
         END SELECT
      CASE (2)
         SELECT CASE (m)
         CASE DEFAULT
            CPABORT("l = 2 and m value out of bounds")
         CASE (2)
            pf = SQRT(15.0_dp/(32.0_dp*pi))
            yp = (r(1)*r(1) - r(2)*r(2))
            ym = 2.0_dp*r(1)*r(2)
            y = pf*CMPLX(yp, ym, KIND=dp)
         CASE (1)
            pf = SQRT(15.0_dp/(8.0_dp*pi))
            yp = r(3)*r(1)
            ym = r(3)*r(2)
            y = pf*CMPLX(yp, ym, KIND=dp)
         CASE (0)
            pf = SQRT(5.0_dp/(16.0_dp*pi))
            y = pf*(3.0_dp*r(3)*r(3) - 1.0_dp)
         CASE (-1)
            pf = SQRT(15.0_dp/(8.0_dp*pi))
            yp = r(3)*r(1)
            ym = r(3)*r(2)
            y = pf*CMPLX(yp, -ym, KIND=dp)
         CASE (-2)
            pf = SQRT(15.0_dp/(32.0_dp*pi))
            yp = (r(1)*r(1) - r(2)*r(2))
            ym = 2.0_dp*r(1)*r(2)
            y = pf*CMPLX(yp, -ym, KIND=dp)
         END SELECT
      CASE (3)
         SELECT CASE (m)
         CASE DEFAULT
            CPABORT("l = 3 and m value out of bounds")
         CASE (3)
            pf = SQRT(35.0_dp/(64.0_dp*pi))
            yp = r(1)*(r(1)**2 - 3.0_dp*r(2)**2)
            ym = r(2)*(3.0_dp*r(1)**2 - r(2)**2)
            y = pf*CMPLX(yp, ym, KIND=dp)
         CASE (2)
            pf = SQRT(105.0_dp/(32.0_dp*pi))
            yp = r(3)*(r(1)**2 - r(2)**2)
            ym = 2.0_dp*r(1)*r(2)*r(3)
            y = pf*CMPLX(yp, ym, KIND=dp)
         CASE (1)
            pf = SQRT(21.0_dp/(64.0_dp*pi))
            yp = r(1)*(5.0_dp*r(3)*r(3) - 1.0_dp)
            ym = r(2)*(5.0_dp*r(3)*r(3) - 1.0_dp)
            y = pf*CMPLX(yp, ym, KIND=dp)
         CASE (0)
            pf = SQRT(7.0_dp/(16.0_dp*pi))
            y = pf*r(3)*(5.0_dp*r(3)*r(3) - 3.0_dp)
         CASE (-1)
            pf = SQRT(21.0_dp/(64.0_dp*pi))
            yp = r(1)*(5.0_dp*r(3)*r(3) - 1.0_dp)
            ym = r(2)*(5.0_dp*r(3)*r(3) - 1.0_dp)
            y = pf*CMPLX(yp, -ym, KIND=dp)
         CASE (-2)
            pf = SQRT(105.0_dp/(32.0_dp*pi))
            yp = r(3)*(r(1)**2 - r(2)**2)
            ym = 2.0_dp*r(1)*r(2)*r(3)
            y = pf*CMPLX(yp, -ym, KIND=dp)
         CASE (-3)
            pf = SQRT(35.0_dp/(64.0_dp*pi))
            yp = r(1)*(r(1)**2 - 3.0_dp*r(2)**2)
            ym = r(2)*(3.0_dp*r(1)**2 - r(2)**2)
            y = pf*CMPLX(yp, -ym, KIND=dp)
         END SELECT
      CASE DEFAULT
         IF (m < -l .OR. m > l) CPABORT("m value out of bounds")
         lpm = fac(l + ABS(m))
         lmm = fac(l - ABS(m))
         t = 4.0_dp*pi
         IF (ABS(lpm) < EPSILON(1.0_dp)) THEN
            pf = REAL(2*l + 1, KIND=dp)/t
         ELSE
            pf = (REAL(2*l + 1, KIND=dp)*lmm)/(t*lpm)
         END IF
         pf = SQRT(pf)
         z = r(3)
         plm = legendre(z, l, m)
         IF (m == 0) THEN
            y = pf*plm
         ELSE
            rxy = SQRT(r(1)**2 + r(2)**2)
            IF (rxy < EPSILON(1.0_dp)) THEN
               y = 0.0_dp
            ELSE
               cp = r(1)/rxy
               sp = r(2)/rxy
               IF (m > 0) THEN
                  yp = cosn(m, cp, sp)
                  ym = sinn(m, cp, sp)
                  y = pf*plm*CMPLX(yp, ym, KIND=dp)
               ELSE
                  yp = cosn(-m, cp, sp)
                  ym = sinn(-m, cp, sp)
                  y = pf*plm*CMPLX(yp, -ym, KIND=dp)
               END IF
            END IF
         END IF
      END SELECT

   END SUBROUTINE ccy_lm

! Calculation of derivatives of Spherical Harmonics

! **************************************************************************************************
!> \brief ...
!> \param c ...
!> \param dy ...
!> \param l ...
!> \param m ...
! **************************************************************************************************
   SUBROUTINE dry_lm(c, dy, l, m)
!
! Real Spherical Harmonics
!                   _                   _
!                  |  [(2l+1)(l-|m|)!]   |1/2 m         cos(m p)   m>=0
!  Y_lm ( t, p ) = |---------------------|   P_l(cos(t))
!                  |[2Pi(1+d_m0)(l+|m|)!]|              sin(|m| p) m<0
!
! Input: c == (t,p)
! Output: dy == (dy/dt, dy/dp)
!
! x == sin(t)*cos(p)
! y == sin(t)*sin(p)
! z == cos(t)
!

      REAL(KIND=dp), DIMENSION(2), INTENT(IN)            :: c
      REAL(KIND=dp), DIMENSION(2), INTENT(OUT)           :: dy
      INTEGER, INTENT(IN)                                :: l, m

      REAL(KIND=dp)                                      :: cp, ct, dplm, lmm, lpm, p, pf, rxy, sp, &
                                                            st, t, tt, y, z
      REAL(KIND=dp), DIMENSION(3)                        :: r

      t = c(1)
      ct = COS(t)
      st = SIN(t)
      p = c(2)
      cp = COS(p)
      sp = SIN(p)
      r(1) = st*cp
      r(2) = st*sp
      r(3) = ct

! dY/dp
      IF (m == 0) THEN
         dy(2) = 0.0_dp
      ELSE
         CALL y_lm(r, y, l, -m)
         dy(2) = -REAL(m, KIND=dp)*y
      END IF

! dY/dt
      SELECT CASE (l)
      CASE (:-1)
         CPABORT("Negative l value")
      CASE (0)
         IF (m /= 0) CPABORT("l = 0 and m value out of bounds")
         dy(1) = 0.0_dp
      CASE (1)
         SELECT CASE (m)
         CASE DEFAULT
            CPABORT("l = 1 and m value out of bounds")
         CASE (1)
            pf = SQRT(3.0_dp/(4.0_dp*pi))
            dy(1) = pf*ct*cp
         CASE (0)
            pf = SQRT(3.0_dp/(4.0_dp*pi))
            dy(1) = -pf*st
         CASE (-1)
            pf = SQRT(3.0_dp/(4.0_dp*pi))
            dy(1) = pf*ct*sp
         END SELECT
      CASE (2)
         SELECT CASE (m)
         CASE DEFAULT
            CPABORT("l = 2 and m value out of bounds")
         CASE (2)
            pf = SQRT(15.0_dp/(16.0_dp*pi))
            dy(1) = pf*2.0_dp*st*ct*COS(2._dp*p)
         CASE (1)
            pf = SQRT(15.0_dp/(4.0_dp*pi))
            dy(1) = pf*cp*(ct*ct - st*st)
         CASE (0)
            pf = SQRT(5.0_dp/(16.0_dp*pi))
            dy(1) = -pf*6.0_dp*ct*st
         CASE (-1)
            pf = SQRT(15.0_dp/(4.0_dp*pi))
            dy(1) = pf*sp*(ct*ct - st*st)
         CASE (-2)
            pf = SQRT(15.0_dp/(16.0_dp*pi))
            dy(1) = pf*2.0_dp*st*ct*SIN(2._dp*p)
         END SELECT
      CASE (3)
         SELECT CASE (m)
         CASE DEFAULT
            CPABORT("l = 3 and m value out of bounds")
         CASE (3)
            pf = SQRT(35.0_dp/(32.0_dp*pi))
            dy(1) = pf*3.0_dp*COS(3._dp*p)*ct*st*st
         CASE (2)
            pf = SQRT(105.0_dp/(16.0_dp*pi))
            dy(1) = pf*2.0_dp*COS(2._dp*p)*ct*st
         CASE (1)
            pf = SQRT(21.0_dp/(32.0_dp*pi))
            dy(1) = pf*cp*(ct*(5.0_dp*ct - 1.0_dp) - 5.0_dp*st*st)
         CASE (0)
            pf = SQRT(7.0_dp/(16.0_dp*pi))
            dy(1) = pf*r(3)*(3.0_dp - 15.0_dp*ct*ct)*st
         CASE (-1)
            pf = SQRT(21.0_dp/(32.0_dp*pi))
            dy(1) = pf*sp*(ct*(5.0_dp*ct - 1.0_dp) - 5.0_dp*st*st)
         CASE (-2)
            pf = SQRT(105.0_dp/(16.0_dp*pi))
            dy(1) = pf*2.0_dp*SIN(2._dp*p)*ct*st
         CASE (-3)
            pf = SQRT(35.0_dp/(32.0_dp*pi))
            dy(1) = pf*3.0_dp*SIN(3._dp*p)*ct*st*st
         END SELECT
      CASE DEFAULT
         IF (m < -l .OR. m > l) CPABORT("m value out of bounds")
         lpm = fac(l + ABS(m))
         lmm = fac(l - ABS(m))
         IF (m == 0) THEN
            tt = 4.0_dp*pi
         ELSE
            tt = 2.0_dp*pi
         END IF
         IF (ABS(lpm) < EPSILON(1.0_dp)) THEN
            pf = REAL(2*l + 1, KIND=dp)/tt
         ELSE
            pf = (REAL(2*l + 1, KIND=dp)*lmm)/(tt*lpm)
         END IF
         pf = SQRT(pf)
         z = ct
         dplm = dlegendre(z, l, m)
         IF (m == 0) THEN
            y = pf*dplm
         ELSE
            rxy = SQRT(r(1)**2 + r(2)**2)
            IF (rxy < EPSILON(1.0_dp)) THEN
               y = 0.0_dp
            ELSE
               IF (m > 0) THEN
                  y = pf*dplm*cosn(m, cp, sp)
               ELSE
                  y = pf*dplm*sinn(-m, cp, sp)
               END IF
            END IF
         END IF
      END SELECT

   END SUBROUTINE dry_lm

! **************************************************************************************************
!> \brief ...
!> \param c ...
!> \param dy ...
!> \param l ...
!> \param m ...
! **************************************************************************************************
   SUBROUTINE dcy_lm(c, dy, l, m)
!
! Complex Spherical Harmonics
!                   _                _
!                  | [(2l+1)(l-|m|)!] |1/2 m
!  Y_lm ( t, p ) = |------------------|   P_l(cos(t)) Exp[ i m p ]
!                  |  [4Pi(l+|m|)!]|  |
!
! Input: r == (x,y,z) : normalised    x^2 + y^2 + z^2 = 1
!
! Input: c == (t,p)
! Output: dy == (dy/dt, dy/dp)
!
! x == sin(t)*cos(p)
! y == sin(t)*sin(p)
! z == cos(t)
!

      REAL(KIND=dp), DIMENSION(2), INTENT(IN)            :: c
      COMPLEX(KIND=dp), DIMENSION(2), INTENT(OUT)        :: dy
      INTEGER, INTENT(IN)                                :: l, m

      REAL(KIND=dp), DIMENSION(2)                        :: dd

      CPABORT("Not implemented")
      CALL dry_lm(c, dd, l, m)
      dy(1) = CMPLX(dd(1), 0.0_dp, KIND=dp)
      dy(2) = CMPLX(dd(2), 0.0_dp, KIND=dp)

   END SUBROUTINE dcy_lm

! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \param l ...
!> \param m ...
!> \return ...
! **************************************************************************************************
   FUNCTION legendre(x, l, m) RESULT(plm)

      REAL(KIND=dp), INTENT(IN)                          :: x
      INTEGER, INTENT(IN)                                :: l, m
      REAL(KIND=dp)                                      :: plm

      INTEGER                                            :: il, im, mm
      REAL(KIND=dp)                                      :: fact, pll, pmm, pmmp1, somx2

      IF (ABS(x) > 1.0_dp) CPABORT("x value > 1")
      SELECT CASE (l)
      CASE (:-1)
         CPABORT("Negative l value")
      CASE (0)
         plm = 1.0_dp
      CASE (1)
         SELECT CASE (ABS(m))
         CASE DEFAULT
            CPABORT("l = 1 and m value out of bounds")
         CASE (1)
            plm = SQRT(1.0_dp - x*x)
         CASE (0)
            plm = x
         END SELECT
      CASE (2)
         SELECT CASE (ABS(m))
         CASE DEFAULT
            CPABORT("l = 2 and m value out of bounds")
         CASE (2)
            plm = 3.0_dp*(1.0_dp - x*x)
         CASE (1)
            plm = 3.0_dp*x*SQRT(1.0_dp - x*x)
         CASE (0)
            plm = 1.5_dp*x*x - 0.5_dp
         END SELECT
      CASE (3)
         SELECT CASE (ABS(m))
         CASE DEFAULT
            CPABORT("l = 3 and m value out of bounds")
         CASE (3)
            plm = 15.0_dp*(1.0_dp - x*x)**1.5_dp
         CASE (2)
            plm = 15.0_dp*x*(1.0_dp - x*x)
         CASE (1)
            plm = (7.5_dp*x*x - 1.5_dp)*SQRT(1.0_dp - x*x)
         CASE (0)
            plm = 2.5_dp*x**3 - 1.5_dp*x
         END SELECT
      CASE (4)
         SELECT CASE (ABS(m))
         CASE DEFAULT
            CPABORT("l = 4 and m value out of bounds")
         CASE (4)
            plm = 105.0_dp*(1.0_dp - x*x)**2
         CASE (3)
            plm = 105.0_dp*x*(1.0_dp - x*x)**1.5_dp
         CASE (2)
            plm = (52.5_dp*x*x - 7.5_dp)*(1.0_dp - x*x)
         CASE (1)
            plm = (17.5_dp*x**3 - 7.5_dp*x)*SQRT(1.0_dp - x*x)
         CASE (0)
            plm = 4.375_dp*x**4 - 3.75_dp*x**2 + 0.375_dp
         END SELECT
      CASE (5)
         SELECT CASE (ABS(m))
         CASE DEFAULT
            CPABORT("l = 5 and m value out of bounds")
         CASE (5)
            plm = 945.0_dp*(1.0_dp - x*x)**2.5_dp
         CASE (4)
            plm = 945.0_dp*x*(1.0_dp - x*x)**2
         CASE (3)
            plm = -(-472.5_dp*x*x + 52.5_dp)*(1.0_dp - x*x)**1.5_dp
         CASE (2)
            plm = (157.5_dp*x**3 - 52.5_dp*x)*(1.0_dp - x*x)
         CASE (1)
            plm = -(-39.375_dp*x**4 + 26.25_dp*x*x - &
                    1.875_dp)*SQRT(1.0_dp - x*x)
         CASE (0)
            plm = 7.875_dp*x**5 - 8.75_dp*x**3 + 1.875_dp*x
         END SELECT
      CASE (6)
         SELECT CASE (ABS(m))
         CASE DEFAULT
            CPABORT("l = 6 and m value out of bounds")
         CASE (6)
            plm = 10395.0_dp*(1.0_dp - x*x)**3
         CASE (5)
            plm = 10395.0_dp*x*(1.0_dp - x*x)**2.5_dp
         CASE (4)
            plm = (5197.5_dp*x*x - 472.5_dp)*(1.0_dp - x*x)**2
         CASE (3)
            plm = -(-1732.5_dp*x**3 + 472.5_dp*x)* &
                  (1.0_dp - x*x)**1.5_dp
         CASE (2)
            plm = (433.125_dp*x**4 - 236.25_dp*x**2 + &
                   13.125_dp)*(1.0_dp - x*x)
         CASE (1)
            plm = -(-86.625_dp*x**5 + 78.75_dp*x**3 - &
                    13.125_dp*x)*SQRT(1.0_dp - x*x)
         CASE (0)
            plm = 14.4375_dp*x**6 - 19.6875_dp*x**4 + &
                  6.5625_dp*x**2 - 0.3125_dp
         END SELECT
      CASE DEFAULT
         mm = ABS(m)
         IF (mm > l) CPABORT("m out of bounds")
! use recurence from numerical recipies
         pmm = 1.0_dp
         IF (mm > 0) THEN
            somx2 = SQRT((1.0_dp - x)*(1.0_dp + x))
            fact = 1.0_dp
            DO im = 1, mm
               pmm = pmm*fact*somx2
               fact = fact + 2.0_dp
            END DO
         END IF
         IF (l == mm) THEN
            plm = pmm
         ELSE
            pmmp1 = x*REAL(2*mm + 1, KIND=dp)*pmm
            IF (l == mm + 1) THEN
               plm = pmmp1
            ELSE
               DO il = mm + 2, l
                  pll = (x*REAL(2*il - 1, KIND=dp)*pmmp1 - &
                         REAL(il + mm - 1, KIND=dp)*pmm)/REAL(il - mm, KIND=dp)
                  pmm = pmmp1
                  pmmp1 = pll
               END DO
               plm = pll
            END IF
         END IF
      END SELECT

   END FUNCTION legendre

! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \param l ...
!> \param m ...
!> \return ...
! **************************************************************************************************
   FUNCTION dlegendre(x, l, m) RESULT(dplm)
      REAL(KIND=dp), INTENT(IN)                          :: x
      INTEGER, INTENT(IN)                                :: l, m
      REAL(KIND=dp)                                      :: dplm

      INTEGER                                            :: mm

      IF (ABS(x) > 1.0_dp) CPABORT("x value > 1")
      SELECT CASE (l)
      CASE (0)
         dplm = 0.0_dp
      CASE (1)
         SELECT CASE (ABS(m))
         CASE DEFAULT
            CPABORT("l = 1 and m value out of bounds")
         CASE (1)
            dplm = -x/SQRT(1.0_dp - x*x)
         CASE (0)
            dplm = 1.0_dp
         END SELECT
      CASE (2)
         SELECT CASE (ABS(m))
         CASE DEFAULT
            CPABORT("l = 2 and m value out of bounds")
         CASE (2)
            dplm = -6.0_dp*x
         CASE (1)
            dplm = 3.0_dp*SQRT(1.0_dp - x*x) - 3.0_dp*x*x/SQRT(1.0_dp - x*x)
         CASE (0)
            dplm = 3.0_dp*x
         END SELECT
      CASE (3)
         SELECT CASE (ABS(m))
         CASE DEFAULT
            CPABORT("l = 3 and m value out of bounds")
         CASE (3)
            dplm = -45.0_dp*SQRT(1.0_dp - x*x)*x
         CASE (2)
            dplm = 15.0_dp*(1.0_dp - x*x) - 30.0_dp*x*x
         CASE (1)
            dplm = 15.0_dp*x*SQRT(1.0_dp - x*x) - (x*(7.5_dp*x*x - 1.5_dp))/SQRT(1.0_dp - x*x)
         CASE (0)
            dplm = 7.5_dp*x*x - 1.5_dp
         END SELECT
      CASE (4)
         SELECT CASE (ABS(m))
         CASE DEFAULT
            CPABORT("l = 4 and m value out of bounds")
         CASE (4)
            dplm = -420*x*(1 - x*x)
         CASE (3)
            dplm = 105.0_dp*((1.0_dp - x*x)**1.5_dp - 3.0_dp*x*x*(1.0_dp - x*x)**0.5_dp)
         CASE (2)
            dplm = 105.0_dp*x*(1.0_dp - x*x) - 2.0_dp*x*(52.5_dp*x*x - 7.5_dp)
         CASE (1)
            IF (ABS(x) - 1.0_dp < EPSILON(1.0_dp)) THEN
               dplm = 0.0_dp
            ELSE
               dplm = (17.5_dp*3.0_dp*x**2 - 7.5_dp)*SQRT(1.0_dp - x*x) - &
                      x*(17.5_dp*x**3 - 7.5_dp*x)/SQRT(1.0_dp - x*x)
            END IF
         CASE (0)
            dplm = 4.375_dp*4.0_dp*x**3 - 3.75_dp*2.0_dp*x
         END SELECT
      CASE (5)
         SELECT CASE (ABS(m))
         CASE DEFAULT
            CPABORT("l = 5 and m value out of bounds")
         CASE (5)
            dplm = -945.0_dp*5.0_dp*x*(1.0_dp - x*x)**1.5_dp
         CASE (4)
            dplm = 945.0_dp*((1.0_dp - x*x)**2 - 2.0_dp*x*x*(1.0_dp - x*x))
         CASE (3)
            dplm = 945.0_dp*x*(1.0_dp - x*x)**1.5_dp - &
                   3.0_dp*x*(472.5_dp*x*x - 52.5_dp)*(1.0_dp - x*x)**0.5_dp
         CASE (2)
            dplm = (3.0_dp*157.5_dp*x**2 - 52.5_dp)*(1.0_dp - x*x) - &
                   (157.5_dp*x**3 - 52.5_dp*x)*(-2.0_dp*x)
         CASE (1)
            IF (ABS(x) - 1.0_dp < EPSILON(1.0_dp)) THEN
               dplm = 0.0_dp
            ELSE
               dplm = -(-39.375_dp*4.0_dp*x*x*x + 2.0_dp*26.25_dp*x)*SQRT(1.0_dp - x*x) + &
                      x*(-39.375_dp*x**4 + 26.25_dp*x*x - 1.875_dp)/SQRT(1.0_dp - x*x)
            END IF
         CASE (0)
            dplm = 5.0_dp*7.875_dp*x**4 - 3.0_dp*8.75_dp*x**2 + 1.875_dp
         END SELECT
      CASE (6)
         SELECT CASE (ABS(m))
         CASE DEFAULT
            CPABORT("l = 6 and m value out of bounds")
         CASE (6)
            dplm = -10395.0_dp*6.0_dp*x*(1.0_dp - x*x)**2
         CASE (5)
            dplm = 10395.0_dp*((1.0_dp - x*x)**2.5_dp - 5.0_dp*x*x*(1.0_dp - x*x)**1.5_dp)
         CASE (4)
            dplm = 2.0_dp*5197.5_dp*x*(1.0_dp - x*x)**2 - &
                   4.0_dp*x*(5197.5_dp*x*x - 472.5_dp)*(1.0_dp - x*x)
         CASE (3)
            dplm = -(-3.0_dp*1732.5_dp*x*x + 472.5_dp)*(1.0_dp - x*x)**1.5_dp + &
                   (-1732.5_dp*x**3 + 472.5_dp*x)*3.0_dp*x*(1.0_dp - x*x)**0.5_dp
         CASE (2)
            dplm = (433.125_dp*4.0_dp*x**3 - 2.0_dp*236.25_dp*x)*(1.0_dp - x*x) - &
                   2.0_dp*x*(433.125_dp*x**4 - 236.25_dp*x**2 + 13.125_dp)
         CASE (1)
            IF (ABS(x) - 1.0_dp < EPSILON(1.0_dp)) THEN
               dplm = 0.0_dp
            ELSE
               dplm = -(-5.0_dp*86.625_dp*x**4 + 3.0_dp*78.75_dp**2 - 13.125_dp)*SQRT(1.0_dp - x*x) + &
                      x*(-86.625_dp*x**5 + 78.75_dp*x**3 - 13.125_dp*x)/SQRT(1.0_dp - x*x)
            END IF
         CASE (0)
            dplm = 14.4375_dp*6.0_dp*x**5 - 19.6875_dp*4.0_dp*x**3 + &
                   6.5625_dp*2.0_dp*x
         END SELECT
      CASE DEFAULT
         mm = ABS(m)
         IF (mm > l) CPABORT("m out of bounds")

         !From Wikipedia: dPlm(x) = -(l+1)x/(x^2-1)*Plm(x) + (l-m+1)/(x^2-1)Pl+1m(x)
         IF (ABS(x) - 1.0_dp < EPSILON(1.0_dp)) THEN
            dplm = 0.0_dp
         ELSE
            dplm = -REAL(l + 1, dp)*x/(x**2 - 1.0_dp)*legendre(x, l, m) &
                   + REAL(l - m + 1, dp)/(x**2 - 1.0_dp)*legendre(x, l + 1, m)
         END IF
      END SELECT

   END FUNCTION dlegendre

! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \param l ...
!> \return ...
! **************************************************************************************************
   FUNCTION dPof1(x, l)

      REAL(KIND=dp), INTENT(IN)                          :: x
      INTEGER, INTENT(IN)                                :: l
      REAL(KIND=dp)                                      :: dPof1

      IF (ABS(x) - 1.0_dp > EPSILON(1.0_dp)) THEN
         CPABORT("|x| is not 1")
      END IF
      IF (x > 0.0_dp) THEN
         SELECT CASE (l)
         CASE (:-1)
            CPABORT("Negative l value")
         CASE (0)
            dPof1 = 0.0_dp
         CASE (1)
            dPof1 = 1.0_dp
         CASE (2)
            dPof1 = 3.0_dp
         CASE (3)
            dPof1 = 6.0_dp
         CASE (4)
            dPof1 = 10.0_dp
         CASE (5)
            dPof1 = 15.0_dp
         CASE (6)
            dPof1 = 21.0_dp
         CASE (7:)
            CPABORT("Not implemented")
         END SELECT
      ELSE
         SELECT CASE (l)
         CASE (:-1)
            CPABORT("Negative l value")
         CASE (0)
            dPof1 = 0.0_dp
         CASE (1)
            dPof1 = 1.0_dp
         CASE (2)
            dPof1 = -3.0_dp
         CASE (3)
            dPof1 = 6.0_dp
         CASE (4)
            dPof1 = -10.0_dp
         CASE (5)
            dPof1 = 15.0_dp
         CASE (6)
            dPof1 = -21.0_dp
         CASE (7:)
            CPABORT("Not implemented")
         END SELECT
      END IF

   END FUNCTION dPof1

! **************************************************************************************************
!> \brief ...
!> \param n ...
!> \param k ...
!> \return ...
! **************************************************************************************************
   FUNCTION choose(n, k)

      INTEGER, INTENT(IN)                                :: n, k
      REAL(KIND=dp)                                      :: choose

      IF (n >= k) THEN
         choose = REAL(NINT(fac(n)/(fac(k)*fac(n - k))), KIND=dp)
      ELSE
         choose = 0.0_dp
      END IF

   END FUNCTION choose

! **************************************************************************************************
!> \brief ...
!> \param n ...
!> \param c ...
!> \param s ...
!> \return ...
! **************************************************************************************************
   FUNCTION cosn(n, c, s)

      INTEGER, INTENT(IN)                                :: n
      REAL(KIND=dp), INTENT(IN)                          :: c, s
      REAL(KIND=dp)                                      :: cosn

      INTEGER                                            :: i, j

      cosn = 0.0_dp
      IF (ABS(c) < EPSILON(1.0_dp) .OR. n == 0) THEN
         IF (MOD(n, 2) == 0) THEN
            cosn = (-1.0_dp)**INT(n/2)
         ELSE
            cosn = 0.0_dp
         END IF
      ELSE
         DO i = n, 0, -2
            IF (i == 0) THEN
               j = n
               IF (j == 0) THEN
                  cosn = cosn + choose(n, j)
               ELSE IF (MOD(j/2, 2) == 0) THEN
                  cosn = cosn + choose(n, j)*s**j
               ELSE
                  cosn = cosn - choose(n, j)*s**j
               END IF
            ELSE
               j = n - i
               IF (j == 0) THEN
                  cosn = cosn + choose(n, j)*c**i
               ELSE IF (MOD(j/2, 2) == 0) THEN
                  cosn = cosn + choose(n, j)*c**i*s**j
               ELSE
                  cosn = cosn - choose(n, j)*c**i*s**j
               END IF
            END IF
         END DO
      END IF

   END FUNCTION cosn

! **************************************************************************************************
!> \brief ...
!> \param n ...
!> \param c ...
!> \param s ...
!> \return ...
! **************************************************************************************************
   FUNCTION dcosn_dcp(n, c, s)

      INTEGER, INTENT(IN)                                :: n
      REAL(KIND=dp), INTENT(IN)                          :: c, s
      REAL(KIND=dp)                                      :: dcosn_dcp

      INTEGER                                            :: i, j

      dcosn_dcp = 0.0_dp

      IF (s < EPSILON(1.0_dp)) THEN
         dcosn_dcp = 0.0_dp
      ELSE
         DO i = n, 0, -2
            IF (i == 0) THEN
               dcosn_dcp = dcosn_dcp
            ELSE
               j = n - i
               IF (MOD(j/2, 2) == 0) THEN
                  dcosn_dcp = dcosn_dcp + choose(n, j)*REAL(i, dp)*c**(i - 1)*s**j
               ELSE
                  dcosn_dcp = dcosn_dcp - choose(n, j)*REAL(i, dp)*c**(i - 1)*s**j
               END IF
            END IF
         END DO
      END IF

   END FUNCTION dcosn_dcp

! **************************************************************************************************
!> \brief ...
!> \param n ...
!> \param c ...
!> \param s ...
!> \return ...
! **************************************************************************************************
   FUNCTION dcosn_dsp(n, c, s)

      INTEGER, INTENT(IN)                                :: n
      REAL(KIND=dp), INTENT(IN)                          :: c, s
      REAL(KIND=dp)                                      :: dcosn_dsp

      INTEGER                                            :: i, j

      dcosn_dsp = 0.0_dp
      IF (c < EPSILON(1.0_dp) .OR. s < EPSILON(1.0_dp)) THEN
         dcosn_dsp = 0.0_dp
      ELSE
         DO i = n, 0, -2
            j = n - i
            IF (j == 0) THEN
               dcosn_dsp = dcosn_dsp
            ELSE
               IF (MOD(j/2, 2) == 0) THEN
                  dcosn_dsp = dcosn_dsp + choose(n, j)*REAL(j, dp)*s**(j - 1)*c**i
               ELSE
                  dcosn_dsp = dcosn_dsp - choose(n, j)*REAL(j, dp)*s**(j - 1)*c**i
               END IF
            END IF
         END DO
      END IF

   END FUNCTION dcosn_dsp

! **************************************************************************************************
!> \brief ...
!> \param n ...
!> \param c ...
!> \param s ...
!> \return ...
! **************************************************************************************************
   FUNCTION sinn(n, c, s)

      INTEGER, INTENT(IN)                                :: n
      REAL(KIND=dp), INTENT(IN)                          :: c, s
      REAL(KIND=dp)                                      :: sinn

      INTEGER                                            :: i, j

      sinn = 0.0_dp

      IF (ABS(s) < EPSILON(1.0_dp) .OR. n == 0) THEN
         sinn = 0.0_dp
      ELSE IF (ABS(c) < EPSILON(1.0_dp)) THEN
         IF (MOD(n, 2) == 0) THEN
            sinn = 0.0_dp
         ELSE
            sinn = s*(-1.0_dp)**INT((n - 1)/2)
         END IF
      ELSE
         DO i = n - 1, 0, -2
            IF (i == 0) THEN
               j = n
               IF (MOD(j/2, 2) == 0) THEN
                  sinn = sinn + choose(n, j)*s**j
               ELSE
                  sinn = sinn - choose(n, j)*s**j
               END IF
            ELSE
               j = n - i
               IF (MOD(j/2, 2) == 0) THEN
                  sinn = sinn + choose(n, j)*c**i*s**j
               ELSE
                  sinn = sinn - choose(n, j)*c**i*s**j
               END IF
            END IF
         END DO
      END IF

   END FUNCTION sinn

! **************************************************************************************************
!> \brief ...
!> \param n ...
!> \param c ...
!> \param s ...
!> \return ...
! **************************************************************************************************
   FUNCTION dsinn_dcp(n, c, s)

      INTEGER, INTENT(IN)                                :: n
      REAL(KIND=dp), INTENT(IN)                          :: c, s
      REAL(KIND=dp)                                      :: dsinn_dcp

      INTEGER                                            :: i, j

      dsinn_dcp = 0.0_dp

      IF (c < EPSILON(1.0_dp) .OR. s < EPSILON(1.0_dp)) THEN
         dsinn_dcp = 0.0_dp
      ELSE
         DO i = n - 1, 0, -2
            IF (i == 0) THEN
               dsinn_dcp = dsinn_dcp
            ELSE
               j = n - i
               IF (MOD(j/2, 2) == 0) THEN
                  dsinn_dcp = dsinn_dcp + choose(n, j)*REAL(i, dp)*c**(i - 1)*s**j
               ELSE
                  dsinn_dcp = dsinn_dcp - choose(n, j)*REAL(i, dp)*c**(i - 1)*s**j
               END IF
            END IF
         END DO
      END IF

   END FUNCTION dsinn_dcp

! **************************************************************************************************
!> \brief ...
!> \param n ...
!> \param c ...
!> \param s ...
!> \return ...
! **************************************************************************************************
   FUNCTION dsinn_dsp(n, c, s)

      INTEGER, INTENT(IN)                                :: n
      REAL(KIND=dp), INTENT(IN)                          :: c, s
      REAL(KIND=dp)                                      :: dsinn_dsp

      INTEGER                                            :: i, j

      dsinn_dsp = 0.0_dp

      IF (c < EPSILON(1.0_dp) .OR. s < EPSILON(1.0_dp)) THEN
         dsinn_dsp = 0.0_dp
      ELSE
         DO i = n - 1, 0, -2
            j = n - i
            IF (j == 0) THEN
               dsinn_dsp = dsinn_dsp
            ELSE
               IF (MOD(j/2, 2) == 0) THEN
                  dsinn_dsp = dsinn_dsp + choose(n, j)*c**i*REAL(j, dp)*s**(j - 1)
               ELSE
                  dsinn_dsp = dsinn_dsp - choose(n, j)*c**i*REAL(j, dp)*s**(j - 1)
               END IF
            END IF
         END DO
      END IF

   END FUNCTION dsinn_dsp

END MODULE spherical_harmonics
