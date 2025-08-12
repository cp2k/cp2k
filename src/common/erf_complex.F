!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! Copyright (c) 2016 Anton Shterenlikht, The University of Bristol, UK
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are
! met:
!
! 1. Redistributions of source code must retain the above copyright
! notice, this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
! TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! Module with error function and related functions.
! Two algorithms are implemented: Poppe and Wijers (faddeyeva_fast, erfz_fast)
! and Zaghloul and Ali (faddeyeva_accurate, erfz_accurate). The second algorithm is
! supposed to be more accurate for some values. However, the first
! algorithm can be over an order of magnitude faster.
!
! This code was written before the author became aware of
! M. R. Zaghloul, Remark on "Algorithm 916: Computing the
! Faddeyeva and Voight functions": efficiency improvements and
! Fortran translation, ACM Trans. Math. Software 42, Article 26, 2016.

! **************************************************************************************************
!> \brief Module to compute the error function of a complex argument
!> \par History
!>      08.2025 Adapted to use CP2K intrinsics and constants
!> \author Stefano Battaglia
! **************************************************************************************************
MODULE erf_complex

   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi

   IMPLICIT NONE

   REAL(kind=dp), PARAMETER :: zero = 0.0_dp, one = 1.0_dp, &
                               two = 2.0_dp, half = 0.5_dp, rmin = TINY(one), &
                               eps0 = EPSILON(one), sqrt_log_rmin = SQRT(-LOG(rmin)), &
                               pi2 = pi*pi, one_sqrt_pi = one/SQRT(pi)
   COMPLEX(kind=dp), PARAMETER :: &
      cmplxj = CMPLX(zero, one, kind=dp), &
      cmplx0 = CMPLX(zero, zero, kind=dp)

   PRIVATE
   PUBLIC :: faddeyeva_fast, faddeyeva_accurate, erfz_fast, erfz_accurate

CONTAINS

! **************************************************************************************************
!> \brief Computes the Faddeyeva function w(z) = exp(-z**2) * erfc(-i*z)
!> \param z complex argument
!> \param err desired accuracy (positive)
!> \return Faddeyeva function w(z)
! **************************************************************************************************
   elemental COMPLEX(kind=dp) FUNCTION faddeyeva_accurate(z, err)

      ! This is the Faddeyeva or the plasma dispersion function,
      ! w(z) = exp(-z**2) * erfc(-i*z). erfc(z) is the complex complementary
      ! error function of z. z is a complex number.
      !
      ! Adapted from the Matlab code implementing TOMS
      ! Algorithm 916: http://www.netlib.org/toms/
      !
      ! file: 916.zip
      ! ref: TOMS 38,2 (Dec 2011) Article: 15
      ! for: Computing the Faddeyeva and Voigt Functions
      ! by: Mofreh R. Zaghloul and Ahmed N. Ali
      !
      ! Most of the code is calculation of equations (13)-(19) of the
      ! above paper.
      !
      ! Inputs:
      !    z - the argument of the function
      !  err - The desired accuracy (positive). err must be (err .le. 1.0e-4)
      !        and  (err .ge. 0.06447*epsilon). For efficiency,
      !        no checks are made!
      !        The lowest accuracy and the fastest calculation are obtained
      !        with err .eq. 1.0e-4. For higher accuracy use smaller err.

      COMPLEX(kind=dp), INTENT(in)                       :: z
      REAL(kind=dp), INTENT(in)                          :: err

      INTEGER                                            :: n, n3, n3_3
      REAL(kind=dp) :: a, a_pi, a_sqr, aux13, cos_2yx, del2_tmp, del3, del3_3_tmp, del3_tmp, del5, &
         delta3, delta5, den1, erfcsy, exp1, exp2, exp3, exp3_3_den, exp3_den, exp_del1, &
         exp_x_sqr, four_a_sqr, half_a, l_old, myerr, sigma1, sigma2, sigma3, sigma4, sigma4_5, &
         sigma5, sin_2yx, two_a, two_a_pi, two_a_sqr, two_a_x, two_exp_x_sqr_ysqr, two_yx, v_old, &
         x, x_sqr, xsign, y, y_sqr, ysign

      x = REAL(z)
      y = AIMAG(z)

      ! For purely imaginary z, use intrinsic scaled complement of
      ! the error function, erfc_scaled (F2008 and beyond).
      ! Return immediately.
      IF (ABS(x) .EQ. zero) THEN
         faddeyeva_accurate = erfc_scaled(y)
         RETURN
      END IF

      myerr = MAX(err, eps0)
      a = SQRT(-pi2/LOG(err/2.0_dp))
      half_a = half*a
      a_sqr = a**2
      two_a = 2*a
      two_a_sqr = 2*a_sqr
      four_a_sqr = 4*a_sqr
      a_pi = a/pi
      two_a_pi = 2*a_pi
      erfcsy = erfc_scaled(ABS(y))
      xsign = SIGN(one, x)
      ysign = SIGN(one, y)
      x = ABS(x)
      y = MAX(rmin, ABS(y))
      x_sqr = x**2
      y_sqr = y**2
      two_yx = 2*y*x
      two_a_x = two_a*x
      exp_x_sqr = EXP(-x_sqr)
      cos_2yx = COS(two_yx)
      sin_2yx = SIN(two_yx)
      v_old = exp_x_sqr* &
              (erfcsy*cos_2yx + two_a_pi*SIN(two_yx/2)**2/y)
      l_old = -erfcsy + a_pi/y
      sigma3 = rmin
      sigma5 = rmin
      sigma4_5 = zero
      delta3 = one
      delta5 = one
      n = 0
      n3 = CEILING(x/a)
      n3_3 = n3 - 1

      outer: IF ((sqrt_log_rmin - x) .GT. 0) THEN
         sigma1 = rmin
         sigma2 = rmin
         sigma4 = rmin
         exp1 = EXP(-two_a_x)
         exp2 = EXP(four_a_sqr*n3 - 2*two_a_x - two_a_sqr)
         exp3 = EXP(-(two_a_sqr*n3 - two_a_x - two_a_sqr))
         del2_tmp = one
         del3_tmp = EXP(-(a_sqr*n3**2 - two_a_x*n3 &
                          - two_a_sqr*n3 + x_sqr + two_a_x + a_sqr))
         del3_3_tmp = EXP(a_sqr - (two_a_sqr*n3 - two_a_x))

         loop1: DO
            IF (delta3 .LT. myerr .AND. delta5 .LT. myerr .AND. &
                n .GT. 50) EXIT loop1
            n = n + 1
            den1 = a_sqr*n**2 + y_sqr
            exp_del1 = EXP(-(a_sqr*n**2))/den1
            del2_tmp = del2_tmp*exp1

            minor: IF (n3_3 .GE. 1) THEN
               del3_tmp = del3_tmp*exp3
               exp3_den = del3_tmp*exp_del1* &
                          (den1/(a_sqr*n3**2 + y_sqr))
               del3_3_tmp = del3_3_tmp*exp2
               exp3_3_den = exp3_den*del3_3_tmp* &
                            ((a_sqr*n3**2 + y_sqr)/ &
                             (a_sqr*n3_3**2 + y_sqr))
               del5 = n3_3*exp3_3_den + n3*exp3_den
               del3 = exp3_3_den + exp3_den
            ELSE
               del3_tmp = del3_tmp*exp3
               del3 = del3_tmp*exp_del1* &
                      (den1/(a_sqr*n3**2 + y_sqr))
               del5 = n3*del3
            END IF minor

            delta3 = del3/sigma3
            delta5 = del5/sigma5
            sigma1 = sigma1 + exp_del1
            sigma2 = sigma2 + del2_tmp*exp_x_sqr*exp_del1
            sigma3 = sigma3 + del3
            sigma4 = sigma4 + n*del2_tmp*exp_x_sqr*exp_del1
            sigma5 = sigma5 + del5

            IF (x .GE. 5.0e-4_dp) THEN
               sigma4_5 = -sigma4 + sigma5
            ELSE
               sigma4_5 = sigma4_5 + 2*n**2*two_a_x*exp_x_sqr &
                          *exp_del1*(one + 1.666666666666667e-1_dp &
                                     *(two_a_x*n)**2 + 8.333333333333333e-3_dp &
                                     *(two_a_x*n)**4)
            END IF

            n3 = n3 + 1
            n3_3 = n3_3 - 1

         END DO loop1

         ! Second line of Eqn (13)
         aux13 = y*two_a_pi* &
                 (-cos_2yx*exp_x_sqr*sigma1 + half*(sigma2 + sigma3))
         mumu: IF (y .LE. 5.0_dp .AND. two_yx .GT. rmin) THEN
            faddeyeva_accurate = v_old + aux13 + cmplxj*xsign &
                                 *(sin_2yx*exp_x_sqr*(l_old + two_a_pi*y*sigma1) &
                                   + two_a_pi*half_a*sigma4_5)
         ELSE IF (y .LE. 5.0_dp .AND. two_yx .LE. rmin) THEN
            faddeyeva_accurate = v_old + aux13 + cmplxj*xsign &
                                 *(two*y*exp_x_sqr*(x*l_old + x*two_a_pi*y &
                                                    *sigma1) + two_a_pi*half_a*sigma4_5)
         ELSE
            faddeyeva_accurate = v_old + aux13 + cmplxj*xsign &
                                 *(sin_2yx*exp_x_sqr*MIN(zero, ABS(l_old &
                                                                   + (two_a_pi*y*sigma1))) + two_a_pi*half_a &
                                   *sigma4_5)
         END IF mumu

      ELSE IF (x .GE. sqrt_log_rmin .AND. x .LT. 1.0e15_dp) THEN

         exp2 = EXP(four_a_sqr*n3 - 2*two_a_x - two_a_sqr)
         del3_3_tmp = EXP(a_sqr + two_a_x - two_a_sqr*n3)

         loop2: DO
            IF (delta3 .LT. myerr .AND. delta5 .LT. myerr .AND. &
                n .GT. 50) EXIT loop2
            n = n + 1
            IF (n3_3 .GE. 1) THEN
               exp3_den = EXP(-(a*n3 - x)*(a*n3 - x)) &
                          /(a_sqr*n3**2 + y_sqr)
               del3_3_tmp = del3_3_tmp*exp2
               exp3_3_den = exp3_den*del3_3_tmp*((a_sqr*n3**2 + y_sqr) &
                                                 /(a_sqr*n3_3**2 + y_sqr))
               del5 = n3_3*exp3_3_den + n3*exp3_den
               del3 = exp3_3_den + exp3_den
            ELSE
               del3 = EXP(-(a*n3 - x)**2)/(a_sqr*n3**2 + y_sqr)
               del5 = n3*del3
            END IF

            delta3 = del3/sigma3
            delta5 = del5/sigma5
            sigma3 = sigma3 + del3
            sigma5 = sigma5 + del5
            n3 = n3 + 1
            n3_3 = n3_3 - 1

         END DO loop2

         faddeyeva_accurate = v_old + y*a_pi*sigma3 + cmplxj*xsign*(sin_2yx &
                                                                    *exp_x_sqr*l_old + two_a_pi*half_a*sigma5)

      ELSE
         faddeyeva_accurate = one_sqrt_pi*((y + cmplxj*xsign*x)/(x_sqr + y_sqr))
      END IF outer

      IF (ysign .LT. zero) THEN
         two_exp_x_sqr_ysqr = two*EXP(-x_sqr + y_sqr)
         faddeyeva_accurate = two_exp_x_sqr_ysqr*cos_2yx - REAL(faddeyeva_accurate) - cmplxj &
                              *(-xsign*two_exp_x_sqr_ysqr*sin_2yx - AIMAG(faddeyeva_accurate))
      END IF

   END FUNCTION faddeyeva_accurate

! **************************************************************************************************
!> \brief Computes the error function of a complex argument using the Zaghloul and Ali algorithm
!> \param z complex argument
!> \param err desired accuracy (positive)
!> \return error function of z
! **************************************************************************************************
   elemental COMPLEX(kind=dp) FUNCTION erfz_accurate(z, err)

      ! This is an error function of a complex argument, which uses faddeyeva_accurate(z).

      COMPLEX(kind=dp), INTENT(in)                       :: z
      REAL(kind=dp), INTENT(in)                          :: err

      erfz_accurate = one - faddeyeva_accurate(cmplxj*z, err)*EXP(-z**2)

   END FUNCTION erfz_accurate

! **************************************************************************************************
!> \brief Computes the Faddeyeva function w(z) = exp(-z**2) * erfc(-i*z)
!> \param z complex argument
!> \return Faddeyeva function w(z)
! **************************************************************************************************
   elemental COMPLEX(kind=dp) FUNCTION faddeyeva_fast(z)

      ! A modified version of algorithm 680, rewritten in Fortran 2008.
      ! G.P.M. Poppe, C.M.J. Wijers, More efficient computation of
      ! the complex error-function, ACM Trans. Math. Software 16:38-46, 1990.
      !  and
      ! G.P.M. Poppe, C.M.J. Wijers, Algorithm 680, Evaluation of the
      ! complex error function, ACM Trans. Math. Software 16:47, 1990.
      !
      ! Given a complex number z, this function computes
      ! the value of the Faddeeva-function w(z) = exp(-z**2)*erfc(-i*z),
      ! where erfc is the complex complementary error-function and i
      ! means sqrt(-1).  The accuracy of the algorithm for z in the 1st
      ! and 2nd quadrant is 14 significant digits; in the 3rd and 4th
      ! it is 13 significant digits outside a circular region with radius
      ! 0.126 around a zero of the function.

      COMPLEX(kind=dp), INTENT(in)                       :: z

      REAL(kind=dp), PARAMETER :: factor = 1.12837916709551257388_dp

      INTEGER                                            :: i, j, kapn, n, np1, nu
      LOGICAL                                            :: a, b
      REAL(kind=dp)                                      :: c, daux, h, h2, qlambda, qrho, rx, ry, &
                                                            sx, sy, tx, ty, u, u1, u2, v, v1, v2, &
                                                            w1, x, xabs, xabsq, xaux, xi, xquad, &
                                                            xsum, y, yabs, yi, yquad, ysum

      !  factor is 2/sqrt(pi)

      ! To avoid the complier uninitialised varning
      h2 = zero

      xi = REAL(z)
      yi = AIMAG(z)
      xabs = ABS(xi)
      yabs = ABS(yi)
      x = xabs/6.3_dp
      y = yabs/4.4_dp
      qrho = x**2 + y**2
      xabsq = xabs**2
      xquad = xabsq - yabs**2
      yquad = 2*xabs*yabs

      a = qrho .LT. 0.085264_dp

      branch1: IF (a) THEN

         ! If ( qrho .lt. 0.085264 ) then the Faddeeva-function is evaluated
         !  using a power-series (abramowitz/stegun, equation (7.1.5), p.297)
         !  n is the minimum number of terms needed to obtain the required
         !  accuracy

         qrho = (one - 0.85_dp*y)*SQRT(qrho)
         n = NINT(6.0_dp + 72.0_dp*qrho)
         j = 2*n + 1
         xsum = one/REAL(j, kind=dp)
         ysum = zero

         DO i = n, 1, -1
            j = j - 2
            xaux = (xsum*xquad - ysum*yquad)/REAL(i, kind=dp)
            ysum = (xsum*yquad + ysum*xquad)/REAL(i, kind=dp)
            xsum = xaux + one/REAL(j, kind=dp)
         END DO

         u1 = -factor*(xsum*yabs + ysum*xabs) + one
         v1 = factor*(xsum*xabs - ysum*yabs)
         daux = EXP(-xquad)
         u2 = daux*COS(yquad)
         v2 = -daux*SIN(yquad)
         u = u1*u2 - v1*v2
         v = u1*v2 + v1*u2
      ELSE

         bran2: IF (qrho .GT. one) THEN

            ! If ( qrho .gt. 1) then w(z) is evaluated using the laplace
            ! continued fraction. nu is the minimum number of terms needed
            ! to obtain the required accuracy.

            h = zero
            kapn = 0
            qrho = SQRT(qrho)
            nu = INT(3.0_dp + (1442.0_dp/(26.0_dp*qrho &
                                          + 77.0_dp)))

         ELSE

            ! If ( qrho .ge. 0.085264 .and. qrho .le. one ) then
            ! w(z) is evaluated by a truncated Taylor expansion,
            ! where the Laplace continued fraction is used to calculate
            ! the derivatives of w(z). KAPN is the minimum number of terms
            ! in the Taylor expansion needed to obtain the required accuracy.
            ! NU is the minimum number of terms of the continued fraction
            ! needed to calculate the derivatives with the required accuracy.

            qrho = (one - y)*SQRT(one - qrho)
            h = 1.88_dp*qrho
            h2 = two*h
            kapn = NINT(7.0_dp + 34.0_dp*qrho)
            nu = NINT(16.0_dp + 26.0_dp*qrho)

         END IF bran2

         b = h .GT. zero

         ! To avoid uninitialise compiler warning. qlambda is used
         ! only if (b), so can define to any value otherwise.
         qlambda = zero
         IF (b) qlambda = h2**kapn

         rx = zero
         ry = zero
         sx = zero
         sy = zero

         DO n = nu, 0, -1
            np1 = n + 1
            tx = yabs + h + np1*rx
            ty = xabs - np1*ry
            c = half/(tx**2 + ty**2)
            rx = c*tx
            ry = c*ty
            IF (b .AND. n .LE. kapn) THEN
               tx = qlambda + sx
               sx = rx*tx - ry*sy
               sy = ry*tx + rx*sy
               qlambda = qlambda/h2
            END IF
         END DO

         IF (h .EQ. zero) THEN
            u = factor*rx
            v = factor*ry
         ELSE
            u = factor*sx
            v = factor*sy
         END IF

         IF (yabs .EQ. zero) u = EXP(-xabs**2)

      END IF branch1

      ! Evaluation of w(z) in the other quadrants

      IF (yi .LT. zero) THEN

         IF (a) THEN
            u2 = two*u2
            v2 = two*v2
         ELSE
            xquad = -xquad
            w1 = two*EXP(xquad)
            u2 = w1*COS(yquad)
            v2 = -w1*SIN(yquad)
         END IF

         u = u2 - u
         v = v2 - v
         IF (xi .GT. zero) v = -v
      ELSE
         IF (xi .LT. zero) v = -v
      END IF

      faddeyeva_fast = CMPLX(u, v, kind=dp)

   END FUNCTION faddeyeva_fast

! **************************************************************************************************
!> \brief Computes the error function of a complex argument using the Poppe and Wijers algorithm
!> \param z complex argument
!> \return error function of z
! **************************************************************************************************
   elemental COMPLEX(kind=dp) FUNCTION erfz_fast(z)

      ! This is an error function of a complex argument, which uses faddeyeva_fast(z).

      COMPLEX(kind=dp), INTENT(in)                       :: z

      erfz_fast = one - faddeyeva_fast(cmplxj*z)*EXP(-z**2)

   END FUNCTION erfz_fast

   !*********************************************************************
END MODULE erf_complex
