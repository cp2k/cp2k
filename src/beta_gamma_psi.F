!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
MODULE beta_gamma_psi
   ! not tested in the case where dp would stand for single precision
   ! Routines for the beta function are untested

   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: psi

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param i ...
!> \return ...
! **************************************************************************************************
   FUNCTION ipmpar(i) RESULT(fn_val)
!-----------------------------------------------------------------------

!     IPMPAR PROVIDES THE INTEGER MACHINE CONSTANTS FOR THE COMPUTER
!     THAT IS USED. IT IS ASSUMED THAT THE ARGUMENT I IS AN INTEGER
!     HAVING ONE OF THE VALUES 1-10. IPMPAR(I) HAS THE VALUE ...

!  INTEGERS.

!     ASSUME INTEGERS ARE REPRESENTED IN THE N-DIGIT, BASE-A FORM

!               SIGN ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) )

!               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,N-1.

!     IPMPAR(1) = A, THE BASE (radix).

!     IPMPAR(2) = N, THE NUMBER OF BASE-A DIGITS (digits).

!     IPMPAR(3) = A**N - 1, THE LARGEST MAGNITUDE (huge).

!  FLOATING-POINT NUMBERS.

!     IT IS ASSUMED THAT THE SINGLE AND DOUBLE PRECISION FLOATING
!     POINT ARITHMETICS HAVE THE SAME BASE, SAY B, AND THAT THE
!     NONZERO NUMBERS ARE REPRESENTED IN THE FORM

!               SIGN (B**E) * (X(1)/B + ... + X(M)/B**M)

!               WHERE X(I) = 0,1,...,B-1 FOR I=1,...,M,
!               X(1) .GE. 1, AND EMIN .LE. E .LE. EMAX.

!     IPMPAR(4) = B, THE BASE.

!  SINGLE-PRECISION

!     IPMPAR(5) = M, THE NUMBER OF BASE-B DIGITS.

!     IPMPAR(6) = EMIN, THE SMALLEST EXPONENT E.

!     IPMPAR(7) = EMAX, THE LARGEST EXPONENT E.

!  DOUBLE-PRECISION

!     IPMPAR(8) = M, THE NUMBER OF BASE-B DIGITS.

!     IPMPAR(9) = EMIN, THE SMALLEST EXPONENT E.

!     IPMPAR(10) = EMAX, THE LARGEST EXPONENT E.

!-----------------------------------------------------------------------

      INTEGER, INTENT(IN)                                :: i
      INTEGER                                            :: fn_val

      SELECT CASE (i)
      CASE (1)
         fn_val = RADIX(i)
      CASE (2)
         fn_val = DIGITS(i)
      CASE (3)
         fn_val = HUGE(i)
      CASE (4)
         fn_val = RADIX(1.0)
      CASE (5)
         fn_val = DIGITS(1.0)
      CASE (6)
         fn_val = MINEXPONENT(1.0)
      CASE (7)
         fn_val = MAXEXPONENT(1.0)
      CASE (8)
         fn_val = DIGITS(1.0e0_dp)
      CASE (9)
         fn_val = MINEXPONENT(1.0e0_dp)
      CASE (10)
         fn_val = MAXEXPONENT(1.0e0_dp)
      CASE DEFAULT
         CPABORT("unknown case")
      END SELECT

   END FUNCTION ipmpar

! **************************************************************************************************
!> \brief ...
!> \param i ...
!> \return ...
! **************************************************************************************************
   FUNCTION dpmpar(i) RESULT(fn_val)
!-----------------------------------------------------------------------

!     DPMPAR PROVIDES THE DOUBLE PRECISION MACHINE CONSTANTS FOR
!     THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT
!     I IS AN INTEGER HAVING ONE OF THE VALUES 1, 2, OR 3. IF THE
!     DOUBLE PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND
!     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX, THEN

!        DPMPAR(1) = B**(1 - M), THE MACHINE PRECISION,

!        DPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,

!        DPMPAR(3) = B**EMAX*(1 - B**(-M)), THE LARGEST MAGNITUDE.
!-----------------------------------------------------------------------

      INTEGER, INTENT(IN)                                :: i
      REAL(dp)                                           :: fn_val

      REAL(dp), PARAMETER                                :: one = 1._dp

! Local variable

      SELECT CASE (i)
      CASE (1)
         fn_val = EPSILON(one)
      CASE (2)
         fn_val = TINY(one)
      CASE (3)
         fn_val = HUGE(one)
      END SELECT

   END FUNCTION dpmpar

! **************************************************************************************************
!> \brief ...
!> \param l ...
!> \return ...
! **************************************************************************************************
   FUNCTION dxparg(l) RESULT(fn_val)
!--------------------------------------------------------------------
!     IF L = 0 THEN  DXPARG(L) = THE LARGEST POSITIVE W FOR WHICH
!     DEXP(W) CAN BE COMPUTED.
!
!     IF L IS NONZERO THEN  DXPARG(L) = THE LARGEST NEGATIVE W FOR
!     WHICH THE COMPUTED VALUE OF DEXP(W) IS NONZERO.
!
!     NOTE... ONLY AN APPROXIMATE VALUE FOR DXPARG(L) IS NEEDED.
!--------------------------------------------------------------------
      INTEGER, INTENT(IN)                                :: l
      REAL(dp)                                           :: fn_val

      REAL(dp), PARAMETER                                :: one = 1._dp

! Local variable

      IF (l == 0) THEN
         fn_val = LOG(HUGE(one))
      ELSE
         fn_val = LOG(TINY(one))
      END IF

   END FUNCTION dxparg

! **************************************************************************************************
!> \brief ...
!> \param a ...
!> \return ...
! **************************************************************************************************
   FUNCTION alnrel(a) RESULT(fn_val)
!-----------------------------------------------------------------------
!            EVALUATION OF THE FUNCTION LN(1 + A)
!-----------------------------------------------------------------------
      REAL(dp), INTENT(IN)                               :: a
      REAL(dp)                                           :: fn_val

      REAL(dp), PARAMETER :: half = 0.5e0_dp, one = 1.e0_dp, p1 = -.129418923021993e+01_dp, &
         p2 = .405303492862024e+00_dp, p3 = -.178874546012214e-01_dp, &
         q1 = -.162752256355323e+01_dp, q2 = .747811014037616e+00_dp, &
         q3 = -.845104217945565e-01_dp, two = 2.e0_dp, zero = 0.e0_dp

      REAL(dp)                                           :: t, t2, w, x

!--------------------------

      IF (ABS(a) <= 0.375e0_dp) THEN
         t = a/(a + two)
         t2 = t*t
         w = (((p3*t2 + p2)*t2 + p1)*t2 + one)/(((q3*t2 + q2)*t2 + q1)*t2 + one)
         fn_val = two*t*w
      ELSE
         x = one + a
         IF (a < zero) x = (a + half) + half
         fn_val = LOG(x)
      END IF

   END FUNCTION alnrel

! **************************************************************************************************
!> \brief ...
!> \param a ...
!> \return ...
! **************************************************************************************************
   FUNCTION gam1(a) RESULT(fn_val)
!-----------------------------------------------------------------------
!     COMPUTATION OF 1/GAMMA(A+1) - 1  FOR -0.5 <= A <= 1.5
!-----------------------------------------------------------------------
      REAL(dp), INTENT(IN)                               :: a
      REAL(dp)                                           :: fn_val

      REAL(dp), PARAMETER :: p(7) = (/.577215664901533e+00_dp, -.409078193005776e+00_dp, &
         -.230975380857675e+00_dp, .597275330452234e-01_dp, .766968181649490e-02_dp, &
         -.514889771323592e-02_dp, .589597428611429e-03_dp/), q(5) = (/.100000000000000e+01_dp, &
         .427569613095214e+00_dp, .158451672430138e+00_dp, .261132021441447e-01_dp, &
         .423244297896961e-02_dp/), r(9) = (/-.422784335098468e+00_dp, -.771330383816272e+00_dp, &
         -.244757765222226e+00_dp, .118378989872749e+00_dp, .930357293360349e-03_dp, &
         -.118290993445146e-01_dp, .223047661158249e-02_dp, .266505979058923e-03_dp, &
         -.132674909766242e-03_dp/)
      REAL(dp), PARAMETER :: s1 = .273076135303957e+00_dp, s2 = .559398236957378e-01_dp

      REAL(dp)                                           :: bot, d, t, top, w

      t = a
      d = a - 0.5e0_dp
      IF (d > 0.0e0_dp) t = d - 0.5e0_dp

      IF (t > 0.e0_dp) THEN
         top = (((((p(7)*t + p(6))*t + p(5))*t + p(4))*t + p(3))*t + p(2))*t + p(1)
         bot = (((q(5)*t + q(4))*t + q(3))*t + q(2))*t + 1.0e0_dp
         w = top/bot
         IF (d > 0.0e0_dp) THEN
            fn_val = (t/a)*((w - 0.5e0_dp) - 0.5e0_dp)
         ELSE
            fn_val = a*w
         END IF
      ELSE IF (t < 0.e0_dp) THEN
         top = (((((((r(9)*t + r(8))*t + r(7))*t + r(6))*t + r(5))*t &
                  + r(4))*t + r(3))*t + r(2))*t + r(1)
         bot = (s2*t + s1)*t + 1.0e0_dp
         w = top/bot
         IF (d > 0.0e0_dp) THEN
            fn_val = t*w/a
         ELSE
            fn_val = a*((w + 0.5e0_dp) + 0.5e0_dp)
         END IF
      ELSE
         fn_val = 0.0e0_dp
      END IF

   END FUNCTION gam1

! **************************************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \return ...
! **************************************************************************************************
   FUNCTION algdiv(a, b) RESULT(fn_val)
!-----------------------------------------------------------------------

!     COMPUTATION OF LN(GAMMA(B)/GAMMA(A+B)) WHEN B >= 8

!                         --------

!     IN THIS ALGORITHM, DEL(X) IS THE FUNCTION DEFINED BY
!     LN(GAMMA(X)) = (X - 0.5)*LN(X) - X + 0.5*LN(2*PI) + DEL(X).

!-----------------------------------------------------------------------
      REAL(dp), INTENT(IN)                               :: a, b
      REAL(dp)                                           :: fn_val

      REAL(dp), PARAMETER :: c0 = .833333333333333e-01_dp, c1 = -.277777777760991e-02_dp, &
         c2 = .793650666825390e-03_dp, c3 = -.595202931351870e-03_dp, &
         c4 = .837308034031215e-03_dp, c5 = -.165322962780713e-02_dp

      REAL(dp)                                           :: c, d, h, s11, s3, s5, s7, s9, t, u, v, &
                                                            w, x, x2

      IF (a > b) THEN
         h = b/a
         c = 1.0e0_dp/(1.0e0_dp + h)
         x = h/(1.0e0_dp + h)
         d = a + (b - 0.5e0_dp)
      ELSE
         h = a/b
         c = h/(1.0e0_dp + h)
         x = 1.0e0_dp/(1.0e0_dp + h)
         d = b + (a - 0.5e0_dp)
      END IF

!                SET SN = (1 - X**N)/(1 - X)

      x2 = x*x
      s3 = 1.0e0_dp + (x + x2)
      s5 = 1.0e0_dp + (x + x2*s3)
      s7 = 1.0e0_dp + (x + x2*s5)
      s9 = 1.0e0_dp + (x + x2*s7)
      s11 = 1.0e0_dp + (x + x2*s9)

!                SET W = DEL(B) - DEL(A + B)

      t = (1.0e0_dp/b)**2
      w = ((((c5*s11*t + c4*s9)*t + c3*s7)*t + c2*s5)*t + c1*s3)*t + c0
      w = w*(c/b)

!                    COMBINE THE RESULTS

      u = d*alnrel(a/b)
      v = a*(LOG(b) - 1.0e0_dp)
      IF (u > v) THEN
         fn_val = (w - v) - u
      ELSE
         fn_val = (w - u) - v
      END IF

   END FUNCTION algdiv

! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \return ...
! **************************************************************************************************
   FUNCTION rexp(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!            EVALUATION OF THE FUNCTION EXP(X) - 1
!-----------------------------------------------------------------------
      REAL(dp), INTENT(IN)                               :: x
      REAL(dp)                                           :: fn_val

      REAL(dp), PARAMETER :: p1 = .914041914819518e-09_dp, p2 = .238082361044469e-01_dp, &
         q1 = -.499999999085958e+00_dp, q2 = .107141568980644e+00_dp, &
         q3 = -.119041179760821e-01_dp, q4 = .595130811860248e-03_dp

      REAL(dp)                                           :: e

      IF (ABS(x) < 0.15e0_dp) THEN
         fn_val = x*(((p2*x + p1)*x + 1.0e0_dp)/((((q4*x + q3)*x + q2)*x + q1)*x + 1.0e0_dp))
         RETURN
      END IF

      IF (x < 0.0e0_dp) THEN
         IF (x > -37.0e0_dp) THEN
            fn_val = (EXP(x) - 0.5e0_dp) - 0.5e0_dp
         ELSE
            fn_val = -1.0e0_dp
         END IF
      ELSE
         e = EXP(x)
         fn_val = e*(0.5e0_dp + (0.5e0_dp - 1.0e0_dp/e))
      END IF

   END FUNCTION rexp

! **************************************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \param x ...
!> \param y ...
!> \param w ...
!> \param eps ...
!> \param ierr ...
! **************************************************************************************************
   SUBROUTINE bgrat(a, b, x, y, w, eps, ierr)
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR IX(A,B) WHEN A IS LARGER THAN B.
!     THE RESULT OF THE EXPANSION IS ADDED TO W. IT IS ASSUMED
!     THAT A <= 15 AND B <= 1.  EPS IS THE TOLERANCE USED.
!     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
!-----------------------------------------------------------------------
      REAL(dp), INTENT(IN)                               :: a, b, x, y
      REAL(dp), INTENT(INOUT)                            :: w
      REAL(dp), INTENT(IN)                               :: eps
      INTEGER, INTENT(OUT)                               :: ierr

      REAL(dp), PARAMETER                                :: half = 0.5e0_dp, one = 1.e0_dp, &
                                                            quarter = 0.25e0_dp, zero = 0.e0_dp

      INTEGER                                            :: i, n
      REAL(dp)                                           :: bm1, bp2n, c(30), cn, coef, d(30), dj, &
                                                            j, l, lnx, n2, nu, q, r, s, sum, t, &
                                                            t2, tol, u, v, z

      bm1 = (b - half) - half
      nu = a + half*bm1
      IF (y > 0.375e0_dp) THEN
         lnx = LOG(x)
      ELSE
         lnx = alnrel(-y)
      END IF
      z = -nu*lnx
      IF (b*z == zero) THEN
         !               THE EXPANSION CANNOT BE COMPUTED
         ierr = 1
         RETURN
      END IF

!                 COMPUTATION OF THE EXPANSION
!                 SET R = EXP(-Z)*Z**B/GAMMA(B)

      r = b*(one + gam1(b))*EXP(b*LOG(z))
      r = r*EXP(a*lnx)*EXP(half*bm1*lnx)
      u = algdiv(b, a) + b*LOG(nu)
      u = r*EXP(-u)
      IF (u == zero) THEN
         !               THE EXPANSION CANNOT BE COMPUTED
         ierr = 1
         RETURN
      END IF
      CALL grat1(b, z, r, q=q, eps=eps)

      tol = 15.0e0_dp*eps
      v = quarter*(one/nu)**2
      t2 = quarter*lnx*lnx
      l = w/u
      j = q/r
      sum = j
      t = one
      cn = one
      n2 = zero
      DO n = 1, 30
         bp2n = b + n2
         j = (bp2n*(bp2n + one)*j + (z + bp2n + one)*t)*v
         n2 = n2 + 2.0e0_dp
         t = t*t2
         cn = cn/(n2*(n2 + one))
         c(n) = cn
         s = zero
         IF (.NOT. (n == 1)) THEN
            coef = b - n
            DO i = 1, n - 1
               s = s + coef*c(i)*d(n - i)
               coef = coef + b
            END DO
         END IF
         d(n) = bm1*cn + s/n
         dj = d(n)*j
         sum = sum + dj
         IF (sum <= zero) THEN
            !               THE EXPANSION CANNOT BE COMPUTED
            ierr = 1
            RETURN
         END IF
         IF (ABS(dj) <= tol*(sum + l)) EXIT
      END DO

!                    ADD THE RESULTS TO W

      ierr = 0
      w = w + u*sum

   END SUBROUTINE bgrat

! **************************************************************************************************
!> \brief ...
!> \param a ...
!> \param x ...
!> \param r ...
!> \param p ...
!> \param q ...
!> \param eps ...
! **************************************************************************************************
   SUBROUTINE grat1(a, x, r, p, q, eps)
!-----------------------------------------------------------------------
!           EVALUATION OF P(A,X) AND Q(A,X) WHERE A <= 1 AND
!        THE INPUT ARGUMENT R HAS THE VALUE E**(-X)*X**A/GAMMA(A)
!-----------------------------------------------------------------------
      REAL(dp), INTENT(IN)                               :: a, x, r
      REAL(dp), INTENT(OUT), OPTIONAL                    :: p, q
      REAL(dp), INTENT(IN)                               :: eps

      REAL(dp), PARAMETER                                :: half = 0.5e0_dp, one = 1.e0_dp, &
                                                            quarter = 0.25e0_dp, three = 3.e0_dp, &
                                                            two = 2.e0_dp, zero = 0.e0_dp

      REAL(dp)                                           :: a2n, a2nm1, an, b2n, b2nm1, c, g, h, j, &
                                                            l, pp, qq, sum, t, tol, w, z

      IF (a*x == zero) THEN
         IF (x <= a) THEN
            pp = zero
            qq = one
         ELSE
            pp = one
            qq = zero
         END IF
         IF (PRESENT(p)) p = pp
         IF (PRESENT(q)) q = qq
         RETURN
      END IF
      IF (a == half) THEN
      IF (x < quarter) THEN
         pp = ERF(SQRT(x))
         qq = half + (half - pp)
      ELSE
         qq = ERFC(SQRT(x))
         pp = half + (half - qq)
      END IF
      IF (PRESENT(p)) p = pp
      IF (PRESENT(q)) q = qq
      RETURN
      END IF
      IF (x < 1.1e0_dp) THEN
         !             TAYLOR SERIES FOR P(A,X)/X**A

         an = three
         c = x
         sum = x/(a + three)
         tol = three*eps/(a + one)
         an = an + one
         c = -c*(x/an)
         t = c/(a + an)
         sum = sum + t
         DO WHILE (ABS(t) > tol)
            an = an + one
            c = -c*(x/an)
            t = c/(a + an)
            sum = sum + t
         END DO
         j = a*x*((sum/6.0e0_dp - half/(a + two))*x + one/(a + one))

         z = a*LOG(x)
         h = gam1(a)
         g = one + h
         IF ((x < quarter .AND. z > -.13394e0_dp) .OR. a < x/2.59e0_dp) THEN
            l = rexp(z)
            qq = ((half + (half + l))*j - l)*g - h
            IF (qq <= zero) THEN
               pp = one
               qq = zero
            ELSE
               pp = half + (half - qq)
            END IF
         ELSE
            w = EXP(z)
            pp = w*g*(half + (half - j))
            qq = half + (half - pp)
         END IF
      ELSE
         !              CONTINUED FRACTION EXPANSION

         tol = 8.0e0_dp*eps
         a2nm1 = one
         a2n = one
         b2nm1 = x
         b2n = x + (one - a)
         c = one
         DO
            a2nm1 = x*a2n + c*a2nm1
            b2nm1 = x*b2n + c*b2nm1
            c = c + one
            a2n = a2nm1 + (c - a)*a2n
            b2n = b2nm1 + (c - a)*b2n
            a2nm1 = a2nm1/b2n
            b2nm1 = b2nm1/b2n
            a2n = a2n/b2n
            b2n = one
            IF (ABS(a2n - a2nm1/b2nm1) < tol*a2n) EXIT
         END DO

         qq = r*a2n
         pp = half + (half - qq)
      END IF

      IF (PRESENT(p)) p = pp
      IF (PRESENT(q)) q = qq

   END SUBROUTINE grat1

! **************************************************************************************************
!> \brief ...
!> \param mu ...
!> \param x ...
!> \return ...
! **************************************************************************************************
   FUNCTION esum(mu, x) RESULT(fn_val)
!-----------------------------------------------------------------------
!                    EVALUATION OF EXP(MU + X)
!-----------------------------------------------------------------------
      INTEGER, INTENT(IN)                                :: mu
      REAL(dp), INTENT(IN)                               :: x
      REAL(dp)                                           :: fn_val

      REAL(dp)                                           :: w

      IF (x > 0.0e0_dp) THEN
         IF (mu > 0 .OR. mu + x < 0.0_dp) THEN
            w = mu
            fn_val = EXP(w)*EXP(x)
         ELSE
            w = mu + x
            fn_val = EXP(w)
         END IF
      ELSE
         IF (mu < 0 .OR. mu + x < 0.0_dp) THEN
            w = mu
            fn_val = EXP(w)*EXP(x)
         ELSE
            w = mu + x
            fn_val = EXP(w)
         END IF
      END IF

   END FUNCTION esum

! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \return ...
! **************************************************************************************************
   FUNCTION rlog1(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!             EVALUATION OF THE FUNCTION X - LN(1 + X)
!-----------------------------------------------------------------------
!     A = RLOG (0.7)
!     B = RLOG (4/3)
!------------------------
      REAL(dp), INTENT(IN)                               :: x
      REAL(dp)                                           :: fn_val

      REAL(dp), PARAMETER :: a = .566749439387324e-01_dp, b = .456512608815524e-01_dp, &
         p0 = .333333333333333e+00_dp, p1 = -.224696413112536e+00_dp, &
         p2 = .620886815375787e-02_dp, q1 = -.127408923933623e+01_dp, q2 = .354508718369557e+00_dp

      REAL(dp)                                           :: r, t, u, up2, w, w1

      IF (x < -0.39e0_dp .OR. x > 0.57e0_dp) THEN
         w = (x + 0.5e0_dp) + 0.5e0_dp
         fn_val = x - LOG(w)
         RETURN
      END IF

!                 ARGUMENT REDUCTION

      IF (x < -0.18e0_dp) THEN
         u = (x + 0.3e0_dp)/0.7e0_dp
         up2 = u + 2.0e0_dp
         w1 = a - u*0.3e0_dp
      ELSEIF (x > 0.18e0_dp) THEN
         t = 0.75e0_dp*x
         u = t - 0.25e0_dp
         up2 = t + 1.75e0_dp
         w1 = b + u/3.0e0_dp
      ELSE
         u = x
         up2 = u + 2.0e0_dp
         w1 = 0.0e0_dp
      END IF

!                  SERIES EXPANSION

      r = u/up2
      t = r*r
      w = ((p2*t + p1)*t + p0)/((q2*t + q1)*t + 1.0e0_dp)
      fn_val = r*(u - 2.0e0_dp*t*w) + w1

   END FUNCTION rlog1

! **************************************************************************************************
!> \brief ...
!> \param a ...
!> \return ...
! **************************************************************************************************
   FUNCTION gamln1(a) RESULT(fn_val)
!-----------------------------------------------------------------------
!     EVALUATION OF LN(GAMMA(1 + A)) FOR -0.2 .LE. A .LE. 1.25
!-----------------------------------------------------------------------
      REAL(dp), INTENT(IN)                               :: a
      REAL(dp)                                           :: fn_val

      REAL(dp), PARAMETER :: p0 = .577215664901533e+00_dp, p1 = .844203922187225e+00_dp, &
         p2 = -.168860593646662e+00_dp, p3 = -.780427615533591e+00_dp, &
         p4 = -.402055799310489e+00_dp, p5 = -.673562214325671e-01_dp, &
         p6 = -.271935708322958e-02_dp, q1 = .288743195473681e+01_dp, &
         q2 = .312755088914843e+01_dp, q3 = .156875193295039e+01_dp, q4 = .361951990101499e+00_dp, &
         q5 = .325038868253937e-01_dp, q6 = .667465618796164e-03_dp, r0 = .422784335098467e+00_dp, &
         r1 = .848044614534529e+00_dp, r2 = .565221050691933e+00_dp, r3 = .156513060486551e+00_dp, &
         r4 = .170502484022650e-01_dp, r5 = .497958207639485e-03_dp
      REAL(dp), PARAMETER :: s1 = .124313399877507e+01_dp, s2 = .548042109832463e+00_dp, &
         s3 = .101552187439830e+00_dp, s4 = .713309612391000e-02_dp, s5 = .116165475989616e-03_dp

      REAL(dp)                                           :: w, x

      IF (a < 0.6e0_dp) THEN
         w = ((((((p6*a + p5)*a + p4)*a + p3)*a + p2)*a + p1)*a + p0)/ &
             ((((((q6*a + q5)*a + q4)*a + q3)*a + q2)*a + q1)*a + 1.0e0_dp)
         fn_val = -a*w
      ELSE
         x = (a - 0.5e0_dp) - 0.5e0_dp
         w = (((((r5*x + r4)*x + r3)*x + r2)*x + r1)*x + r0)/ &
             (((((s5*x + s4)*x + s3)*x + s2)*x + s1)*x + 1.0e0_dp)
         fn_val = x*w
      END IF

   END FUNCTION gamln1

! **************************************************************************************************
!> \brief ...
!> \param xx ...
!> \return ...
! **************************************************************************************************
   FUNCTION psi(xx) RESULT(fn_val)
!---------------------------------------------------------------------

!                 EVALUATION OF THE DIGAMMA FUNCTION

!                           -----------

!     PSI(XX) IS ASSIGNED THE VALUE 0 WHEN THE DIGAMMA FUNCTION CANNOT
!     BE COMPUTED.

!     THE MAIN COMPUTATION INVOLVES EVALUATION OF RATIONAL CHEBYSHEV
!     APPROXIMATIONS PUBLISHED IN MATH. COMP. 27, 123-127(1973) BY
!     CODY, STRECOK AND THACHER.

!---------------------------------------------------------------------
!     PSI WAS WRITTEN AT ARGONNE NATIONAL LABORATORY FOR THE FUNPACK
!     PACKAGE OF SPECIAL FUNCTION SUBROUTINES. PSI WAS MODIFIED BY
!     A.H. MORRIS (NSWC).
!---------------------------------------------------------------------
      REAL(dp), INTENT(IN)                               :: xx
      REAL(dp)                                           :: fn_val

      REAL(dp), PARAMETER :: dx0 = 1.461632144968362341262659542325721325e0_dp, p1(7) = (/ &
         .895385022981970e-02_dp, .477762828042627e+01_dp, .142441585084029e+03_dp, &
         .118645200713425e+04_dp, .363351846806499e+04_dp, .413810161269013e+04_dp, &
         .130560269827897e+04_dp/), p2(4) = (/-.212940445131011e+01_dp, -.701677227766759e+01_dp, &
         -.448616543918019e+01_dp, -.648157123766197e+00_dp/), piov4 = .785398163397448e0_dp, q1(6)&
         = (/.448452573429826e+02_dp, .520752771467162e+03_dp, .221000799247830e+04_dp, &
         .364127349079381e+04_dp, .190831076596300e+04_dp, .691091682714533e-05_dp/)
      REAL(dp), PARAMETER :: q2(4) = (/.322703493791143e+02_dp, .892920700481861e+02_dp, &
         .546117738103215e+02_dp, .777788548522962e+01_dp/)

      INTEGER                                            :: i, m, n, nq
      REAL(dp)                                           :: aug, den, sgn, upper, w, x, xmax1, xmx0, &
                                                            xsmall, z

!---------------------------------------------------------------------
!     PIOV4 = PI/4
!     DX0 = ZERO OF PSI TO EXTENDED PRECISION
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!     COEFFICIENTS FOR RATIONAL APPROXIMATION OF
!     PSI(X) / (X - X0),  0.5 <= X <= 3.0
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!     COEFFICIENTS FOR RATIONAL APPROXIMATION OF
!     PSI(X) - LN(X) + 1 / (2*X),  X > 3.0
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!     MACHINE DEPENDENT CONSTANTS ...
!        XMAX1  = THE SMALLEST POSITIVE FLOATING POINT CONSTANT
!                 WITH ENTIRELY INTEGER REPRESENTATION.  ALSO USED
!                 AS NEGATIVE OF LOWER BOUND ON ACCEPTABLE NEGATIVE
!                 ARGUMENTS AND AS THE POSITIVE ARGUMENT BEYOND WHICH
!                 PSI MAY BE REPRESENTED AS ALOG(X).
!        XSMALL = ABSOLUTE ARGUMENT BELOW WHICH PI*COTAN(PI*X)
!                 MAY BE REPRESENTED BY 1/X.
!---------------------------------------------------------------------

      xmax1 = ipmpar(3)
      xmax1 = MIN(xmax1, 1.0e0_dp/dpmpar(1))
      xsmall = 1.e-9_dp
!---------------------------------------------------------------------
      x = xx
      aug = 0.0e0_dp
      IF (x < 0.5e0_dp) THEN
!---------------------------------------------------------------------
!     X .LT. 0.5,  USE REFLECTION FORMULA
!     PSI(1-X) = PSI(X) + PI * COTAN(PI*X)
!---------------------------------------------------------------------
         IF (ABS(x) <= xsmall) THEN
            IF (x == 0.0e0_dp) THEN
               !     ERROR RETURN
               fn_val = 0.0e0_dp
               RETURN
            END IF
!---------------------------------------------------------------------
!     0 .LT. ABS(X) .LE. XSMALL.  USE 1/X AS A SUBSTITUTE
!     FOR  PI*COTAN(PI*X)
!---------------------------------------------------------------------
            aug = -1.0e0_dp/x
            x = 1.0e0_dp - x
         ELSE
!---------------------------------------------------------------------
!     REDUCTION OF ARGUMENT FOR COTAN
!---------------------------------------------------------------------
            w = -x
            sgn = piov4
            IF (w <= 0.0e0_dp) THEN
               w = -w
               sgn = -sgn
            END IF
!---------------------------------------------------------------------
!     MAKE AN ERROR EXIT IF X .LE. -XMAX1
!---------------------------------------------------------------------
            IF (w >= xmax1) THEN
               !     ERROR RETURN
               fn_val = 0.0e0_dp
               RETURN
            END IF
            nq = INT(w)
            w = w - nq
            nq = INT(w*4.0e0_dp)
            w = 4.0e0_dp*(w - nq*.25e0_dp)
!---------------------------------------------------------------------
!     W IS NOW RELATED TO THE FRACTIONAL PART OF  4.0 * X.
!     ADJUST ARGUMENT TO CORRESPOND TO VALUES IN FIRST
!     QUADRANT AND DETERMINE SIGN
!---------------------------------------------------------------------
            n = nq/2
            IF ((n + n) /= nq) w = 1.0e0_dp - w
            z = piov4*w
            m = n/2
            IF ((m + m) /= n) sgn = -sgn
!---------------------------------------------------------------------
!     DETERMINE FINAL VALUE FOR  -PI*COTAN(PI*X)
!---------------------------------------------------------------------
            n = (nq + 1)/2
            m = n/2
            m = m + m
            IF (m /= n) THEN
               aug = sgn*((SIN(z)/COS(z))*4.0e0_dp)
            ELSE
               !---------------------------------------------------------------------
               !     CHECK FOR SINGULARITY
               !---------------------------------------------------------------------
               IF (z == 0.0e0_dp) THEN
                  !     ERROR RETURN
                  fn_val = 0.0e0_dp
                  RETURN
               END IF
!---------------------------------------------------------------------
!     USE COS/SIN AS A SUBSTITUTE FOR COTAN, AND
!     SIN/COS AS A SUBSTITUTE FOR TAN
!---------------------------------------------------------------------
               aug = sgn*((COS(z)/SIN(z))*4.0e0_dp)
            END IF
            x = 1.0e0_dp - x
         END IF
      END IF
      IF (x <= 3.0e0_dp) THEN
!---------------------------------------------------------------------
!     0.5 .LE. X .LE. 3.0
!---------------------------------------------------------------------
         den = x
         upper = p1(1)*x

         DO i = 1, 5
            den = (den + q1(i))*x
            upper = (upper + p1(i + 1))*x
         END DO

         den = (upper + p1(7))/(den + q1(6))
         xmx0 = x - dx0
         fn_val = den*xmx0 + aug
         RETURN
      END IF
!---------------------------------------------------------------------
!     IF X .GE. XMAX1, PSI = LN(X)
!---------------------------------------------------------------------
      IF (x < xmax1) THEN
!---------------------------------------------------------------------
!     3.0 .LT. X .LT. XMAX1
!---------------------------------------------------------------------
         w = 1.0e0_dp/(x*x)
         den = w
         upper = p2(1)*w

         DO i = 1, 3
            den = (den + q2(i))*w
            upper = (upper + p2(i + 1))*w
         END DO

         aug = upper/(den + q2(4)) - 0.5e0_dp/x + aug
      END IF
      fn_val = aug + LOG(x)

   END FUNCTION psi

! **************************************************************************************************
!> \brief ...
!> \param a0 ...
!> \param b0 ...
!> \return ...
! **************************************************************************************************
   FUNCTION betaln(a0, b0) RESULT(fn_val)
!-----------------------------------------------------------------------
!     EVALUATION OF THE LOGARITHM OF THE BETA FUNCTION
!-----------------------------------------------------------------------
!     E = 0.5*LN(2*PI)
!--------------------------
      REAL(dp), INTENT(IN)                               :: a0, b0
      REAL(dp)                                           :: fn_val

      REAL(dp), PARAMETER                                :: e = .918938533204673e0_dp

      INTEGER                                            :: i, n
      REAL(dp)                                           :: a, b, c, h, u, v, w, z

!--------------------------

      a = MIN(a0, b0)
      b = MAX(a0, b0)
!-----------------------------------------------------------------------
!                   PROCEDURE WHEN A .GE. 8
!-----------------------------------------------------------------------
      IF (a >= 8.0e0_dp) THEN
         w = bcorr(a, b)
         h = a/b
         c = h/(1.0e0_dp + h)
         u = -(a - 0.5e0_dp)*LOG(c)
         v = b*alnrel(h)
         IF (u > v) THEN
            fn_val = (((-0.5e0_dp*LOG(b) + e) + w) - v) - u
         ELSE
            fn_val = (((-0.5e0_dp*LOG(b) + e) + w) - u) - v
         END IF
!-----------------------------------------------------------------------
!                   PROCEDURE WHEN A .LT. 1
!-----------------------------------------------------------------------
      ELSEIF (a < 1.0e0_dp) THEN
         IF (b < 8.0e0_dp) THEN
            fn_val = LOG_GAMMA(a) + (LOG_GAMMA(b) - LOG_GAMMA(a + b))
         ELSE
            fn_val = LOG_GAMMA(a) + algdiv(a, b)
         END IF
!-----------------------------------------------------------------------
!                PROCEDURE WHEN 1 .LE. A .LT. 8
!-----------------------------------------------------------------------
      ELSEIF (a <= 2.0e0_dp) THEN
         IF (b <= 2.0e0_dp) THEN
            fn_val = LOG_GAMMA(a) + LOG_GAMMA(b) - gsumln(a, b)
            RETURN
         END IF
         w = 0.0e0_dp
         IF (b < 8.0e0_dp) THEN
            !                 REDUCTION OF B WHEN B .LT. 8

            n = INT(b - 1.0e0_dp)
            z = 1.0e0_dp
            DO i = 1, n
               b = b - 1.0e0_dp
               z = z*(b/(a + b))
            END DO
            fn_val = w + LOG(z) + (LOG_GAMMA(a) + (LOG_GAMMA(b) - gsumln(a, b)))
            RETURN
         END IF
         fn_val = LOG_GAMMA(a) + algdiv(a, b)
!                REDUCTION OF A WHEN B .LE. 1000
      ELSEIF (b <= 1000.0e0_dp) THEN
         n = INT(a - 1.0e0_dp)
         w = 1.0e0_dp
         DO i = 1, n
            a = a - 1.0e0_dp
            h = a/b
            w = w*(h/(1.0e0_dp + h))
         END DO
         w = LOG(w)
         IF (b >= 8.0e0_dp) THEN
            fn_val = w + LOG_GAMMA(a) + algdiv(a, b)
            RETURN
         END IF

         !                 REDUCTION OF B WHEN B .LT. 8

         n = INT(b - 1.0e0_dp)
         z = 1.0e0_dp
         DO i = 1, n
            b = b - 1.0e0_dp
            z = z*(b/(a + b))
         END DO
         fn_val = w + LOG(z) + (LOG_GAMMA(a) + (LOG_GAMMA(b) - gsumln(a, b)))
      ELSE
!                REDUCTION OF A WHEN B .GT. 1000

         n = INT(a - 1.0e0_dp)
         w = 1.0e0_dp
         DO i = 1, n
            a = a - 1.0e0_dp
            w = w*(a/(1.0e0_dp + a/b))
         END DO
         fn_val = (LOG(w) - n*LOG(b)) + (LOG_GAMMA(a) + algdiv(a, b))
      END IF

   END FUNCTION betaln

! **************************************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \return ...
! **************************************************************************************************
   FUNCTION gsumln(a, b) RESULT(fn_val)
!-----------------------------------------------------------------------
!          EVALUATION OF THE FUNCTION LN(GAMMA(A + B))
!          FOR 1 .LE. A .LE. 2  AND  1 .LE. B .LE. 2
!-----------------------------------------------------------------------
      REAL(dp), INTENT(IN)                               :: a, b
      REAL(dp)                                           :: fn_val

      REAL(dp)                                           :: x

      x = a + b - 2.e0_dp
      IF (x <= 0.25e0_dp) THEN
         fn_val = gamln1(1.0e0_dp + x)
      ELSEIF (x <= 1.25e0_dp) THEN
         fn_val = gamln1(x) + alnrel(x)
      ELSE
         fn_val = gamln1(x - 1.0e0_dp) + LOG(x*(1.0e0_dp + x))
      END IF
   END FUNCTION gsumln

! **************************************************************************************************
!> \brief ...
!> \param a0 ...
!> \param b0 ...
!> \return ...
! **************************************************************************************************
   FUNCTION bcorr(a0, b0) RESULT(fn_val)
!-----------------------------------------------------------------------

!     EVALUATION OF  DEL(A0) + DEL(B0) - DEL(A0 + B0)  WHERE
!     LN(GAMMA(A)) = (A - 0.5)*LN(A) - A + 0.5*LN(2*PI) + DEL(A).
!     IT IS ASSUMED THAT A0 .GE. 8 AND B0 .GE. 8.

!-----------------------------------------------------------------------
      REAL(dp), INTENT(IN)                               :: a0, b0
      REAL(dp)                                           :: fn_val

      REAL(dp), PARAMETER :: c0 = .833333333333333e-01_dp, c1 = -.277777777760991e-02_dp, &
         c2 = .793650666825390e-03_dp, c3 = -.595202931351870e-03_dp, &
         c4 = .837308034031215e-03_dp, c5 = -.165322962780713e-02_dp

      REAL(dp)                                           :: a, b, c, h, s11, s3, s5, s7, s9, t, w, &
                                                            x, x2

      a = MIN(a0, b0)
      b = MAX(a0, b0)

      h = a/b
      c = h/(1.0e0_dp + h)
      x = 1.0e0_dp/(1.0e0_dp + h)
      x2 = x*x

!                SET SN = (1 - X**N)/(1 - X)

      s3 = 1.0e0_dp + (x + x2)
      s5 = 1.0e0_dp + (x + x2*s3)
      s7 = 1.0e0_dp + (x + x2*s5)
      s9 = 1.0e0_dp + (x + x2*s7)
      s11 = 1.0e0_dp + (x + x2*s9)

!                SET W = DEL(B) - DEL(A + B)

      t = (1.0e0_dp/b)**2
      w = ((((c5*s11*t + c4*s9)*t + c3*s7)*t + c2*s5)*t + c1*s3)*t + c0
      w = w*(c/b)

!                   COMPUTE  DEL(A) + W

      t = (1.0e0_dp/a)**2
      fn_val = (((((c5*t + c4)*t + c3)*t + c2)*t + c1)*t + c0)/a + w
      RETURN
   END FUNCTION bcorr

! **************************************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \param x ...
!> \param y ...
!> \param w ...
!> \param w1 ...
!> \param ierr ...
! **************************************************************************************************
   SUBROUTINE bratio(a, b, x, y, w, w1, ierr)
!-----------------------------------------------------------------------

!            EVALUATION OF THE INCOMPLETE BETA FUNCTION IX(A,B)

!                     --------------------

!     IT IS ASSUMED THAT A AND B ARE NONNEGATIVE, AND THAT X <= 1
!     AND Y = 1 - X.  BRATIO ASSIGNS W AND W1 THE VALUES

!                      W  = IX(A,B)
!                      W1 = 1 - IX(A,B)

!     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
!     IF NO INPUT ERRORS ARE DETECTED THEN IERR IS SET TO 0 AND
!     W AND W1 ARE COMPUTED. OTHERWISE, IF AN ERROR IS DETECTED,
!     THEN W AND W1 ARE ASSIGNED THE VALUE 0 AND IERR IS SET TO
!     ONE OF THE FOLLOWING VALUES ...

!        IERR = 1  IF A OR B IS NEGATIVE
!        IERR = 2  IF A = B = 0
!        IERR = 3  IF X .LT. 0 OR X .GT. 1
!        IERR = 4  IF Y .LT. 0 OR Y .GT. 1
!        IERR = 5  IF X + Y .NE. 1
!        IERR = 6  IF X = A = 0
!        IERR = 7  IF Y = B = 0

!--------------------
!     WRITTEN BY ALFRED H. MORRIS, JR.
!        NAVAL SURFACE WARFARE CENTER
!        DAHLGREN, VIRGINIA
!     REVISED ... APRIL 1993
!-----------------------------------------------------------------------
      REAL(dp), INTENT(IN)                               :: a, b, x, y
      REAL(dp), INTENT(OUT)                              :: w, w1
      INTEGER, INTENT(OUT)                               :: ierr

      INTEGER                                            :: ierr1, ind, n
      REAL(dp)                                           :: a0, b0, eps, lambda, t, x0, y0, z

!-----------------------------------------------------------------------
!     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE SMALLEST
!            FLOATING POINT NUMBER FOR WHICH 1.0 + EPS .GT. 1.0

      eps = dpmpar(1)

      w = 0.0e0_dp
      w1 = 0.0e0_dp
      IF (a < 0.0e0_dp .OR. b < 0.0e0_dp) THEN
         ierr = 1
         RETURN
      END IF
      IF (a == 0.0e0_dp .AND. b == 0.0e0_dp) THEN
         ierr = 2
         RETURN
      END IF
      IF (x < 0.0e0_dp .OR. x > 1.0e0_dp) THEN
         ierr = 3
         RETURN
      END IF
      IF (y < 0.0e0_dp .OR. y > 1.0e0_dp) THEN
         ierr = 4
         RETURN
      END IF
      z = ((x + y) - 0.5e0_dp) - 0.5e0_dp
      IF (ABS(z) > 3.0e0_dp*eps) THEN
         ierr = 5
         RETURN
      END IF

      ierr = 0

      IF (x == 0.0e0_dp .OR. a == 0.0e0_dp) THEN
         IF (x == 0.0e0_dp .AND. a == 0.0e0_dp) THEN
            ierr = 6
         ELSE
            w = 0.0e0_dp
            w1 = 1.0e0_dp
         END IF
         RETURN
      END IF

      IF (y == 0.0e0_dp .OR. b == 0.0e0_dp) THEN
         IF (y == 0.0e0_dp .AND. b == 0.0e0_dp) THEN
            ierr = 7
         ELSE
            w = 1.0e0_dp
            w1 = 0.0e0_dp
         END IF
         RETURN
      END IF

      eps = MAX(eps, 1.e-15_dp)
      IF (MAX(a, b) < 1.e-3_dp*eps) THEN
!           PROCEDURE FOR A AND B .LT. 1.E-3*EPS
         w = b/(a + b)
         w1 = a/(a + b)
         RETURN
      END IF

      ind = 0
      a0 = a
      b0 = b
      x0 = x
      y0 = y
      IF (MIN(a0, b0) > 1.0e0_dp) GO TO 30

!             PROCEDURE FOR A0 .LE. 1 OR B0 .LE. 1

      IF (x > 0.5e0_dp) THEN
         ind = 1
         a0 = b
         b0 = a
         x0 = y
         y0 = x
      END IF

      IF (b0 < MIN(eps, eps*a0)) GO TO 80
      IF (a0 < MIN(eps, eps*b0) .AND. b0*x0 <= 1.0e0_dp) GO TO 90
      IF (MAX(a0, b0) > 1.0e0_dp) GO TO 20
      IF (a0 >= MIN(0.2e0_dp, b0)) GO TO 100
      IF (x0**a0 <= 0.9e0_dp) GO TO 100
      IF (x0 >= 0.3e0_dp) GO TO 110
      n = 20
      GO TO 130

20    IF (b0 <= 1.0e0_dp) GO TO 100
      IF (x0 >= 0.3e0_dp) GO TO 110
      IF (x0 < 0.1e0_dp .AND. (x0*b0)**a0 <= 0.7e0_dp) GO TO 100
      IF (b0 > 15.0e0_dp) GO TO 131
      n = 20
      GO TO 130

!             PROCEDURE FOR A0 .GT. 1 AND B0 .GT. 1

30    IF (a > b) THEN
         lambda = (a + b)*y - b
      ELSE
         lambda = a - (a + b)*x
      END IF

      IF (lambda < 0.0e0_dp) THEN
         ind = 1
         a0 = b
         b0 = a
         x0 = y
         y0 = x
         lambda = ABS(lambda)
      END IF

      IF (b0 < 40.0e0_dp .AND. b0*x0 <= 0.7e0_dp) GO TO 100
      IF (b0 < 40.0e0_dp) GO TO 140
      IF (a0 > b0) GO TO 50
      IF (a0 <= 100.0e0_dp) GO TO 120
      IF (lambda > 0.03e0_dp*a0) GO TO 120
      GO TO 180
50    IF (b0 <= 100.0e0_dp) GO TO 120
      IF (lambda > 0.03e0_dp*b0) GO TO 120
      GO TO 180

!            EVALUATION OF THE APPROPRIATE ALGORITHM

80    w = fpser(a0, b0, x0, eps)
      w1 = 0.5e0_dp + (0.5e0_dp - w)
      GO TO 220

90    w1 = apser(a0, b0, x0, eps)
      w = 0.5e0_dp + (0.5e0_dp - w1)
      GO TO 220

100   w = bpser(a0, b0, x0, eps)
      w1 = 0.5e0_dp + (0.5e0_dp - w)
      GO TO 220

110   w1 = bpser(b0, a0, y0, eps)
      w = 0.5e0_dp + (0.5e0_dp - w1)
      GO TO 220

120   w = bfrac(a0, b0, x0, y0, lambda, 15.0e0_dp*eps)
      w1 = 0.5e0_dp + (0.5e0_dp - w)
      GO TO 220

130   w1 = bup(b0, a0, y0, x0, n, eps)
      b0 = b0 + n
131   CALL bgrat(b0, a0, y0, x0, w1, eps, ierr1)
      IF (ierr1 > 0) CPABORT("Error in BGRAT")
      w = 0.5e0_dp + (0.5e0_dp - w1)
      GO TO 220

140   n = INT(b0)
      b0 = b0 - n
      IF (b0 /= 0.0e0_dp) GO TO 141
      n = n - 1
      b0 = 1.0e0_dp
141   w = bup(b0, a0, y0, x0, n, eps)
      IF (x0 > 0.7e0_dp) GO TO 150
      w = w + bpser(a0, b0, x0, eps)
      w1 = 0.5e0_dp + (0.5e0_dp - w)
      GO TO 220

150   IF (a0 > 15.0e0_dp) GO TO 151
      n = 20
      w = w + bup(a0, b0, x0, y0, n, eps)
      a0 = a0 + n
151   CALL bgrat(a0, b0, x0, y0, w, eps, ierr1)
      w1 = 0.5e0_dp + (0.5e0_dp - w)
      GO TO 220

180   w = basym(a0, b0, lambda, 100.0e0_dp*eps)
      w1 = 0.5e0_dp + (0.5e0_dp - w)
      GO TO 220

!               TERMINATION OF THE PROCEDURE

220   IF (ind == 0) RETURN
      t = w
      w = w1
      w1 = t

   END SUBROUTINE bratio

! **************************************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \param x ...
!> \param eps ...
!> \return ...
! **************************************************************************************************
   FUNCTION fpser(a, b, x, eps) RESULT(fn_val)
!-----------------------------------------------------------------------

!                 EVALUATION OF I (A,B)
!                                X

!          FOR B .LT. MIN(EPS,EPS*A) AND X .LE. 0.5.

!-----------------------------------------------------------------------
      REAL(dp), INTENT(IN)                               :: a, b, x, eps
      REAL(dp)                                           :: fn_val

      REAL(dp)                                           :: an, c, s, t, tol

!                  SET  FPSER = X**A

      fn_val = 1.0e0_dp
      IF (a > 1.e-3_dp*eps) THEN
         fn_val = 0.0e0_dp
         t = a*LOG(x)
         IF (t < dxparg(1)) RETURN
         fn_val = EXP(t)
      END IF

!                NOTE THAT 1/B(A,B) = B

      fn_val = (b/a)*fn_val
      tol = eps/a
      an = a + 1.0e0_dp
      t = x
      s = t/an
      an = an + 1.0e0_dp
      t = x*t
      c = t/an
      s = s + c
      DO WHILE (ABS(c) > tol)
         an = an + 1.0e0_dp
         t = x*t
         c = t/an
         s = s + c
      END DO

      fn_val = fn_val*(1.0e0_dp + a*s)

   END FUNCTION fpser

! **************************************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \param x ...
!> \param eps ...
!> \return ...
! **************************************************************************************************
   FUNCTION apser(a, b, x, eps) RESULT(fn_val)
!-----------------------------------------------------------------------
!     APSER YIELDS THE INCOMPLETE BETA RATIO I(SUB(1-X))(B,A) FOR
!     A .LE. MIN(EPS,EPS*B), B*X .LE. 1, AND X .LE. 0.5. USED WHEN
!     A IS VERY SMALL. USE ONLY IF ABOVE INEQUALITIES ARE SATISFIED.
!-----------------------------------------------------------------------
      REAL(dp), INTENT(IN)                               :: a, b, x, eps
      REAL(dp)                                           :: fn_val

      REAL(dp), PARAMETER                                :: g = .577215664901533e0_dp

      REAL(dp)                                           :: aj, bx, c, j, s, t, tol

      bx = b*x
      t = x - bx
      IF (b*eps > 2.e-2_dp) THEN
         c = LOG(bx) + g + t
      ELSE
         c = LOG(x) + psi(b) + g + t
      END IF

      tol = 5.0e0_dp*eps*ABS(c)
      j = 1.0e0_dp
      s = 0.0e0_dp
      j = j + 1.0e0_dp
      t = t*(x - bx/j)
      aj = t/j
      s = s + aj
      DO WHILE (ABS(aj) > tol)
         t = t*(x - bx/j)
         aj = t/j
         s = s + aj
      END DO

      fn_val = -a*(c + s)

   END FUNCTION apser

! **************************************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \param x ...
!> \param eps ...
!> \return ...
! **************************************************************************************************
   FUNCTION bpser(a, b, x, eps) RESULT(fn_val)
!-----------------------------------------------------------------------
!     POWER SERIES EXPANSION FOR EVALUATING IX(A,B) WHEN B .LE. 1
!     OR B*X .LE. 0.7.  EPS IS THE TOLERANCE USED.
!-----------------------------------------------------------------------
      REAL(dp), INTENT(IN)                               :: a, b, x, eps
      REAL(dp)                                           :: fn_val

      INTEGER                                            :: i, m
      REAL(dp)                                           :: a0, apb, b0, c, n, sum, t, tol, u, w, z

      fn_val = 0.0e0_dp
      IF (x == 0.0e0_dp) RETURN
!-----------------------------------------------------------------------
!            COMPUTE THE FACTOR X**A/(A*BETA(A,B))
!-----------------------------------------------------------------------
      a0 = MIN(a, b)
      b0 = MAX(a, b)
      IF (a0 >= 1.0e0_dp) THEN
         z = a*LOG(x) - betaln(a, b)
         fn_val = EXP(z)/a
      ELSEIF (b0 >= 8.0e0_dp) THEN
         u = gamln1(a0) + algdiv(a0, b0)
         z = a*LOG(x) - u
         fn_val = (a0/a)*EXP(z)
      ELSEIF (b0 > 1.0e0_dp) THEN
!         PROCEDURE FOR A0 .LT. 1 AND 1 .LT. B0 .LT. 8

         u = gamln1(a0)
         m = INT(b0 - 1.0e0_dp)
         IF (m >= 1) THEN
            c = 1.0e0_dp
            DO i = 1, m
               b0 = b0 - 1.0e0_dp
               c = c*(b0/(a0 + b0))
            END DO
            u = LOG(c) + u
         END IF

         z = a*LOG(x) - u
         b0 = b0 - 1.0e0_dp
         apb = a0 + b0
         IF (apb > 1.0e0_dp) THEN
            u = a0 + b0 - 1.e0_dp
            t = (1.0e0_dp + gam1(u))/apb
         ELSE
            t = 1.0e0_dp + gam1(apb)
         END IF
         fn_val = EXP(z)*(a0/a)*(1.0e0_dp + gam1(b0))/t
      ELSE

         !            PROCEDURE FOR A0 .LT. 1 AND B0 .LE. 1

         fn_val = x**a
         IF (fn_val == 0.0e0_dp) RETURN

         apb = a + b
         IF (apb > 1.0e0_dp) THEN
            u = a + b - 1.e0_dp
            z = (1.0e0_dp + gam1(u))/apb
         ELSE
            z = 1.0e0_dp + gam1(apb)
         END IF

         c = (1.0e0_dp + gam1(a))*(1.0e0_dp + gam1(b))/z
         fn_val = fn_val*c*(b/apb)
      END IF

!            PROCEDURE FOR A0 .LT. 1 AND B0 .GE. 8

      IF (fn_val == 0.0e0_dp .OR. a <= 0.1e0_dp*eps) RETURN
!-----------------------------------------------------------------------
!                     COMPUTE THE SERIES
!-----------------------------------------------------------------------
      sum = 0.0e0_dp
      n = 0.0e0_dp
      c = 1.0e0_dp
      tol = eps/a
      n = n + 1.0e0_dp
      c = c*(0.5e0_dp + (0.5e0_dp - b/n))*x
      w = c/(a + n)
      sum = sum + w
      DO WHILE (ABS(w) > tol)
         n = n + 1.0e0_dp
         c = c*(0.5e0_dp + (0.5e0_dp - b/n))*x
         w = c/(a + n)
         sum = sum + w
      END DO
      fn_val = fn_val*(1.0e0_dp + a*sum)

   END FUNCTION bpser

! **************************************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \param x ...
!> \param y ...
!> \param n ...
!> \param eps ...
!> \return ...
! **************************************************************************************************
   FUNCTION bup(a, b, x, y, n, eps) RESULT(fn_val)
!-----------------------------------------------------------------------
!     EVALUATION OF IX(A,B) - IX(A+N,B) WHERE N IS A POSITIVE INTEGER.
!     EPS IS THE TOLERANCE USED.
!-----------------------------------------------------------------------
      REAL(dp), INTENT(IN)                               :: a, b, x, y
      INTEGER, INTENT(IN)                                :: n
      REAL(dp), INTENT(IN)                               :: eps
      REAL(dp)                                           :: fn_val

      INTEGER                                            :: i, k, kp1, mu, nm1
      REAL(dp)                                           :: ap1, apb, d, l, r, t, w

!          OBTAIN THE SCALING FACTOR EXP(-MU) AND
!             EXP(MU)*(X**A*Y**B/BETA(A,B))/A

      apb = a + b
      ap1 = a + 1.0e0_dp
      mu = 0
      d = 1.0e0_dp
      IF (.NOT. (n == 1 .OR. a < 1.0e0_dp .OR. apb < 1.1e0_dp*ap1)) THEN
         mu = INT(ABS(dxparg(1)))
         k = INT(dxparg(0))
         IF (k < mu) mu = k
         t = mu
         d = EXP(-t)
      END IF

      fn_val = brcmp1(mu, a, b, x, y)/a
      IF (n == 1 .OR. fn_val == 0.0e0_dp) RETURN
      nm1 = n - 1
      w = d

!          LET K BE THE INDEX OF THE MAXIMUM TERM

      k = 0
      IF (b > 1.0e0_dp) THEN
         IF (y <= 1.e-4_dp) THEN
            k = nm1
            DO i = 1, k
               l = i - 1
               d = ((apb + l)/(ap1 + l))*x*d
               w = w + d
            END DO
            IF (k == nm1) THEN
               fn_val = fn_val*w
               RETURN
            END IF
         ELSE
            r = (b - 1.0e0_dp)*x/y - a
            IF (r >= 1.0e0_dp) THEN
               k = nm1
               t = nm1
               IF (r < t) k = INT(r)

!          ADD THE INCREASING TERMS OF THE SERIES

               DO i = 1, k
                  l = i - 1
                  d = ((apb + l)/(ap1 + l))*x*d
                  w = w + d
               END DO
               IF (k == nm1) THEN
                  fn_val = fn_val*w
                  RETURN
               END IF
            END IF
         END IF
      END IF

!          ADD THE REMAINING TERMS OF THE SERIES

      kp1 = k + 1
      DO i = kp1, nm1
         l = i - 1
         d = ((apb + l)/(ap1 + l))*x*d
         w = w + d
         IF (d <= eps*w) EXIT
      END DO

!               TERMINATE THE PROCEDURE

      fn_val = fn_val*w

   END FUNCTION bup

! **************************************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \param x ...
!> \param y ...
!> \param lambda ...
!> \param eps ...
!> \return ...
! **************************************************************************************************
   FUNCTION bfrac(a, b, x, y, lambda, eps) RESULT(fn_val)
!-----------------------------------------------------------------------
!     CONTINUED FRACTION EXPANSION FOR IX(A,B) WHEN A,B .GT. 1.
!     IT IS ASSUMED THAT  LAMBDA = (A + B)*Y - B.
!-----------------------------------------------------------------------
      REAL(dp), INTENT(IN)                               :: a, b, x, y, lambda, eps
      REAL(dp)                                           :: fn_val

      REAL(dp)                                           :: alpha, an, anp1, beta, bn, bnp1, c, c0, &
                                                            c1, e, n, p, r, r0, s, t, w, yp1

      fn_val = brcomp(a, b, x, y)
      IF (fn_val == 0.0e0_dp) RETURN

      c = 1.0e0_dp + lambda
      c0 = b/a
      c1 = 1.0e0_dp + 1.0e0_dp/a
      yp1 = y + 1.0e0_dp

      n = 0.0e0_dp
      p = 1.0e0_dp
      s = a + 1.0e0_dp
      an = 0.0e0_dp
      bn = 1.0e0_dp
      anp1 = 1.0e0_dp
      bnp1 = c/c1
      r = c1/c

!        CONTINUED FRACTION CALCULATION

      DO WHILE (.TRUE.)
         n = n + 1.0e0_dp
         t = n/a
         w = n*(b - n)*x
         e = a/s
         alpha = (p*(p + c0)*e*e)*(w*x)
         IF (alpha <= 0.0e0_dp) THEN
!                 TERMINATION
            fn_val = fn_val*r
            RETURN
         END IF
         e = (1.0e0_dp + t)/(c1 + t + t)
         beta = n + w/s + e*(c + n*yp1)
         p = 1.0e0_dp + t
         s = s + 2.0e0_dp

!        UPDATE AN, BN, ANP1, AND BNP1

         t = alpha*an + beta*anp1
         an = anp1
         anp1 = t
         t = alpha*bn + beta*bnp1
         bn = bnp1
         bnp1 = t
         r0 = r
         r = anp1/bnp1
         IF (ABS(r - r0) <= eps*r) THEN
!                 TERMINATION
            fn_val = fn_val*r
            RETURN
         END IF

!        RESCALE AN, BN, ANP1, AND BNP1

         an = an/bnp1
         bn = bn/bnp1
         anp1 = r
         bnp1 = 1.0e0_dp
      END DO

   END FUNCTION bfrac

! **************************************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \param x ...
!> \param y ...
!> \return ...
! **************************************************************************************************
   FUNCTION brcomp(a, b, x, y) RESULT(fn_val)
!-----------------------------------------------------------------------
!               EVALUATION OF X**A*Y**B/BETA(A,B)
!-----------------------------------------------------------------------
      REAL(dp), INTENT(IN)                               :: a, b, x, y
      REAL(dp)                                           :: fn_val

      REAL(dp), PARAMETER                                :: const = 0.398942280401433e0_dp

      INTEGER                                            :: i, n
      REAL(dp)                                           :: a0, apb, b0, c, e, h, lambda, lnx, lny, &
                                                            t, u, v, x0, y0, z

!-----------------
!     CONST = 1/SQRT(2*PI)
!-----------------

      fn_val = 0.0e0_dp
      IF (x == 0.0e0_dp .OR. y == 0.0e0_dp) RETURN
      a0 = MIN(a, b)
      IF (a0 >= 8.0e0_dp) THEN
!-----------------------------------------------------------------------
!              PROCEDURE FOR A .GE. 8 AND B .GE. 8
!-----------------------------------------------------------------------
         IF (a > b) THEN
            h = b/a
            x0 = 1.0e0_dp/(1.0e0_dp + h)
            y0 = h/(1.0e0_dp + h)
            lambda = (a + b)*y - b
         ELSE
            h = a/b
            x0 = h/(1.0e0_dp + h)
            y0 = 1.0e0_dp/(1.0e0_dp + h)
            lambda = a - (a + b)*x
         END IF

         e = -lambda/a
         IF (ABS(e) > 0.6e0_dp) THEN
            u = e - LOG(x/x0)
         ELSE
            u = rlog1(e)
         END IF

         e = lambda/b
         IF (ABS(e) > 0.6e0_dp) THEN
            v = e - LOG(y/y0)
         ELSE
            v = rlog1(e)
         END IF

         z = EXP(-(a*u + b*v))
         fn_val = const*SQRT(b*x0)*z*EXP(-bcorr(a, b))
         RETURN
      END IF

      IF (x > 0.375e0_dp) THEN
         IF (y > 0.375e0_dp) THEN
            lnx = LOG(x)
            lny = LOG(y)
         ELSE
            lnx = alnrel(-y)
            lny = LOG(y)
         END IF
      ELSE
         lnx = LOG(x)
         lny = alnrel(-x)
      END IF

      z = a*lnx + b*lny
      IF (a0 >= 1.0e0_dp) THEN
         z = z - betaln(a, b)
         fn_val = EXP(z)
         RETURN
      END IF
!-----------------------------------------------------------------------
!              PROCEDURE FOR A .LT. 1 OR B .LT. 1
!-----------------------------------------------------------------------
      b0 = MAX(a, b)
      IF (b0 >= 8.0e0_dp) THEN
!                   ALGORITHM FOR B0 .GE. 8
         u = gamln1(a0) + algdiv(a0, b0)
         fn_val = a0*EXP(z - u)
      END IF
      IF (b0 <= 1.0e0_dp) THEN

!                   ALGORITHM FOR B0 .LE. 1

         fn_val = EXP(z)
         IF (fn_val == 0.0e0_dp) RETURN

         apb = a + b
         IF (apb > 1.0e0_dp) THEN
            u = a + b - 1.e0_dp
            z = (1.0e0_dp + gam1(u))/apb
         ELSE
            z = 1.0e0_dp + gam1(apb)
         END IF

         c = (1.0e0_dp + gam1(a))*(1.0e0_dp + gam1(b))/z
         fn_val = fn_val*(a0*c)/(1.0e0_dp + a0/b0)
         RETURN
      END IF

!                ALGORITHM FOR 1 .LT. B0 .LT. 8

      u = gamln1(a0)
      n = INT(b0 - 1.0e0_dp)
      IF (n >= 1) THEN
         c = 1.0e0_dp
         DO i = 1, n
            b0 = b0 - 1.0e0_dp
            c = c*(b0/(a0 + b0))
         END DO
         u = LOG(c) + u
      END IF

      z = z - u
      b0 = b0 - 1.0e0_dp
      apb = a0 + b0
      IF (apb > 1.0e0_dp) THEN
         u = a0 + b0 - 1.e0_dp
         t = (1.0e0_dp + gam1(u))/apb
      ELSE
         t = 1.0e0_dp + gam1(apb)
      END IF
      fn_val = a0*EXP(z)*(1.0e0_dp + gam1(b0))/t

   END FUNCTION brcomp

! **************************************************************************************************
!> \brief ...
!> \param mu ...
!> \param a ...
!> \param b ...
!> \param x ...
!> \param y ...
!> \return ...
! **************************************************************************************************
   FUNCTION brcmp1(mu, a, b, x, y) RESULT(fn_val)
!-----------------------------------------------------------------------
!          EVALUATION OF  EXP(MU) * (X**A*Y**B/BETA(A,B))
!-----------------------------------------------------------------------
      INTEGER, INTENT(IN)                                :: mu
      REAL(dp), INTENT(IN)                               :: a, b, x, y
      REAL(dp)                                           :: fn_val

      REAL(dp), PARAMETER                                :: const = 0.398942280401433e0_dp

      INTEGER                                            :: i, n
      REAL(dp)                                           :: a0, apb, b0, c, e, h, lambda, lnx, lny, &
                                                            t, u, v, x0, y0, z

!-----------------
!     CONST = 1/SQRT(2*PI)
!-----------------

      a0 = MIN(a, b)
      IF (a0 >= 8.0e0_dp) THEN
!-----------------------------------------------------------------------
!              PROCEDURE FOR A .GE. 8 AND B .GE. 8
!-----------------------------------------------------------------------
         IF (a > b) THEN
            h = b/a
            x0 = 1.0e0_dp/(1.0e0_dp + h)
            y0 = h/(1.0e0_dp + h)
            lambda = (a + b)*y - b
         END IF
         h = a/b
         x0 = h/(1.0e0_dp + h)
         y0 = 1.0e0_dp/(1.0e0_dp + h)
         lambda = a - (a + b)*x

         e = -lambda/a
         IF (ABS(e) > 0.6e0_dp) THEN
            u = e - LOG(x/x0)
         ELSE
            u = rlog1(e)
         END IF

         e = lambda/b
         IF (ABS(e) <= 0.6e0_dp) THEN
            v = rlog1(e)
         ELSE
            v = e - LOG(y/y0)
         END IF

         z = esum(mu, -(a*u + b*v))
         fn_val = const*SQRT(b*x0)*z*EXP(-bcorr(a, b))
      END IF

      IF (x > 0.375e0_dp) THEN
         IF (y > 0.375e0_dp) THEN
            lnx = LOG(x)
            lny = LOG(y)
         ELSE
            lnx = alnrel(-y)
            lny = LOG(y)
         END IF
      ELSE
         lnx = LOG(x)
         lny = alnrel(-x)
      END IF
      z = a*lnx + b*lny
      IF (a0 >= 1.0e0_dp) THEN
         z = z - betaln(a, b)
         fn_val = esum(mu, z)
         RETURN
      END IF
!-----------------------------------------------------------------------
!              PROCEDURE FOR A .LT. 1 OR B .LT. 1
!-----------------------------------------------------------------------
      b0 = MAX(a, b)
      IF (b0 >= 8.0e0_dp) THEN
!                   ALGORITHM FOR B0 .GE. 8

         u = gamln1(a0) + algdiv(a0, b0)
         fn_val = a0*esum(mu, z - u)
         RETURN
      END IF
      IF (b0 <= 1.0e0_dp) THEN

!                   ALGORITHM FOR B0 .LE. 1

         fn_val = esum(mu, z)
         IF (fn_val == 0.0e0_dp) RETURN

         apb = a + b
         IF (apb > 1.0e0_dp) THEN
            u = a + b - 1.e0_dp
            z = (1.0e0_dp + gam1(u))/apb
         ELSE
            z = 1.0e0_dp + gam1(apb)
         END IF

         c = (1.0e0_dp + gam1(a))*(1.0e0_dp + gam1(b))/z
         fn_val = fn_val*(a0*c)/(1.0e0_dp + a0/b0)
         RETURN
      END IF

!                ALGORITHM FOR 1 .LT. B0 .LT. 8

      u = gamln1(a0)
      n = INT(b0 - 1.0e0_dp)
      IF (n >= 1) THEN
         c = 1.0e0_dp
         DO i = 1, n
            b0 = b0 - 1.0e0_dp
            c = c*(b0/(a0 + b0))
         END DO
         u = LOG(c) + u
      END IF

      z = z - u
      b0 = b0 - 1.0e0_dp
      apb = a0 + b0
      IF (apb > 1.0e0_dp) THEN
         u = a0 + b0 - 1.e0_dp
         t = (1.0e0_dp + gam1(u))/apb
      ELSE
         t = 1.0e0_dp + gam1(apb)
      END IF
      fn_val = a0*esum(mu, z)*(1.0e0_dp + gam1(b0))/t

   END FUNCTION brcmp1

! **************************************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \param lambda ...
!> \param eps ...
!> \return ...
! **************************************************************************************************
   FUNCTION basym(a, b, lambda, eps) RESULT(fn_val)
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR IX(A,B) FOR LARGE A AND B.
!     LAMBDA = (A + B)*Y - B  AND EPS IS THE TOLERANCE USED.
!     IT IS ASSUMED THAT LAMBDA IS NONNEGATIVE AND THAT
!     A AND B ARE GREATER THAN OR EQUAL TO 15.
!-----------------------------------------------------------------------
      REAL(dp), INTENT(IN)                               :: a, b, lambda, eps
      REAL(dp)                                           :: fn_val

      INTEGER, PARAMETER                                 :: num = 20
      REAL(dp), PARAMETER                                :: e0 = 1.12837916709551e0_dp, &
                                                            e1 = .353553390593274e0_dp

      INTEGER                                            :: i, im1, imj, j, m, mm1, mmj, n, np1
      REAL(dp)                                           :: a0(21), b0(21), bsum, c(21), d(21), &
                                                            dsum, f, h, h2, hn, j0, j1, r, r0, r1, &
                                                            s, sum, t, t0, t1, u, w, w0, z, z0, &
                                                            z2, zn, znm1

!------------------------
!     ****** NUM IS THE MAXIMUM VALUE THAT N CAN TAKE IN THE DO LOOP
!            ENDING AT STATEMENT 50. IT IS REQUIRED THAT NUM BE EVEN.
!            THE ARRAYS A0, B0, C, D HAVE DIMENSION NUM + 1.
!------------------------
!     E0 = 2/SQRT(PI)
!     E1 = 2**(-3/2)
!------------------------

      fn_val = 0.0e0_dp
      IF (a >= b) THEN
         h = b/a
         r0 = 1.0e0_dp/(1.0e0_dp + h)
         r1 = (b - a)/a
         w0 = 1.0e0_dp/SQRT(b*(1.0e0_dp + h))
      ELSE
         h = a/b
         r0 = 1.0e0_dp/(1.0e0_dp + h)
         r1 = (b - a)/b
         w0 = 1.0e0_dp/SQRT(a*(1.0e0_dp + h))
      END IF

20    f = a*rlog1(-lambda/a) + b*rlog1(lambda/b)
      t = EXP(-f)
      IF (t == 0.0e0_dp) RETURN
      z0 = SQRT(f)
      z = 0.5e0_dp*(z0/e1)
      z2 = f + f

      a0(1) = (2.0e0_dp/3.0e0_dp)*r1
      c(1) = -0.5e0_dp*a0(1)
      d(1) = -c(1)
      j0 = (0.5e0_dp/e0)*ERFC_SCALED(z0)
      j1 = e1
      sum = j0 + d(1)*w0*j1

      s = 1.0e0_dp
      h2 = h*h
      hn = 1.0e0_dp
      w = w0
      znm1 = z
      zn = z2
      DO n = 2, num, 2
         hn = h2*hn
         a0(n) = 2.0e0_dp*r0*(1.0e0_dp + h*hn)/(n + 2.0e0_dp)
         np1 = n + 1
         s = s + hn
         a0(np1) = 2.0e0_dp*r1*s/(n + 3.0e0_dp)

         DO i = n, np1
            r = -0.5e0_dp*(i + 1.0e0_dp)
            b0(1) = r*a0(1)
            DO m = 2, i
               bsum = 0.0e0_dp
               mm1 = m - 1
               DO j = 1, mm1
                  mmj = m - j
                  bsum = bsum + (j*r - mmj)*a0(j)*b0(mmj)
               END DO
               b0(m) = r*a0(m) + bsum/m
            END DO
            c(i) = b0(i)/(i + 1.0e0_dp)

            dsum = 0.0e0_dp
            im1 = i - 1
            DO j = 1, im1
               imj = i - j
               dsum = dsum + d(imj)*c(j)
            END DO
            d(i) = -(dsum + c(i))
         END DO

         j0 = e1*znm1 + (n - 1.0e0_dp)*j0
         j1 = e1*zn + n*j1
         znm1 = z2*znm1
         zn = z2*zn
         w = w0*w
         t0 = d(n)*w*j0
         w = w0*w
         t1 = d(np1)*w*j1
         sum = sum + (t0 + t1)
         IF ((ABS(t0) + ABS(t1)) <= eps*sum) EXIT
      END DO

      u = EXP(-bcorr(a, b))
      fn_val = e0*t*u*sum

   END FUNCTION basym

END MODULE beta_gamma_psi
