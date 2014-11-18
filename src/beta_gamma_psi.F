!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
MODULE beta_gamma_psi
 ! not tested in the case where dp would stand for single precision
  
  USE kinds,                           ONLY: dp
#include "./common/cp_common_uses.f90"

 IMPLICIT NONE
 PRIVATE

 PUBLIC :: gamln,erf,psi

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param i ...
!> \retval fn_val ...
! *****************************************************************************
FUNCTION ipmpar (i) RESULT(fn_val)
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

    INTEGER, INTENT(IN)                      :: i
    INTEGER                                  :: fn_val

SELECT CASE(i)
  CASE( 1)
    fn_val = RADIX(i)
  CASE( 2)
    fn_val = DIGITS(i)
  CASE( 3)
    fn_val = HUGE(i)
  CASE( 4)
    fn_val = RADIX(1.0)
  CASE( 5)
    fn_val = DIGITS(1.0)
  CASE( 6)
    fn_val = MINEXPONENT(1.0)
  CASE( 7)
    fn_val = MAXEXPONENT(1.0)
  CASE( 8)
    fn_val = DIGITS(1.0e0_dp)
  CASE( 9)
    fn_val = MINEXPONENT(1.0e0_dp)
  CASE(10)
    fn_val = MAXEXPONENT(1.0e0_dp)
  CASE DEFAULT
    STOP "unknown case"
END SELECT

RETURN
END FUNCTION ipmpar

! *****************************************************************************
!> \brief ...
!> \param i ...
!> \retval fn_val ...
! *****************************************************************************
FUNCTION dpmpar (i) RESULT(fn_val)
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

    INTEGER, INTENT(IN)                      :: i
    REAL(dp)                                 :: fn_val

    REAL(dp), PARAMETER                      :: one = 1._dp

! Local variable

SELECT CASE (i)
  CASE (1)
    fn_val = EPSILON(one)
  CASE (2)
    fn_val = TINY(one)
  CASE (3)
    fn_val = HUGE(one)
END SELECT

RETURN
END FUNCTION dpmpar

! *****************************************************************************
!> \brief ...
!> \retval fn_val ...
! *****************************************************************************
FUNCTION epsln () RESULT(fn_val)
!--------------------------------------------------------------------
!     THE EVALUATION OF LN(EPS) WHERE EPS IS THE SMALLEST NUMBER
!     SUCH THAT 1.0 + EPS .GT. 1.0 .  L IS A DUMMY ARGUMENT.
!--------------------------------------------------------------------
    REAL                                     :: fn_val

    REAL, PARAMETER                          :: one = 1.0

! Local variable

fn_val = LOG( EPSILON(one) )
RETURN
END FUNCTION epsln

! *****************************************************************************
!> \brief ...
!> \param l ...
!> \retval fn_val ...
! *****************************************************************************
FUNCTION exparg (l) RESULT(fn_val)
!--------------------------------------------------------------------
!     IF L = 0 THEN  EXPARG(L) = THE LARGEST POSITIVE W FOR WHICH
!     EXP(W) CAN BE COMPUTED.
!
!     IF L IS NONZERO THEN  EXPARG(L) = THE LARGEST NEGATIVE W FOR
!     WHICH THE COMPUTED VALUE OF EXP(W) IS NONZERO.
!
!     NOTE... ONLY AN APPROXIMATE VALUE FOR EXPARG(L) IS NEEDED.
!--------------------------------------------------------------------
    INTEGER, INTENT(IN)                      :: l
    REAL                                     :: fn_val

    REAL, PARAMETER                          :: one = 1.0

! Local variable

IF (l == 0) THEN
  fn_val = LOG( HUGE(one) )
ELSE
  fn_val = LOG( TINY(one) )
END IF
RETURN
END FUNCTION exparg

! *****************************************************************************
!> \brief ...
!> \retval fn_val ...
! *****************************************************************************
FUNCTION depsln () RESULT(fn_val)
!--------------------------------------------------------------------
!     THE EVALUATION OF LN(EPS) WHERE EPS IS THE SMALLEST NUMBER
!     SUCH THAT 1.e0_dp + EPS .GT. 1.e0_dp .  L IS A DUMMY ARGUMENT.
!--------------------------------------------------------------------
    REAL(dp)                                 :: fn_val

    REAL(dp), PARAMETER                      :: one = 1._dp

! Local variable

fn_val = LOG( EPSILON(one) )
RETURN
END FUNCTION depsln

! *****************************************************************************
!> \brief ...
!> \param l ...
!> \retval fn_val ...
! *****************************************************************************
FUNCTION dxparg (l) RESULT(fn_val)
!--------------------------------------------------------------------
!     IF L = 0 THEN  DXPARG(L) = THE LARGEST POSITIVE W FOR WHICH
!     DEXP(W) CAN BE COMPUTED.
!
!     IF L IS NONZERO THEN  DXPARG(L) = THE LARGEST NEGATIVE W FOR
!     WHICH THE COMPUTED VALUE OF DEXP(W) IS NONZERO.
!
!     NOTE... ONLY AN APPROXIMATE VALUE FOR DXPARG(L) IS NEEDED.
!--------------------------------------------------------------------
    INTEGER, INTENT(IN)                      :: l
    REAL(dp)                                 :: fn_val

    REAL(dp), PARAMETER                      :: one = 1._dp

! Local variable

IF (l == 0) THEN
  fn_val = LOG( HUGE(one) )
ELSE
  fn_val = LOG( TINY(one) )
END IF
RETURN
END FUNCTION dxparg

! *****************************************************************************
!> \brief ...
!> \param a ...
!> \retval fn_val ...
! *****************************************************************************
FUNCTION alnrel(a) RESULT(fn_val)
!-----------------------------------------------------------------------
!            EVALUATION OF THE FUNCTION LN(1 + A)
!-----------------------------------------------------------------------
    REAL(dp), INTENT(IN)                     :: a
    REAL(dp)                                 :: fn_val

    REAL(dp), PARAMETER :: half = 0.5e0_dp, one = 1.e0_dp, &
      p1 = -.129418923021993e+01_dp, p2 = .405303492862024e+00_dp, &
      p3 = -.178874546012214e-01_dp, q1 = -.162752256355323e+01_dp, &
      q2 = .747811014037616e+00_dp, q3 = -.845104217945565e-01_dp, &
      two = 2.e0_dp, zero = 0.e0_dp

    REAL(dp)                                 :: t, t2, w, x

!--------------------------

IF (ABS(a) <= 0.375e0_dp) THEN
  t = a/(a + two)
  t2 = t*t
  w = (((p3*t2 + p2)*t2 + p1)*t2 + one)/ (((q3*t2 + q2)*t2 + q1)*t2 + one)
  fn_val = two*t*w
ELSE
  x = one + a
  IF (a < zero) x = (a + half) + half
  fn_val = LOG(x)
END IF

RETURN
END FUNCTION alnrel

! *****************************************************************************
!> \brief ...
!> \param x ...
!> \retval fn_val ...
! *****************************************************************************
FUNCTION erf (x) RESULT(fn_val)
!-----------------------------------------------------------------------
!             EVALUATION OF THE REAL ERROR FUNCTION
!-----------------------------------------------------------------------
    REAL(dp), INTENT(IN)                     :: x
    REAL(dp)                                 :: fn_val

    REAL(dp), PARAMETER :: a(5) = (/ .771058495001320e-04_dp, &
      -.133733772997339e-02_dp,  .323076579225834e-01_dp,  &
      .479137145607681e-01_dp, .128379167095513e+00_dp /), b(3) = (/ &
      .301048631703895e-02_dp,  .538971687740286e-01_dp,&
      .375795757275549e+00_dp /), c = .564189583547756e0_dp, p(8) = (/&
      -1.36864857382717e-07_dp,  5.64195517478974e-01_dp, &
      7.21175825088309e+00_dp,  4.31622272220567e+01_dp, &
      1.52989285046940e+02_dp,3.39320816734344e+02_dp, 4.51918953711873e+02_dp&
      ,  3.00459261020162e+02_dp /)
    REAL(dp), PARAMETER :: q(8) = (/ 1.00000000000000e+00_dp,  &
      1.27827273196294e+01_dp, 7.70001529352295e+01_dp,  &
      2.77585444743988e+02_dp,6.38980264465631e+02_dp,  &
      9.31354094850610e+02_dp, 7.90950925327898e+02_dp,  &
      3.00459260956983e+02_dp /), r(5) = (/ 2.10144126479064e+00_dp,  &
      2.62370141675169e+01_dp, 2.13688200555087e+01_dp,  &
      4.65807828718470e+00_dp, 2.82094791773523e-01_dp /), s(4) = (/ &
      9.41537750555460e+01_dp,  1.87114811799590e+02_dp, &
      9.90191814623914e+01_dp, 1.80124575948747e+01_dp /)

    REAL(dp)                                 :: ax, bot, t, top, x2

ax = ABS(x)
IF (ax < 0.5e0_dp) THEN
  t = x*x
  top = ((((a(1)*t + a(2))*t + a(3))*t + a(4))*t + a(5)) + 1.0e0_dp
  bot = ((b(1)*t + b(2))*t + b(3))*t + 1.0e0_dp
  fn_val = x*(top/bot)
  RETURN

ELSE IF (ax < 4.0e0_dp) THEN
  top = ((((((p(1)*ax + p(2))*ax + p(3))*ax + p(4))*ax + p(5))*ax  &
        + p(6))*ax + p(7))*ax + p(8)
  bot = ((((((q(1)*ax + q(2))*ax + q(3))*ax + q(4))*ax + q(5))*ax  &
        + q(6))*ax + q(7))*ax + q(8)
  fn_val = 0.5e0_dp + (0.5e0_dp - EXP(-x*x)*top/bot)
  IF (x < 0.0e0_dp) fn_val = -fn_val
  RETURN

ELSE IF (ax < 5.8e0_dp) THEN
  x2 = x*x
  t = 1.0e0_dp/x2
  top = (((r(1)*t + r(2))*t + r(3))*t + r(4))*t + r(5)
  bot = (((s(1)*t + s(2))*t + s(3))*t + s(4))*t + 1.0e0_dp
  fn_val = (c - top/(x2*bot)) / ax
  fn_val = 0.5e0_dp + (0.5e0_dp - EXP(-x2)*fn_val)
  IF (x < 0.0e0_dp) fn_val = -fn_val
  RETURN

ELSE
  fn_val = SIGN(1.0e0_dp,x)
END IF

RETURN
END FUNCTION erf

! *****************************************************************************
!> \brief ...
!> \param ind ...
!> \param x ...
!> \retval fn_val ...
! *****************************************************************************
FUNCTION erfc1 (ind, x) RESULT(fn_val)
!-----------------------------------------------------------------------
!         EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION

!          ERFC1(IND,X) = ERFC(X)            IF IND = 0
!          ERFC1(IND,X) = EXP(X*X)*ERFC(X)   OTHERWISE
!-----------------------------------------------------------------------
    INTEGER, INTENT(IN)                      :: ind
    REAL(dp), INTENT(IN)                     :: x
    REAL(dp)                                 :: fn_val

    REAL(dp), PARAMETER :: a(5) = (/ .771058495001320e-04_dp, &
      -.133733772997339e-02_dp,  .323076579225834e-01_dp,  &
      .479137145607681e-01_dp,  .128379167095513e+00_dp /), b(3) = (/ &
      .301048631703895e-02_dp,  .538971687740286e-01_dp,&
      .375795757275549e+00_dp /), c = .564189583547756e0_dp, p(8) = (/&
      -1.36864857382717e-07_dp,  5.64195517478974e-01_dp,  &
      7.21175825088309e+00_dp,  4.31622272220567e+01_dp,  &
      1.52989285046940e+02_dp,3.39320816734344e+02_dp,  &
      4.51918953711873e+02_dp,  3.00459261020162e+02_dp /)
    REAL(dp), PARAMETER :: q(8) = (/ 1.00000000000000e+00_dp,  &
      1.27827273196294e+01_dp,  7.70001529352295e+01_dp,  &
      2.77585444743988e+02_dp,6.38980264465631e+02_dp,  &
      9.31354094850610e+02_dp,  7.90950925327898e+02_dp,  &
      3.00459260956983e+02_dp /), r(5) = (/ 2.10144126479064e+00_dp,  &
      2.62370141675169e+01_dp,  2.13688200555087e+01_dp,  &
      4.65807828718470e+00_dp,  2.82094791773523e-01_dp /), s(4) = (/ &
      9.41537750555460e+01_dp,  1.87114811799590e+02_dp,  &
      9.90191814623914e+01_dp, 1.80124575948747e+01_dp /)

    REAL(dp)                                 :: ax, bot, e, t, top, w

ax = ABS(x)
IF (ax > 0.5e0_dp) GO TO 10
t = x*x
top = ((((a(1)*t + a(2))*t + a(3))*t + a(4))*t + a(5)) + 1.0e0_dp
bot = ((b(1)*t + b(2))*t + b(3))*t + 1.0e0_dp
fn_val = 0.5e0_dp + (0.5e0_dp - x*(top/bot))
IF (ind /= 0) fn_val = EXP(t) * fn_val
RETURN

!                  0.5 < ABS(X) <= 4

10 IF (ax > 4.0e0_dp) GO TO 20
top = ((((((p(1)*ax + p(2))*ax + p(3))*ax + p(4))*ax + p(5))*ax  &
      + p(6))*ax + p(7))*ax + p(8)
bot = ((((((q(1)*ax + q(2))*ax + q(3))*ax + q(4))*ax + q(5))*ax  &
      + q(6))*ax + q(7))*ax + q(8)
fn_val = top/bot
GO TO 40

!                      ABS(X) > 4

20 IF (x <= -5.6e0_dp) GO TO 50
IF (ind /= 0) GO TO 30
IF (x > 100.0e0_dp) GO TO 60
IF (x*x > -dxparg(1)) GO TO 60

30 t = (1.0e0_dp/x)**2
top = (((r(1)*t + r(2))*t + r(3))*t + r(4))*t + r(5)
bot = (((s(1)*t + s(2))*t + s(3))*t + s(4))*t + 1.0e0_dp
fn_val = (c - t*top/bot)/ax

!                      FINAL ASSEMBLY

40 IF (ind /= 0) THEN
  IF (x < 0.0e0_dp) fn_val = 2.0e0_dp*EXP(x*x) - fn_val
  RETURN
END IF
w = x * x
t = w
e = w - t
fn_val = ((0.5e0_dp + (0.5e0_dp - e)) * EXP(-t)) * fn_val
IF (x < 0.0e0_dp) fn_val = 2.0e0_dp - fn_val
RETURN

!             LIMIT VALUE FOR LARGE NEGATIVE X

50 fn_val = 2.0e0_dp
IF (ind /= 0) fn_val = 2.0e0_dp*EXP(x*x)
RETURN

!             LIMIT VALUE FOR LARGE POSITIVE X
!                       WHEN IND = 0

60 fn_val = 0.0e0_dp
RETURN
END FUNCTION erfc1

! *****************************************************************************
!> \brief ...
!> \param a ...
!> \retval fn_val ...
! *****************************************************************************
FUNCTION gam1(a) RESULT(fn_val)
!-----------------------------------------------------------------------
!     COMPUTATION OF 1/GAMMA(A+1) - 1  FOR -0.5 <= A <= 1.5
!-----------------------------------------------------------------------
    REAL(dp), INTENT(IN)                     :: a
    REAL(dp)                                 :: fn_val

    REAL(dp), PARAMETER :: p(7) = (/  .577215664901533e+00_dp, &
      -.409078193005776e+00_dp,-.230975380857675e+00_dp,  &
      .597275330452234e-01_dp, .766968181649490e-02_dp, &
      -.514889771323592e-02_dp, .589597428611429e-03_dp /), q(5) = (/  &
      .100000000000000e+01_dp,  .427569613095214e+00_dp, &
      .158451672430138e+00_dp,  .261132021441447e-01_dp, &
      .423244297896961e-02_dp /), r(9) = (/ -.422784335098468e+00_dp, &
      -.771330383816272e+00_dp, -.244757765222226e+00_dp,  &
      .118378989872749e+00_dp, .930357293360349e-03_dp,&
      -.118290993445146e-01_dp, .223047661158249e-02_dp,  &
      .266505979058923e-03_dp, -.132674909766242e-03_dp /)
    REAL(dp), PARAMETER :: s1 = .273076135303957e+00_dp, &
      s2 = .559398236957378e-01_dp

    REAL(dp)                                 :: bot, d, t, top, w

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
  top = (((((((r(9)*t + r(8))*t + r(7))*t + r(6))*t + r(5))*t  &
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

RETURN
END FUNCTION gam1

! *****************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \retval fn_val ...
! *****************************************************************************
FUNCTION algdiv (a, b) RESULT(fn_val)
!-----------------------------------------------------------------------

!     COMPUTATION OF LN(GAMMA(B)/GAMMA(A+B)) WHEN B >= 8

!                         --------

!     IN THIS ALGORITHM, DEL(X) IS THE FUNCTION DEFINED BY
!     LN(GAMMA(X)) = (X - 0.5)*LN(X) - X + 0.5*LN(2*PI) + DEL(X).

!-----------------------------------------------------------------------
    REAL(dp), INTENT(IN)                     :: a, b
    REAL(dp)                                 :: fn_val

    REAL(dp), PARAMETER :: c0 = .833333333333333e-01_dp, &
      c1 = -.277777777760991e-02_dp, c2 = .793650666825390e-03_dp, &
      c3 = -.595202931351870e-03_dp, c4 = .837308034031215e-03_dp, &
      c5 = -.165322962780713e-02_dp

    REAL(dp)                                 :: c, d, h, s11, s3, s5, s7, s9, &
                                                t, u, v, w, x, x2

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

RETURN
END FUNCTION algdiv

! *****************************************************************************
!> \brief ...
!> \param x ...
!> \retval fn_val ...
! *****************************************************************************
FUNCTION rexp (x) RESULT(fn_val)
!-----------------------------------------------------------------------
!            EVALUATION OF THE FUNCTION EXP(X) - 1
!-----------------------------------------------------------------------
    REAL(dp), INTENT(IN)                     :: x
    REAL(dp)                                 :: fn_val

    REAL(dp), PARAMETER :: p1 = .914041914819518e-09_dp, &
      p2 = .238082361044469e-01_dp, q1 = -.499999999085958e+00_dp, &
      q2 = .107141568980644e+00_dp, q3 = -.119041179760821e-01_dp, &
      q4 = .595130811860248e-03_dp

    REAL(dp)                                 :: e

IF (ABS(x) < 0.15e0_dp) THEN
  fn_val = x*(((p2*x + p1)*x + 1.0e0_dp)/((((q4*x + q3)*x + q2)*x + q1)*x + 1.0e0_dp))
  RETURN
END IF

IF (x < 0.0e0_dp) GO TO 20
e = EXP(x)
fn_val = e*(0.5e0_dp + (0.5e0_dp - 1.0e0_dp/e))
RETURN

20 IF (x > -37.0e0_dp) THEN
  fn_val = (EXP(x) - 0.5e0_dp) - 0.5e0_dp
  RETURN
END IF

fn_val = -1.0e0_dp

RETURN
END FUNCTION rexp

! *****************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \param x ...
!> \param y ...
!> \param w ...
!> \param eps ...
!> \param ierr ...
! *****************************************************************************
SUBROUTINE bgrat (a, b, x, y, w, eps, ierr)
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR IX(A,B) WHEN A IS LARGER THAN B.
!     THE RESULT OF THE EXPANSION IS ADDED TO W. IT IS ASSUMED
!     THAT A <= 15 AND B <= 1.  EPS IS THE TOLERANCE USED.
!     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
!-----------------------------------------------------------------------
    REAL(dp), INTENT(IN)                     :: a, b, x, y
    REAL(dp), INTENT(INOUT)                  :: w
    REAL(dp), INTENT(IN)                     :: eps
    INTEGER, INTENT(OUT)                     :: ierr

    REAL(dp), PARAMETER                      :: half = 0.5e0_dp, &
                                                one = 1.e0_dp, &
                                                quarter = 0.25e0_dp, &
                                                zero = 0.e0_dp

    INTEGER                                  :: i, n
    REAL(dp)                                 :: bm1, bp2n, c(30), cn, coef, &
                                                d(30), dj, j, l, lnx, n2, nu, &
                                                q, r, s, sum, t, t2, tol, u, &
                                                v, z

bm1 = (b - half) - half
nu = a + half*bm1
IF (y > 0.375e0_dp) GO TO 10
lnx = alnrel(-y)
GO TO 11
10 lnx = LOG(x)
11 z = -nu*lnx
IF (b*z == zero) GO TO 100

!                 COMPUTATION OF THE EXPANSION
!                 SET R = EXP(-Z)*Z**B/GAMMA(B)

r = b*(one + gam1(b))*EXP(b*LOG(z))
r = r*EXP(a*lnx)*EXP(half*bm1*lnx)
u = algdiv(b,a) + b*LOG(nu)
u = r*EXP(-u)
IF (u == zero) GO TO 100
CALL grat1 (b, z, r, q=q, eps=eps)

tol = 15.0e0_dp*eps
v = quarter*(one/nu)**2
t2 = quarter*lnx*lnx
l = w/u
j = q/r
sum = j
t = one
cn = one
n2 = zero
DO n = 1,30
  bp2n = b + n2
  j = (bp2n*(bp2n + one)*j + (z + bp2n + one)*t)*v
  n2 = n2 + 2.0e0_dp
  t = t*t2
  cn = cn/(n2*(n2 + one))
  c(n) = cn
  s = zero
  IF (n == 1) GO TO 21
  coef = b - n
  DO i = 1, n-1
    s = s + coef*c(i)*d(n-i)
    coef = coef + b
  END DO
  21 d(n) = bm1*cn + s/n
  dj = d(n)*j
  sum = sum + dj
  IF (sum <= zero) GO TO 100
  IF (ABS(dj) <= tol*(sum + l)) GO TO 30
END DO

!                    ADD THE RESULTS TO W

30 ierr = 0
w = w + u*sum
RETURN

!               THE EXPANSION CANNOT BE COMPUTED

100 ierr = 1
RETURN
END SUBROUTINE bgrat

! *****************************************************************************
!> \brief ...
!> \param a ...
!> \param x ...
!> \param r ...
!> \param p ...
!> \param q ...
!> \param eps ...
! *****************************************************************************
SUBROUTINE grat1 (a, x, r, p, q, eps)
!-----------------------------------------------------------------------
!           EVALUATION OF P(A,X) AND Q(A,X) WHERE A <= 1 AND
!        THE INPUT ARGUMENT R HAS THE VALUE E**(-X)*X**A/GAMMA(A)
!-----------------------------------------------------------------------
    REAL(dp), INTENT(IN)                     :: a, x, r
    REAL(dp), INTENT(OUT), OPTIONAL          :: p, q
    REAL(dp), INTENT(IN)                     :: eps

    REAL(dp), PARAMETER :: half = 0.5e0_dp, one = 1.e0_dp, &
      quarter = 0.25e0_dp, three = 3.e0_dp, two = 2.e0_dp, zero = 0.e0_dp

    REAL(dp)                                 :: a2n, a2nm1, an, b2n, b2nm1, &
                                                c, g, h, j, l, pp, qq, sum, &
                                                t, tol, w, z

IF (a*x == zero) GO TO 130
IF (a == half) GO TO 120
IF (x < 1.1e0_dp) GO TO 10
GO TO 50

!             TAYLOR SERIES FOR P(A,X)/X**A

10 an = three
c = x
sum = x/(a + three)
tol = three*eps/(a + one)
11 an = an + one
c = -c*(x/an)
t = c/(a + an)
sum = sum + t
IF (ABS(t) > tol) GO TO 11
j = a*x*((sum/6.0e0_dp - half/(a + two))*x + one/(a + one))

z = a*LOG(x)
h = gam1(a)
g = one + h
IF (x < quarter) GO TO 20
IF (a < x/2.59e0_dp) GO TO 40
GO TO 30
20 IF (z > -.13394e0_dp) GO TO 40

30 w = EXP(z)
pp = w*g*(half + (half - j))
qq = half + (half - pp)
GO TO 500

40 l = rexp(z)
qq = ((half + (half + l))*j - l)*g - h
IF (qq <= zero) GO TO 110
pp = half + (half - qq)
GO TO 500

!              CONTINUED FRACTION EXPANSION

50 tol = 8.0e0_dp*eps
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
GO TO 500

!                SPECIAL CASES

100 pp = zero
qq = one
GO TO 500

110 pp = one
qq = zero
GO TO 500

120 IF (x < quarter) THEN
  pp = erf(SQRT(x))
  qq = half + (half - pp)
  GO TO 500
ELSE
  qq = erfc1(0,SQRT(x))
  pp = half + (half - qq)
  GO TO 500
END IF

130 IF (x <= a) GO TO 100
GO TO 110

500 IF (PRESENT(p)) p = pp
IF (PRESENT(q)) q = qq

RETURN
END SUBROUTINE grat1

! *****************************************************************************
!> \brief ...
!> \param mu ...
!> \param x ...
!> \retval fn_val ...
! *****************************************************************************
FUNCTION esum (mu, x) RESULT(fn_val)
!-----------------------------------------------------------------------
!                    EVALUATION OF EXP(MU + X)
!-----------------------------------------------------------------------
    INTEGER, INTENT(IN)                      :: mu
    REAL(dp), INTENT(IN)                     :: x
    REAL(dp)                                 :: fn_val

    REAL(dp)                                 :: w

IF (x > 0.0e0_dp) GO TO 10

IF (mu < 0) GO TO 20
w = mu + x
IF (w > 0.0e0_dp) GO TO 20
fn_val = EXP(w)
RETURN

10 IF (mu > 0) GO TO 20
w = mu + x
IF (w < 0.0e0_dp) GO TO 20
fn_val = EXP(w)
RETURN

20 w = mu
fn_val = EXP(w)*EXP(x)

RETURN
END FUNCTION esum

! *****************************************************************************
!> \brief ...
!> \param x ...
!> \retval fn_val ...
! *****************************************************************************
FUNCTION rlog1(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!             EVALUATION OF THE FUNCTION X - LN(1 + X)
!-----------------------------------------------------------------------
!     A = RLOG (0.7)
!     B = RLOG (4/3)
!------------------------
    REAL(dp), INTENT(IN)                     :: x
    REAL(dp)                                 :: fn_val

    REAL(dp), PARAMETER :: a = .566749439387324e-01_dp, &
      b = .456512608815524e-01_dp, p0 = .333333333333333e+00_dp, &
      p1 = -.224696413112536e+00_dp, p2 = .620886815375787e-02_dp, &
      q1 = -.127408923933623e+01_dp, q2 = .354508718369557e+00_dp

    REAL(dp)                                 :: r, t, u, up2, w, w1

IF (x < -0.39e0_dp .OR. x > 0.57e0_dp) GO TO 100
IF (x < -0.18e0_dp) GO TO 10
IF (x >  0.18e0_dp) GO TO 20

!                 ARGUMENT REDUCTION

u = x
up2 = u + 2.0e0_dp
w1 = 0.0e0_dp
GO TO 30

10 u = (x + 0.3e0_dp)/0.7e0_dp
up2 = u + 2.0e0_dp
w1 = a - u*0.3e0_dp
GO TO 30

20 t = 0.75e0_dp*x
u = t - 0.25e0_dp
up2 = t + 1.75e0_dp
w1 = b + u/3.0e0_dp

!                  SERIES EXPANSION

30 r = u/up2
t = r*r
w = ((p2*t + p1)*t + p0)/((q2*t + q1)*t + 1.0e0_dp)
fn_val = r*(u - 2.0e0_dp*t*w) + w1
RETURN

100 w = (x + 0.5e0_dp) + 0.5e0_dp
fn_val = x - LOG(w)
RETURN
END FUNCTION rlog1

! *****************************************************************************
!> \brief ...
!> \param a ...
!> \retval fn_val ...
! *****************************************************************************
FUNCTION gamln (a) RESULT(fn_val)
!-----------------------------------------------------------------------
!            EVALUATION OF LN(GAMMA(A)) FOR POSITIVE A
!-----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS
!          NAVAL SURFACE WARFARE CENTER
!          DAHLGREN, VIRGINIA
!--------------------------
!     D = 0.5*(LN(2*PI) - 1)
!--------------------------
    REAL(dp), INTENT(IN)                     :: a
    REAL(dp)                                 :: fn_val

    REAL(dp), PARAMETER :: c0 = .833333333333333e-01_dp, &
      c1 = -.277777777760991e-02_dp, c2 = .793650666825390e-03_dp, &
      c3 = -.595202931351870e-03_dp, c4 = .837308034031215e-03_dp, &
      c5 = -.165322962780713e-02_dp, d = .418938533204673e0_dp

    INTEGER                                  :: i, n
    REAL(dp)                                 :: t, w

IF (a > 0.8e0_dp) GO TO 10
fn_val = gamln1(a) - LOG(a)
RETURN
10 IF (a > 2.25e0_dp) GO TO 20
t = (a - 0.5e0_dp) - 0.5e0_dp
fn_val = gamln1(t)
RETURN

20 IF (a >= 10.0e0_dp) GO TO 30
n = a - 1.25e0_dp
t = a
w = 1.0e0_dp
DO i = 1, n
  t = t - 1.0e0_dp
  w = t*w
END DO
fn_val = gamln1(t - 1.0e0_dp) + LOG(w)
RETURN

30 t = (1.0e0_dp/a)**2
w = (((((c5*t + c4)*t + c3)*t + c2)*t + c1)*t + c0)/a
fn_val = (d + w) + (a - 0.5e0_dp)*(LOG(a) - 1.0e0_dp)

RETURN
END FUNCTION gamln

! *****************************************************************************
!> \brief ...
!> \param a ...
!> \retval fn_val ...
! *****************************************************************************
FUNCTION gamln1 (a) RESULT(fn_val)
!-----------------------------------------------------------------------
!     EVALUATION OF LN(GAMMA(1 + A)) FOR -0.2 .LE. A .LE. 1.25
!-----------------------------------------------------------------------
    REAL(dp), INTENT(IN)                     :: a
    REAL(dp)                                 :: fn_val

    REAL(dp), PARAMETER :: p0 = .577215664901533e+00_dp, &
      p1 = .844203922187225e+00_dp, p2 = -.168860593646662e+00_dp, &
      p3 = -.780427615533591e+00_dp, p4 = -.402055799310489e+00_dp, &
      p5 = -.673562214325671e-01_dp, p6 = -.271935708322958e-02_dp, &
      q1 = .288743195473681e+01_dp, q2 = .312755088914843e+01_dp, &
      q3 = .156875193295039e+01_dp, q4 = .361951990101499e+00_dp, &
      q5 = .325038868253937e-01_dp, q6 = .667465618796164e-03_dp, &
      r0 = .422784335098467e+00_dp, r1 = .848044614534529e+00_dp, &
      r2 = .565221050691933e+00_dp, r3 = .156513060486551e+00_dp, &
      r4 = .170502484022650e-01_dp, r5 = .497958207639485e-03_dp
    REAL(dp), PARAMETER :: s1 = .124313399877507e+01_dp, &
      s2 = .548042109832463e+00_dp, s3 = .101552187439830e+00_dp, &
      s4 = .713309612391000e-02_dp, s5 = .116165475989616e-03_dp

    REAL(dp)                                 :: w, x

IF (a >= 0.6e0_dp) GO TO 10
w = ((((((p6*a + p5)*a + p4)*a + p3)*a + p2)*a + p1)*a + p0)/  &
    ((((((q6*a + q5)*a + q4)*a + q3)*a + q2)*a + q1)*a + 1.0e0_dp)
fn_val = -a*w
RETURN

10 x = (a - 0.5e0_dp) - 0.5e0_dp
w = (((((r5*x + r4)*x + r3)*x + r2)*x + r1)*x + r0)/  &
    (((((s5*x + s4)*x + s3)*x + s2)*x + s1)*x + 1.0e0_dp)
fn_val = x*w
RETURN
END FUNCTION gamln1

! *****************************************************************************
!> \brief ...
!> \param xx ...
!> \retval fn_val ...
! *****************************************************************************
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
    REAL(dp), INTENT(IN)                     :: xx
    REAL(dp)                                 :: fn_val

    REAL(dp), PARAMETER :: dx0 = 1.461632144968362341262659542325721325e0_dp, &
      p1(7) = (/ .895385022981970e-02_dp, .477762828042627e+01_dp,  &
      .142441585084029e+03_dp, .118645200713425e+04_dp,  &
      .363351846806499e+04_dp,.413810161269013e+04_dp,  &
      .130560269827897e+04_dp /), p2(4) = (/ -.212940445131011e+01_dp, &
      -.701677227766759e+01_dp,  -.448616543918019e+01_dp,&
      -.648157123766197e+00_dp /), piov4 = .785398163397448e0_dp, q1(6) = (/ &
      .448452573429826e+02_dp, .520752771467162e+03_dp,  &
      .221000799247830e+04_dp, .364127349079381e+04_dp,  &
      .190831076596300e+04_dp,.691091682714533e-05_dp /)
    REAL(dp), PARAMETER :: q2(4) = (/  .322703493791143e+02_dp, &
      .892920700481861e+02_dp,  .546117738103215e+02_dp,  &
      .777788548522962e+01_dp /)

    INTEGER                                  :: i, m, n, nq
    REAL(dp)                                 :: aug, den, sgn, upper, w, x, &
                                                xmax1, xmx0, xsmall, z

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
IF (x >= 0.5e0_dp) GO TO 200
!---------------------------------------------------------------------
!     X .LT. 0.5,  USE REFLECTION FORMULA
!     PSI(1-X) = PSI(X) + PI * COTAN(PI*X)
!---------------------------------------------------------------------
IF (ABS(x) > xsmall) GO TO 100
IF (x == 0.0e0_dp) GO TO 400
!---------------------------------------------------------------------
!     0 .LT. ABS(X) .LE. XSMALL.  USE 1/X AS A SUBSTITUTE
!     FOR  PI*COTAN(PI*X)
!---------------------------------------------------------------------
aug = -1.0e0_dp / x
GO TO 150
!---------------------------------------------------------------------
!     REDUCTION OF ARGUMENT FOR COTAN
!---------------------------------------------------------------------
100 w = - x
sgn = piov4
IF (w > 0.0e0_dp) GO TO 120
w = - w
sgn = -sgn
!---------------------------------------------------------------------
!     MAKE AN ERROR EXIT IF X .LE. -XMAX1
!---------------------------------------------------------------------
120 IF (w >= xmax1) GO TO 400
nq = INT(w)
w = w - nq
nq = INT(w*4.0e0_dp)
w = 4.0e0_dp * (w - nq * .25e0_dp)
!---------------------------------------------------------------------
!     W IS NOW RELATED TO THE FRACTIONAL PART OF  4.0 * X.
!     ADJUST ARGUMENT TO CORRESPOND TO VALUES IN FIRST
!     QUADRANT AND DETERMINE SIGN
!---------------------------------------------------------------------
n = nq / 2
IF ((n+n) /= nq) w = 1.0e0_dp - w
z = piov4 * w
m = n / 2
IF ((m+m) /= n) sgn = - sgn
!---------------------------------------------------------------------
!     DETERMINE FINAL VALUE FOR  -PI*COTAN(PI*X)
!---------------------------------------------------------------------
n = (nq + 1) / 2
m = n / 2
m = m + m
IF (m /= n) GO TO 140
!---------------------------------------------------------------------
!     CHECK FOR SINGULARITY
!---------------------------------------------------------------------
IF (z == 0.0e0_dp) GO TO 400
!---------------------------------------------------------------------
!     USE COS/SIN AS A SUBSTITUTE FOR COTAN, AND
!     SIN/COS AS A SUBSTITUTE FOR TAN
!---------------------------------------------------------------------
aug = sgn * ((COS(z) / SIN(z)) * 4.0e0_dp)
GO TO 150
140 aug = sgn * ((SIN(z) / COS(z)) * 4.0e0_dp)
150 x = 1.0e0_dp - x
200 IF (x > 3.0e0_dp) GO TO 300
!---------------------------------------------------------------------
!     0.5 .LE. X .LE. 3.0
!---------------------------------------------------------------------
den = x
upper = p1(1) * x

DO i = 1, 5
  den = (den + q1(i)) * x
  upper = (upper + p1(i+1)) * x
END DO

den = (upper + p1(7)) / (den + q1(6))
xmx0 = x - dx0
fn_val = den * xmx0 + aug
RETURN
!---------------------------------------------------------------------
!     IF X .GE. XMAX1, PSI = LN(X)
!---------------------------------------------------------------------
300 IF (x >= xmax1) GO TO 350
!---------------------------------------------------------------------
!     3.0 .LT. X .LT. XMAX1
!---------------------------------------------------------------------
w = 1.0e0_dp / (x * x)
den = w
upper = p2(1) * w

DO i = 1, 3
  den = (den + q2(i)) * w
  upper = (upper + p2(i+1)) * w
END DO

aug = upper / (den + q2(4)) - 0.5e0_dp / x + aug
350 fn_val = aug + LOG(x)
RETURN
!---------------------------------------------------------------------
!     ERROR RETURN
!---------------------------------------------------------------------
400 fn_val = 0.0e0_dp
RETURN
END FUNCTION psi

! *****************************************************************************
!> \brief ...
!> \param a0 ...
!> \param b0 ...
!> \retval fn_val ...
! *****************************************************************************
FUNCTION betaln (a0, b0) RESULT(fn_val)
!-----------------------------------------------------------------------
!     EVALUATION OF THE LOGARITHM OF THE BETA FUNCTION
!-----------------------------------------------------------------------
!     E = 0.5*LN(2*PI)
!--------------------------
    REAL(dp), INTENT(IN)                     :: a0, b0
    REAL(dp)                                 :: fn_val

    REAL(dp), PARAMETER                      :: e = .918938533204673e0_dp

    INTEGER                                  :: i, n
    REAL(dp)                                 :: a, b, c, h, u, v, w, z

!--------------------------

a = MIN(a0,b0)
b = MAX(a0,b0)
IF (a >= 8.0e0_dp) GO TO 60
IF (a >= 1.0e0_dp) GO TO 20
!-----------------------------------------------------------------------
!                   PROCEDURE WHEN A .LT. 1
!-----------------------------------------------------------------------
IF (b >= 8.0e0_dp) GO TO 10
fn_val = gamln(a) + (gamln(b) - gamln(a + b))
RETURN
10 fn_val = gamln(a) + algdiv(a,b)
RETURN
!-----------------------------------------------------------------------
!                PROCEDURE WHEN 1 .LE. A .LT. 8
!-----------------------------------------------------------------------
20 IF (a > 2.0e0_dp) GO TO 30
IF (b > 2.0e0_dp) GO TO 21
fn_val = gamln(a) + gamln(b) - gsumln(a,b)
RETURN
21 w = 0.0e0_dp
IF (b < 8.0e0_dp) GO TO 40
fn_val = gamln(a) + algdiv(a,b)
RETURN

!                REDUCTION OF A WHEN B .LE. 1000

30 IF (b > 1000.0e0_dp) GO TO 50
n = a - 1.0e0_dp
w = 1.0e0_dp
DO i = 1, n
  a = a - 1.0e0_dp
  h = a/b
  w = w * (h/(1.0e0_dp + h))
END DO
w = LOG(w)
IF (b < 8.0e0_dp) GO TO 40
fn_val = w + gamln(a) + algdiv(a,b)
RETURN

!                 REDUCTION OF B WHEN B .LT. 8

40 n = b - 1.0e0_dp
z = 1.0e0_dp
DO i = 1,n
  b = b - 1.0e0_dp
  z = z * (b/(a + b))
END DO
fn_val = w + LOG(z) + (gamln(a) + (gamln(b) - gsumln(a,b)))
RETURN

!                REDUCTION OF A WHEN B .GT. 1000

50 n = a - 1.0e0_dp
w = 1.0e0_dp
DO i = 1,n
  a = a - 1.0e0_dp
  w = w * (a/(1.0e0_dp + a/b))
END DO
fn_val = (LOG(w) - n*LOG(b)) + (gamln(a) + algdiv(a,b))
RETURN
!-----------------------------------------------------------------------
!                   PROCEDURE WHEN A .GE. 8
!-----------------------------------------------------------------------
60 w = bcorr(a,b)
h = a/b
c = h/(1.0e0_dp + h)
u = -(a - 0.5e0_dp)*LOG(c)
v = b*alnrel(h)
IF (u <= v) GO TO 61
fn_val = (((-0.5e0_dp*LOG(b) + e) + w) - v) - u
RETURN
61 fn_val = (((-0.5e0_dp*LOG(b) + e) + w) - u) - v
RETURN
END FUNCTION betaln

! *****************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \retval fn_val ...
! *****************************************************************************
FUNCTION gsumln (a, b) RESULT(fn_val)
!-----------------------------------------------------------------------
!          EVALUATION OF THE FUNCTION LN(GAMMA(A + B))
!          FOR 1 .LE. A .LE. 2  AND  1 .LE. B .LE. 2
!-----------------------------------------------------------------------
    REAL(dp), INTENT(IN)                     :: a, b
    REAL(dp)                                 :: fn_val

    REAL(dp)                                 :: x

x = a + b - 2.e0_dp
IF (x > 0.25e0_dp) GO TO 10
fn_val = gamln1(1.0e0_dp + x)
RETURN
10 IF (x > 1.25e0_dp) GO TO 20
fn_val = gamln1(x) + alnrel(x)
RETURN
20 fn_val = gamln1(x - 1.0e0_dp) + LOG(x*(1.0e0_dp + x))
RETURN
END FUNCTION gsumln

! *****************************************************************************
!> \brief ...
!> \param a0 ...
!> \param b0 ...
!> \retval fn_val ...
! *****************************************************************************
FUNCTION bcorr (a0, b0) RESULT(fn_val)
!-----------------------------------------------------------------------

!     EVALUATION OF  DEL(A0) + DEL(B0) - DEL(A0 + B0)  WHERE
!     LN(GAMMA(A)) = (A - 0.5)*LN(A) - A + 0.5*LN(2*PI) + DEL(A).
!     IT IS ASSUMED THAT A0 .GE. 8 AND B0 .GE. 8.

!-----------------------------------------------------------------------
    REAL(dp), INTENT(IN)                     :: a0, b0
    REAL(dp)                                 :: fn_val

    REAL(dp), PARAMETER :: c0 = .833333333333333e-01_dp, &
      c1 = -.277777777760991e-02_dp, c2 = .793650666825390e-03_dp, &
      c3 = -.595202931351870e-03_dp, c4 = .837308034031215e-03_dp, &
      c5 = -.165322962780713e-02_dp

    REAL(dp)                                 :: a, b, c, h, s11, s3, s5, s7, &
                                                s9, t, w, x, x2

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

! *****************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \param x ...
!> \param y ...
!> \param w ...
!> \param w1 ...
!> \param ierr ...
! *****************************************************************************
SUBROUTINE bratio (a, b, x, y, w, w1, ierr)
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
    REAL(dp), INTENT(IN)                     :: a, b, x, y
    REAL(dp), INTENT(OUT)                    :: w, w1
    INTEGER, INTENT(OUT)                     :: ierr

    INTEGER                                  :: ierr1, ind, n
    REAL(dp)                                 :: a0, b0, eps, lambda, t, x0, &
                                                y0, z

!-----------------------------------------------------------------------
!     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE SMALLEST
!            FLOATING POINT NUMBER FOR WHICH 1.0 + EPS .GT. 1.0

eps = dpmpar(1)

w = 0.0e0_dp
w1 = 0.0e0_dp
IF (a < 0.0e0_dp .OR. b < 0.0e0_dp) GO TO 300
IF (a == 0.0e0_dp .AND. b == 0.0e0_dp) GO TO 310
IF (x < 0.0e0_dp .OR. x > 1.0e0_dp) GO TO 320
IF (y < 0.0e0_dp .OR. y > 1.0e0_dp) GO TO 330
z = ((x + y) - 0.5e0_dp) - 0.5e0_dp
IF (ABS(z) > 3.0e0_dp*eps) GO TO 340

ierr = 0
IF (x == 0.0e0_dp) GO TO 200
IF (y == 0.0e0_dp) GO TO 210
IF (a == 0.0e0_dp) GO TO 211
IF (b == 0.0e0_dp) GO TO 201

eps = MAX(eps, 1.e-15_dp)
IF (MAX(a,b) < 1.e-3_dp*eps) GO TO 230

ind = 0
a0 = a
b0 = b
x0 = x
y0 = y
IF (MIN(a0, b0) > 1.0e0_dp) GO TO 30

!             PROCEDURE FOR A0 .LE. 1 OR B0 .LE. 1

IF (x <= 0.5e0_dp) GO TO 10
ind = 1
a0 = b
b0 = a
x0 = y
y0 = x

10 IF (b0 < MIN(eps,eps*a0)) GO TO 80
IF (a0 < MIN(eps,eps*b0) .AND. b0*x0 <= 1.0e0_dp) GO TO 90
IF (MAX(a0, b0) > 1.0e0_dp) GO TO 20
IF (a0 >= MIN(0.2e0_dp, b0)) GO TO 100
IF (x0**a0 <= 0.9e0_dp) GO TO 100
IF (x0 >= 0.3e0_dp) GO TO 110
n = 20
GO TO 130

20 IF (b0 <= 1.0e0_dp) GO TO 100
IF (x0 >= 0.3e0_dp) GO TO 110
IF (x0 >= 0.1e0_dp) GO TO 21
IF ((x0*b0)**a0 <= 0.7e0_dp) GO TO 100
21 IF (b0 > 15.0e0_dp) GO TO 131
n = 20
GO TO 130

!             PROCEDURE FOR A0 .GT. 1 AND B0 .GT. 1

30 IF (a > b) GO TO 31
lambda = a - (a + b)*x
GO TO 32
31 lambda = (a + b)*y - b
32 IF (lambda >= 0.0e0_dp) GO TO 40
ind = 1
a0 = b
b0 = a
x0 = y
y0 = x
lambda = ABS(lambda)

40 IF (b0 < 40.0e0_dp .AND. b0*x0 <= 0.7e0_dp) GO TO 100
IF (b0 < 40.0e0_dp) GO TO 140
IF (a0 > b0) GO TO 50
IF (a0 <= 100.0e0_dp) GO TO 120
IF (lambda > 0.03e0_dp*a0) GO TO 120
GO TO 180
50 IF (b0 <= 100.0e0_dp) GO TO 120
IF (lambda > 0.03e0_dp*b0) GO TO 120
GO TO 180

!            EVALUATION OF THE APPROPRIATE ALGORITHM

80 w = fpser(a0, b0, x0, eps)
w1 = 0.5e0_dp + (0.5e0_dp - w)
GO TO 220

90 w1 = apser(a0, b0, x0, eps)
w = 0.5e0_dp + (0.5e0_dp - w1)
GO TO 220

100 w = bpser(a0, b0, x0, eps)
w1 = 0.5e0_dp + (0.5e0_dp - w)
GO TO 220

110 w1 = bpser(b0, a0, y0, eps)
w = 0.5e0_dp + (0.5e0_dp - w1)
GO TO 220

120 w = bfrac(a0, b0, x0, y0, lambda, 15.0e0_dp*eps)
w1 = 0.5e0_dp + (0.5e0_dp - w)
GO TO 220

130 w1 = bup(b0, a0, y0, x0, n, eps)
b0 = b0 + n
131 CALL bgrat (b0, a0, y0, x0, w1, eps, ierr1)
IF (ierr1 > 0) STOP "Error in BGRAT"
w = 0.5e0_dp + (0.5e0_dp - w1)
GO TO 220

140 n = b0
b0 = b0 - n
IF (b0 /= 0.0e0_dp) GO TO 141
n = n - 1
b0 = 1.0e0_dp
141 w = bup(b0, a0, y0, x0, n, eps)
IF (x0 > 0.7e0_dp) GO TO 150
w = w + bpser(a0, b0, x0, eps)
w1 = 0.5e0_dp + (0.5e0_dp - w)
GO TO 220

150 IF (a0 > 15.0e0_dp) GO TO 151
n = 20
w = w + bup(a0, b0, x0, y0, n, eps)
a0 = a0 + n
151 CALL bgrat (a0, b0, x0, y0, w, eps, ierr1)
w1 = 0.5e0_dp + (0.5e0_dp - w)
GO TO 220

180 w = basym(a0, b0, lambda, 100.0e0_dp*eps)
w1 = 0.5e0_dp + (0.5e0_dp - w)
GO TO 220

!               TERMINATION OF THE PROCEDURE

200 IF (a == 0.0e0_dp) GO TO 350
201 w = 0.0e0_dp
w1 = 1.0e0_dp
RETURN

210 IF (b == 0.0e0_dp) GO TO 360
211 w = 1.0e0_dp
w1 = 0.0e0_dp
RETURN

220 IF (ind == 0) RETURN
t = w
w = w1
w1 = t
RETURN

!           PROCEDURE FOR A AND B .LT. 1.E-3*EPS

230 w = b/(a + b)
w1 = a/(a + b)
RETURN

!                       ERROR RETURN

300 ierr = 1
RETURN
310 ierr = 2
RETURN
320 ierr = 3
RETURN
330 ierr = 4
RETURN
340 ierr = 5
RETURN
350 ierr = 6
RETURN
360 ierr = 7
RETURN
END SUBROUTINE bratio

! *****************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \param x ...
!> \param eps ...
!> \retval fn_val ...
! *****************************************************************************
FUNCTION fpser (a, b, x, eps) RESULT(fn_val)
!-----------------------------------------------------------------------

!                 EVALUATION OF I (A,B)
!                                X

!          FOR B .LT. MIN(EPS,EPS*A) AND X .LE. 0.5.

!-----------------------------------------------------------------------
    REAL(dp), INTENT(IN)                     :: a, b, x, eps
    REAL(dp)                                 :: fn_val

    REAL(dp)                                 :: an, c, s, t, tol

!                  SET  FPSER = X**A

fn_val = 1.0e0_dp
IF (a <= 1.e-3_dp*eps) GO TO 10
fn_val = 0.0e0_dp
t = a*LOG(x)
IF (t < dxparg(1)) RETURN
fn_val = EXP(t)

!                NOTE THAT 1/B(A,B) = B

10 fn_val = (b/a)*fn_val
tol = eps/a
an = a + 1.0e0_dp
t = x
s = t/an
20    an = an + 1.0e0_dp
t = x*t
c = t/an
s = s + c
IF (ABS(c) > tol) GO TO 20

fn_val = fn_val*(1.0e0_dp + a*s)
RETURN
END FUNCTION fpser

! *****************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \param x ...
!> \param eps ...
!> \retval fn_val ...
! *****************************************************************************
FUNCTION apser (a, b, x, eps) RESULT(fn_val)
!-----------------------------------------------------------------------
!     APSER YIELDS THE INCOMPLETE BETA RATIO I(SUB(1-X))(B,A) FOR
!     A .LE. MIN(EPS,EPS*B), B*X .LE. 1, AND X .LE. 0.5. USED WHEN
!     A IS VERY SMALL. USE ONLY IF ABOVE INEQUALITIES ARE SATISFIED.
!-----------------------------------------------------------------------
    REAL(dp), INTENT(IN)                     :: a, b, x, eps
    REAL(dp)                                 :: fn_val

    REAL(dp), PARAMETER                      :: g = .577215664901533e0_dp

    REAL(dp)                                 :: aj, bx, c, j, s, t, tol

bx = b*x
t = x - bx
IF (b*eps > 2.e-2_dp) GO TO 10
c = LOG(x) + psi(b) + g + t
GO TO 20
10 c = LOG(bx) + g + t

20 tol = 5.0e0_dp*eps*ABS(c)
j = 1.0e0_dp
s = 0.0e0_dp
30    j = j + 1.0e0_dp
t = t*(x - bx/j)
aj = t/j
s = s + aj
IF (ABS(aj) > tol) GO TO 30

fn_val = -a*(c + s)
RETURN
END FUNCTION apser

! *****************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \param x ...
!> \param eps ...
!> \retval fn_val ...
! *****************************************************************************
FUNCTION bpser (a, b, x, eps) RESULT(fn_val)
!-----------------------------------------------------------------------
!     POWER SERIES EXPANSION FOR EVALUATING IX(A,B) WHEN B .LE. 1
!     OR B*X .LE. 0.7.  EPS IS THE TOLERANCE USED.
!-----------------------------------------------------------------------
    REAL(dp), INTENT(IN)                     :: a, b, x, eps
    REAL(dp)                                 :: fn_val

    INTEGER                                  :: i, m
    REAL(dp)                                 :: a0, apb, b0, c, n, sum, t, &
                                                tol, u, w, z

fn_val = 0.0e0_dp
IF (x == 0.0e0_dp) RETURN
!-----------------------------------------------------------------------
!            COMPUTE THE FACTOR X**A/(A*BETA(A,B))
!-----------------------------------------------------------------------
a0 = MIN(a,b)
IF (a0 < 1.0e0_dp) GO TO 10
z = a*LOG(x) - betaln(a,b)
fn_val = EXP(z)/a
GO TO 70
10 b0 = MAX(a,b)
IF (b0 >= 8.0e0_dp) GO TO 60
IF (b0 > 1.0e0_dp) GO TO 40

!            PROCEDURE FOR A0 .LT. 1 AND B0 .LE. 1

fn_val = x**a
IF (fn_val == 0.0e0_dp) RETURN

apb = a + b
IF (apb > 1.0e0_dp) GO TO 20
z = 1.0e0_dp + gam1(apb)
GO TO 30
20 u = a + b - 1.e0_dp
z = (1.0e0_dp + gam1(u))/apb

30 c = (1.0e0_dp + gam1(a))*(1.0e0_dp + gam1(b))/z
fn_val = fn_val*c*(b/apb)
GO TO 70

!         PROCEDURE FOR A0 .LT. 1 AND 1 .LT. B0 .LT. 8

40 u = gamln1(a0)
m = b0 - 1.0e0_dp
IF (m < 1) GO TO 50
c = 1.0e0_dp
DO i = 1, m
  b0 = b0 - 1.0e0_dp
  c = c*(b0/(a0 + b0))
END DO
u = LOG(c) + u

50 z = a*LOG(x) - u
b0 = b0 - 1.0e0_dp
apb = a0 + b0
IF (apb > 1.0e0_dp) GO TO 51
t = 1.0e0_dp + gam1(apb)
GO TO 52
51 u = a0 + b0 - 1.e0_dp
t = (1.0e0_dp + gam1(u))/apb
52 fn_val = EXP(z)*(a0/a)*(1.0e0_dp + gam1(b0))/t
GO TO 70

!            PROCEDURE FOR A0 .LT. 1 AND B0 .GE. 8

60 u = gamln1(a0) + algdiv(a0,b0)
z = a*LOG(x) - u
fn_val = (a0/a)*EXP(z)
70 IF (fn_val == 0.0e0_dp .OR. a <= 0.1e0_dp*eps) RETURN
!-----------------------------------------------------------------------
!                     COMPUTE THE SERIES
!-----------------------------------------------------------------------
sum = 0.0e0_dp
n = 0.0e0_dp
c = 1.0e0_dp
tol = eps/a
100    n = n + 1.0e0_dp
c = c*(0.5e0_dp + (0.5e0_dp - b/n))*x
w = c/(a + n)
sum = sum + w
IF (ABS(w) > tol) GO TO 100
fn_val = fn_val*(1.0e0_dp + a*sum)
RETURN
END FUNCTION bpser

! *****************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \param x ...
!> \param y ...
!> \param n ...
!> \param eps ...
!> \retval fn_val ...
! *****************************************************************************
FUNCTION bup (a, b, x, y, n, eps) RESULT(fn_val)
!-----------------------------------------------------------------------
!     EVALUATION OF IX(A,B) - IX(A+N,B) WHERE N IS A POSITIVE INTEGER.
!     EPS IS THE TOLERANCE USED.
!-----------------------------------------------------------------------
    REAL(dp), INTENT(IN)                     :: a, b, x, y
    INTEGER, INTENT(IN)                      :: n
    REAL(dp), INTENT(IN)                     :: eps
    REAL(dp)                                 :: fn_val

    INTEGER                                  :: i, k, kp1, mu, nm1
    REAL(dp)                                 :: ap1, apb, d, l, r, t, w

!          OBTAIN THE SCALING FACTOR EXP(-MU) AND
!             EXP(MU)*(X**A*Y**B/BETA(A,B))/A

apb = a + b
ap1 = a + 1.0e0_dp
mu = 0
d = 1.0e0_dp
IF (n == 1 .OR. a < 1.0e0_dp) GO TO 10
IF (apb < 1.1e0_dp*ap1) GO TO 10
mu = ABS(dxparg(1))
k = dxparg(0)
IF (k < mu) mu = k
t = mu
d = EXP(-t)

10 fn_val = brcmp1(mu,a,b,x,y)/a
IF (n == 1 .OR. fn_val == 0.0e0_dp) RETURN
nm1 = n - 1
w = d

!          LET K BE THE INDEX OF THE MAXIMUM TERM

k = 0
IF (b <= 1.0e0_dp) GO TO 40
IF (y > 1.e-4_dp) GO TO 20
k = nm1
GO TO 30
20 r = (b - 1.0e0_dp)*x/y - a
IF (r < 1.0e0_dp) GO TO 40
k = nm1
t = nm1
IF (r < t) k = r

!          ADD THE INCREASING TERMS OF THE SERIES

30 DO i = 1,k
  l = i - 1
  d = ((apb + l)/(ap1 + l))*x*d
  w = w + d
END DO
IF (k == nm1) GO TO 50

!          ADD THE REMAINING TERMS OF THE SERIES

40 kp1 = k + 1
DO i = kp1,nm1
  l = i - 1
  d = ((apb + l)/(ap1 + l))*x*d
  w = w + d
  IF (d <= eps*w) GO TO 50
END DO

!               TERMINATE THE PROCEDURE

50 fn_val = fn_val*w
RETURN
END FUNCTION bup

! *****************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \param x ...
!> \param y ...
!> \param lambda ...
!> \param eps ...
!> \retval fn_val ...
! *****************************************************************************
FUNCTION bfrac (a, b, x, y, lambda, eps) RESULT(fn_val)
!-----------------------------------------------------------------------
!     CONTINUED FRACTION EXPANSION FOR IX(A,B) WHEN A,B .GT. 1.
!     IT IS ASSUMED THAT  LAMBDA = (A + B)*Y - B.
!-----------------------------------------------------------------------
    REAL(dp), INTENT(IN)                     :: a, b, x, y, lambda, eps
    REAL(dp)                                 :: fn_val

    REAL(dp)                                 :: alpha, an, anp1, beta, bn, &
                                                bnp1, c, c0, c1, e, n, p, r, &
                                                r0, s, t, w, yp1

fn_val = brcomp(a,b,x,y)
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

10    n = n + 1.0e0_dp
t = n/a
w = n*(b - n)*x
e = a/s
alpha = (p*(p + c0)*e*e)*(w*x)
IF (alpha <= 0.0e0_dp) GO TO 20
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
IF (ABS(r - r0) <= eps*r) GO TO 20

!        RESCALE AN, BN, ANP1, AND BNP1

an = an/bnp1
bn = bn/bnp1
anp1 = r
bnp1 = 1.0e0_dp
GO TO 10

!                 TERMINATION

20 fn_val = fn_val*r
RETURN
END FUNCTION bfrac

! *****************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \param x ...
!> \param y ...
!> \retval fn_val ...
! *****************************************************************************
FUNCTION brcomp (a, b, x, y) RESULT(fn_val)
!-----------------------------------------------------------------------
!               EVALUATION OF X**A*Y**B/BETA(A,B)
!-----------------------------------------------------------------------
    REAL(dp), INTENT(IN)                     :: a, b, x, y
    REAL(dp)                                 :: fn_val

    REAL(dp), PARAMETER                      :: const = 0.398942280401433e0_dp

    INTEGER                                  :: i, n
    REAL(dp)                                 :: a0, apb, b0, c, e, h, lambda, &
                                                lnx, lny, t, u, v, x0, y0, z

!-----------------
!     CONST = 1/SQRT(2*PI)
!-----------------

fn_val = 0.0e0_dp
IF (x == 0.0e0_dp .OR. y == 0.0e0_dp) RETURN
a0 = MIN(a,b)
IF (a0 >= 8.0e0_dp) GO TO 100

IF (x > 0.375e0_dp) GO TO 10
lnx = LOG(x)
lny = alnrel(-x)
GO TO 20
10 IF (y > 0.375e0_dp) GO TO 11
lnx = alnrel(-y)
lny = LOG(y)
GO TO 20
11 lnx = LOG(x)
lny = LOG(y)

20 z = a*lnx + b*lny
IF (a0 < 1.0e0_dp) GO TO 30
z = z - betaln(a,b)
fn_val = EXP(z)
RETURN
!-----------------------------------------------------------------------
!              PROCEDURE FOR A .LT. 1 OR B .LT. 1
!-----------------------------------------------------------------------
30 b0 = MAX(a,b)
IF (b0 >= 8.0e0_dp) GO TO 80
IF (b0 > 1.0e0_dp) GO TO 60

!                   ALGORITHM FOR B0 .LE. 1

fn_val = EXP(z)
IF (fn_val == 0.0e0_dp) RETURN

apb = a + b
IF (apb > 1.0e0_dp) GO TO 40
z = 1.0e0_dp + gam1(apb)
GO TO 50
40 u = a + b - 1.e0_dp
z = (1.0e0_dp + gam1(u))/apb

50 c = (1.0e0_dp + gam1(a))*(1.0e0_dp + gam1(b))/z
fn_val = fn_val*(a0*c)/(1.0e0_dp + a0/b0)
RETURN

!                ALGORITHM FOR 1 .LT. B0 .LT. 8

60 u = gamln1(a0)
n = b0 - 1.0e0_dp
IF (n < 1) GO TO 70
c = 1.0e0_dp
DO i = 1, n
  b0 = b0 - 1.0e0_dp
  c = c*(b0/(a0 + b0))
END DO
u = LOG(c) + u

70 z = z - u
b0 = b0 - 1.0e0_dp
apb = a0 + b0
IF (apb > 1.0e0_dp) GO TO 71
t = 1.0e0_dp + gam1(apb)
GO TO 72
71 u = a0 + b0 - 1.e0_dp
t = (1.0e0_dp + gam1(u))/apb
72 fn_val = a0*EXP(z)*(1.0e0_dp + gam1(b0))/t
RETURN

!                   ALGORITHM FOR B0 .GE. 8

80 u = gamln1(a0) + algdiv(a0,b0)
fn_val = a0*EXP(z - u)
RETURN
!-----------------------------------------------------------------------
!              PROCEDURE FOR A .GE. 8 AND B .GE. 8
!-----------------------------------------------------------------------
100 IF (a > b) GO TO 101
h = a/b
x0 = h/(1.0e0_dp + h)
y0 = 1.0e0_dp/(1.0e0_dp + h)
lambda = a - (a + b)*x
GO TO 110
101 h = b/a
x0 = 1.0e0_dp/(1.0e0_dp + h)
y0 = h/(1.0e0_dp + h)
lambda = (a + b)*y - b

110 e = -lambda/a
IF (ABS(e) > 0.6e0_dp) GO TO 111
u = rlog1(e)
GO TO 120
111 u = e - LOG(x/x0)

120 e = lambda/b
IF (ABS(e) > 0.6e0_dp) GO TO 121
v = rlog1(e)
GO TO 130
121 v = e - LOG(y/y0)

130 z = EXP(-(a*u + b*v))
fn_val = const*SQRT(b*x0)*z*EXP(-bcorr(a,b))
RETURN
END FUNCTION brcomp

! *****************************************************************************
!> \brief ...
!> \param mu ...
!> \param a ...
!> \param b ...
!> \param x ...
!> \param y ...
!> \retval fn_val ...
! *****************************************************************************
FUNCTION brcmp1 (mu, a, b, x, y) RESULT(fn_val)
!-----------------------------------------------------------------------
!          EVALUATION OF  EXP(MU) * (X**A*Y**B/BETA(A,B))
!-----------------------------------------------------------------------
    INTEGER, INTENT(IN)                      :: mu
    REAL(dp), INTENT(IN)                     :: a, b, x, y
    REAL(dp)                                 :: fn_val

    REAL(dp), PARAMETER                      :: const = 0.398942280401433e0_dp

    INTEGER                                  :: i, n
    REAL(dp)                                 :: a0, apb, b0, c, e, h, lambda, &
                                                lnx, lny, t, u, v, x0, y0, z

!-----------------
!     CONST = 1/SQRT(2*PI)
!-----------------

a0 = MIN(a,b)
IF (a0 >= 8.0e0_dp) GO TO 100

IF (x > 0.375e0_dp) GO TO 10
lnx = LOG(x)
lny = alnrel(-x)
GO TO 20
10 IF (y > 0.375e0_dp) GO TO 11
lnx = alnrel(-y)
lny = LOG(y)
GO TO 20
11 lnx = LOG(x)
lny = LOG(y)

20 z = a*lnx + b*lny
IF (a0 < 1.0e0_dp) GO TO 30
z = z - betaln(a,b)
fn_val = esum(mu,z)
RETURN
!-----------------------------------------------------------------------
!              PROCEDURE FOR A .LT. 1 OR B .LT. 1
!-----------------------------------------------------------------------
30 b0 = MAX(a,b)
IF (b0 >= 8.0e0_dp) GO TO 80
IF (b0 > 1.0e0_dp) GO TO 60

!                   ALGORITHM FOR B0 .LE. 1

fn_val = esum(mu,z)
IF (fn_val == 0.0e0_dp) RETURN

apb = a + b
IF (apb > 1.0e0_dp) GO TO 40
z = 1.0e0_dp + gam1(apb)
GO TO 50
40 u = a + b - 1.e0_dp
z = (1.0e0_dp + gam1(u))/apb

50 c = (1.0e0_dp + gam1(a))*(1.0e0_dp + gam1(b))/z
fn_val = fn_val*(a0*c)/(1.0e0_dp + a0/b0)
RETURN

!                ALGORITHM FOR 1 .LT. B0 .LT. 8

60 u = gamln1(a0)
n = b0 - 1.0e0_dp
IF (n < 1) GO TO 70
c = 1.0e0_dp
DO i = 1, n
  b0 = b0 - 1.0e0_dp
  c = c*(b0/(a0 + b0))
END DO
u = LOG(c) + u

70 z = z - u
b0 = b0 - 1.0e0_dp
apb = a0 + b0
IF (apb > 1.0e0_dp) GO TO 71
t = 1.0e0_dp + gam1(apb)
GO TO 72
71 u = a0 + b0 - 1.e0_dp
t = (1.0e0_dp + gam1(u))/apb
72 fn_val = a0*esum(mu,z)*(1.0e0_dp + gam1(b0))/t
RETURN

!                   ALGORITHM FOR B0 .GE. 8

80 u = gamln1(a0) + algdiv(a0,b0)
fn_val = a0*esum(mu,z - u)
RETURN
!-----------------------------------------------------------------------
!              PROCEDURE FOR A .GE. 8 AND B .GE. 8
!-----------------------------------------------------------------------
100 IF (a > b) GO TO 101
h = a/b
x0 = h/(1.0e0_dp + h)
y0 = 1.0e0_dp/(1.0e0_dp + h)
lambda = a - (a + b)*x
GO TO 110
101 h = b/a
x0 = 1.0e0_dp/(1.0e0_dp + h)
y0 = h/(1.0e0_dp + h)
lambda = (a + b)*y - b

110 e = -lambda/a
IF (ABS(e) > 0.6e0_dp) GO TO 111
u = rlog1(e)
GO TO 120
111 u = e - LOG(x/x0)

120 e = lambda/b
IF (ABS(e) > 0.6e0_dp) GO TO 121
v = rlog1(e)
GO TO 130
121 v = e - LOG(y/y0)

130 z = esum(mu,-(a*u + b*v))
fn_val = const*SQRT(b*x0)*z*EXP(-bcorr(a,b))
RETURN
END FUNCTION brcmp1

! *****************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \param lambda ...
!> \param eps ...
!> \retval fn_val ...
! *****************************************************************************
FUNCTION basym (a, b, lambda, eps) RESULT(fn_val)
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR IX(A,B) FOR LARGE A AND B.
!     LAMBDA = (A + B)*Y - B  AND EPS IS THE TOLERANCE USED.
!     IT IS ASSUMED THAT LAMBDA IS NONNEGATIVE AND THAT
!     A AND B ARE GREATER THAN OR EQUAL TO 15.
!-----------------------------------------------------------------------
    REAL(dp), INTENT(IN)                     :: a, b, lambda, eps
    REAL(dp)                                 :: fn_val

    INTEGER, PARAMETER                       :: num = 20
    REAL(dp), PARAMETER                      :: e0 = 1.12837916709551e0_dp, &
                                                e1 = .353553390593274e0_dp

    INTEGER                                  :: i, im1, imj, j, m, mm1, mmj, &
                                                n, np1
    REAL(dp)                                 :: a0(21), b0(21), bsum, c(21), &
                                                d(21), dsum, f, h, h2, hn, &
                                                j0, j1, r, r0, r1, s, sum, t, &
                                                t0, t1, u, w, w0, z, z0, z2, &
                                                zn, znm1

!------------------------
!     ****** NUM IS THE MAXIMUM VALUE THAT N CAN TAKE IN THE DO LOOP
!            ENDING AT STATEMENT 50. IT IS REQUIRED THAT NUM BE EVEN.
!            THE ARRAYS A0, B0, C, D HAVE DIMENSION NUM + 1.
!------------------------
!     E0 = 2/SQRT(PI)
!     E1 = 2**(-3/2)
!------------------------

fn_val = 0.0e0_dp
IF (a >= b) GO TO 10
h = a/b
r0 = 1.0e0_dp/(1.0e0_dp + h)
r1 = (b - a)/b
w0 = 1.0e0_dp/SQRT(a*(1.0e0_dp + h))
GO TO 20
10 h = b/a
r0 = 1.0e0_dp/(1.0e0_dp + h)
r1 = (b - a)/a
w0 = 1.0e0_dp/SQRT(b*(1.0e0_dp + h))

20 f = a*rlog1(-lambda/a) + b*rlog1(lambda/b)
t = EXP(-f)
IF (t == 0.0e0_dp) RETURN
z0 = SQRT(f)
z = 0.5e0_dp*(z0/e1)
z2 = f + f

a0(1) = (2.0e0_dp/3.0e0_dp)*r1
c(1) = - 0.5e0_dp*a0(1)
d(1) = - c(1)
j0 = (0.5e0_dp/e0)*erfc1(1,z0)
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
  IF ((ABS(t0) + ABS(t1)) <= eps*sum) GO TO 60
END DO

60 u = EXP(-bcorr(a,b))
fn_val = e0*t*u*sum
RETURN
END FUNCTION basym

END MODULE beta_gamma_psi
