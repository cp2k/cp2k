MODULE beta_gamma_psi
 ! not tested in the case where dp would stand for single precision
  USE kinds,                           ONLY: dp

 IMPLICIT NONE
 PRIVATE

 PUBLIC :: gamln,erf,psi

CONTAINS

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
    fn_val = DIGITS(1.0D0)
  CASE( 9)
    fn_val = MINEXPONENT(1.0D0)
  CASE(10)
    fn_val = MAXEXPONENT(1.0D0)
  CASE DEFAULT
    RETURN
END SELECT

RETURN
END FUNCTION ipmpar



FUNCTION spmpar (i) RESULT(fn_val)
!-----------------------------------------------------------------------

!     SPMPAR PROVIDES THE SINGLE PRECISION MACHINE CONSTANTS FOR
!     THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT
!     I IS AN INTEGER HAVING ONE OF THE VALUES 1, 2, OR 3. IF THE
!     SINGLE PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND
!     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX, THEN

!        SPMPAR(1) = B**(1 - M), THE MACHINE PRECISION,

!        SPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,

!        SPMPAR(3) = B**EMAX*(1 - B**(-M)), THE LARGEST MAGNITUDE.
!-----------------------------------------------------------------------

    INTEGER, INTENT(IN)                      :: i
    REAL                                     :: fn_val

    REAL                                     :: one = 1.0

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
END FUNCTION spmpar



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

    REAL(dp)                                 :: one = 1._dp

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


FUNCTION epsln () RESULT(fn_val)
!--------------------------------------------------------------------
!     THE EVALUATION OF LN(EPS) WHERE EPS IS THE SMALLEST NUMBER
!     SUCH THAT 1.0 + EPS .GT. 1.0 .  L IS A DUMMY ARGUMENT.
!--------------------------------------------------------------------
    REAL                                     :: fn_val

    REAL                                     :: one = 1.0

! Local variable

fn_val = LOG( EPSILON(one) )
RETURN
END FUNCTION epsln


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

    REAL                                     :: one = 1.0

! Local variable

IF (l == 0) THEN
  fn_val = LOG( HUGE(one) )
ELSE
  fn_val = LOG( TINY(one) )
END IF
RETURN
END FUNCTION exparg


FUNCTION depsln () RESULT(fn_val)
!--------------------------------------------------------------------
!     THE EVALUATION OF LN(EPS) WHERE EPS IS THE SMALLEST NUMBER
!     SUCH THAT 1.D0 + EPS .GT. 1.D0 .  L IS A DUMMY ARGUMENT.
!--------------------------------------------------------------------
    REAL(dp)                                 :: fn_val

    REAL(dp)                                 :: one = 1._dp

! Local variable

fn_val = LOG( EPSILON(one) )
RETURN
END FUNCTION depsln


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

    REAL(dp)                                 :: one = 1._dp

! Local variable

IF (l == 0) THEN
  fn_val = LOG( HUGE(one) )
ELSE
  fn_val = LOG( TINY(one) )
END IF
RETURN
END FUNCTION dxparg

FUNCTION alnrel(a) RESULT(fn_val)
!-----------------------------------------------------------------------
!            EVALUATION OF THE FUNCTION LN(1 + A)
!-----------------------------------------------------------------------
    REAL(dp), INTENT(IN)                     :: a
    REAL(dp)                                 :: fn_val

    REAL(dp) :: half = 0.5D0, one = 1.D0, p1 = -.129418923021993D+01, &
      p2 = .405303492862024D+00, p3 = -.178874546012214D-01, &
      q1 = -.162752256355323D+01, q2 = .747811014037616D+00, &
      q3 = -.845104217945565D-01, t, t2, two = 2.D0, w, x, zero = 0.D0

!--------------------------

IF (ABS(a) <= 0.375D0) THEN
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



FUNCTION erf (x) RESULT(fn_val)
!-----------------------------------------------------------------------
!             EVALUATION OF THE REAL ERROR FUNCTION
!-----------------------------------------------------------------------
    REAL(dp), INTENT(IN)                     :: x
    REAL(dp)                                 :: fn_val

    REAL(dp) :: a(5) = (/ .771058495001320D-04, -.133733772997339D-02,   &
      .323076579225834D-01,  .479137145607681D-01,   .128379167095513D+00 /), &
      ax, b(3) = (/ .301048631703895D-02,  .538971687740286D-01,   &
      .375795757275549D+00 /), bot, c = .564189583547756D0, p(8) = (/-&
      1.36864857382717D-07,  5.64195517478974D-01, 7.21175825088309D+00,  &
      4.31622272220567D+01, 1.52989285046940D+02,  3.39320816734344D+02, &
      4.51918953711873D+02,  3.00459261020162D+02 /), q(8) = (/ &
      1.00000000000000D+00,  1.27827273196294D+01, 7.70001529352295D+01,  &
      2.77585444743988D+02, 6.38980264465631D+02,  9.31354094850610D+02, &
      7.90950925327898D+02,  3.00459260956983D+02 /)
    REAL(dp) :: r(5) = (/ 2.10144126479064D+00,  2.62370141675169D+01, &
      2.13688200555087D+01,  4.65807828718470D+00, 2.82094791773523D-01 /), s(4) = &
      (/ 9.41537750555460D+01,  1.87114811799590D+02, 9.90191814623914D+01, &
      1.80124575948747D+01 /), t, top, x2

!-------------------------
!-------------------------

ax = ABS(x)
IF (ax < 0.5D0) THEN
  t = x*x
  top = ((((a(1)*t + a(2))*t + a(3))*t + a(4))*t + a(5)) + 1.0D0
  bot = ((b(1)*t + b(2))*t + b(3))*t + 1.0D0
  fn_val = x*(top/bot)
  RETURN

ELSE IF (ax < 4.0D0) THEN
  top = ((((((p(1)*ax + p(2))*ax + p(3))*ax + p(4))*ax + p(5))*ax  &
        + p(6))*ax + p(7))*ax + p(8)
  bot = ((((((q(1)*ax + q(2))*ax + q(3))*ax + q(4))*ax + q(5))*ax  &
        + q(6))*ax + q(7))*ax + q(8)
  fn_val = 0.5D0 + (0.5D0 - EXP(-x*x)*top/bot)
  IF (x < 0.0D0) fn_val = -fn_val
  RETURN

ELSE IF (ax < 5.8D0) THEN
  x2 = x*x
  t = 1.0D0/x2
  top = (((r(1)*t + r(2))*t + r(3))*t + r(4))*t + r(5)
  bot = (((s(1)*t + s(2))*t + s(3))*t + s(4))*t + 1.0D0
  fn_val = (c - top/(x2*bot)) / ax
  fn_val = 0.5D0 + (0.5D0 - EXP(-x2)*fn_val)
  IF (x < 0.0D0) fn_val = -fn_val
  RETURN

ELSE
  fn_val = SIGN(1.0D0,x)
END IF

RETURN
END FUNCTION erf



FUNCTION erfc1 (ind, x) RESULT(fn_val)
!-----------------------------------------------------------------------
!         EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION

!          ERFC1(IND,X) = ERFC(X)            IF IND = 0
!          ERFC1(IND,X) = EXP(X*X)*ERFC(X)   OTHERWISE
!-----------------------------------------------------------------------
    INTEGER, INTENT(IN)                      :: ind
    REAL(dp), INTENT(IN)                     :: x
    REAL(dp)                                 :: fn_val

    REAL(dp) :: a(5) = (/ .771058495001320D-04, -.133733772997339D-02,  &
      .323076579225834D-01,  .479137145607681D-01,  .128379167095513D+00 /), &
      ax, b(3) = (/ .301048631703895D-02,  .538971687740286D-01,  &
      .375795757275549D+00 /), bot, c = .564189583547756D0, e, p(8) = (/-&
      1.36864857382717D-07,  5.64195517478974D-01,  7.21175825088309D+00,  &
      4.31622272220567D+01,  1.52989285046940D+02,  3.39320816734344D+02,  &
      4.51918953711873D+02,  3.00459261020162D+02 /), q(8) = (/ &
      1.00000000000000D+00,  1.27827273196294D+01,  7.70001529352295D+01,  &
      2.77585444743988D+02,  6.38980264465631D+02,  9.31354094850610D+02,  &
      7.90950925327898D+02,  3.00459260956983D+02 /)
    REAL(dp) :: r(5) = (/ 2.10144126479064D+00,  2.62370141675169D+01,  &
      2.13688200555087D+01,  4.65807828718470D+00,  2.82094791773523D-01 /), s(4) = &
         (/ 9.41537750555460D+01,  1.87114811799590D+02,  9.90191814623914D+01, &
               1.80124575948747D+01 /), t, top, w

!-----------------------------------------------------------------------
!-------------------------
!-------------------------
!                     ABS(X) <= 0.5

ax = ABS(x)
IF (ax > 0.5D0) GO TO 10
t = x*x
top = ((((a(1)*t + a(2))*t + a(3))*t + a(4))*t + a(5)) + 1.0D0
bot = ((b(1)*t + b(2))*t + b(3))*t + 1.0D0
fn_val = 0.5D0 + (0.5D0 - x*(top/bot))
IF (ind /= 0) fn_val = EXP(t) * fn_val
RETURN

!                  0.5 < ABS(X) <= 4

10 IF (ax > 4.0D0) GO TO 20
top = ((((((p(1)*ax + p(2))*ax + p(3))*ax + p(4))*ax + p(5))*ax  &
      + p(6))*ax + p(7))*ax + p(8)
bot = ((((((q(1)*ax + q(2))*ax + q(3))*ax + q(4))*ax + q(5))*ax  &
      + q(6))*ax + q(7))*ax + q(8)
fn_val = top/bot
GO TO 40

!                      ABS(X) > 4

20 IF (x <= -5.6D0) GO TO 50
IF (ind /= 0) GO TO 30
IF (x > 100.0D0) GO TO 60
IF (x*x > -dxparg(1)) GO TO 60

30 t = (1.0D0/x)**2
top = (((r(1)*t + r(2))*t + r(3))*t + r(4))*t + r(5)
bot = (((s(1)*t + s(2))*t + s(3))*t + s(4))*t + 1.0D0
fn_val = (c - t*top/bot)/ax

!                      FINAL ASSEMBLY

40 IF (ind /= 0) THEN
  IF (x < 0.0D0) fn_val = 2.0D0*EXP(x*x) - fn_val
  RETURN
END IF
w = x * x
t = w
e = w - t
fn_val = ((0.5D0 + (0.5D0 - e)) * EXP(-t)) * fn_val
IF (x < 0.0D0) fn_val = 2.0D0 - fn_val
RETURN

!             LIMIT VALUE FOR LARGE NEGATIVE X

50 fn_val = 2.0D0
IF (ind /= 0) fn_val = 2.0D0*EXP(x*x)
RETURN

!             LIMIT VALUE FOR LARGE POSITIVE X
!                       WHEN IND = 0

60 fn_val = 0.0D0
RETURN
END FUNCTION erfc1



FUNCTION gam1(a) RESULT(fn_val)
!-----------------------------------------------------------------------
!     COMPUTATION OF 1/GAMMA(A+1) - 1  FOR -0.5 <= A <= 1.5
!-----------------------------------------------------------------------
    REAL(dp), INTENT(IN)                     :: a
    REAL(dp)                                 :: fn_val

    REAL(dp) :: bot, d, p(7) = (/  .577215664901533D+00, -.409078193005776D+00&
      , -.230975380857675D+00,  .597275330452234D-01, .766968181649490D-02, -&
      .514889771323592D-02, .589597428611429D-03 /), q(5) = (/  &
      .100000000000000D+01,  .427569613095214D+00, .158451672430138D+00,  &
      .261132021441447D-01, .423244297896961D-02 /), r(9) = (/ -&
      .422784335098468D+00, -.771330383816272D+00, -.244757765222226D+00,  &
      .118378989872749D+00, .930357293360349D-03, -.118290993445146D-01, &
      .223047661158249D-02,  .266505979058923D-03, -.132674909766242D-03 /), &
      s1 = .273076135303957D+00, s2 = .559398236957378D-01, t, top, w

!------------------------
!------------------------

t = a
d = a - 0.5D0
IF (d > 0.0D0) t = d - 0.5D0

IF (t > 0.D0) THEN
  top = (((((p(7)*t + p(6))*t + p(5))*t + p(4))*t + p(3))*t + p(2))*t + p(1)
  bot = (((q(5)*t + q(4))*t + q(3))*t + q(2))*t + 1.0D0
  w = top/bot
  IF (d > 0.0D0) THEN
    fn_val = (t/a)*((w - 0.5D0) - 0.5D0)
  ELSE
    fn_val = a*w
  END IF
ELSE IF (t < 0.D0) THEN
  top = (((((((r(9)*t + r(8))*t + r(7))*t + r(6))*t + r(5))*t  &
           + r(4))*t + r(3))*t + r(2))*t + r(1)
  bot = (s2*t + s1)*t + 1.0D0
  w = top/bot
  IF (d > 0.0D0) THEN
    fn_val = t*w/a
  ELSE
    fn_val = a*((w + 0.5D0) + 0.5D0)
  END IF
ELSE
  fn_val = 0.0D0
END IF

RETURN
END FUNCTION gam1


FUNCTION algdiv (a, b) RESULT(fn_val)
!-----------------------------------------------------------------------

!     COMPUTATION OF LN(GAMMA(B)/GAMMA(A+B)) WHEN B >= 8

!                         --------

!     IN THIS ALGORITHM, DEL(X) IS THE FUNCTION DEFINED BY
!     LN(GAMMA(X)) = (X - 0.5)*LN(X) - X + 0.5*LN(2*PI) + DEL(X).

!-----------------------------------------------------------------------
    REAL(dp), INTENT(IN)                     :: a, b
    REAL(dp)                                 :: fn_val

    REAL(dp) :: c, c0 = .833333333333333D-01, c1 = -.277777777760991D-02, &
      c2 = .793650666825390D-03, c3 = -.595202931351870D-03, &
      c4 = .837308034031215D-03, c5 = -.165322962780713D-02, d, h, s11, s3, &
      s5, s7, s9, t, u, v, w, x, x2

! EXTERNAL   alnrel
!------------------------

IF (a > b) THEN
  h = b/a
  c = 1.0D0/(1.0D0 + h)
  x = h/(1.0D0 + h)
  d = a + (b - 0.5D0)
ELSE
  h = a/b
  c = h/(1.0D0 + h)
  x = 1.0D0/(1.0D0 + h)
  d = b + (a - 0.5D0)
END IF

!                SET SN = (1 - X**N)/(1 - X)

x2 = x*x
s3 = 1.0D0 + (x + x2)
s5 = 1.0D0 + (x + x2*s3)
s7 = 1.0D0 + (x + x2*s5)
s9 = 1.0D0 + (x + x2*s7)
s11 = 1.0D0 + (x + x2*s9)

!                SET W = DEL(B) - DEL(A + B)

t = (1.0D0/b)**2
w = ((((c5*s11*t + c4*s9)*t + c3*s7)*t + c2*s5)*t + c1*s3)*t + c0
w = w*(c/b)

!                    COMBINE THE RESULTS

u = d*alnrel(a/b)
v = a*(LOG(b) - 1.0D0)
IF (u > v) THEN
  fn_val = (w - v) - u
ELSE
  fn_val = (w - u) - v
END IF

RETURN
END FUNCTION algdiv



FUNCTION rexp (x) RESULT(fn_val)
!-----------------------------------------------------------------------
!            EVALUATION OF THE FUNCTION EXP(X) - 1
!-----------------------------------------------------------------------
    REAL(dp), INTENT(IN)                     :: x
    REAL(dp)                                 :: fn_val

    REAL(dp) :: e, p1 = .914041914819518D-09, p2 = .238082361044469D-01, &
      q1 = -.499999999085958D+00, q2 = .107141568980644D+00, &
      q3 = -.119041179760821D-01, q4 = .595130811860248D-03

!-----------------------

IF (ABS(x) < 0.15D0) THEN
  fn_val = x*(((p2*x + p1)*x + 1.0D0)/((((q4*x + q3)*x + q2)*x + q1)*x + 1.0D0))
  RETURN
END IF

IF (x < 0.0D0) GO TO 20
e = EXP(x)
fn_val = e*(0.5D0 + (0.5D0 - 1.0D0/e))
RETURN

20 IF (x > -37.0D0) THEN
  fn_val = (EXP(x) - 0.5D0) - 0.5D0
  RETURN
END IF

fn_val = -1.0D0

RETURN
END FUNCTION rexp



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

    INTEGER                                  :: i, n
    REAL(dp) :: bm1, bp2n, c(30), cn, coef, d(30), dj, half = 0.5D0, j, l, &
      lnx, n2, nu, one = 1.D0, q, quarter = 0.25D0, r, s, sum, t, t2, tol, u, &
      v, z, zero = 0.D0

bm1 = (b - half) - half
nu = a + half*bm1
IF (y > 0.375D0) GO TO 10
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

tol = 15.0D0*eps
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
  n2 = n2 + 2.0D0
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



SUBROUTINE grat1 (a, x, r, p, q, eps)
!-----------------------------------------------------------------------
!           EVALUATION OF P(A,X) AND Q(A,X) WHERE A <= 1 AND
!        THE INPUT ARGUMENT R HAS THE VALUE E**(-X)*X**A/GAMMA(A)
!-----------------------------------------------------------------------
    REAL(dp), INTENT(IN)                     :: a, x, r
    REAL(dp), INTENT(OUT), OPTIONAL          :: p, q
    REAL(dp), INTENT(IN)                     :: eps

    REAL(dp) :: a2n, a2nm1, an, b2n, b2nm1, c, g, h, half = 0.5D0, j, l, &
      one = 1.D0, pp, qq, quarter = 0.25D0, sum, t, three = 3.D0, tol, &
      two = 2.D0, w, z, zero = 0.D0

IF (a*x == zero) GO TO 130
IF (a == half) GO TO 120
IF (x < 1.1D0) GO TO 10
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
j = a*x*((sum/6.0D0 - half/(a + two))*x + one/(a + one))

z = a*LOG(x)
h = gam1(a)
g = one + h
IF (x < quarter) GO TO 20
IF (a < x/2.59D0) GO TO 40
GO TO 30
20 IF (z > -.13394D0) GO TO 40

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

50 tol = 8.0D0*eps
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



FUNCTION esum (mu, x) RESULT(fn_val)
!-----------------------------------------------------------------------
!                    EVALUATION OF EXP(MU + X)
!-----------------------------------------------------------------------
    INTEGER, INTENT(IN)                      :: mu
    REAL(dp), INTENT(IN)                     :: x
    REAL(dp)                                 :: fn_val

    REAL(dp)                                 :: w

IF (x > 0.0D0) GO TO 10

IF (mu < 0) GO TO 20
w = mu + x
IF (w > 0.0D0) GO TO 20
fn_val = EXP(w)
RETURN

10 IF (mu > 0) GO TO 20
w = mu + x
IF (w < 0.0D0) GO TO 20
fn_val = EXP(w)
RETURN

20 w = mu
fn_val = EXP(w)*EXP(x)

RETURN
END FUNCTION esum



FUNCTION rlog1(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!             EVALUATION OF THE FUNCTION X - LN(1 + X)
!-----------------------------------------------------------------------
!     A = RLOG (0.7)
!     B = RLOG (4/3)
!------------------------
    REAL(dp), INTENT(IN)                     :: x
    REAL(dp)                                 :: fn_val

    REAL(dp) :: a = .566749439387324D-01, b = .456512608815524D-01, &
      p0 = .333333333333333D+00, p1 = -.224696413112536D+00, &
      p2 = .620886815375787D-02, q1 = -.127408923933623D+01, &
      q2 = .354508718369557D+00, r, t, u, up2, w, w1

!------------------------

IF (x < -0.39D0 .OR. x > 0.57D0) GO TO 100
IF (x < -0.18D0) GO TO 10
IF (x >  0.18D0) GO TO 20

!                 ARGUMENT REDUCTION

u = x
up2 = u + 2.0D0
w1 = 0.0D0
GO TO 30

10 u = (x + 0.3D0)/0.7D0
up2 = u + 2.0D0
w1 = a - u*0.3D0
GO TO 30

20 t = 0.75D0*x
u = t - 0.25D0
up2 = t + 1.75D0
w1 = b + u/3.0D0

!                  SERIES EXPANSION

30 r = u/up2
t = r*r
w = ((p2*t + p1)*t + p0)/((q2*t + q1)*t + 1.0D0)
fn_val = r*(u - 2.0D0*t*w) + w1
RETURN


100 w = (x + 0.5D0) + 0.5D0
fn_val = x - LOG(w)
RETURN
END FUNCTION rlog1



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

    INTEGER                                  :: i, n
    REAL(dp) :: c0 = .833333333333333D-01, c1 = -.277777777760991D-02, &
      c2 = .793650666825390D-03, c3 = -.595202931351870D-03, &
      c4 = .837308034031215D-03, c5 = -.165322962780713D-02, &
      d = .418938533204673D0, t, w

!--------------------------

IF (a > 0.8D0) GO TO 10
fn_val = gamln1(a) - LOG(a)
RETURN
10 IF (a > 2.25D0) GO TO 20
t = (a - 0.5D0) - 0.5D0
fn_val = gamln1(t)
RETURN

20 IF (a >= 10.0D0) GO TO 30
n = a - 1.25D0
t = a
w = 1.0D0
DO i = 1, n
  t = t - 1.0D0
  w = t*w
END DO
fn_val = gamln1(t - 1.0D0) + LOG(w)
RETURN

30 t = (1.0D0/a)**2
w = (((((c5*t + c4)*t + c3)*t + c2)*t + c1)*t + c0)/a
fn_val = (d + w) + (a - 0.5D0)*(LOG(a) - 1.0D0)

RETURN
END FUNCTION gamln



FUNCTION gamln1 (a) RESULT(fn_val)
!-----------------------------------------------------------------------
!     EVALUATION OF LN(GAMMA(1 + A)) FOR -0.2 .LE. A .LE. 1.25
!-----------------------------------------------------------------------
    REAL(dp), INTENT(IN)                     :: a
    REAL(dp)                                 :: fn_val

    REAL(dp) :: p0 = .577215664901533D+00, p1 = .844203922187225D+00, &
      p2 = -.168860593646662D+00, p3 = -.780427615533591D+00, &
      p4 = -.402055799310489D+00, p5 = -.673562214325671D-01, &
      p6 = -.271935708322958D-02, q1 = .288743195473681D+01, &
      q2 = .312755088914843D+01, q3 = .156875193295039D+01, &
      q4 = .361951990101499D+00, q5 = .325038868253937D-01, &
      q6 = .667465618796164D-03, r0 = .422784335098467D+00, &
      r1 = .848044614534529D+00, r2 = .565221050691933D+00, &
      r3 = .156513060486551D+00, r4 = .170502484022650D-01, &
      r5 = .497958207639485D-03, s1 = .124313399877507D+01, &
      s2 = .548042109832463D+00, s3 = .101552187439830D+00, &
      s4 = .713309612391000D-02
    REAL(dp) :: s5 = .116165475989616D-03, w, x

!----------------------

IF (a >= 0.6D0) GO TO 10
w = ((((((p6*a + p5)*a + p4)*a + p3)*a + p2)*a + p1)*a + p0)/  &
    ((((((q6*a + q5)*a + q4)*a + q3)*a + q2)*a + q1)*a + 1.0D0)
fn_val = -a*w
RETURN

10 x = (a - 0.5D0) - 0.5D0
w = (((((r5*x + r4)*x + r3)*x + r2)*x + r1)*x + r0)/  &
    (((((s5*x + s4)*x + s3)*x + s2)*x + s1)*x + 1.0D0)
fn_val = x*w
RETURN
END FUNCTION gamln1



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

    INTEGER                                  :: i, m, n, nq
    REAL(dp) :: aug, den, dx0 = 1.461632144968362341262659542325721325D0, p1(7&
      ) = (/ .895385022981970D-02, .477762828042627D+01,  .142441585084029D+03&
      , .118645200713425D+04,  .363351846806499D+04, .413810161269013D+04,  &
      .130560269827897D+04 /), p2(4) = (/ -.212940445131011D+01, -&
      .701677227766759D+01,  -.448616543918019D+01, -.648157123766197D+00 /), &
      piov4 = .785398163397448D0, q1(6) = (/ .448452573429826D+02, &
      .520752771467162D+03,  .221000799247830D+04, .364127349079381D+04,  &
      .190831076596300D+04, .691091682714533D-05 /), q2(4) = (/  &
      .322703493791143D+02,  .892920700481861D+02,  .546117738103215D+02,  &
      .777788548522962D+01 /)
    REAL(dp) :: sgn, upper, w, x, xmax1, xmx0, xsmall, z

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
xmax1 = MIN(xmax1, 1.0D0/dpmpar(1))
xsmall = 1.d-9
!---------------------------------------------------------------------
x = xx
aug = 0.0D0
IF (x >= 0.5D0) GO TO 200
!---------------------------------------------------------------------
!     X .LT. 0.5,  USE REFLECTION FORMULA
!     PSI(1-X) = PSI(X) + PI * COTAN(PI*X)
!---------------------------------------------------------------------
IF (ABS(x) > xsmall) GO TO 100
IF (x == 0.0D0) GO TO 400
!---------------------------------------------------------------------
!     0 .LT. ABS(X) .LE. XSMALL.  USE 1/X AS A SUBSTITUTE
!     FOR  PI*COTAN(PI*X)
!---------------------------------------------------------------------
aug = -1.0D0 / x
GO TO 150
!---------------------------------------------------------------------
!     REDUCTION OF ARGUMENT FOR COTAN
!---------------------------------------------------------------------
100 w = - x
sgn = piov4
IF (w > 0.0D0) GO TO 120
w = - w
sgn = -sgn
!---------------------------------------------------------------------
!     MAKE AN ERROR EXIT IF X .LE. -XMAX1
!---------------------------------------------------------------------
120 IF (w >= xmax1) GO TO 400
nq = INT(w)
w = w - nq
nq = INT(w*4.0D0)
w = 4.0D0 * (w - nq * .25D0)
!---------------------------------------------------------------------
!     W IS NOW RELATED TO THE FRACTIONAL PART OF  4.0 * X.
!     ADJUST ARGUMENT TO CORRESPOND TO VALUES IN FIRST
!     QUADRANT AND DETERMINE SIGN
!---------------------------------------------------------------------
n = nq / 2
IF ((n+n) /= nq) w = 1.0D0 - w
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
IF (z == 0.0D0) GO TO 400
!---------------------------------------------------------------------
!     USE COS/SIN AS A SUBSTITUTE FOR COTAN, AND
!     SIN/COS AS A SUBSTITUTE FOR TAN
!---------------------------------------------------------------------
aug = sgn * ((COS(z) / SIN(z)) * 4.0D0)
GO TO 150
140 aug = sgn * ((SIN(z) / COS(z)) * 4.0D0)
150 x = 1.0D0 - x
200 IF (x > 3.0D0) GO TO 300
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
w = 1.0D0 / (x * x)
den = w
upper = p2(1) * w

DO i = 1, 3
  den = (den + q2(i)) * w
  upper = (upper + p2(i+1)) * w
END DO

aug = upper / (den + q2(4)) - 0.5D0 / x + aug
350 fn_val = aug + LOG(x)
RETURN
!---------------------------------------------------------------------
!     ERROR RETURN
!---------------------------------------------------------------------
400 fn_val = 0.0D0
RETURN
END FUNCTION psi



FUNCTION betaln (a0, b0) RESULT(fn_val)
!-----------------------------------------------------------------------
!     EVALUATION OF THE LOGARITHM OF THE BETA FUNCTION
!-----------------------------------------------------------------------
!     E = 0.5*LN(2*PI)
!--------------------------
    REAL(dp), INTENT(IN)                     :: a0, b0
    REAL(dp)                                 :: fn_val

    INTEGER                                  :: i, n
    REAL(dp)                                 :: a, b, c, &
                                                e = .918938533204673D0, h, u, &
                                                v, w, z

!--------------------------

a = MIN(a0,b0)
b = MAX(a0,b0)
IF (a >= 8.0D0) GO TO 60
IF (a >= 1.0D0) GO TO 20
!-----------------------------------------------------------------------
!                   PROCEDURE WHEN A .LT. 1
!-----------------------------------------------------------------------
IF (b >= 8.0D0) GO TO 10
fn_val = gamln(a) + (gamln(b) - gamln(a + b))
RETURN
10 fn_val = gamln(a) + algdiv(a,b)
RETURN
!-----------------------------------------------------------------------
!                PROCEDURE WHEN 1 .LE. A .LT. 8
!-----------------------------------------------------------------------
20 IF (a > 2.0D0) GO TO 30
IF (b > 2.0D0) GO TO 21
fn_val = gamln(a) + gamln(b) - gsumln(a,b)
RETURN
21 w = 0.0D0
IF (b < 8.0D0) GO TO 40
fn_val = gamln(a) + algdiv(a,b)
RETURN

!                REDUCTION OF A WHEN B .LE. 1000

30 IF (b > 1000.0D0) GO TO 50
n = a - 1.0D0
w = 1.0D0
DO i = 1, n
  a = a - 1.0D0
  h = a/b
  w = w * (h/(1.0D0 + h))
END DO
w = LOG(w)
IF (b < 8.0D0) GO TO 40
fn_val = w + gamln(a) + algdiv(a,b)
RETURN

!                 REDUCTION OF B WHEN B .LT. 8

40 n = b - 1.0D0
z = 1.0D0
DO i = 1,n
  b = b - 1.0D0
  z = z * (b/(a + b))
END DO
fn_val = w + LOG(z) + (gamln(a) + (gamln(b) - gsumln(a,b)))
RETURN

!                REDUCTION OF A WHEN B .GT. 1000

50 n = a - 1.0D0
w = 1.0D0
DO i = 1,n
  a = a - 1.0D0
  w = w * (a/(1.0D0 + a/b))
END DO
fn_val = (LOG(w) - n*LOG(b)) + (gamln(a) + algdiv(a,b))
RETURN
!-----------------------------------------------------------------------
!                   PROCEDURE WHEN A .GE. 8
!-----------------------------------------------------------------------
60 w = bcorr(a,b)
h = a/b
c = h/(1.0D0 + h)
u = -(a - 0.5D0)*LOG(c)
v = b*alnrel(h)
IF (u <= v) GO TO 61
fn_val = (((-0.5D0*LOG(b) + e) + w) - v) - u
RETURN
61 fn_val = (((-0.5D0*LOG(b) + e) + w) - u) - v
RETURN
END FUNCTION betaln



FUNCTION gsumln (a, b) RESULT(fn_val)
!-----------------------------------------------------------------------
!          EVALUATION OF THE FUNCTION LN(GAMMA(A + B))
!          FOR 1 .LE. A .LE. 2  AND  1 .LE. B .LE. 2
!-----------------------------------------------------------------------
    REAL(dp), INTENT(IN)                     :: a, b
    REAL(dp)                                 :: fn_val

    REAL(dp)                                 :: x

x = a + b - 2.d0
IF (x > 0.25D0) GO TO 10
fn_val = gamln1(1.0D0 + x)
RETURN
10 IF (x > 1.25D0) GO TO 20
fn_val = gamln1(x) + alnrel(x)
RETURN
20 fn_val = gamln1(x - 1.0D0) + LOG(x*(1.0D0 + x))
RETURN
END FUNCTION gsumln



FUNCTION bcorr (a0, b0) RESULT(fn_val)
!-----------------------------------------------------------------------

!     EVALUATION OF  DEL(A0) + DEL(B0) - DEL(A0 + B0)  WHERE
!     LN(GAMMA(A)) = (A - 0.5)*LN(A) - A + 0.5*LN(2*PI) + DEL(A).
!     IT IS ASSUMED THAT A0 .GE. 8 AND B0 .GE. 8.

!-----------------------------------------------------------------------
    REAL(dp), INTENT(IN)                     :: a0, b0
    REAL(dp)                                 :: fn_val

    REAL(dp) :: a, b, c, c0 = .833333333333333D-01, &
      c1 = -.277777777760991D-02, c2 = .793650666825390D-03, &
      c3 = -.595202931351870D-03, c4 = .837308034031215D-03, &
      c5 = -.165322962780713D-02, h, s11, s3, s5, s7, s9, t, w, x, x2

!------------------------

a = MIN(a0, b0)
b = MAX(a0, b0)

h = a/b
c = h/(1.0D0 + h)
x = 1.0D0/(1.0D0 + h)
x2 = x*x

!                SET SN = (1 - X**N)/(1 - X)

s3 = 1.0D0 + (x + x2)
s5 = 1.0D0 + (x + x2*s3)
s7 = 1.0D0 + (x + x2*s5)
s9 = 1.0D0 + (x + x2*s7)
s11 = 1.0D0 + (x + x2*s9)

!                SET W = DEL(B) - DEL(A + B)

t = (1.0D0/b)**2
w = ((((c5*s11*t + c4*s9)*t + c3*s7)*t + c2*s5)*t + c1*s3)*t + c0
w = w*(c/b)

!                   COMPUTE  DEL(A) + W

t = (1.0D0/a)**2
fn_val = (((((c5*t + c4)*t + c3)*t + c2)*t + c1)*t + c0)/a + w
RETURN
END FUNCTION bcorr



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

!-----------------------------------------------------------------------
w = 0.0D0
w1 = 0.0D0
IF (a < 0.0D0 .OR. b < 0.0D0) GO TO 300
IF (a == 0.0D0 .AND. b == 0.0D0) GO TO 310
IF (x < 0.0D0 .OR. x > 1.0D0) GO TO 320
IF (y < 0.0D0 .OR. y > 1.0D0) GO TO 330
z = ((x + y) - 0.5D0) - 0.5D0
IF (ABS(z) > 3.0D0*eps) GO TO 340

ierr = 0
IF (x == 0.0D0) GO TO 200
IF (y == 0.0D0) GO TO 210
IF (a == 0.0D0) GO TO 211
IF (b == 0.0D0) GO TO 201

eps = MAX(eps, 1.d-15)
IF (MAX(a,b) < 1.d-3*eps) GO TO 230

ind = 0
a0 = a
b0 = b
x0 = x
y0 = y
IF (MIN(a0, b0) > 1.0D0) GO TO 30

!             PROCEDURE FOR A0 .LE. 1 OR B0 .LE. 1

IF (x <= 0.5D0) GO TO 10
ind = 1
a0 = b
b0 = a
x0 = y
y0 = x

10 IF (b0 < MIN(eps,eps*a0)) GO TO 80
IF (a0 < MIN(eps,eps*b0) .AND. b0*x0 <= 1.0D0) GO TO 90
IF (MAX(a0, b0) > 1.0D0) GO TO 20
IF (a0 >= MIN(0.2D0, b0)) GO TO 100
IF (x0**a0 <= 0.9D0) GO TO 100
IF (x0 >= 0.3D0) GO TO 110
n = 20
GO TO 130

20 IF (b0 <= 1.0D0) GO TO 100
IF (x0 >= 0.3D0) GO TO 110
IF (x0 >= 0.1D0) GO TO 21
IF ((x0*b0)**a0 <= 0.7D0) GO TO 100
21 IF (b0 > 15.0D0) GO TO 131
n = 20
GO TO 130

!             PROCEDURE FOR A0 .GT. 1 AND B0 .GT. 1

30 IF (a > b) GO TO 31
lambda = a - (a + b)*x
GO TO 32
31 lambda = (a + b)*y - b
32 IF (lambda >= 0.0D0) GO TO 40
ind = 1
a0 = b
b0 = a
x0 = y
y0 = x
lambda = ABS(lambda)

40 IF (b0 < 40.0D0 .AND. b0*x0 <= 0.7D0) GO TO 100
IF (b0 < 40.0D0) GO TO 140
IF (a0 > b0) GO TO 50
IF (a0 <= 100.0D0) GO TO 120
IF (lambda > 0.03D0*a0) GO TO 120
GO TO 180
50 IF (b0 <= 100.0D0) GO TO 120
IF (lambda > 0.03D0*b0) GO TO 120
GO TO 180

!            EVALUATION OF THE APPROPRIATE ALGORITHM

80 w = fpser(a0, b0, x0, eps)
w1 = 0.5D0 + (0.5D0 - w)
GO TO 220

90 w1 = apser(a0, b0, x0, eps)
w = 0.5D0 + (0.5D0 - w1)
GO TO 220

100 w = bpser(a0, b0, x0, eps)
w1 = 0.5D0 + (0.5D0 - w)
GO TO 220

110 w1 = bpser(b0, a0, y0, eps)
w = 0.5D0 + (0.5D0 - w1)
GO TO 220

120 w = bfrac(a0, b0, x0, y0, lambda, 15.0D0*eps)
w1 = 0.5D0 + (0.5D0 - w)
GO TO 220

130 w1 = bup(b0, a0, y0, x0, n, eps)
b0 = b0 + n
131 CALL bgrat (b0, a0, y0, x0, w1, eps, ierr1)
IF (ierr1 > 0) STOP "Error in BGRAT"
w = 0.5D0 + (0.5D0 - w1)
GO TO 220

140 n = b0
b0 = b0 - n
IF (b0 /= 0.0D0) GO TO 141
n = n - 1
b0 = 1.0D0
141 w = bup(b0, a0, y0, x0, n, eps)
IF (x0 > 0.7D0) GO TO 150
w = w + bpser(a0, b0, x0, eps)
w1 = 0.5D0 + (0.5D0 - w)
GO TO 220

150 IF (a0 > 15.0D0) GO TO 151
n = 20
w = w + bup(a0, b0, x0, y0, n, eps)
a0 = a0 + n
151 CALL bgrat (a0, b0, x0, y0, w, eps, ierr1)
w1 = 0.5D0 + (0.5D0 - w)
GO TO 220

180 w = basym(a0, b0, lambda, 100.0D0*eps)
w1 = 0.5D0 + (0.5D0 - w)
GO TO 220

!               TERMINATION OF THE PROCEDURE

200 IF (a == 0.0D0) GO TO 350
201 w = 0.0D0
w1 = 1.0D0
RETURN

210 IF (b == 0.0D0) GO TO 360
211 w = 1.0D0
w1 = 0.0D0
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

fn_val = 1.0D0
IF (a <= 1.d-3*eps) GO TO 10
fn_val = 0.0D0
t = a*LOG(x)
IF (t < dxparg(1)) RETURN
fn_val = EXP(t)

!                NOTE THAT 1/B(A,B) = B

10 fn_val = (b/a)*fn_val
tol = eps/a
an = a + 1.0D0
t = x
s = t/an
20    an = an + 1.0D0
t = x*t
c = t/an
s = s + c
IF (ABS(c) > tol) GO TO 20

fn_val = fn_val*(1.0D0 + a*s)
RETURN
END FUNCTION fpser



FUNCTION apser (a, b, x, eps) RESULT(fn_val)
!-----------------------------------------------------------------------
!     APSER YIELDS THE INCOMPLETE BETA RATIO I(SUB(1-X))(B,A) FOR
!     A .LE. MIN(EPS,EPS*B), B*X .LE. 1, AND X .LE. 0.5. USED WHEN
!     A IS VERY SMALL. USE ONLY IF ABOVE INEQUALITIES ARE SATISFIED.
!-----------------------------------------------------------------------
    REAL(dp), INTENT(IN)                     :: a, b, x, eps
    REAL(dp)                                 :: fn_val

    REAL(dp)                                 :: aj, bx, c, &
                                                g = .577215664901533D0, j, s, &
                                                t, tol

!--------------------
!--------------------

bx = b*x
t = x - bx
IF (b*eps > 2.d-2) GO TO 10
c = LOG(x) + psi(b) + g + t
GO TO 20
10 c = LOG(bx) + g + t

20 tol = 5.0D0*eps*ABS(c)
j = 1.0D0
s = 0.0D0
30    j = j + 1.0D0
t = t*(x - bx/j)
aj = t/j
s = s + aj
IF (ABS(aj) > tol) GO TO 30

fn_val = -a*(c + s)
RETURN
END FUNCTION apser



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

fn_val = 0.0D0
IF (x == 0.0D0) RETURN
!-----------------------------------------------------------------------
!            COMPUTE THE FACTOR X**A/(A*BETA(A,B))
!-----------------------------------------------------------------------
a0 = MIN(a,b)
IF (a0 < 1.0D0) GO TO 10
z = a*LOG(x) - betaln(a,b)
fn_val = EXP(z)/a
GO TO 70
10 b0 = MAX(a,b)
IF (b0 >= 8.0D0) GO TO 60
IF (b0 > 1.0D0) GO TO 40

!            PROCEDURE FOR A0 .LT. 1 AND B0 .LE. 1

fn_val = x**a
IF (fn_val == 0.0D0) RETURN

apb = a + b
IF (apb > 1.0D0) GO TO 20
z = 1.0D0 + gam1(apb)
GO TO 30
20 u = a + b - 1.d0
z = (1.0D0 + gam1(u))/apb

30 c = (1.0D0 + gam1(a))*(1.0D0 + gam1(b))/z
fn_val = fn_val*c*(b/apb)
GO TO 70

!         PROCEDURE FOR A0 .LT. 1 AND 1 .LT. B0 .LT. 8

40 u = gamln1(a0)
m = b0 - 1.0D0
IF (m < 1) GO TO 50
c = 1.0D0
DO i = 1, m
  b0 = b0 - 1.0D0
  c = c*(b0/(a0 + b0))
END DO
u = LOG(c) + u

50 z = a*LOG(x) - u
b0 = b0 - 1.0D0
apb = a0 + b0
IF (apb > 1.0D0) GO TO 51
t = 1.0D0 + gam1(apb)
GO TO 52
51 u = a0 + b0 - 1.d0
t = (1.0D0 + gam1(u))/apb
52 fn_val = EXP(z)*(a0/a)*(1.0D0 + gam1(b0))/t
GO TO 70

!            PROCEDURE FOR A0 .LT. 1 AND B0 .GE. 8

60 u = gamln1(a0) + algdiv(a0,b0)
z = a*LOG(x) - u
fn_val = (a0/a)*EXP(z)
70 IF (fn_val == 0.0D0 .OR. a <= 0.1D0*eps) RETURN
!-----------------------------------------------------------------------
!                     COMPUTE THE SERIES
!-----------------------------------------------------------------------
sum = 0.0D0
n = 0.0D0
c = 1.0D0
tol = eps/a
100    n = n + 1.0D0
c = c*(0.5D0 + (0.5D0 - b/n))*x
w = c/(a + n)
sum = sum + w
IF (ABS(w) > tol) GO TO 100
fn_val = fn_val*(1.0D0 + a*sum)
RETURN
END FUNCTION bpser



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
ap1 = a + 1.0D0
mu = 0
d = 1.0D0
IF (n == 1 .OR. a < 1.0D0) GO TO 10
IF (apb < 1.1D0*ap1) GO TO 10
mu = ABS(dxparg(1))
k = dxparg(0)
IF (k < mu) mu = k
t = mu
d = EXP(-t)

10 fn_val = brcmp1(mu,a,b,x,y)/a
IF (n == 1 .OR. fn_val == 0.0D0) RETURN
nm1 = n - 1
w = d

!          LET K BE THE INDEX OF THE MAXIMUM TERM

k = 0
IF (b <= 1.0D0) GO TO 40
IF (y > 1.d-4) GO TO 20
k = nm1
GO TO 30
20 r = (b - 1.0D0)*x/y - a
IF (r < 1.0D0) GO TO 40
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
IF (fn_val == 0.0D0) RETURN

c = 1.0D0 + lambda
c0 = b/a
c1 = 1.0D0 + 1.0D0/a
yp1 = y + 1.0D0

n = 0.0D0
p = 1.0D0
s = a + 1.0D0
an = 0.0D0
bn = 1.0D0
anp1 = 1.0D0
bnp1 = c/c1
r = c1/c

!        CONTINUED FRACTION CALCULATION

10    n = n + 1.0D0
t = n/a
w = n*(b - n)*x
e = a/s
alpha = (p*(p + c0)*e*e)*(w*x)
IF (alpha <= 0.0D0) GO TO 20
e = (1.0D0 + t)/(c1 + t + t)
beta = n + w/s + e*(c + n*yp1)
p = 1.0D0 + t
s = s + 2.0D0

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
bnp1 = 1.0D0
GO TO 10

!                 TERMINATION

20 fn_val = fn_val*r
RETURN
END FUNCTION bfrac



FUNCTION brcomp (a, b, x, y) RESULT(fn_val)
!-----------------------------------------------------------------------
!               EVALUATION OF X**A*Y**B/BETA(A,B)
!-----------------------------------------------------------------------
    REAL(dp), INTENT(IN)                     :: a, b, x, y
    REAL(dp)                                 :: fn_val

    INTEGER                                  :: i, n
    REAL(dp)                                 :: a0, apb, b0, c, &
                                                const = .398942280401433D0, &
                                                e, h, lambda, lnx, lny, t, u, &
                                                v, x0, y0, z

!-----------------
!     CONST = 1/SQRT(2*PI)
!-----------------

fn_val = 0.0D0
IF (x == 0.0D0 .OR. y == 0.0D0) RETURN
a0 = MIN(a,b)
IF (a0 >= 8.0D0) GO TO 100

IF (x > 0.375D0) GO TO 10
lnx = LOG(x)
lny = alnrel(-x)
GO TO 20
10 IF (y > 0.375D0) GO TO 11
lnx = alnrel(-y)
lny = LOG(y)
GO TO 20
11 lnx = LOG(x)
lny = LOG(y)

20 z = a*lnx + b*lny
IF (a0 < 1.0D0) GO TO 30
z = z - betaln(a,b)
fn_val = EXP(z)
RETURN
!-----------------------------------------------------------------------
!              PROCEDURE FOR A .LT. 1 OR B .LT. 1
!-----------------------------------------------------------------------
30 b0 = MAX(a,b)
IF (b0 >= 8.0D0) GO TO 80
IF (b0 > 1.0D0) GO TO 60

!                   ALGORITHM FOR B0 .LE. 1

fn_val = EXP(z)
IF (fn_val == 0.0D0) RETURN

apb = a + b
IF (apb > 1.0D0) GO TO 40
z = 1.0D0 + gam1(apb)
GO TO 50
40 u = a + b - 1.d0
z = (1.0D0 + gam1(u))/apb

50 c = (1.0D0 + gam1(a))*(1.0D0 + gam1(b))/z
fn_val = fn_val*(a0*c)/(1.0D0 + a0/b0)
RETURN

!                ALGORITHM FOR 1 .LT. B0 .LT. 8

60 u = gamln1(a0)
n = b0 - 1.0D0
IF (n < 1) GO TO 70
c = 1.0D0
DO i = 1, n
  b0 = b0 - 1.0D0
  c = c*(b0/(a0 + b0))
END DO
u = LOG(c) + u

70 z = z - u
b0 = b0 - 1.0D0
apb = a0 + b0
IF (apb > 1.0D0) GO TO 71
t = 1.0D0 + gam1(apb)
GO TO 72
71 u = a0 + b0 - 1.d0
t = (1.0D0 + gam1(u))/apb
72 fn_val = a0*EXP(z)*(1.0D0 + gam1(b0))/t
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
x0 = h/(1.0D0 + h)
y0 = 1.0D0/(1.0D0 + h)
lambda = a - (a + b)*x
GO TO 110
101 h = b/a
x0 = 1.0D0/(1.0D0 + h)
y0 = h/(1.0D0 + h)
lambda = (a + b)*y - b

110 e = -lambda/a
IF (ABS(e) > 0.6D0) GO TO 111
u = rlog1(e)
GO TO 120
111 u = e - LOG(x/x0)

120 e = lambda/b
IF (ABS(e) > 0.6D0) GO TO 121
v = rlog1(e)
GO TO 130
121 v = e - LOG(y/y0)

130 z = EXP(-(a*u + b*v))
fn_val = const*SQRT(b*x0)*z*EXP(-bcorr(a,b))
RETURN
END FUNCTION brcomp



FUNCTION brcmp1 (mu, a, b, x, y) RESULT(fn_val)
!-----------------------------------------------------------------------
!          EVALUATION OF  EXP(MU) * (X**A*Y**B/BETA(A,B))
!-----------------------------------------------------------------------
    INTEGER, INTENT(IN)                      :: mu
    REAL(dp), INTENT(IN)                     :: a, b, x, y
    REAL(dp)                                 :: fn_val

    INTEGER                                  :: i, n
    REAL(dp)                                 :: a0, apb, b0, c, &
                                                const = .398942280401433D0, &
                                                e, h, lambda, lnx, lny, t, u, &
                                                v, x0, y0, z

!-----------------
!     CONST = 1/SQRT(2*PI)
!-----------------

a0 = MIN(a,b)
IF (a0 >= 8.0D0) GO TO 100

IF (x > 0.375D0) GO TO 10
lnx = LOG(x)
lny = alnrel(-x)
GO TO 20
10 IF (y > 0.375D0) GO TO 11
lnx = alnrel(-y)
lny = LOG(y)
GO TO 20
11 lnx = LOG(x)
lny = LOG(y)

20 z = a*lnx + b*lny
IF (a0 < 1.0D0) GO TO 30
z = z - betaln(a,b)
fn_val = esum(mu,z)
RETURN
!-----------------------------------------------------------------------
!              PROCEDURE FOR A .LT. 1 OR B .LT. 1
!-----------------------------------------------------------------------
30 b0 = MAX(a,b)
IF (b0 >= 8.0D0) GO TO 80
IF (b0 > 1.0D0) GO TO 60

!                   ALGORITHM FOR B0 .LE. 1

fn_val = esum(mu,z)
IF (fn_val == 0.0D0) RETURN

apb = a + b
IF (apb > 1.0D0) GO TO 40
z = 1.0D0 + gam1(apb)
GO TO 50
40 u = a + b - 1.d0
z = (1.0D0 + gam1(u))/apb

50 c = (1.0D0 + gam1(a))*(1.0D0 + gam1(b))/z
fn_val = fn_val*(a0*c)/(1.0D0 + a0/b0)
RETURN

!                ALGORITHM FOR 1 .LT. B0 .LT. 8

60 u = gamln1(a0)
n = b0 - 1.0D0
IF (n < 1) GO TO 70
c = 1.0D0
DO i = 1, n
  b0 = b0 - 1.0D0
  c = c*(b0/(a0 + b0))
END DO
u = LOG(c) + u

70 z = z - u
b0 = b0 - 1.0D0
apb = a0 + b0
IF (apb > 1.0D0) GO TO 71
t = 1.0D0 + gam1(apb)
GO TO 72
71 u = a0 + b0 - 1.d0
t = (1.0D0 + gam1(u))/apb
72 fn_val = a0*esum(mu,z)*(1.0D0 + gam1(b0))/t
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
x0 = h/(1.0D0 + h)
y0 = 1.0D0/(1.0D0 + h)
lambda = a - (a + b)*x
GO TO 110
101 h = b/a
x0 = 1.0D0/(1.0D0 + h)
y0 = h/(1.0D0 + h)
lambda = (a + b)*y - b

110 e = -lambda/a
IF (ABS(e) > 0.6D0) GO TO 111
u = rlog1(e)
GO TO 120
111 u = e - LOG(x/x0)

120 e = lambda/b
IF (ABS(e) > 0.6D0) GO TO 121
v = rlog1(e)
GO TO 130
121 v = e - LOG(y/y0)

130 z = esum(mu,-(a*u + b*v))
fn_val = const*SQRT(b*x0)*z*EXP(-bcorr(a,b))
RETURN
END FUNCTION brcmp1



FUNCTION basym (a, b, lambda, eps) RESULT(fn_val)
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR IX(A,B) FOR LARGE A AND B.
!     LAMBDA = (A + B)*Y - B  AND EPS IS THE TOLERANCE USED.
!     IT IS ASSUMED THAT LAMBDA IS NONNEGATIVE AND THAT
!     A AND B ARE GREATER THAN OR EQUAL TO 15.
!-----------------------------------------------------------------------
    REAL(dp), INTENT(IN)                     :: a, b, lambda, eps
    REAL(dp)                                 :: fn_val

    INTEGER                                  :: i, im1, imj, j, m, mm1, mmj, &
                                                n, np1, num = 20
    REAL(dp) :: a0(21), b0(21), bsum, c(21), d(21), dsum, &
      e0 = 1.12837916709551D0, e1 = .353553390593274D0, f, h, h2, hn, j0, j1, &
      r, r0, r1, s, sum, t, t0, t1, u, w, w0, z, z0, z2, zn, znm1

!------------------------
!     ****** NUM IS THE MAXIMUM VALUE THAT N CAN TAKE IN THE DO LOOP
!            ENDING AT STATEMENT 50. IT IS REQUIRED THAT NUM BE EVEN.
!            THE ARRAYS A0, B0, C, D HAVE DIMENSION NUM + 1.
!------------------------
!     E0 = 2/SQRT(PI)
!     E1 = 2**(-3/2)
!------------------------

fn_val = 0.0D0
IF (a >= b) GO TO 10
h = a/b
r0 = 1.0D0/(1.0D0 + h)
r1 = (b - a)/b
w0 = 1.0D0/SQRT(a*(1.0D0 + h))
GO TO 20
10 h = b/a
r0 = 1.0D0/(1.0D0 + h)
r1 = (b - a)/a
w0 = 1.0D0/SQRT(b*(1.0D0 + h))

20 f = a*rlog1(-lambda/a) + b*rlog1(lambda/b)
t = EXP(-f)
IF (t == 0.0D0) RETURN
z0 = SQRT(f)
z = 0.5D0*(z0/e1)
z2 = f + f

a0(1) = (2.0D0/3.0D0)*r1
c(1) = - 0.5D0*a0(1)
d(1) = - c(1)
j0 = (0.5D0/e0)*erfc1(1,z0)
j1 = e1
sum = j0 + d(1)*w0*j1

s = 1.0D0
h2 = h*h
hn = 1.0D0
w = w0
znm1 = z
zn = z2
DO n = 2, num, 2
  hn = h2*hn
  a0(n) = 2.0D0*r0*(1.0D0 + h*hn)/(n + 2.0D0)
  np1 = n + 1
  s = s + hn
  a0(np1) = 2.0D0*r1*s/(n + 3.0D0)
  
  DO i = n, np1
    r = -0.5D0*(i + 1.0D0)
    b0(1) = r*a0(1)
    DO m = 2, i
      bsum = 0.0D0
      mm1 = m - 1
      DO j = 1, mm1
        mmj = m - j
        bsum = bsum + (j*r - mmj)*a0(j)*b0(mmj)
      END DO
      b0(m) = r*a0(m) + bsum/m
    END DO
    c(i) = b0(i)/(i + 1.0D0)
    
    dsum = 0.0D0
    im1 = i - 1
    DO j = 1, im1
      imj = i - j
      dsum = dsum + d(imj)*c(j)
    END DO
    d(i) = -(dsum + c(i))
  END DO
  
  j0 = e1*znm1 + (n - 1.0D0)*j0
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
