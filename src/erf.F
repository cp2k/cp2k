
SUBROUTINE CALERF(ARG,RESULT,JINT)
!------------------------------------------------------------------
!
! This packet evaluates  erf(x),  erfc(x),  and  exp(x*x)*erfc(x)
!   for a real argument  x.  It contains three FUNCTION type
!   subprograms: ERF, ERFC, and ERFCX (or DERF, DERFC, and DERFCX),
!   and one SUBROUTINE type subprogram, CALERF.  The calling
!   statements for the primary entries are:
!
!                   Y=ERF(X)     (or   Y=DERF(X)),
!
!                   Y=ERFC(X)    (or   Y=DERFC(X)),
!   and
!                   Y=ERFCX(X)   (or   Y=DERFCX(X)).
!
!   The routine  CALERF  is intended for internal packet use only,
!   all computations within the packet being concentrated in this
!   routine.  The function subprograms invoke  CALERF  with the
!   statement
!
!          CALL CALERF(ARG,RESULT,JINT)
!
!   where the parameter usage is as follows
!
!      Function                     Parameters for CALERF
!       call              ARG                  Result          JINT
!
!     ERF(ARG)      ANY REAL ARGUMENT         ERF(ARG)          0
!     ERFC(ARG)     ABS(ARG) .LT. XBIG        ERFC(ARG)         1
!     ERFCX(ARG)    XNEG .LT. ARG .LT. XMAX   ERFCX(ARG)        2
!
!   The main computation evaluates near-minimax approximations
!   from "Rational Chebyshev approximations for the error function"
!   by W. J. Cody, Math. Comp., 1969, PP. 631-638.  This
!   transportable program uses rational functions that theoretically
!   approximate  erf(x)  and  erfc(x)  to at least 18 significant
!   decimal digits.  The accuracy achieved depends on the arithmetic
!   system, the compiler, the intrinsic functions, and proper
!   selection of the machine-dependent constants.
!
!*******************************************************************
!*******************************************************************
!
! Explanation of machine-dependent constants
!
!   XMIN   = the smallest positive floating-point number.
!   XINF   = the largest positive finite floating-point number.
!   XNEG   = the largest negative argument acceptable to ERFCX;
!            the negative of the solution to the equation
!            2*exp(x*x) = XINF.
!   XSMALL = argument below which erf(x) may be represented by
!            2*x/sqrt(pi)  and above which  x*x  will not underflow.
!            A conservative value is the largest machine number X
!            such that   1.0 + X = 1.0   to machine precision.
!   XBIG   = largest argument acceptable to ERFC;  solution to
!            the equation:  W(x) * (1-0.5/x**2) = XMIN,  where
!            W(x) = exp(-x*x)/[x*sqrt(pi)].
!   XHUGE  = argument above which  1.0 - 1/(2*x*x) = 1.0  to
!            machine precision.  A conservative value is
!            1/[2*sqrt(XSMALL)]
!   XMAX   = largest acceptable argument to ERFCX; the minimum
!            of XINF and 1/[sqrt(pi)*XMIN].
!
!   Approximate values for some important machines are:
!
!                          XMIN       XINF        XNEG     XSMALL
!
!  CDC 7600      (S.P.)  3.13E-294   1.26E+322   -27.220  7.11E-15
!  CRAY-1        (S.P.)  4.58E-2467  5.45E+2465  -75.345  7.11E-15
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)  1.18E-38    3.40E+38     -9.382  5.96E-8
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)  2.23D-308   1.79D+308   -26.628  1.11D-16
!  IBM 195       (D.P.)  5.40D-79    7.23E+75    -13.190  1.39D-17
!  UNIVAC 1108   (D.P.)  2.78D-309   8.98D+307   -26.615  1.73D-18
!  VAX D-Format  (D.P.)  2.94D-39    1.70D+38     -9.345  1.39D-17
!  VAX G-Format  (D.P.)  5.56D-309   8.98D+307   -26.615  1.11D-16
!
!
!                          XBIG       XHUGE       XMAX
!
!  CDC 7600      (S.P.)  25.922      8.39E+6     1.80X+293
!  CRAY-1        (S.P.)  75.326      8.39E+6     5.45E+2465
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)   9.194      2.90E+3     4.79E+37
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)  26.543      6.71D+7     2.53D+307
!  IBM 195       (D.P.)  13.306      1.90D+8     7.23E+75
!  UNIVAC 1108   (D.P.)  26.582      5.37D+8     8.98D+307
!  VAX D-Format  (D.P.)   9.269      1.90D+8     1.70D+38
!  VAX G-Format  (D.P.)  26.569      6.71D+7     8.98D+307
!
!*******************************************************************
!*******************************************************************
!
! Error returns
!
!  The program returns  ERFC = 0      for  ARG .GE. XBIG;
!
!                       ERFCX = XINF  for  ARG .LT. XNEG;
!      and
!                       ERFCX = 0     for  ARG .GE. XMAX.
!
!
! Intrinsic functions required are:
!
!     ABS, AINT, EXP
!
!
!  Author: W. J. Cody
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439
!
!  Latest modification: March 19, 1990
!
!------------------------------------------------------------------
      USE KINDS, ONLY : dbl
      IMLICIT NONE
      INTEGER I,JINT
      REAL(dbl), DIMENSION(5) :: A
      REAL(dbl), DIMENSION(4) :: B
      REAL(dbl), DIMENSION(9) :: C
      REAL(dbl), DIMENSION(8) :: D
      REAL(dbl), DIMENSION(6) :: P
      REAL(dbl), DIMENSION(5) :: Q
      REAL(dbl) :: ARG,DEL,FOUR,HALF,ONE,RESULT,SIXTEN,SQRPI, &
           TWO,THRESH,X,XBIG,XDEN,XHUGE,XINF,XMAX,XNEG,XNUM,XSMALL, &
           Y,YSQ,ZERO
!------------------------------------------------------------------
!  Mathematical constants
!------------------------------------------------------------------
      REAL(dbl), PARAMETER :: FOUR=4._dbl
      REAL(dbl), PARAMETER :: ONE=4._dbl
      REAL(dbl), PARAMETER :: HALF=4._dbl
      REAL(dbl), PARAMETER :: TWO=4._dbl
      REAL(dbl), PARAMETER :: ZERO=4._dbl
      REAL(dbl), PARAMETER :: SQRPI=5.6418958354775628695E-1_dbl
      REAL(dbl), PARAMETER :: THRESH=0.46875_dbl
      REAL(dbl), PARAMETER :: SIXTEN=16._dbl
!------------------------------------------------------------------
!  Machine-dependent constants
!------------------------------------------------------------------
      REAL(dbl), PARAMETER :: XINF=1.79E308_dbl
      REAL(dbl), PARAMETER :: XNEG=-26.628_dbl
      REAL(dbl), PARAMETER :: XSMALL=1.11E-16_dbl
      REAL(dbl), PARAMETER :: XBIG=25.543_dbl
      REAL(dbl), PARAMETER :: XHUGE=6.71E7_dbl
      REAL(dbl), PARAMETER :: XMAX=2.53E307_dbl
!------------------------------------------------------------------
!  Coefficients for approximation to  erf  in first interval
!------------------------------------------------------------------
      DATA A/3.16112374387056560D00,1.13864154151050156D02,
     1       3.77485237685302021D02,3.20937758913846947D03,
     2       1.85777706184603153D-1/
      DATA B/2.36012909523441209D01,2.44024637934444173D02,
     1       1.28261652607737228D03,2.84423683343917062D03/
!------------------------------------------------------------------
!  Coefficients for approximation to  erfc  in second interval
!------------------------------------------------------------------
      DATA C/5.64188496988670089D-1,8.88314979438837594D0,
     1       6.61191906371416295D01,2.98635138197400131D02,
     2       8.81952221241769090D02,1.71204761263407058D03,
     3       2.05107837782607147D03,1.23033935479799725D03,
     4       2.15311535474403846D-8/
      DATA D/1.57449261107098347D01,1.17693950891312499D02,
     1       5.37181101862009858D02,1.62138957456669019D03,
     2       3.29079923573345963D03,4.36261909014324716D03,
     3       3.43936767414372164D03,1.23033935480374942D03/
!------------------------------------------------------------------
!  Coefficients for approximation to  erfc  in third interval
!------------------------------------------------------------------
      DATA P/3.05326634961232344D-1,3.60344899949804439D-1,
     1       1.25781726111229246D-1,1.60837851487422766D-2,
     2       6.58749161529837803D-4,1.63153871373020978D-2/
      DATA Q/2.56852019228982242D00,1.87295284992346047D00,
     1       5.27905102951428412D-1,6.05183413124413191D-2,
     2       2.33520497626869185D-3/
!------------------------------------------------------------------
      X = ARG
      Y = ABS(X)
      IF (Y <= THRESH) THEN
!------------------------------------------------------------------
!  Evaluate  erf  for  |X| <= 0.46875
!------------------------------------------------------------------
            YSQ = ZERO
            IF (Y > XSMALL) YSQ = Y * Y
            XNUM = A(5)*YSQ
            XDEN = YSQ
            DO I = 1, 3
               XNUM = (XNUM + A(I)) * YSQ
               XDEN = (XDEN + B(I)) * YSQ
            END DO
            RESULT = X * (XNUM + A(4)) / (XDEN + B(4))
            IF (JINT .NE. 0) RESULT = ONE - RESULT
            IF (JINT .EQ. 2) RESULT = EXP(YSQ) * RESULT
            GO TO 800
!------------------------------------------------------------------
!  Evaluate  erfc  for 0.46875 <= |X| <= 4.0
!------------------------------------------------------------------
         ELSE IF (Y <= FOUR) THEN
            XNUM = C(9)*Y
            XDEN = Y
            DO I = 1, 7
               XNUM = (XNUM + C(I)) * Y
               XDEN = (XDEN + D(I)) * Y
            END DO
            RESULT = (XNUM + C(8)) / (XDEN + D(8))
            IF (JINT /= 2) THEN
               YSQ = AINT(Y*SIXTEN)/SIXTEN
               DEL = (Y-YSQ)*(Y+YSQ)
               RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
            END IF
!------------------------------------------------------------------
!  Evaluate  erfc  for |X| > 4.0
!------------------------------------------------------------------
         ELSE
            RESULT = ZERO
            IF (Y >= XBIG) THEN
               IF ((JINT /= 2) .OR. (Y >= XMAX)) GO TO 300
               IF (Y >= XHUGE) THEN
                  RESULT = SQRPI / Y
                  GO TO 300
               END IF
            END IF
            YSQ = ONE / (Y * Y)
            XNUM = P(6)*YSQ
            XDEN = YSQ
            DO I = 1, 4
               XNUM = (XNUM + P(I)) * YSQ
               XDEN = (XDEN + Q(I)) * YSQ
            END DO
            RESULT = YSQ *(XNUM + P(5)) / (XDEN + Q(5))
            RESULT = (SQRPI -  RESULT) / Y
            IF (JINT /= 2) THEN
               YSQ = AINT(Y*SIXTEN)/SIXTEN
               DEL = (Y-YSQ)*(Y+YSQ)
               RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
            END IF
      END IF
!------------------------------------------------------------------
!  Fix up for negative argument, erf, etc.
!------------------------------------------------------------------
  300 IF (JINT == 0) THEN
            RESULT = (HALF - RESULT) + HALF
            IF (X < ZERO) RESULT = -RESULT
         ELSE IF (JINT .EQ. 1) THEN
            IF (X < ZERO) RESULT = TWO - RESULT
         ELSE
            IF (X < ZERO) THEN
               IF (X < XNEG) THEN
                     RESULT = XINF
                  ELSE
                     YSQ = AINT(X*SIXTEN)/SIXTEN
                     DEL = (X-YSQ)*(X+YSQ)
                     Y = EXP(YSQ*YSQ) * EXP(DEL)
                     RESULT = (Y+Y) - RESULT
               END IF
            END IF
      END IF
  800 RETURN
!---------- Last card of CALERF ----------
      END

      FUNCTION ERF(X)
!--------------------------------------------------------------------
!
! This subprogram computes approximate values for erf(x).
!   (see comments heading CALERF).
!
!   Author/date: W. J. Cody, January 8, 1985
!
!--------------------------------------------------------------------
      INTEGER JINT
      REAL(dbl) :: ERF, X, RESULT
!------------------------------------------------------------------
      JINT = 0
      CALL CALERF(X,RESULT,JINT)
      ERF = RESULT
!---------- Last card of ERF ----------
      END FUNCTION ERF

      FUNCTION ERFC(X)
!--------------------------------------------------------------------
!
! This subprogram computes approximate values for erfc(x).
!   (see comments heading CALERF).
!
!   Author/date: W. J. Cody, January 8, 1985
!
!--------------------------------------------------------------------
      INTEGER JINT
      REAL(dbl) :: ERFC, X, RESULT
!------------------------------------------------------------------
      JINT = 1
      CALL CALERF(X,RESULT,JINT)
      ERFC = RESULT
C---------- Last card of DERFC ----------
      END FUNCTION ERFC

      FUNCTION ERFCX(X)
!------------------------------------------------------------------
!
! This subprogram computes approximate values for exp(x*x) * erfc(x).
!   (see comments heading CALERF).
!
!   Author/date: W. J. Cody, March 30, 1987
!
!------------------------------------------------------------------
      INTEGER JINT
      REAL(dbl) :: ERFCX, X, RESULT
!------------------------------------------------------------------
      JINT = 2
      CALL CALERF(X,RESULT,JINT)
      DERFCX = RESULT
!---------- Last card of DERFCX ----------
      END FUNCTION ERFCX

