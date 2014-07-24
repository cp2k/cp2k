!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculates error function and complementary error function
!> \par History
!>      adapted to F95
! *****************************************************************************
MODULE erf_fn
  USE f77_blas
  USE kinds,                           ONLY: dp
#include "cp_common_uses.f90"

      IMPLICIT NONE
      PRIVATE
  PUBLIC :: erf, erfc

    CONTAINS

! *****************************************************************************
!> \brief ...
!> \param ARG ...
! *****************************************************************************
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
!     ERFC(ARG)     ABS(ARG) < XBIG           ERFC(ARG)         1
!     ERFCX(ARG)    XNEG < ARG < XMAX         ERFCX(ARG)        2
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
!                          XMIN            XINF        XNEG     XSMALL
!
!  CDC 7600      (S.P.)  3.13E-294      1.26E+322     -27.220  7.11E-15
!  CRAY-1        (S.P.)  4.58E-2467     5.45E+2465    -75.345  7.11E-15
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)  1.18E-38       3.40E+38       -9.382  5.96E-8
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)  2.23e-308_dp   1.79e+308_dp  -26.628  1.11e-16_dp
!  IBM 195       (D.P.)  5.40e-79_dp    7.23E+75      -13.190  1.39e-17_dp
!  UNIVAC 1108   (D.P.)  2.78e-309_dp   8.98e+307_dp  -26.615  1.73e-18_dp
!  VAX D-Format  (D.P.)  2.94e-39_dp    1.70e+38_dp    -9.345  1.39e-17_dp
!  VAX G-Format  (D.P.)  5.56e-309_dp   8.98e+307_dp  -26.615  1.11e-16_dp
!
!
!                          XBIG       XHUGE       XMAX
!
!  CDC 7600      (S.P.)  25.922      8.39E+6     1.80X+293
!  CRAY-1        (S.P.)  75.326      8.39E+6     5.45E+2465
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)   9.194      2.90E+3     4.79E+37
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)  26.543      6.71e+7_dp  2.53e+307_dp
!  IBM 195       (D.P.)  13.306      1.90e+8_dp  7.23E+75
!  UNIVAC 1108   (D.P.)  26.582      5.37e+8_dp  8.98e+307_dp
!  VAX D-Format  (D.P.)   9.269      1.90e+8_dp  1.70e+38_dp
!  VAX G-Format  (D.P.)  26.569      6.71e+7_dp  8.98e+307_dp
!
!
! Error returns
!
!  The program returns  ERFC = 0      for  ARG >= XBIG;
!
!                       ERFCX = XINF  for  ARG < XNEG;
!      and
!                       ERFCX = 0     for  ARG >= XMAX.
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

    REAL(KIND=dp)                            :: ARG, RESULT
    INTEGER                                  :: JINT

    INTEGER                                  :: I
    REAL(KIND=dp) :: A, B, C, D, DEL, FOUR, HALF, ONE, P, Q, SIXTEN, SQRPI, &
      THRESH, TWO, X, XBIG, XDEN, XHUGE, XINF, XMAX, XNEG, XNUM, XSMALL, Y, &
      YSQ, ZERO

  DIMENSION A(5),B(4),C(9),D(8),P(6),Q(5)
!------------------------------------------------------------------
!  Mathematical constants
!------------------------------------------------------------------
  DATA FOUR,ONE,HALF,TWO,ZERO/4.0_dp,1.0_dp,0.5_dp,2.0_dp,0.0_dp/, &
       SQRPI/5.6418958354775628695E-1_dp/,THRESH/0.46875_dp/, &
       SIXTEN/16.0_dp/
!------------------------------------------------------------------
!  Machine-dependent constants
!------------------------------------------------------------------
  DATA XINF,XNEG,XSMALL/1.79E+308_dp,-26.628_dp,1.11E-16_dp/, &
       XBIG,XHUGE,XMAX/26.543_dp,6.71E+7_dp,2.53E+307_dp/
!------------------------------------------------------------------
!  Coefficients for approximation to  erf  in first interval
!------------------------------------------------------------------
  DATA A/3.16112374387056560E+00_dp,1.13864154151050156E+02_dp, &
       3.77485237685302021E+02_dp,3.20937758913846947E+03_dp, &
       1.85777706184603153E-1_dp/
  DATA B/2.36012909523441209E+01_dp,2.44024637934444173E+02_dp, &
       1.28261652607737228E+03_dp,2.84423683343917062E+03_dp/
!------------------------------------------------------------------
!  Coefficients for approximation to  erfc  in second interval
!------------------------------------------------------------------
  DATA C/5.64188496988670089E-1_dp,8.88314979438837594E+0_dp, &
       6.61191906371416295E+01_dp,2.98635138197400131E+02_dp, &
       8.81952221241769090E+02_dp,1.71204761263407058E+03_dp, &
       2.05107837782607147E+03_dp,1.23033935479799725E+03_dp, &
       2.15311535474403846E-8_dp/
  DATA D/1.57449261107098347E+01_dp,1.17693950891312499E+02_dp, &
       5.37181101862009858E+02_dp,1.62138957456669019E+03_dp, &
       3.29079923573345963E+03_dp,4.36261909014324716E+03_dp, &
       3.43936767414372164E+03_dp,1.23033935480374942E+03_dp/
!------------------------------------------------------------------
!  Coefficients for approximation to  erfc  in third interval
!------------------------------------------------------------------
  DATA P /3.05326634961232344E-1_dp,3.60344899949804439E-1_dp, &
       1.25781726111229246E-1_dp,1.60837851487422766E-2_dp, &
       6.58749161529837803E-4_dp,1.63153871373020978E-2_dp/
  DATA Q /2.56852019228982242E+00_dp,1.87295284992346047E+00_dp, &
       5.27905102951428412E-1_dp,6.05183413124413191E-2_dp, &
       2.33520497626869185E-3_dp/

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
     IF (JINT /= 0) RESULT = ONE - RESULT
     IF (JINT == 2) RESULT = EXP(YSQ) * RESULT
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
  ELSE IF (JINT == 1) THEN
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
END SUBROUTINE CALERF

! *****************************************************************************
!> \brief ...
!> \param X ...
! *****************************************************************************
FUNCTION ERF(X)
!--------------------------------------------------------------------
!
! This subprogram computes approximate values for erf(x).
!   (see comments heading CALERF).
!
!   Author/date: W. J. Cody, January 8, 1985
!
!--------------------------------------------------------------------
    REAL(KIND=dp)                            :: X, ERF

    INTEGER                                  :: JINT
    REAL(KIND=dp)                            :: RESULT

  JINT = 0
  CALL CALERF(X,RESULT,JINT)
  ERF = RESULT
  RETURN
!---------- Last card of DERF ----------
END FUNCTION erf

! *****************************************************************************
!> \brief ...
!> \param X ...
! *****************************************************************************
FUNCTION ERFC(X)
!--------------------------------------------------------------------
!
! This subprogram computes approximate values for erfc(x).
!   (see comments heading CALERF).
!
!   Author/date: W. J. Cody, January 8, 1985
!
!--------------------------------------------------------------------
    REAL(KIND=dp)                            :: X, ERFC

    INTEGER                                  :: JINT
    REAL(KIND=dp)                            :: RESULT

  JINT = 1
  CALL CALERF(X,RESULT,JINT)
  ERFC = RESULT
  RETURN
!---------- Last card of DERFC ----------
END FUNCTION erfc

! *****************************************************************************
!> \brief ...
!> \param X ...
! *****************************************************************************
FUNCTION ERFCX(X)
!------------------------------------------------------------------
!
! This subprogram computes approximate values for exp(x*x) * erfc(x).
!   (see comments heading CALERF).
!
!   Author/date: W. J. Cody, March 30, 1987
!
!------------------------------------------------------------------
    REAL(KIND=dp)                            :: X, ERFCX

    INTEGER                                  :: JINT
    REAL(KIND=dp)                            :: RESULT

  JINT = 2
  CALL CALERF(X,RESULT,JINT)
  ERFCX = RESULT
  RETURN
!---------- Last card of DERFCX ----------
END FUNCTION erfcx
END  MODULE erf_fn
