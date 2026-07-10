!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
MODULE beta_gamma_psi
   ! not tested in the case where dp would stand for single precision
   ! Routine for the digamma function.

   USE kinds,                           ONLY: dp

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: psi

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param xx ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION psi(xx) RESULT(fn_val)
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

      REAL(dp), PARAMETER :: dx0 = 1.461632144968362341262659542325721325e0_dp, p1(7) = [ &
         .895385022981970e-02_dp, .477762828042627e+01_dp, .142441585084029e+03_dp, &
         .118645200713425e+04_dp, .363351846806499e+04_dp, .413810161269013e+04_dp, &
         .130560269827897e+04_dp], p2(4) = [-.212940445131011e+01_dp, -.701677227766759e+01_dp, &
         -.448616543918019e+01_dp, -.648157123766197e+00_dp], piov4 = .785398163397448e0_dp, q1(6) &
         = [.448452573429826e+02_dp, .520752771467162e+03_dp, .221000799247830e+04_dp, &
         .364127349079381e+04_dp, .190831076596300e+04_dp, .691091682714533e-05_dp]
      REAL(dp), PARAMETER :: q2(4) = [.322703493791143e+02_dp, .892920700481861e+02_dp, &
         .546117738103215e+02_dp, .777788548522962e+01_dp]

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

      xmax1 = MIN(REAL(HUGE(0), KIND=dp), 1.0e0_dp/EPSILON(1.0e0_dp))
      xsmall = 1.e-9_dp
!---------------------------------------------------------------------
      x = xx
      aug = 0.0e0_dp
      IF (x < 0.5e0_dp) THEN
!---------------------------------------------------------------------
!     X < 0.5,  USE REFLECTION FORMULA
!     PSI(1-X) = PSI(X) + PI * COTAN(PI*X)
!---------------------------------------------------------------------
         IF (ABS(x) <= xsmall) THEN
            IF (x == 0.0e0_dp) THEN
               !     ERROR RETURN
               fn_val = 0.0e0_dp
               RETURN
            END IF
!---------------------------------------------------------------------
!     0 < ABS(X) <= XSMALL.  USE 1/X AS A SUBSTITUTE
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
!     MAKE AN ERROR EXIT IF X <= -XMAX1
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
!     0.5 <= X <= 3.0
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
!     IF X >= XMAX1, PSI = LN(X)
!---------------------------------------------------------------------
      IF (x < xmax1) THEN
!---------------------------------------------------------------------
!     3.0 < X < XMAX1
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

END MODULE beta_gamma_psi
