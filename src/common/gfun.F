!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of the G function G_n(t) for 1/R^2 operators
!
!         (1) J. 0. Jensen, A. H. Cameri, C. P. Vlahacos, D. Zeroka, H. F. Hameka, C. N. Merrow,
!             Evaluation of one-electron integrals for arbitrary operators V(r) over cartesian Gaussians:
!             Application to inverse-square distance and Yukawa operators.
!             J. Comput. Chem. 14(8), 986 (1993).
!             doi: 10.1002/jcc.540140814
!         (2) B. Gao, A. J. Thorvaldsen, K. Ruud,
!             GEN1INT: A unified procedure for the evaluation of one-electron integrals over Gaussian
!             basis functions and their geometric derivatives.
!             Int. J. Quantum Chem. 111(4), 858 (2011).
!             doi: 10.1002/qua.22886
!         (3) libgrpp : specfun_gfun.c
!         (4) William Cody, Kathleen Paciorek, Henry Thacher,
!             Chebyshev Approximations for Dawson's Integral,
!             Mathematics of Computation,
!             Volume 24, Number 109, January 1970, pages 171-178.
!
!> \author JHU
! **************************************************************************************************
MODULE gfun

   USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gfun'

   PUBLIC :: gfun_values

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param nmax ...
!> \param t ...
!> \param g ...
! **************************************************************************************************
   SUBROUTINE gfun_values(nmax, t, g)

      INTEGER, INTENT(IN)                                :: nmax
      REAL(KIND=dp), INTENT(IN)                          :: t
      REAL(KIND=dp), DIMENSION(0:nmax), INTENT(OUT)      :: g

      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: st

      g = 0.0_dp

      IF (t <= 12.0_dp) THEN
         ! downward recursion
         g(nmax) = gfun_taylor(nmax, t); 
         DO i = nmax, 1, -1
            g(i - 1) = (1.0_dp - 2.0_dp*t*g(i))/(2.0_dp*i - 1.0_dp)
         END DO
      ELSE
         ! upward recursion
         st = SQRT(t)
         g(0) = daw(st)/st
         DO i = 0, nmax - 1
            g(i + 1) = (1.0_dp - (2.0_dp*i + 1.0_dp)*g(i))/(2.0_dp*t)
         END DO
      END IF

   END SUBROUTINE gfun_values

! **************************************************************************************************
!> \brief ...
!> \param n ...
!> \param x ...
!> \return ...
! **************************************************************************************************
   FUNCTION gfun_taylor(n, x) RESULT(g)
      INTEGER, INTENT(IN)                                :: n
      REAL(KIND=dp), INTENT(IN)                          :: x
      REAL(KIND=dp)                                      :: g

      REAL(KIND=dp), PARAMETER                           :: eps = 1.E-15_dp

      INTEGER                                            :: k
      REAL(KIND=dp)                                      :: ex, gk, tkk

      ex = EXP(-x)
      tkk = 1.0_dp
      g = ex/REAL(2*n + 1, KIND=dp)
      DO k = 1, 100
         tkk = tkk*x/REAL(k, KIND=dp)
         gk = ex*tkk/REAL(2*n + 2*k + 1, KIND=dp)
         g = g + gk
         IF (gk < eps) EXIT
      END DO
      IF (gk > eps) THEN
         CPWARN("gfun_taylor did not converge")
      END IF

   END FUNCTION gfun_taylor

!*****************************************************************************80
!
!  DAW evaluates Dawson's integral function.
!
!  Discussion:
!
!    This routine evaluates Dawson's integral,
!
!      F(x) = exp ( - x * x ) * Integral ( 0 <= t <= x ) exp ( t * t ) dt
!
!    for a real argument x.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Kathleen Paciorek, Henry Thacher,
!    Chebyshev Approximations for Dawson's Integral,
!    Mathematics of Computation,
!    Volume 24, Number 109, January 1970, pages 171-178.
!
!  Parameters:
!
!    Input, real ( kind = dp ) XX, the argument of the function.
!
!    Output, real ( kind = dp ) DAW, the value of the function.
!
! **************************************************************************************************
!> \brief ...
!> \param xx ...
!> \return ...
! **************************************************************************************************
   FUNCTION daw(xx)

      REAL(kind=dp)                                      :: xx, daw

      INTEGER                                            :: i
      REAL(kind=dp)                                      :: frac, one225, p1(10), p2(10), p3(10), &
                                                            p4(10), q1(10), q2(9), q3(9), q4(9), &
                                                            six25, sump, sumq, two5, w2, x, &
                                                            xlarge, xmax, xsmall, y

!
!  Mathematical constants.
!
      DATA six25/6.25D+00/
      DATA one225/12.25d0/
      DATA two5/25.0d0/
!
!  Machine-dependent constants
!
      DATA xsmall/1.05d-08/
      DATA xlarge/9.49d+07/
      DATA xmax/2.24d+307/
!
!  Coefficients for R(9,9) approximation for  |x| < 2.5
!
      DATA p1/-2.69020398788704782410d-12, 4.18572065374337710778d-10, &
         -1.34848304455939419963d-08, 9.28264872583444852976d-07, &
         -1.23877783329049120592d-05, 4.07205792429155826266d-04, &
         -2.84388121441008500446d-03, 4.70139022887204722217d-02, &
         -1.38868086253931995101d-01, 1.00000000000000000004d+00/
      DATA q1/1.71257170854690554214d-10, 1.19266846372297253797d-08, &
         4.32287827678631772231d-07, 1.03867633767414421898d-05, &
         1.78910965284246249340d-04, 2.26061077235076703171d-03, &
         2.07422774641447644725d-02, 1.32212955897210128811d-01, &
         5.27798580412734677256d-01, 1.00000000000000000000d+00/
!
!  Coefficients for R(9,9) approximation in J-fraction form
!  for  x in [2.5, 3.5)
!
      DATA p2/-1.70953804700855494930d+00, -3.79258977271042880786d+01, &
         2.61935631268825992835d+01, 1.25808703738951251885d+01, &
         -2.27571829525075891337d+01, 4.56604250725163310122d+00, &
         -7.33080089896402870750d+00, 4.65842087940015295573d+01, &
         -1.73717177843672791149d+01, 5.00260183622027967838d-01/
      DATA q2/1.82180093313514478378d+00, 1.10067081034515532891d+03, &
         -7.08465686676573000364d+00, 4.53642111102577727153d+02, &
         4.06209742218935689922d+01, 3.02890110610122663923d+02, &
         1.70641269745236227356d+02, 9.51190923960381458747d+02, &
         2.06522691539642105009d-01/
!
!  Coefficients for R(9,9) approximation in J-fraction form
!  for  x in [3.5, 5.0]
!
      DATA p3/-4.55169503255094815112d+00, -1.86647123338493852582d+01, &
         -7.36315669126830526754d+00, -6.68407240337696756838d+01, &
         4.84507265081491452130d+01, 2.69790586735467649969d+01, &
         -3.35044149820592449072d+01, 7.50964459838919612289d+00, &
         -1.48432341823343965307d+00, 4.99999810924858824981d-01/
      DATA q3/4.47820908025971749852d+01, 9.98607198039452081913d+01, &
         1.40238373126149385228d+01, 3.48817758822286353588d+03, &
         -9.18871385293215873406d+00, 1.24018500009917163023d+03, &
         -6.88024952504512254535d+01, -2.31251575385145143070d+00, &
         2.50041492369922381761d-01/
!
!  Coefficients for R(9,9) approximation in J-fraction form
!  for 5.0 < |x|.
!
      DATA p4/-8.11753647558432685797d+00, -3.84043882477454453430d+01, &
         -2.23787669028751886675d+01, -2.88301992467056105854d+01, &
         -5.99085540418222002197d+00, -1.13867365736066102577d+01, &
         -6.52828727526980741590d+00, -4.50002293000355585708d+00, &
         -2.50000000088955834952d+00, 5.00000000000000488400d-01/
      DATA q4/2.69382300417238816428d+02, 5.04198958742465752861d+01, &
         6.11539671480115846173d+01, 2.08210246935564547889d+02, &
         1.97325365692316183531d+01, -1.22097010558934838708d+01, &
         -6.99732735041547247161d+00, -2.49999970104184464568d+00, &
         7.49999999999027092188d-01/

      x = xx

      IF (xlarge < ABS(x)) THEN

         IF (ABS(x) <= xmax) THEN
            daw = 0.5D+00/x
         ELSE
            daw = 0.0D+00
         END IF

      ELSE IF (ABS(x) < xsmall) THEN

         daw = x

      ELSE

         y = x*x
!
!  ABS(X) < 2.5.
!
         IF (y < six25) THEN

            sump = p1(1)
            sumq = q1(1)
            DO i = 2, 10
               sump = sump*y + p1(i)
               sumq = sumq*y + q1(i)
            END DO

            daw = x*sump/sumq
!
!  2.5 <= ABS(X) < 3.5.
!
         ELSE IF (y < one225) THEN

            frac = 0.0D+00
            DO i = 1, 9
               frac = q2(i)/(p2(i) + y + frac)
            END DO

            daw = (p2(10) + frac)/x
!
!  3.5 <= ABS(X) < 5.0.
!
         ELSE IF (y < two5) THEN

            frac = 0.0D+00
            DO i = 1, 9
               frac = q3(i)/(p3(i) + y + frac)
            END DO

            daw = (p3(10) + frac)/x

         ELSE
!
!  5.0 <= ABS(X) <= XLARGE.
!
            w2 = 1.0D+00/x/x

            frac = 0.0D+00
            DO i = 1, 9
               frac = q4(i)/(p4(i) + y + frac)
            END DO
            frac = p4(10) + frac

            daw = (0.5D+00 + 0.5D+00*w2*frac)/x

         END IF

      END IF

   END FUNCTION daw

END MODULE gfun
