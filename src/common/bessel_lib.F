!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculates Bessel functions
!> \note
!>      Functions adapted from netlib
!> \par History
!>      March-2006: Bessel Transform (JGH)
!> \author JGH (10-02-2001)
! **************************************************************************************************
MODULE bessel_lib

   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fac,&
                                              pi
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: bessel0

CONTAINS

! **************************************************************************************************
!> \brief Calculates spherical Bessel functions
!> \param x ...
!> \param l ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL IMPURE FUNCTION bessel0(x, l)
      !
      ! Abramowitz & Stegun using Formulas 10.1.2, 10.1.8, 10.1.9
      ! Adapted from P. Bloechl
      !
      REAL(KIND=dp), INTENT(IN)                          :: x
      INTEGER, INTENT(IN)                                :: l
      REAL(KIND=dp)                                      :: bessel0

      REAL(KIND=dp), PARAMETER                           :: tol = 1.e-12_dp

      INTEGER                                            :: i, ii, il, isvar, k
      REAL(KIND=dp)                                      :: arg, fact, xsq
      REAL(KIND=dp), DIMENSION(4)                        :: trig

      IF (x > SQRT(REAL(l, KIND=dp) + 0.5_dp)) THEN
         arg = x - 0.5_dp*REAL(l, KIND=dp)*pi
         trig(1) = SIN(arg)/x
         trig(2) = COS(arg)/x
         trig(3) = -trig(1)
         trig(4) = -trig(2)
         bessel0 = trig(1)
         IF (l /= 0) THEN
            xsq = 0.5_dp/x
            fact = 1._dp
            DO k = 1, l
               ii = MOD(k, 4) + 1
               fact = fac(k + l)/fac(k)/fac(l - k)*xsq**k
               bessel0 = bessel0 + fact*trig(ii)
            END DO
         END IF
      ELSE
         ! Taylor expansion for small arguments
         isvar = 1
         DO il = 1, l
            isvar = isvar*(2*il + 1)
         END DO
         IF (l /= 0._dp) THEN
            fact = x**l/REAL(isvar, KIND=dp)
         ELSE
            fact = 1._dp/REAL(isvar, KIND=dp)
         END IF
         bessel0 = fact
         xsq = -0.5_dp*x*x
         isvar = 2*l + 1
         DO i = 1, 1000
            isvar = isvar + 2
            fact = fact*xsq/REAL(i*isvar, KIND=dp)
            bessel0 = bessel0 + fact
            IF (ABS(fact) < tol) EXIT
         END DO
         IF (ABS(fact) > tol) CPABORT("BESSEL0 NOT CONVERGED")
      END IF

   END FUNCTION bessel0

END MODULE bessel_lib
