!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
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
   PUBLIC :: bessk0, bessk1, bessel0

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param x must be positive
!> \return ...
! **************************************************************************************************
   ELEMENTAL FUNCTION bessk0(x)

      REAL(KIND=dp), INTENT(IN)                          :: x
      REAL(KIND=dp)                                      :: bessk0

      REAL(KIND=dp), PARAMETER :: p1 = -0.57721566_dp, p2 = 0.42278420_dp, p3 = 0.23069756_dp, &
         p4 = 0.3488590e-1_dp, p5 = 0.262698e-2_dp, p6 = 0.10750e-3_dp, p7 = 0.74e-5_dp, &
         q1 = 1.25331414_dp, q2 = -0.7832358e-1_dp, q3 = 0.2189568e-1_dp, q4 = -0.1062446e-1_dp, &
         q5 = 0.587872e-2_dp, q6 = -0.251540e-2_dp, q7 = 0.53208e-3_dp

      REAL(KIND=dp)                                      :: y

      IF (x < 2.0_dp) THEN
         y = x*x/4.0_dp
         bessk0 = (-LOG(x/2.0_dp)*bessi0(x)) + (p1 + y* &
                                                (p2 + y*(p3 + y*(p4 + y*(p5 + y*(p6 + y*p7))))))
      ELSE
         y = (2.0_dp/x)
         bessk0 = (EXP(-x)/SQRT(x))*(q1 + y*(q2 + y* &
                                             (q3 + y*(q4 + y*(q5 + y*(q6 + y*q7))))))
      END IF

   END FUNCTION bessk0

! **************************************************************************************************
!> \brief ...
!> \param x must be positive
!> \return ...
! **************************************************************************************************
   ELEMENTAL FUNCTION bessk1(x)

      REAL(KIND=dp), INTENT(IN)                          :: x
      REAL(KIND=dp)                                      :: bessk1

      REAL(KIND=dp), PARAMETER :: p1 = 1.0_dp, p2 = 0.15443144_dp, p3 = -0.67278579_dp, &
         p4 = -0.18156897_dp, p5 = -0.1919402e-1_dp, p6 = -0.110404e-2_dp, p7 = -0.4686e-4_dp, &
         q1 = 1.25331414_dp, q2 = 0.23498619_dp, q3 = -0.3655620e-1_dp, q4 = 0.1504268e-1_dp, &
         q5 = -0.780353e-2_dp, q6 = 0.325614e-2_dp, q7 = -0.68245e-3_dp

      REAL(KIND=dp)                                      :: y

      IF (x < 2.0_dp) THEN
         y = x*x/4.0_dp
         bessk1 = (LOG(x/2.0_dp)*bessi1(x)) + (1.0_dp/x)* &
                  (p1 + y*(p2 + y*(p3 + y*(p4 + y*(p5 + y* &
                                                   (p6 + y*p7))))))
      ELSE
         y = 2.0_dp/x
         bessk1 = (EXP(-x)/SQRT(x))*(q1 + y*(q2 + y* &
                                             (q3 + y*(q4 + y*(q5 + y*(q6 + y*q7))))))
      END IF

   END FUNCTION bessk1

! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL FUNCTION bessi0(x)

      REAL(KIND=dp), INTENT(IN)                          :: x
      REAL(KIND=dp)                                      :: bessi0

      REAL(KIND=dp), PARAMETER :: p1 = 1.0_dp, p2 = 3.5156229_dp, p3 = 3.0899424_dp, &
         p4 = 1.2067492_dp, p5 = 0.2659732_dp, p6 = 0.360768e-1_dp, p7 = 0.45813e-2_dp, &
         q1 = 0.39894228_dp, q2 = 0.1328592e-1_dp, q3 = 0.225319e-2_dp, q4 = -0.157565e-2_dp, &
         q5 = 0.916281e-2_dp, q6 = -0.2057706e-1_dp, q7 = 0.2635537e-1_dp, q8 = -0.1647633e-1_dp, &
         q9 = 0.392377e-2_dp

      REAL(KIND=dp)                                      :: ax, y

      IF (ABS(x) < 3.75_dp) THEN
         y = (x/3.75_dp)**2
         bessi0 = p1 + y*(p2 + y*(p3 + y*(p4 + y* &
                                          (p5 + y*(p6 + y*p7)))))
      ELSE
         ax = ABS(x)
         y = 3.75_dp/ax
         bessi0 = (EXP(ax)/SQRT(ax))*(q1 + y*(q2 + y* &
                                              (q3 + y*(q4 + y*(q5 + y*(q6 + y*(q7 + y* &
                                                                               (q8 + y*q9))))))))
      END IF

   END FUNCTION bessi0

! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL FUNCTION bessi1(x)

      REAL(KIND=dp), INTENT(IN)                          :: x
      REAL(KIND=dp)                                      :: bessi1

      REAL(KIND=dp), PARAMETER :: p1 = 0.5_dp, p2 = 0.87890594_dp, p3 = 0.51498869_dp, &
         p4 = 0.15084934e0_dp, p5 = 0.2658733e-1_dp, p6 = 0.301532e-2_dp, p7 = 0.32411e-3_dp, &
         q1 = 0.39894228_dp, q2 = -0.3988024e-1_dp, q3 = -0.362018e-2_dp, q4 = 0.163801e-2_dp, &
         q5 = -0.1031555e-1_dp, q6 = 0.2282967e-1_dp, q7 = -0.2895312e-1_dp, q8 = 0.1787654e-1_dp, &
         q9 = -0.420059e-2_dp

      REAL(KIND=dp)                                      :: ax, y

      IF (ABS(x) < 3.75_dp) THEN
         y = (x/3.75_dp)**2
         bessi1 = p1 + y*(p2 + y*(p3 + y*(p4 + y* &
                                          (p5 + y*(p6 + y*p7)))))
      ELSE
         ax = ABS(x)
         y = 3.75_dp/ax
         bessi1 = (EXP(ax)/SQRT(ax))*(q1 + y*(q2 + y* &
                                              (q3 + y*(q4 + y*(q5 + y*(q6 + y*(q7 + y* &
                                                                               (q8 + y*q9))))))))
         IF (x < 0.0_dp) bessi1 = -bessi1
      END IF

   END FUNCTION bessi1

! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \param l ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL IMPURE FUNCTION bessel0(x, l)
      !
      ! Calculates spherical Bessel functions
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
         ENDDO
         IF (ABS(fact) > tol) CPABORT("BESSEL0 NOT CONVERGED")
      END IF

   END FUNCTION bessel0

END MODULE bessel_lib

