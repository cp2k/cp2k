!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief All kind of helpful little routines
!> \par History
!>      - Cleaned (22.04.2021, MK)
!> \author CJM & JGH
! **************************************************************************************************
MODULE ao_util

   USE kinds,                           ONLY: dp,&
                                              sp
   USE mathconstants,                   ONLY: dfac,&
                                              fac,&
                                              rootpi
   USE orbital_pointers,                ONLY: coset
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ao_util'

   ! Public subroutines

   PUBLIC :: exp_radius, &
             exp_radius_very_extended, &
             gauss_exponent, &
             gaussint_sph, &
             trace_r_AxB

CONTAINS

! **************************************************************************************************
!> \brief  The exponent of a primitive Gaussian function for a given radius
!>         and threshold is calculated.
!> \param l ...
!> \param radius ...
!> \param threshold ...
!> \param prefactor ...
!> \return ...
!> \date   07.03.1999
!> \par Variables
!>      - exponent : Exponent of the primitive Gaussian function.
!>      - l        : Angular momentum quantum number l.
!>      - prefactor: Prefactor of the Gaussian function (e.g. a contraction
!>                   coefficient).
!>      - radius   : Calculated radius of the Gaussian function.
!>      - threshold: Threshold for radius.
!> \author MK
!> \version 1.0
! **************************************************************************************************
   FUNCTION gauss_exponent(l, radius, threshold, prefactor) RESULT(exponent)
      INTEGER, INTENT(IN)                                :: l
      REAL(KIND=dp), INTENT(IN)                          :: radius, threshold, prefactor
      REAL(KIND=dp)                                      :: exponent

      exponent = 0.0_dp

      IF (radius < 1.0E-6_dp) RETURN
      IF (threshold < 1.0E-12_dp) RETURN

      exponent = LOG(ABS(prefactor)*radius**l/threshold)/radius**2

   END FUNCTION gauss_exponent

! **************************************************************************************************
!> \brief  The radius of a primitive Gaussian function for a given threshold
!>         is calculated.
!>               g(r) = prefactor*r**l*exp(-alpha*r**2) - threshold = 0
!> \param l Angular momentum quantum number l.
!> \param alpha Exponent of the primitive Gaussian function.
!> \param threshold Threshold for function g(r).
!> \param prefactor Prefactor of the Gaussian function (e.g. a contraction
!>                     coefficient).
!> \param epsabs Absolute tolerance (radius)
!> \param epsrel Relative tolerance (radius)
!> \param rlow Optional lower bound radius, e.g., when searching maximum radius
!> \return Calculated radius of the Gaussian function
!> \par Variables
!>        - g       : function g(r)
!>        - prec_exp: single/double precision exponential function
!>        - itermax : Maximum number of iterations
!>        - contract: Golden Section Search (GSS): [0.38, 0.62]
!>                    Equidistant sampling: [0.2, 0.4, 0.6, 0.8]
!>                    Bisection (original approach): [0.5]
!>                    Array size may trigger vectorization
! **************************************************************************************************
   FUNCTION exp_radius(l, alpha, threshold, prefactor, epsabs, epsrel, rlow) RESULT(radius)
      INTEGER, INTENT(IN)                                :: l
      REAL(KIND=dp), INTENT(IN)                          :: alpha, threshold, prefactor
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: epsabs, epsrel, rlow
      REAL(KIND=dp)                                      :: radius

      INTEGER, PARAMETER                                 :: itermax = 5000, prec_exp = sp
      REAL(KIND=dp), PARAMETER                           :: contract(*) = (/0.38, 0.62/), &
                                                            mineps = 1.0E-12, next = 2.0, &
                                                            step = 1.0

      INTEGER                                            :: i, j
      REAL(KIND=dp)                                      :: a, d, dr, eps, r, rd, t
      REAL(KIND=dp), DIMENSION(SIZE(contract))           :: g, s

      IF (l .LT. 0) THEN
         CPABORT("The angular momentum quantum number is negative")
      END IF

      IF (alpha .EQ. 0.0_dp) THEN
         CPABORT("The Gaussian function exponent is zero")
      ELSE
         a = ABS(alpha)
      END IF

      IF (threshold .NE. 0.0_dp) THEN
         t = ABS(threshold)
      ELSE
         CPABORT("The requested threshold is zero")
      END IF

      radius = 0.0_dp
      IF (PRESENT(rlow)) radius = rlow
      IF (prefactor .EQ. 0.0_dp) RETURN

      ! MAX: facilitate early exit
      r = MAX(SQRT(0.5_dp*REAL(l, dp)/a), radius)

      d = ABS(prefactor); g(1) = d
      IF (l .NE. 0) THEN
         g(1) = g(1)*EXP(REAL(-a*r*r, KIND=prec_exp))*r**l
      END IF
      ! original approach may return radius=0
      ! with g(r) != g(radius)
      !radius = r
      IF (g(1) .LT. t) RETURN ! early exit

      radius = r*next + step
      DO i = 1, itermax
         g(1) = d*EXP(REAL(-a*radius*radius, KIND=prec_exp))*radius**l
         IF (g(1) .LT. t) EXIT
         r = radius; radius = r*next + step
      END DO

      ! consider absolute and relative accuracy (interval width)
      IF (PRESENT(epsabs)) THEN
         eps = epsabs
      ELSE IF (.NOT. PRESENT(epsrel)) THEN
         eps = mineps
      ELSE
         eps = HUGE(eps)
      END IF
      IF (PRESENT(epsrel)) eps = MIN(eps, epsrel*r)

      dr = 0.0_dp
      DO i = i + 1, itermax
         rd = radius - r
         ! check if finished or no further progress
         IF ((rd .LT. eps) .OR. (rd .EQ. dr)) RETURN
         s = r + rd*contract ! interval contraction
         g = d*EXP(REAL(-a*s*s, KIND=prec_exp))*s**l
         DO j = 1, SIZE(contract)
            IF (g(j) .LT. t) THEN
               radius = s(j) ! interval [r, sj)
               EXIT
            ELSE
               r = s(j) ! interval [sj, radius)
            END IF
         END DO
         dr = rd
      END DO
      IF (i .GE. itermax) THEN
         CPABORT("Maximum number of iterations reached")
      END IF

   END FUNCTION exp_radius

! **************************************************************************************************
!> \brief computes the radius of the Gaussian outside of which it is smaller
!>      than eps
!> \param la_min ...
!> \param la_max ...
!> \param lb_min ...
!> \param lb_max ...
!> \param pab ...
!> \param o1 ...
!> \param o2 ...
!> \param ra ...
!> \param rb ...
!> \param rp ...
!> \param zetp ...
!> \param eps ...
!> \param prefactor ...
!> \param cutoff ...
!> \param epsabs ...
!> \return ...
!> \par History
!>      03.2007 new version that assumes that the Gaussians origante from spherical
!>              Gaussians
!> \note
!>      can optionally screen by the maximum element of the pab block
! **************************************************************************************************
   FUNCTION exp_radius_very_extended(la_min, la_max, lb_min, lb_max, pab, o1, o2, ra, rb, rp, &
                                     zetp, eps, prefactor, cutoff, epsabs) RESULT(radius)

      INTEGER, INTENT(IN)                                :: la_min, la_max, lb_min, lb_max
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: pab
      INTEGER, OPTIONAL                                  :: o1, o2
      REAL(KIND=dp), INTENT(IN)                          :: ra(3), rb(3), rp(3), zetp, eps, &
                                                            prefactor, cutoff
      REAL(KIND=dp), OPTIONAL                            :: epsabs
      REAL(KIND=dp)                                      :: radius

      INTEGER                                            :: i, ico, j, jco, la(3), lb(3), lxa, lxb, &
                                                            lya, lyb, lza, lzb
      REAL(KIND=dp)                                      :: bini, binj, coef(0:20), epsin_local, &
                                                            polycoef(0:60), prefactor_local, &
                                                            rad_a, rad_b, s1, s2

! get the local prefactor, we'll now use the largest density matrix element of the block to screen

      epsin_local = 1.0E-2_dp
      IF (PRESENT(epsabs)) epsin_local = epsabs

      IF (PRESENT(pab)) THEN
         prefactor_local = cutoff
         DO lxa = 0, la_max
         DO lxb = 0, lb_max
            DO lya = 0, la_max - lxa
            DO lyb = 0, lb_max - lxb
               DO lza = MAX(la_min - lxa - lya, 0), la_max - lxa - lya
               DO lzb = MAX(lb_min - lxb - lyb, 0), lb_max - lxb - lyb
                  la = (/lxa, lya, lza/)
                  lb = (/lxb, lyb, lzb/)
                  ico = coset(lxa, lya, lza)
                  jco = coset(lxb, lyb, lzb)
                  prefactor_local = MAX(ABS(pab(o1 + ico, o2 + jco)), prefactor_local)
               END DO
               END DO
            END DO
            END DO
         END DO
         END DO
         prefactor_local = prefactor*prefactor_local
      ELSE
         prefactor_local = prefactor*MAX(1.0_dp, cutoff)
      END IF

      !
      ! assumes that we can compute the radius for the case where
      ! the Gaussians a and b are both on the z - axis, but at the same
      ! distance as the original a and b
      !
      rad_a = SQRT(SUM((ra - rp)**2))
      rad_b = SQRT(SUM((rb - rp)**2))

      polycoef(0:la_max + lb_max) = 0.0_dp
      DO lxa = 0, la_max
      DO lxb = 0, lb_max
         coef(0:la_max + lb_max) = 0.0_dp
         bini = 1.0_dp
         s1 = 1.0_dp
         DO i = 0, lxa
            binj = 1.0_dp
            s2 = 1.0_dp
            DO j = 0, lxb
               coef(lxa + lxb - i - j) = coef(lxa + lxb - i - j) + bini*binj*s1*s2
               binj = (binj*(lxb - j))/(j + 1)
               s2 = s2*(rad_b)
            END DO
            bini = (bini*(lxa - i))/(i + 1)
            s1 = s1*(rad_a)
         END DO
         DO i = 0, lxa + lxb
            polycoef(i) = MAX(polycoef(i), coef(i))
         END DO
      END DO
      END DO

      polycoef(0:la_max + lb_max) = polycoef(0:la_max + lb_max)*prefactor_local
      radius = 0.0_dp
      DO i = 0, la_max + lb_max
         radius = MAX(radius, exp_radius(i, zetp, eps, polycoef(i), epsin_local, rlow=radius))
      END DO

   END FUNCTION exp_radius_very_extended

! **************************************************************************************************
!> \brief ...
!> \param alpha ...
!> \param l ...
!> \return ...
! **************************************************************************************************
   FUNCTION gaussint_sph(alpha, l)

      !  calculates the radial integral over a spherical Gaussian
      !  of the form
      !     r**(2+l) * exp(-alpha * r**2)
      !

      REAL(dp), INTENT(IN)                               :: alpha
      INTEGER, INTENT(IN)                                :: l
      REAL(dp)                                           :: gaussint_sph

      IF ((l/2)*2 == l) THEN
         !even l:
         gaussint_sph = ROOTPI*0.5_dp**(l/2 + 2)*dfac(l + 1) &
                        /SQRT(alpha)**(l + 3)
      ELSE
         !odd l:
         gaussint_sph = 0.5_dp*fac((l + 1)/2)/SQRT(alpha)**(l + 3)
      END IF

   END FUNCTION gaussint_sph

! **************************************************************************************************
!> \brief ...
!> \param A ...
!> \param lda ...
!> \param B ...
!> \param ldb ...
!> \param m ...
!> \param n ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION trace_r_AxB(A, lda, B, ldb, m, n)

      INTEGER, INTENT(in)                                :: lda
      REAL(dp), INTENT(in)                               :: A(lda, *)
      INTEGER, INTENT(in)                                :: ldb
      REAL(dp), INTENT(in)                               :: B(ldb, *)
      INTEGER, INTENT(in)                                :: m, n
      REAL(dp)                                           :: trace_r_AxB

      INTEGER                                            :: i1, i2, imod, mminus3
      REAL(dp)                                           :: t1, t2, t3, t4

      t1 = 0._dp
      t2 = 0._dp
      t3 = 0._dp
      t4 = 0._dp
      imod = MODULO(m, 4)
      SELECT CASE (imod)
      CASE (0)
         DO i2 = 1, n
            DO i1 = 1, m, 4
               t1 = t1 + A(i1, i2)*B(i1, i2)
               t2 = t2 + A(i1 + 1, i2)*B(i1 + 1, i2)
               t3 = t3 + A(i1 + 2, i2)*B(i1 + 2, i2)
               t4 = t4 + A(i1 + 3, i2)*B(i1 + 3, i2)
            END DO
         END DO
      CASE (1)
         mminus3 = m - 3
         DO i2 = 1, n
            DO i1 = 1, mminus3, 4
               t1 = t1 + A(i1, i2)*B(i1, i2)
               t2 = t2 + A(i1 + 1, i2)*B(i1 + 1, i2)
               t3 = t3 + A(i1 + 2, i2)*B(i1 + 2, i2)
               t4 = t4 + A(i1 + 3, i2)*B(i1 + 3, i2)
            END DO
            t1 = t1 + A(m, i2)*B(m, i2)
         END DO
      CASE (2)
         mminus3 = m - 3
         DO i2 = 1, n
            DO i1 = 1, mminus3, 4
               t1 = t1 + A(i1, i2)*B(i1, i2)
               t2 = t2 + A(i1 + 1, i2)*B(i1 + 1, i2)
               t3 = t3 + A(i1 + 2, i2)*B(i1 + 2, i2)
               t4 = t4 + A(i1 + 3, i2)*B(i1 + 3, i2)
            END DO
            t1 = t1 + A(m - 1, i2)*B(m - 1, i2)
            t2 = t2 + A(m, i2)*B(m, i2)
         END DO
      CASE (3)
         mminus3 = m - 3
         DO i2 = 1, n
            DO i1 = 1, mminus3, 4
               t1 = t1 + A(i1, i2)*B(i1, i2)
               t2 = t2 + A(i1 + 1, i2)*B(i1 + 1, i2)
               t3 = t3 + A(i1 + 2, i2)*B(i1 + 2, i2)
               t4 = t4 + A(i1 + 3, i2)*B(i1 + 3, i2)
            END DO
            t1 = t1 + A(m - 2, i2)*B(m - 2, i2)
            t2 = t2 + A(m - 1, i2)*B(m - 1, i2)
            t3 = t3 + A(m, i2)*B(m, i2)
         END DO
      END SELECT
      trace_r_AxB = t1 + t2 + t3 + t4

   END FUNCTION trace_r_AxB

END MODULE ao_util

