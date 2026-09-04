!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Test of the Gaunt coefficients
!> \note
!>      Checks that the Gaunt coefficients returned by libwignernj agree with a
!>      Lebedev quadrature of the corresponding triple product of the spherical
!>      harmonics of spherical_harmonics::y_lm. This guards the phase and
!>      normalisation conventions on both sides; the values themselves are
!>      exact by construction in libwignernj.
!> \par History
!>      JGH  1-Mar-2002 : Clebsch-Gordon coefficients
!>      09.2026 rewritten against libwignernj [Susi Lehtola]
!> \author JGH  1-Mar-2002
! **************************************************************************************************
MODULE cg_test

   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE kinds,                           ONLY: dp
   USE lebedev,                         ONLY: deallocate_lebedev_grids,&
                                              get_number_of_lebedev_grid,&
                                              init_lebedev_grids,&
                                              lebedev_grid
   USE mathconstants,                   ONLY: pi
   USE spherical_harmonics,             ONLY: y_lm
   USE wignernj_interface,              ONLY: complex_gaunt,&
                                              real_gaunt
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cg_test'
   PUBLIC :: gaunt_coefficient_test

CONTAINS

! **************************************************************************************************
!> \brief Compare the Gaunt coefficients against a Lebedev quadrature of the
!>        triple product of spherical harmonics
! **************************************************************************************************
   SUBROUTINE gaunt_coefficient_test()

      INTEGER, PARAMETER                                 :: lmax = 7
      REAL(KIND=dp), PARAMETER                           :: eps = 1.0E-10_dp

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: yc
      INTEGER                                            :: i1, i2, i3, iw, l, ll, m, n, na, nfun, &
                                                            nfun1
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: lof, mof
      REAL(KIND=dp)                                      :: dmax_complex, dmax_real, ref
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: wa
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: yr

      iw = cp_logger_get_default_io_unit()

      IF (iw <= 0) RETURN

      WRITE (iw, '(/,A,/)') " Test of the Gaunt coefficients"
      WRITE (iw, '(T40,A,T77,I4)') " Maximum l value tested:", lmax

      na = 500
      CALL init_lebedev_grids()
      ll = get_number_of_lebedev_grid(n=na)
      na = lebedev_grid(ll)%n
      ALLOCATE (wa(na))
      wa(1:na) = 4.0_dp*pi*lebedev_grid(ll)%w(1:na)

      ! The two harmonics that are multiplied run up to lmax, their product
      ! is expanded in harmonics up to 2*lmax
      nfun1 = (lmax + 1)**2
      nfun = (2*lmax + 1)**2
      ALLOCATE (lof(nfun), mof(nfun))
      ALLOCATE (yr(na, nfun), yc(na, nfun))
      n = 0
      DO l = 0, 2*lmax
         DO m = -l, l
            n = n + 1
            lof(n) = l
            mof(n) = m
            CALL y_lm(lebedev_grid(ll)%r, yr(:, n), l, m)
            CALL y_lm(lebedev_grid(ll)%r, yc(:, n), l, m)
         END DO
      END DO
      CPASSERT(n == nfun)

      dmax_real = 0.0_dp
      dmax_complex = 0.0_dp
      DO i1 = 1, nfun1
         DO i2 = 1, nfun1
            DO i3 = 1, nfun
               ! Real spherical harmonics
               ref = real_gaunt(lof(i1), mof(i1), lof(i2), mof(i2), lof(i3), mof(i3))
               dmax_real = MAX(dmax_real, &
                               ABS(ref - SUM(wa*yr(:, i1)*yr(:, i2)*yr(:, i3))))
               ! Complex spherical harmonics. Unlike the real ones, the complex
               ! harmonics of y_lm are defined without the Condon-Shortley phase
               ! (see cvy_lm: P_l^|m| with no (-1)^m), so they are related to the
               ! Condon-Shortley harmonics that libwignernj uses by
               !     Y_cp2k(l,m) = (-1)^max(m,0) Y_cs(l,m).
               ! The Gaunt coefficient vanishes unless m1 + m2 + m3 = 0, in which
               ! case the integral is real.
               ref = complex_gaunt(lof(i1), mof(i1), lof(i2), mof(i2), lof(i3), mof(i3))
               ref = ref*(-1.0_dp)**(MAX(mof(i1), 0) + MAX(mof(i2), 0) + MAX(mof(i3), 0))
               dmax_complex = MAX(dmax_complex, &
                                  ABS(CMPLX(ref, 0.0_dp, KIND=dp) - &
                                      SUM(wa*yc(:, i1)*yc(:, i2)*yc(:, i3))))
            END DO
         END DO
         WRITE (iw, '(A,I2,A)') " Test for l = ", lof(i1), " done"
      END DO

      WRITE (iw, '(/,T30,A,T71,ES10.3)') " Real spherical harmonics, max deviation ", dmax_real
      WRITE (iw, '(T30,A,T71,ES10.3)') " Complex spherical harmonics, max deviation ", dmax_complex

      DEALLOCATE (wa, lof, mof, yr, yc)
      CALL deallocate_lebedev_grids()

      IF (MAX(dmax_real, dmax_complex) > eps) THEN
         CPABORT("The Gaunt coefficients disagree with the numerical integration")
      END IF

   END SUBROUTINE gaunt_coefficient_test

END MODULE cg_test
