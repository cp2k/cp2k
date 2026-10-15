!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE ai_os_rr

   USE gamma, ONLY: fgamma => fgamma_0
   USE kinds, ONLY: dp
   USE orbital_pointers, ONLY: coset, indco, ncoset
#include "../base/base_uses.f90"

   IMPLICIT NONE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_os_rr'
   PRIVATE

   ! *** Public subroutines ***

   PUBLIC :: os_rr_ovlp, os_rr_ovlp3, os_rr_coul, os_rr_radial, os_rr_coulomb3

   INTERFACE os_rr_ovlp
      MODULE PROCEDURE os_rr_ovlp_real, os_rr_ovlp_complex
   END INTERFACE os_rr_ovlp

   INTERFACE os_rr_ovlp3
      MODULE PROCEDURE os_rr_ovlp3_real, os_rr_ovlp3_complex
   END INTERFACE os_rr_ovlp3

   INTERFACE os_rr_radial
      MODULE PROCEDURE os_rr_radial_abm
   END INTERFACE os_rr_radial

CONTAINS

   #:for suffix, integral_type in [('real', 'REAL'), ('complex', 'COMPLEX')]
! **************************************************************************************************
!> \brief   Calculation of the basic Obara-Saika recurrence relation
!>          for real or complex product-center displacements. Each Cartesian
!>          table is normalized to a unit zeroth-degree integral.
!> \param rap product-center displacement relative to center A
!> \param la_max ...
!> \param rbp product-center displacement relative to center B
!> \param lb_max ...
!> \param zet real sum of the Gaussian exponents
!> \param ldrr ...
!> \param rr ...
!> \date    02.03.2009
!> \author  VW
!> \version 1.0
! **************************************************************************************************
      SUBROUTINE os_rr_ovlp_${suffix}$ (rap, la_max, rbp, lb_max, zet, ldrr, rr)
         ${integral_type}$ (dp), DIMENSION(3), INTENT(IN)     :: rap
         INTEGER, INTENT(IN)                                :: la_max
         ${integral_type}$ (dp), DIMENSION(3), INTENT(IN)     :: rbp
         INTEGER, INTENT(IN)                                :: lb_max
         REAL(dp), INTENT(IN)                               :: zet
         INTEGER, INTENT(IN)                                :: ldrr
         ${integral_type}$ (dp), DIMENSION(0:ldrr - 1, 0:ldrr - 1, 3) :: rr

         INTEGER                                            :: idir

         DO idir = 1, 3
            CALL os_rr_ovlp_1d_${suffix}$ (rap(idir), la_max, rbp(idir), lb_max, zet, rr(:, :, idir))
         END DO
      END SUBROUTINE os_rr_ovlp_${suffix}$

! **************************************************************************************************
!> \brief Normalized one-dimensional two-center Gaussian moments.
!> \param rap product-center displacement relative to A
!> \param la_max highest A degree
!> \param rbp product-center displacement relative to B
!> \param lb_max highest B degree
!> \param zet sum of Gaussian exponents
!> \param rr normalized moments with zeroth-degree value one
! **************************************************************************************************
      PURE SUBROUTINE os_rr_ovlp_1d_${suffix}$ (rap, la_max, rbp, lb_max, zet, rr)
         ${integral_type}$ (dp), INTENT(IN) :: rap, rbp
         INTEGER, INTENT(IN) :: la_max, lb_max
         REAL(dp), INTENT(IN) :: zet
         ${integral_type}$ (dp), INTENT(INOUT) :: rr(0:, 0:)

         INTEGER :: la, lam1, lap1, lb, lbm1, lbp1
         REAL(dp)                                           :: g

         g = 0.5_dp/zet
         rr(0, 0) = 1.0_dp
         ! recursion along la for lb=0
         IF (la_max > 0) THEN
            rr(1, 0) = rap
            DO la = 1, la_max - 1
               lap1 = la + 1
               lam1 = la - 1
               rr(lap1, 0) = REAL(la, dp)*g*rr(lam1, 0) + rap*rr(la, 0)
            END DO
         END IF
         ! recursion along lb for all la
         IF (lb_max > 0) THEN
            rr(0, 1) = rbp
            DO la = 1, la_max
               lam1 = la - 1
               rr(la, 1) = REAL(la, dp)*g*rr(lam1, 0) + rbp*rr(la, 0)
            END DO
            DO lb = 1, lb_max - 1
               lbp1 = lb + 1
               lbm1 = lb - 1
               rr(0, lbp1) = REAL(lb, dp)*g*rr(0, lbm1) + rbp*rr(0, lb)
               DO la = 1, la_max
                  lam1 = la - 1
                  rr(la, lbp1) = g*(REAL(la, dp)*rr(lam1, lb) + REAL(lb, dp)*rr(la, lbm1)) + &
                                 rbp*rr(la, lb)
               END DO
            END DO
         END IF
      END SUBROUTINE os_rr_ovlp_1d_${suffix}$

! **************************************************************************************************
!> \brief Normalized separable three-center Gaussian overlap tables.
!>        The C recurrence uses the directly supplied G-C displacement.
!> \param rap product-center displacement G-A
!> \param la_max highest A degree
!> \param rbp product-center displacement G-B
!> \param lb_max highest B degree
!> \param rcp product-center displacement G-C
!> \param lc_max highest C degree
!> \param zet sum of all three Gaussian exponents
!> \param rr normalized tables indexed by A/B/C degrees and Cartesian axis
! **************************************************************************************************
      SUBROUTINE os_rr_ovlp3_${suffix}$ (rap, la_max, rbp, lb_max, rcp, lc_max, zet, rr)
         ${integral_type}$ (dp), INTENT(IN) :: rap(3), rbp(3), rcp(3)
         INTEGER, INTENT(IN) :: la_max, lb_max, lc_max
         REAL(dp), INTENT(IN) :: zet
         ${integral_type}$ (dp), INTENT(INOUT) :: rr(0:, 0:, 0:, :)

         INTEGER :: a, b, c, idir
         REAL(dp) :: g

         g = 0.5_dp/zet
         DO idir = 1, 3
            CALL os_rr_ovlp_1d_${suffix}$ (rap(idir), la_max, rbp(idir), lb_max, zet, rr(:, :, 0, idir))
            DO c = 0, lc_max - 1
               DO b = 0, lb_max
                  DO a = 0, la_max
                     rr(a, b, c + 1, idir) = rcp(idir)*rr(a, b, c, idir) + &
                                             g*(REAL(a, dp)*rr(MAX(0, a - 1), b, c, idir) + &
                                                REAL(b, dp)*rr(a, MAX(0, b - 1), c, idir) + &
                                                REAL(c, dp)*rr(a, b, MAX(0, c - 1), idir))
                  END DO
               END DO
            END DO
         END DO
      END SUBROUTINE os_rr_ovlp3_${suffix}$
   #:endfor

   #:for layout in ['abm', 'mab']
      #:set idx = (lambda a, b, m: a + ', ' + b + ', ' + m) if layout == 'abm' else (lambda a, b, m: m + ', ' + a + ', ' + b)
      #:set bounds = ':, :, 0:' if layout == 'abm' else '0:, :, :'
      #:set outer, inner = ('m', 'ia') if layout == 'abm' else ('ia', 'm')
      #:set starts = {'m': '0', 'ia': 'ncoset(la - 1) + 1'}
      #:set stops = {'m': 'mmax - la', 'ia': 'ncoset(la)'}
! **************************************************************************************************
!> \brief Propagate caller-supplied radial seeds through the A/B angular recurrence.
!>        Both supported layouts use the same algebra and require no transpose.
!> \param rap displacement P-A
!> \param la_max highest A shell
!> \param rbp displacement P-B
!> \param lb_max highest B shell
!> \param rpw displacement multiplying the next radial order
!> \param zetp reciprocal exponent sum on A/B
!> \param eta reduction factor of the next radial order (one for point Coulomb)
!> \param mmax highest supplied seed order, including subsequent C recurrence
!> \param rr caller-initialized s/s seeds and resulting angular table
! **************************************************************************************************
      SUBROUTINE os_rr_radial_${layout}$ (rap, la_max, rbp, lb_max, rpw, zetp, eta, mmax, rr)
         REAL(dp), INTENT(IN) :: rap(3), rbp(3), rpw(3), zetp, eta
         INTEGER, INTENT(IN) :: la_max, lb_max, mmax
         REAL(dp), INTENT(INOUT) :: rr(${bounds}$)

         INTEGER :: axis, degree, ia, ia1, ia2, ib, ib1, ib2, ico, iinner, iouter, la, lb, m, power(3)
         INTEGER :: lowered(3, ncoset(MAX(la_max, lb_max)))
         REAL(dp) :: g

         g = 0.5_dp*zetp
         DO ico = 1, ncoset(MAX(la_max, lb_max))
            power = indco(:, ico)
            lowered(:, ico) = [coset(MAX(0, power(1) - 1), power(2), power(3)), &
                               coset(power(1), MAX(0, power(2) - 1), power(3)), &
                               coset(power(1), power(2), MAX(0, power(3) - 1))]
         END DO
         DO la = 1, la_max
            DO iouter = ${starts[outer]}$, ${stops[outer]}$
            DO iinner = ${starts[inner]}$, ${stops[inner]}$
               ${outer}$ = iouter
               ${inner}$ = iinner
               axis = 3
               IF (indco(3, ia) == 0) axis = 2
               IF (indco(axis, ia) == 0) axis = 1
               ia1 = lowered(axis, ia)
               ia2 = lowered(axis, ia1)
               degree = indco(axis, ia) - 1
               rr(${idx('ia', '1', 'm')}$) = rap(axis)*rr(${idx('ia1', '1', 'm')}$) + &
                                             rpw(axis)*rr(${idx('ia1', '1', 'm + 1')}$)
               IF (degree > 0) THEN
                  rr(${idx('ia', '1', 'm')}$) = rr(${idx('ia', '1', 'm')}$) + g*REAL(degree, dp)* &
                                                (rr(${idx('ia2', '1', 'm')}$) - eta*rr(${idx('ia2', '1', 'm + 1')}$))
               END IF
            END DO
            END DO
         END DO
         DO lb = 1, lb_max
            DO ib = ncoset(lb - 1) + 1, ncoset(lb)
               axis = 3
               IF (indco(3, ib) == 0) axis = 2
               IF (indco(axis, ib) == 0) axis = 1
               ib1 = lowered(axis, ib)
               ib2 = lowered(axis, ib1)
               degree = indco(axis, ib) - 1
               DO la = 0, la_max
                  #:set stops = {'m': 'mmax - la - lb', 'ia': 'ncoset(la)'}
                  DO iouter = ${starts[outer]}$, ${stops[outer]}$
                  DO iinner = ${starts[inner]}$, ${stops[inner]}$
                     ${outer}$ = iouter
                     ${inner}$ = iinner
                     ia1 = lowered(axis, ia)
                     rr(${idx('ia', 'ib', 'm')}$) = rbp(axis)*rr(${idx('ia', 'ib1', 'm')}$) + &
                                                    rpw(axis)*rr(${idx('ia', 'ib1', 'm + 1')}$)
                     IF (degree > 0) THEN
                        rr(${idx('ia', 'ib', 'm')}$) = rr(${idx('ia', 'ib', 'm')}$) + g*REAL(degree, dp)* &
                                                       (rr(${idx('ia', 'ib2', 'm')}$) - eta*rr(${idx('ia', 'ib2', 'm + 1')}$))
                     END IF
                     IF (indco(axis, ia) > 0) THEN
                        rr(${idx('ia', 'ib', 'm')}$) = rr(${idx('ia', 'ib', 'm')}$) + g*REAL(indco(axis, ia), dp)* &
                                                       (rr(${idx('ia1', 'ib1', 'm')}$) - eta*rr(${idx('ia1', 'ib1', 'm + 1')}$))
                     END IF
                  END DO
                  END DO
               END DO
            END DO
         END DO
      END SUBROUTINE os_rr_radial_${layout}$
   #:endfor

! **************************************************************************************************
!> \brief Angular three-center Coulomb recurrence from caller-supplied s/s/s seeds.
!> \param rap displacement P-A
!> \param la_max highest A shell
!> \param rbp displacement P-B
!> \param lb_max highest B shell
!> \param rcp displacement P-C
!> \param lc_max highest C shell
!> \param zetp reciprocal exponent sum on A/B
!> \param zetq reciprocal exponent on C
!> \param zetw reciprocal exponent sum on A/B/C
!> \param rr caller-initialized seeds and angular table, with radial order last
! **************************************************************************************************
   SUBROUTINE os_rr_coulomb3(rap, la_max, rbp, lb_max, rcp, lc_max, zetp, zetq, zetw, rr)
      REAL(dp), INTENT(IN) :: rap(3), rbp(3), rcp(3), zetp, zetq, zetw
      INTEGER, INTENT(IN) :: la_max, lb_max, lc_max
      REAL(dp), INTENT(INOUT) :: rr(:, :, :, 0:)

      INTEGER :: axis, ia, ia1, ib, ib1, ic, ic1, ic2, ico, la, lb, lc, m, mmax, power(3)
      INTEGER :: lowered(3, ncoset(MAX(la_max, lb_max, lc_max)))
      REAL(dp) :: eta, eta_c, g, g_c, rcw(3), rpw(3)

      mmax = la_max + lb_max + lc_max
      eta = zetw/zetq
      eta_c = zetw/zetp
      rpw = -eta*rcp
      ! Direct evaluation avoids subtracting P-C and its nearly equal Coulomb shift.
      rcw = eta_c*rcp
      CALL os_rr_radial(rap, la_max, rbp, lb_max, rpw, zetp, eta, mmax, rr(:, :, 1, :))
      IF (lc_max == 0) RETURN
      g = 0.5_dp*zetw
      g_c = 0.5_dp*zetq
      DO ico = 1, ncoset(MAX(la_max, lb_max, lc_max))
         power = indco(:, ico)
         lowered(:, ico) = [coset(MAX(0, power(1) - 1), power(2), power(3)), &
                            coset(power(1), MAX(0, power(2) - 1), power(3)), &
                            coset(power(1), power(2), MAX(0, power(3) - 1))]
      END DO
      DO lc = 1, lc_max
         DO ic = ncoset(lc - 1) + 1, ncoset(lc)
            axis = 3
            IF (indco(3, ic) == 0) axis = 2
            IF (indco(axis, ic) == 0) axis = 1
            ic1 = lowered(axis, ic)
            ic2 = lowered(axis, ic1)
            DO lb = 0, lb_max
               DO la = 0, la_max
                  DO m = 0, mmax - la - lb - lc
                     DO ib = ncoset(lb - 1) + 1, ncoset(lb)
                        ib1 = lowered(axis, ib)
                        DO ia = ncoset(la - 1) + 1, ncoset(la)
                           ia1 = lowered(axis, ia)
                           rr(ia, ib, ic, m) = rcw(axis)*rr(ia, ib, ic1, m + 1)
                           IF (indco(axis, ic) > 1) THEN
                              rr(ia, ib, ic, m) = rr(ia, ib, ic, m) + g_c*REAL(indco(axis, ic) - 1, dp)* &
                                                  (rr(ia, ib, ic2, m) - eta_c*rr(ia, ib, ic2, m + 1))
                           END IF
                           IF (indco(axis, ia) > 0) THEN
                              rr(ia, ib, ic, m) = rr(ia, ib, ic, m) + &
                                                  g*REAL(indco(axis, ia), dp)*rr(ia1, ib, ic1, m + 1)
                           END IF
                           IF (indco(axis, ib) > 0) THEN
                              rr(ia, ib, ic, m) = rr(ia, ib, ic, m) + &
                                                  g*REAL(indco(axis, ib), dp)*rr(ia, ib1, ic1, m + 1)
                           END IF
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE os_rr_coulomb3

! **************************************************************************************************
!> \brief   Calculation of the Obara-Saika recurrence relation for 1/r_C
!> \param rap ...
!> \param la_max ...
!> \param rbp ...
!> \param lb_max ...
!> \param rcp ...
!> \param zet ...
!> \param ldrr1 ...
!> \param ldrr2 ...
!> \param rr ...
!> \date    02.03.2009
!> \author  VW
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE os_rr_coul(rap, la_max, rbp, lb_max, rcp, zet, ldrr1, ldrr2, rr)
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: rap
      INTEGER, INTENT(IN)                                :: la_max
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: rbp
      INTEGER, INTENT(IN)                                :: lb_max
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: rcp
      REAL(dp), INTENT(IN)                               :: zet
      INTEGER, INTENT(IN)                                :: ldrr1, ldrr2
      REAL(dp), DIMENSION(0:ldrr1 - 1, ldrr2, *), &
         INTENT(INOUT)                                   :: rr

      INTEGER :: mmax
      REAL(dp) :: t

      mmax = la_max + lb_max
      t = zet*DOT_PRODUCT(rcp, rcp)
      CALL fgamma(mmax, t, rr(0:mmax, 1, 1))
      CALL os_rr_radial_mab(rap, la_max, rbp, lb_max, -rcp, 1.0_dp/zet, 1.0_dp, mmax, &
                            rr(0:mmax, 1:ncoset(la_max), 1:ncoset(lb_max)))

   END SUBROUTINE os_rr_coul

END MODULE ai_os_rr
