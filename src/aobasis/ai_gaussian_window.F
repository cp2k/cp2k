!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Analytic overlaps weighted by a Cartesian Gaussian spatial window.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE ai_gaussian_window
   USE ai_os_rr,                        ONLY: os_rr_ovlp
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
   USE orbital_pointers,                ONLY: coset,&
                                              ncoset
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_gaussian_window'
   PUBLIC :: gaussian_window_overlap
CONTAINS
! **************************************************************************************************
!> \brief Integrate (a|exp(-sum(alpha*(r-center)**2))|b) using the existing overlap recurrence.
!> \param la_max Maximum angular momentum on A
!> \param la_min Minimum angular momentum on A
!> \param zeta Primitive exponents on A
!> \param lb_max Maximum angular momentum on B
!> \param lb_min Minimum angular momentum on B
!> \param zetb Primitive exponents on B
!> \param ra Cartesian center of A
!> \param rb Cartesian center of B
!> \param center Window center in the same frame
!> \param alpha Nonnegative Cartesian window exponents; zero means no restriction in that direction
!> \param sab Compact primitive integrals, in the same ordering as overlap_ab and scon
! **************************************************************************************************
   SUBROUTINE gaussian_window_overlap(la_max, la_min, zeta, lb_max, lb_min, zetb, ra, rb, center, alpha, sab)
      INTEGER, INTENT(IN)                                :: la_max, la_min
      REAL(dp), INTENT(IN)                               :: zeta(:)
      INTEGER, INTENT(IN)                                :: lb_max, lb_min
      REAL(dp), INTENT(IN)                               :: zetb(:), ra(3), rb(3), center(3), &
                                                            alpha(3)
      REAL(dp), INTENT(OUT)                              :: sab(:, :)

      INTEGER                                            :: ax, axis, ay, az, bx, by, bz, ia, ib, &
                                                            ip, jp, la, lb, ld, na, nb, oa, ob
      REAL(dp)                                           :: a, b, f0, p(3), pair_center(3), rap(3), &
                                                            rbp(3)
      REAL(dp), ALLOCATABLE                              :: rr(:, :, :), work(:, :, :)

      CPASSERT(ALL(alpha >= 0.0_dp) .AND. ALL(zeta > 0.0_dp) .AND. ALL(zetb > 0.0_dp))
      oa = ncoset(la_min - 1)
      ob = ncoset(lb_min - 1)
      na = ncoset(la_max) - oa
      nb = ncoset(lb_max) - ob
      CPASSERT(SIZE(sab, 1) >= na*SIZE(zeta) .AND. SIZE(sab, 2) >= nb*SIZE(zetb))
      ld = MAX(la_max, lb_max) + 1
      ALLOCATE (rr(0:ld - 1, 0:ld - 1, 3), work(0:ld - 1, 0:ld - 1, 3))
      sab = 0.0_dp
      DO jp = 1, SIZE(zetb)
         b = zetb(jp)
         DO ip = 1, SIZE(zeta)
            a = zeta(ip)
            ! Complete the square separately in each direction; the window is anisotropic.
            pair_center = ra + b*(rb - ra)/(a + b)
            p = a + b + alpha
            rap = b*(rb - ra)/(a + b) + alpha*(center - pair_center)/p
            rbp = rap + ra - rb
            f0 = PRODUCT(SQRT(pi/p))*EXP(-a*b/(a + b)*SUM((ra - rb)**2) - &
                                         SUM((a + b)*alpha/p*(pair_center - center)**2))
            DO axis = 1, 3
               CALL os_rr_ovlp(rap, la_max, rbp, lb_max, p(axis), ld, work)
               rr(:, :, axis) = work(:, :, axis)
            END DO
            DO lb = lb_min, lb_max
               DO bx = 0, lb
                  DO by = 0, lb - bx
                     bz = lb - bx - by
                     ib = (jp - 1)*nb + coset(bx, by, bz) - ob
                     DO la = la_min, la_max
                        DO ax = 0, la
                           DO ay = 0, la - ax
                              az = la - ax - ay
                              ia = (ip - 1)*na + coset(ax, ay, az) - oa
                              sab(ia, ib) = f0*rr(ax, bx, 1)*rr(ay, by, 2)*rr(az, bz, 3)
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE gaussian_window_overlap
END MODULE ai_gaussian_window
