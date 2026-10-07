!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of the angular momentum integrals over
!>      Cartesian Gaussian-type functions.
!> \par Literature
!>      S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!> \par History
!>      none
!> \author JGH (20.02.2005)
! **************************************************************************************************
MODULE ai_angmom

   USE ai_cartesian_operators,          ONLY: apply_cartesian_1d
   USE ai_os_rr,                        ONLY: os_rr_ovlp
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
   USE orbital_pointers,                ONLY: indco,&
                                              ncoset
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Public subroutines ***

   PUBLIC :: angmom

CONTAINS

! **************************************************************************************************
!> \brief Real angular-momentum kernel (r-C) cross grad, equal to i times L.
!> \param la_max ...
!> \param npgfa ...
!> \param zeta ...
!> \param rpgfa ...
!> \param la_min lowest occupied shell; the output retains the full cumulative Cartesian layout
!> \param lb_max ...
!> \param npgfb ...
!> \param zetb ...
!> \param rpgfb ...
!> \param rac ...
!> \param rbc ...
!> \param angab ...
! **************************************************************************************************
   SUBROUTINE angmom(la_max, npgfa, zeta, rpgfa, la_min, &
                     lb_max, npgfb, zetb, rpgfb, rac, rbc, angab)
      INTEGER, INTENT(IN)                                :: la_max, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta, rpgfa
      INTEGER, INTENT(IN)                                :: la_min, lb_max, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb, rpgfb
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rac, rbc
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)     :: angab

      INTEGER                                            :: avec(3), bvec(3), coa, cob, ia, ib, &
                                                            idir, ipgf, jdir, jpgf, kdir, ldrr, &
                                                            ma, mb
      REAL(KIND=dp)                                      :: a, b, dab, f0, rab(3), rab2, rap(3), &
                                                            rbp(3), zet
      REAL(KIND=dp), ALLOCATABLE                         :: drr(:, :, :), position(:, :, :), &
                                                            rr(:, :, :)

      MARK_USED(la_min)
      rab = rbc - rac
      dab = NORM2(rab)
      rab2 = dab**2
      ldrr = MAX(la_max, lb_max) + 2
      ALLOCATE (rr(0:ldrr - 1, 0:ldrr - 1, 3), position(0:ldrr - 1, 0:ldrr - 1, 3), &
                drr(0:ldrr - 1, 0:ldrr - 1, 3))
      angab = 0.0_dp
      DO ipgf = 1, npgfa
         ma = (ipgf - 1)*ncoset(la_max)
         DO jpgf = 1, npgfb
            IF (rpgfa(ipgf) + rpgfb(jpgf) < dab) CYCLE
            mb = (jpgf - 1)*ncoset(lb_max)
            a = zeta(ipgf)
            b = zetb(jpgf)
            zet = a + b
            rap = b*rab/zet
            rbp = -a*rab/zet
            f0 = (pi/zet)**1.5_dp*EXP(-a*b*rab2/zet)
            CALL os_rr_ovlp(rap, la_max + 1, rbp, lb_max + 1, zet, ldrr, rr)
            CALL apply_cartesian_1d(rr, position, la_max, lb_max, .FALSE., .FALSE., 0.0_dp, rac)
            CALL apply_cartesian_1d(rr, drr, la_max, lb_max, .TRUE., .TRUE., b, [0.0_dp, 0.0_dp, 0.0_dp])

            ! The ket electronic derivative is the negative of its center derivative.
            ! Form only the three cross-product components before contraction.
            DO cob = 1, ncoset(lb_max)
               bvec = indco(:, cob)
               ib = mb + cob
               DO coa = 1, ncoset(la_max)
                  avec = indco(:, coa)
                  ia = ma + coa
                  DO idir = 1, 3
                     jdir = MOD(idir, 3) + 1
                     kdir = MOD(idir + 1, 3) + 1
                     angab(ia, ib, idir) = f0*rr(avec(idir), bvec(idir), idir)* &
                                           (position(avec(kdir), bvec(kdir), kdir)*drr(avec(jdir), bvec(jdir), jdir) - &
                                            position(avec(jdir), bvec(jdir), jdir)*drr(avec(kdir), bvec(kdir), kdir))
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE angmom

END MODULE ai_angmom
