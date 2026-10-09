!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of the kinetic energy integrals over Cartesian
!>      Gaussian-type functions.
!>
!>      [a|T|b] = [a|-nabla**2/2|b]
!> \par Literature
!>      S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!> \par History
!>      - Derivatives added (10.05.2002,MK)
!>      - Fully refactored (07.07.2014,JGH)
!> \author Matthias Krack (31.07.2000)
! **************************************************************************************************
MODULE ai_kinetic
   USE ai_cartesian_operators,          ONLY: apply_cartesian_1d
   USE ai_os_rr,                        ONLY: os_rr_ovlp
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
   USE orbital_pointers,                ONLY: indco,&
                                              ncoset
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_kinetic'

! *** Public subroutines ***

   PUBLIC :: kinetic

CONTAINS

! **************************************************************************************************
!> \brief   Calculation of the two-center kinetic energy integrals [a|T|b] over
!>          Cartesian Gaussian-type functions.
!> \param la_max Maximum L of basis on A
!> \param la_min Minimum L of basis on A
!> \param npgfa  Number of primitive functions in set of basis on A
!> \param rpgfa  Range of functions on A (used for prescreening)
!> \param zeta   Exponents of basis on center A
!> \param lb_max Maximum L of basis on A
!> \param lb_min Minimum L of basis on A
!> \param npgfb  Number of primitive functions in set of basis on B
!> \param rpgfb  Range of functions on B (used for prescreening)
!> \param zetb   Exponents of basis on center B
!> \param rab    Distance vector between centers A and B
!> \param kab    Kinetic energy integrals, optional
!> \param dab    First derivatives of Kinetic energy integrals, optional
!> \date    07.07.2014
!> \author  JGH
! **************************************************************************************************
   SUBROUTINE kinetic(la_max, la_min, npgfa, rpgfa, zeta, &
                      lb_max, lb_min, npgfb, rpgfb, zetb, &
                      rab, kab, dab)
      INTEGER, INTENT(IN)                                :: la_max, la_min, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rpgfa, zeta
      INTEGER, INTENT(IN)                                :: lb_max, lb_min, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rpgfb, zetb
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT), &
         OPTIONAL                                        :: kab
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), &
         OPTIONAL                                        :: dab

      INTEGER                                            :: avec(3), bvec(3), coa, cob, ia, ib, &
                                                            idir, ipgf, jdir, jpgf, kdir, ldrr, &
                                                            lma, lmb, ma, mb, na, nb, ofa, ofb
      REAL(KIND=dp)                                      :: a, b, f0, rab2, tab, xhi, zet
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: drr, dtt, rr, tt
      REAL(KIND=dp), DIMENSION(3)                        :: rap, rbp

      CPASSERT(PRESENT(kab) .OR. PRESENT(dab))

! Distance of the centers a and b

      rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
      tab = SQRT(rab2)

      ! Maximum l for auxiliary integrals
      IF (PRESENT(kab)) THEN
         lma = la_max + 1
         lmb = lb_max + 1
      END IF
      IF (PRESENT(dab)) THEN
         lma = la_max + 2
         lmb = lb_max + 1
      END IF
      ldrr = MAX(lma, lmb) + 1

      ! Allocate space for auxiliary integrals
      ALLOCATE (rr(0:ldrr - 1, 0:ldrr - 1, 3), tt(0:ldrr - 1, 0:ldrr - 1, 3), &
                drr(0:ldrr - 1, 0:ldrr - 1, 3))
      IF (PRESENT(dab)) ALLOCATE (dtt(0:ldrr - 1, 0:ldrr - 1, 3))

      ! Number of integrals, check size of arrays
      ofa = ncoset(la_min - 1)
      ofb = ncoset(lb_min - 1)
      na = ncoset(la_max) - ofa
      nb = ncoset(lb_max) - ofb
      IF (PRESENT(kab)) THEN
         CPASSERT((SIZE(kab, 1) >= na*npgfa))
         CPASSERT((SIZE(kab, 2) >= nb*npgfb))
      END IF
      IF (PRESENT(dab)) THEN
         CPASSERT((SIZE(dab, 1) >= na*npgfa))
         CPASSERT((SIZE(dab, 2) >= nb*npgfb))
         CPASSERT((SIZE(dab, 3) >= 3))
      END IF

      ! Loops over all pairs of primitive Gaussian-type functions
      ma = 0
      DO ipgf = 1, npgfa
         mb = 0
         DO jpgf = 1, npgfb
            ! Distance Screening
            IF (rpgfa(ipgf) + rpgfb(jpgf) < tab) THEN
               IF (PRESENT(kab)) kab(ma + 1:ma + na, mb + 1:mb + nb) = 0.0_dp
               IF (PRESENT(dab)) dab(ma + 1:ma + na, mb + 1:mb + nb, 1:3) = 0.0_dp
               mb = mb + nb
               CYCLE
            END IF

            ! Calculate some prefactors
            a = zeta(ipgf)
            b = zetb(jpgf)
            zet = a + b
            xhi = a*b/zet
            rap = b*rab/zet
            rbp = -a*rab/zet

            ! [s|s] integral
            f0 = 0.5_dp*(pi/zet)**(1.5_dp)*EXP(-xhi*rab2)

            ! Calculate the recurrence relation, overlap
            CALL os_rr_ovlp(rap, lma, rbp, lmb, zet, ldrr, rr)

            ! Mixed center derivatives give the overlap of the two electronic gradients.
            CALL apply_cartesian_1d(rr, drr, lma - 1, lmb, .FALSE., .TRUE., a, [0.0_dp, 0.0_dp, 0.0_dp])
            CALL apply_cartesian_1d(drr, tt, lma - 1, lmb - 1, .TRUE., .TRUE., b, [0.0_dp, 0.0_dp, 0.0_dp])
            IF (PRESENT(dab)) THEN
               CALL apply_cartesian_1d(tt, dtt, la_max, lb_max, .FALSE., .TRUE., a, [0.0_dp, 0.0_dp, 0.0_dp])
            END IF

            DO cob = ofb + 1, ncoset(lb_max)
               bvec = indco(:, cob)
               ib = mb + cob - ofb
               DO coa = ofa + 1, ncoset(la_max)
                  avec = indco(:, coa)
                  ia = ma + coa - ofa
                  IF (PRESENT(kab)) kab(ia, ib) = 0.0_dp
                  DO idir = 1, 3
                     jdir = MOD(idir, 3) + 1
                     kdir = MOD(idir + 1, 3) + 1
                     IF (PRESENT(kab)) THEN
                        kab(ia, ib) = kab(ia, ib) + tt(avec(idir), bvec(idir), idir)* &
                                      rr(avec(jdir), bvec(jdir), jdir)*rr(avec(kdir), bvec(kdir), kdir)
                     END IF
                     IF (PRESENT(dab)) THEN
                        dab(ia, ib, idir) = dtt(avec(idir), bvec(idir), idir)* &
                                            rr(avec(jdir), bvec(jdir), jdir)*rr(avec(kdir), bvec(kdir), kdir) + &
                                            drr(avec(idir), bvec(idir), idir)* &
                                            (tt(avec(jdir), bvec(jdir), jdir)*rr(avec(kdir), bvec(kdir), kdir) + &
                                             rr(avec(jdir), bvec(jdir), jdir)*tt(avec(kdir), bvec(kdir), kdir))
                     END IF
                  END DO
                  IF (PRESENT(kab)) kab(ia, ib) = f0*kab(ia, ib)
                  IF (PRESENT(dab)) dab(ia, ib, 1:3) = f0*dab(ia, ib, 1:3)
               END DO
            END DO

            mb = mb + nb
         END DO
         ma = ma + na
      END DO

      DEALLOCATE (rr, tt)

   END SUBROUTINE kinetic

END MODULE ai_kinetic
