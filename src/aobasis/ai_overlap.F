!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of the overlap integrals over Cartesian Gaussian-type
!>      functions.
!> \par Literature
!>      S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!> \par History
!>      - Derivatives added (02.05.2002,MK)
!>      - New OS routine with simpler logic (11.07.2014, JGH)
!> \author Matthias Krack (08.10.1999)
! **************************************************************************************************
MODULE ai_overlap
   USE ai_os_rr,                        ONLY: os_rr_ovlp
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi,&
                                              twopi,&
                                              z_one
   USE orbital_pointers,                ONLY: coset,&
                                              nco,&
                                              ncoset,&
                                              nso
   USE orbital_transformation_matrices, ONLY: orbtramat
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_overlap'

! *** Public subroutines ***
   PUBLIC :: overlap, overlap_ab, overlap_aab, overlap_ab_s, overlap_ab_sp, &
             overlap_abb

CONTAINS

! **************************************************************************************************
!> \brief   Purpose: Calculation of the two-center overlap integrals [a|b] over
!>          Cartesian Gaussian-type functions.
!> \param la_max_set Max L on center A
!> \param la_min_set Min L on center A
!> \param npgfa      Number of primitives on center A
!> \param rpgfa      Range of functions on A, used for screening
!> \param zeta       Exponents on center A
!> \param lb_max_set Max L on center B
!> \param lb_min_set Min L on center B
!> \param npgfb      Number of primitives on center B
!> \param rpgfb      Range of functions on B, used for screening
!> \param zetb       Exponents on center B
!> \param rab        Distance vector A-B
!> \param dab        Distance A-B
!> \param sab        Final Integrals, basic and derivatives
!> \param da_max_set Some additional derivative information
!> \param return_derivatives   Return integral derivatives
!> \param s          Work space
!> \param lds        Leading dimension of s
!> \date    19.09.2000
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE overlap(la_max_set, la_min_set, npgfa, rpgfa, zeta, &
                      lb_max_set, lb_min_set, npgfb, rpgfb, zetb, &
                      rab, dab, sab, da_max_set, return_derivatives, s, lds)
      INTEGER, INTENT(IN)                                :: la_max_set, la_min_set, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rpgfa, zeta
      INTEGER, INTENT(IN)                                :: lb_max_set, lb_min_set, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rpgfb, zetb
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), INTENT(IN)                          :: dab
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: sab
      INTEGER, INTENT(IN)                                :: da_max_set
      LOGICAL, INTENT(IN)                                :: return_derivatives
      INTEGER, INTENT(IN)                                :: lds
      REAL(KIND=dp), DIMENSION(lds, lds, *), &
         INTENT(INOUT)                                   :: s

      INTEGER :: ax, ay, az, bx, by, bz, cda, cdax, cday, cdaz, coa, coamx, coamy, coamz, coapx, &
         coapy, coapz, cob, da, da_max, dax, day, daz, i, ipgf, j, jk, jpgf, jstart, k, la, &
         la_max, la_start, lb, lb_max, lb_start, ldrr, na, nb
      REAL(KIND=dp)                                      :: f0, fax, fay, faz, ftz, zetp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: rr
      REAL(KIND=dp), DIMENSION(3)                        :: rap, rbp

      da_max = da_max_set
      la_max = la_max_set + da_max_set

      lb_max = lb_max_set
      ldrr = MAX(la_max, lb_max) + 1
      ALLOCATE (rr(0:ldrr - 1, 0:ldrr - 1, 3))

!   *** Loop over all pairs of primitive Gaussian-type functions ***

      na = 0
      DO ipgf = 1, npgfa

         nb = 0

         DO jpgf = 1, npgfb

!       *** Screening ***

            IF (rpgfa(ipgf) + rpgfb(jpgf) < dab) THEN
               DO j = nb + 1, nb + ncoset(lb_max_set)
                  DO i = na + 1, na + ncoset(la_max_set)
                     sab(i, j) = 0.0_dp
                  END DO
               END DO
               IF (return_derivatives) THEN
                  DO k = 2, ncoset(da_max_set)
                     jstart = (k - 1)*SIZE(sab, 1)
                     DO j = jstart + nb + 1, jstart + nb + ncoset(lb_max_set)
                        DO i = na + 1, na + ncoset(la_max_set)
                           sab(i, j) = 0.0_dp
                        END DO
                     END DO
                  END DO
               END IF
               nb = nb + ncoset(lb_max_set)
               CYCLE
            END IF

!       *** Calculate some prefactors ***

            zetp = 1.0_dp/(zeta(ipgf) + zetb(jpgf))

            f0 = SQRT((pi*zetp)**3)*EXP(-zeta(ipgf)*zetb(jpgf)*zetp*dab*dab)
            rap(:) = zetb(jpgf)*zetp*rab(:)
            rbp(:) = -zeta(ipgf)*zetp*rab(:)

            CALL os_rr_ovlp(rap, la_max, rbp, lb_max, 1.0_dp/zetp, ldrr, rr)

            DO lb = 0, lb_max
               DO bx = 0, lb
                  DO by = 0, lb - bx
                     bz = lb - bx - by
                     cob = coset(bx, by, bz)
                     DO la = 0, la_max
                        DO ax = 0, la
                           DO ay = 0, la - ax
                              az = la - ax - ay
                              coa = coset(ax, ay, az)
                              s(coa, cob, 1) = f0*rr(ax, bx, 1)*rr(ay, by, 2)*rr(az, bz, 3)
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO

!       *** Store the primitive overlap integrals ***

            DO j = 1, ncoset(lb_max_set)
               DO i = 1, ncoset(la_max_set)
                  sab(na + i, nb + j) = s(i, j, 1)
               END DO
            END DO

!       *** Calculate the requested derivatives with respect  ***
!       *** to the nuclear coordinates of the atomic center a ***

            IF (return_derivatives) THEN
               la_start = 0
               lb_start = 0
            ELSE
               la_start = la_min_set
               lb_start = lb_min_set
            END IF

            DO da = 0, da_max - 1
               ftz = 2.0_dp*zeta(ipgf)
               DO dax = 0, da
                  DO day = 0, da - dax
                     daz = da - dax - day
                     cda = coset(dax, day, daz)
                     cdax = coset(dax + 1, day, daz)
                     cday = coset(dax, day + 1, daz)
                     cdaz = coset(dax, day, daz + 1)

!             *** [da/dAi|b] = 2*zeta*[a+1i|b] - Ni(a)[a-1i|b] ***

                     DO la = la_start, la_max - da - 1
                        DO ax = 0, la
                           fax = REAL(ax, dp)
                           DO ay = 0, la - ax
                              fay = REAL(ay, dp)
                              az = la - ax - ay
                              faz = REAL(az, dp)
                              coa = coset(ax, ay, az)
                              coamx = coset(ax - 1, ay, az)
                              coamy = coset(ax, ay - 1, az)
                              coamz = coset(ax, ay, az - 1)
                              coapx = coset(ax + 1, ay, az)
                              coapy = coset(ax, ay + 1, az)
                              coapz = coset(ax, ay, az + 1)
                              DO lb = lb_start, lb_max_set
                                 DO bx = 0, lb
                                    DO by = 0, lb - bx
                                       bz = lb - bx - by
                                       cob = coset(bx, by, bz)
                                       s(coa, cob, cdax) = ftz*s(coapx, cob, cda) - &
                                                           fax*s(coamx, cob, cda)
                                       s(coa, cob, cday) = ftz*s(coapy, cob, cda) - &
                                                           fay*s(coamy, cob, cda)
                                       s(coa, cob, cdaz) = ftz*s(coapz, cob, cda) - &
                                                           faz*s(coamz, cob, cda)
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO

                  END DO
               END DO
            END DO

!       *** Return all the calculated derivatives of the ***
!       *** primitive overlap integrals, if requested    ***

            IF (return_derivatives) THEN
               DO k = 2, ncoset(da_max_set)
                  jstart = (k - 1)*SIZE(sab, 1)
                  DO j = 1, ncoset(lb_max_set)
                     jk = jstart + j
                     DO i = 1, ncoset(la_max_set)
                        sab(na + i, nb + jk) = s(i, j, k)
                     END DO
                  END DO
               END DO
            END IF

            nb = nb + ncoset(lb_max_set)

         END DO

         na = na + ncoset(la_max_set)
      END DO

      DEALLOCATE (rr)

   END SUBROUTINE overlap

! **************************************************************************************************
!> \brief   Calculation of the two-center overlap integrals [a|b] over
!>          Cartesian Gaussian-type functions. First and second derivatives
!> \param la_max     Max L on center A
!> \param la_min     Min L on center A
!> \param npgfa      Number of primitives on center A
!> \param rpgfa      Range of functions on A, used for screening
!> \param zeta       Exponents on center A
!> \param lb_max     Max L on center B
!> \param lb_min     Min L on center B
!> \param npgfb      Number of primitives on center B
!> \param rpgfb      Range of functions on B, used for screening
!> \param zetb       Exponents on center B
!> \param rab        Distance vector A-B
!> \param sab        Final overlap integrals
!> \param dab        First derivative overlap integrals
!> \param ddab       Second derivative overlap integrals
!> \param rr_work    Optional caller-provided one-dimensional overlap recurrence workspace
!> \date    01.07.2014
!> \author  JGH
! **************************************************************************************************
   SUBROUTINE overlap_ab(la_max, la_min, npgfa, rpgfa, zeta, &
                         lb_max, lb_min, npgfb, rpgfb, zetb, &
                         rab, sab, dab, ddab, rr_work)
      INTEGER, INTENT(IN)                                :: la_max, la_min, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rpgfa, zeta
      INTEGER, INTENT(IN)                                :: lb_max, lb_min, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rpgfb, zetb
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT), &
         OPTIONAL                                        :: sab
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), &
         OPTIONAL                                        :: dab, ddab
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT), &
         OPTIONAL, TARGET                                :: rr_work

      INTEGER                                            :: ax, ay, az, bx, by, bz, coa, cob, ia, &
                                                            ib, ipgf, jpgf, la, lb, ldrr, lma, &
                                                            lmb, ma, mb, na, nb, ofa, ofb
      REAL(KIND=dp)                                      :: a, ambm, ambp, apbm, apbp, b, dumx, &
                                                            dumy, dumz, f0, rab2, tab, xhi, zet
      REAL(KIND=dp), DIMENSION(3)                        :: rap, rbp
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: rr

      ! Distance of the centers a and b

      rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
      tab = SQRT(rab2)

      ! Maximum l for auxiliary integrals
      CPASSERT(PRESENT(sab) .OR. PRESENT(dab) .OR. PRESENT(ddab))
      IF (PRESENT(sab)) THEN
         lma = la_max
         lmb = lb_max
      END IF
      IF (PRESENT(dab)) THEN
         lma = la_max + 1
         lmb = lb_max
      END IF
      IF (PRESENT(ddab)) THEN
         lma = la_max + 1
         lmb = lb_max + 1
      END IF
      ldrr = MAX(lma, lmb) + 1

      ! Allocate or attach the workspace for auxiliary integrals
      NULLIFY (rr)
      IF (PRESENT(rr_work)) THEN
         CPASSERT(SIZE(rr_work) >= ldrr*ldrr*3)
         rr(0:ldrr - 1, 0:ldrr - 1, 1:3) => rr_work(1:ldrr*ldrr*3)
      ELSE
         ALLOCATE (rr(0:ldrr - 1, 0:ldrr - 1, 3))
      END IF

      ! Number of integrals, check size of arrays
      ofa = ncoset(la_min - 1)
      ofb = ncoset(lb_min - 1)
      na = ncoset(la_max) - ofa
      nb = ncoset(lb_max) - ofb
      IF (PRESENT(sab)) THEN
         CPASSERT((SIZE(sab, 1) >= na*npgfa))
         CPASSERT((SIZE(sab, 2) >= nb*npgfb))
      END IF
      IF (PRESENT(dab)) THEN
         CPASSERT((SIZE(dab, 1) >= na*npgfa))
         CPASSERT((SIZE(dab, 2) >= nb*npgfb))
         CPASSERT((SIZE(dab, 3) >= 3))
      END IF
      IF (PRESENT(ddab)) THEN
         CPASSERT((SIZE(ddab, 1) >= na*npgfa))
         CPASSERT((SIZE(ddab, 2) >= nb*npgfb))
         CPASSERT((SIZE(ddab, 3) >= 6))
      END IF

      ! Loops over all pairs of primitive Gaussian-type functions
      ma = 0
      DO ipgf = 1, npgfa
         mb = 0
         DO jpgf = 1, npgfb
            ! Distance Screening
            IF (rpgfa(ipgf) + rpgfb(jpgf) < tab) THEN
               IF (PRESENT(sab)) sab(ma + 1:ma + na, mb + 1:mb + nb) = 0.0_dp
               IF (PRESENT(dab)) dab(ma + 1:ma + na, mb + 1:mb + nb, 1:3) = 0.0_dp
               IF (PRESENT(ddab)) ddab(ma + 1:ma + na, mb + 1:mb + nb, 1:6) = 0.0_dp
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
            f0 = (pi/zet)**(1.5_dp)*EXP(-xhi*rab2)

            ! Calculate the recurrence relation
            CALL os_rr_ovlp(rap, lma, rbp, lmb, zet, ldrr, rr)

            DO lb = lb_min, lb_max
            DO bx = 0, lb
            DO by = 0, lb - bx
               bz = lb - bx - by
               cob = coset(bx, by, bz) - ofb
               ib = mb + cob
               DO la = la_min, la_max
               DO ax = 0, la
               DO ay = 0, la - ax
                  az = la - ax - ay
                  coa = coset(ax, ay, az) - ofa
                  ia = ma + coa
                  ! integrals
                  IF (PRESENT(sab)) THEN
                     sab(ia, ib) = f0*rr(ax, bx, 1)*rr(ay, by, 2)*rr(az, bz, 3)
                  END IF
                  ! first derivatives
                  IF (PRESENT(dab)) THEN
                     ! (da|b) = 2*a*(a+1|b) - N(a)*(a-1|b)
                     ! dx
                     dumx = 2.0_dp*a*rr(ax + 1, bx, 1)
                     IF (ax > 0) dumx = dumx - REAL(ax, dp)*rr(ax - 1, bx, 1)
                     dab(ia, ib, 1) = f0*dumx*rr(ay, by, 2)*rr(az, bz, 3)
                     ! dy
                     dumy = 2.0_dp*a*rr(ay + 1, by, 2)
                     IF (ay > 0) dumy = dumy - REAL(ay, dp)*rr(ay - 1, by, 2)
                     dab(ia, ib, 2) = f0*rr(ax, bx, 1)*dumy*rr(az, bz, 3)
                     ! dz
                     dumz = 2.0_dp*a*rr(az + 1, bz, 3)
                     IF (az > 0) dumz = dumz - REAL(az, dp)*rr(az - 1, bz, 3)
                     dab(ia, ib, 3) = f0*rr(ax, bx, 1)*rr(ay, by, 2)*dumz
                  END IF
                  ! 2nd derivatives
                  IF (PRESENT(ddab)) THEN
                     ! (dda|b) = -4*a*b*(a+1|b+1) + 2*a*N(b)*(a+1|b-1)
                     !           + 2*b*N(a)*(a-1|b+1) - N(a)*N(b)*(a-1|b-1)
                     ! dx dx
                     apbp = f0*rr(ax + 1, bx + 1, 1)*rr(ay, by, 2)*rr(az, bz, 3)
                     IF (bx > 0) THEN
                        apbm = f0*rr(ax + 1, bx - 1, 1)*rr(ay, by, 2)*rr(az, bz, 3)
                     ELSE
                        apbm = 0.0_dp
                     END IF
                     IF (ax > 0) THEN
                        ambp = f0*rr(ax - 1, bx + 1, 1)*rr(ay, by, 2)*rr(az, bz, 3)
                     ELSE
                        ambp = 0.0_dp
                     END IF
                     IF (ax > 0 .AND. bx > 0) THEN
                        ambm = f0*rr(ax - 1, bx - 1, 1)*rr(ay, by, 2)*rr(az, bz, 3)
                     ELSE
                        ambm = 0.0_dp
                     END IF
                     ddab(ia, ib, 1) = -4.0_dp*a*b*apbp + 2.0_dp*a*REAL(bx, dp)*apbm &
                                       + 2.0_dp*b*REAL(ax, dp)*ambp - REAL(ax, dp)*REAL(bx, dp)*ambm
                     ! dx dy
                     apbp = f0*rr(ax + 1, bx, 1)*rr(ay, by + 1, 2)*rr(az, bz, 3)
                     IF (by > 0) THEN
                        apbm = f0*rr(ax + 1, bx, 1)*rr(ay, by - 1, 2)*rr(az, bz, 3)
                     ELSE
                        apbm = 0.0_dp
                     END IF
                     IF (ax > 0) THEN
                        ambp = f0*rr(ax - 1, bx, 1)*rr(ay, by + 1, 2)*rr(az, bz, 3)
                     ELSE
                        ambp = 0.0_dp
                     END IF
                     IF (ax > 0 .AND. by > 0) THEN
                        ambm = f0*rr(ax - 1, bx, 1)*rr(ay, by - 1, 2)*rr(az, bz, 3)
                     ELSE
                        ambm = 0.0_dp
                     END IF
                     ddab(ia, ib, 2) = -4.0_dp*a*b*apbp + 2.0_dp*a*REAL(by, dp)*apbm &
                                       + 2.0_dp*b*REAL(ax, dp)*ambp - REAL(ax, dp)*REAL(by, dp)*ambm
                     ! dx dz
                     apbp = f0*rr(ax + 1, bx, 1)*rr(ay, by, 2)*rr(az, bz + 1, 3)
                     IF (bz > 0) THEN
                        apbm = f0*rr(ax + 1, bx, 1)*rr(ay, by, 2)*rr(az, bz - 1, 3)
                     ELSE
                        apbm = 0.0_dp
                     END IF
                     IF (ax > 0) THEN
                        ambp = f0*rr(ax - 1, bx, 1)*rr(ay, by, 2)*rr(az, bz + 1, 3)
                     ELSE
                        ambp = 0.0_dp
                     END IF
                     IF (ax > 0 .AND. bz > 0) THEN
                        ambm = f0*rr(ax - 1, bx, 1)*rr(ay, by, 2)*rr(az, bz - 1, 3)
                     ELSE
                        ambm = 0.0_dp
                     END IF
                     ddab(ia, ib, 3) = -4.0_dp*a*b*apbp + 2.0_dp*a*REAL(bz, dp)*apbm &
                                       + 2.0_dp*b*REAL(ax, dp)*ambp - REAL(ax, dp)*REAL(bz, dp)*ambm
                     ! dy dy
                     apbp = f0*rr(ax, bx, 1)*rr(ay + 1, by + 1, 2)*rr(az, bz, 3)
                     IF (by > 0) THEN
                        apbm = f0*rr(ax, bx, 1)*rr(ay + 1, by - 1, 2)*rr(az, bz, 3)
                     ELSE
                        apbm = 0.0_dp
                     END IF
                     IF (ay > 0) THEN
                        ambp = f0*rr(ax, bx, 1)*rr(ay - 1, by + 1, 2)*rr(az, bz, 3)
                     ELSE
                        ambp = 0.0_dp
                     END IF
                     IF (ay > 0 .AND. by > 0) THEN
                        ambm = f0*rr(ax, bx, 1)*rr(ay - 1, by - 1, 2)*rr(az, bz, 3)
                     ELSE
                        ambm = 0.0_dp
                     END IF
                     ddab(ia, ib, 4) = -4.0_dp*a*b*apbp + 2.0_dp*a*REAL(by, dp)*apbm &
                                       + 2.0_dp*b*REAL(ay, dp)*ambp - REAL(ay, dp)*REAL(by, dp)*ambm
                     ! dy dz
                     apbp = f0*rr(ax, bx, 1)*rr(ay + 1, by, 2)*rr(az, bz + 1, 3)
                     IF (bz > 0) THEN
                        apbm = f0*rr(ax, bx, 1)*rr(ay + 1, by, 2)*rr(az, bz - 1, 3)
                     ELSE
                        apbm = 0.0_dp
                     END IF
                     IF (ay > 0) THEN
                        ambp = f0*rr(ax, bx, 1)*rr(ay - 1, by, 2)*rr(az, bz + 1, 3)
                     ELSE
                        ambp = 0.0_dp
                     END IF
                     IF (ay > 0 .AND. bz > 0) THEN
                        ambm = f0*rr(ax, bx, 1)*rr(ay - 1, by, 2)*rr(az, bz - 1, 3)
                     ELSE
                        ambm = 0.0_dp
                     END IF
                     ddab(ia, ib, 5) = -4.0_dp*a*b*apbp + 2.0_dp*a*REAL(bz, dp)*apbm &
                                       + 2.0_dp*b*REAL(ay, dp)*ambp - REAL(ay, dp)*REAL(bz, dp)*ambm
                     ! dz dz
                     apbp = f0*rr(ax, bx, 1)*rr(ay, by, 2)*rr(az + 1, bz + 1, 3)
                     IF (bz > 0) THEN
                        apbm = f0*rr(ax, bx, 1)*rr(ay, by, 2)*rr(az + 1, bz - 1, 3)
                     ELSE
                        apbm = 0.0_dp
                     END IF
                     IF (az > 0) THEN
                        ambp = f0*rr(ax, bx, 1)*rr(ay, by, 2)*rr(az - 1, bz + 1, 3)
                     ELSE
                        ambp = 0.0_dp
                     END IF
                     IF (az > 0 .AND. bz > 0) THEN
                        ambm = f0*rr(ax, bx, 1)*rr(ay, by, 2)*rr(az - 1, bz - 1, 3)
                     ELSE
                        ambm = 0.0_dp
                     END IF
                     ddab(ia, ib, 6) = -4.0_dp*a*b*apbp + 2.0_dp*a*REAL(bz, dp)*apbm &
                                       + 2.0_dp*b*REAL(az, dp)*ambp - REAL(az, dp)*REAL(bz, dp)*ambm
                  END IF
                  !
               END DO
               END DO
               END DO !la
            END DO
            END DO
            END DO !lb

            mb = mb + nb
         END DO
         ma = ma + na
      END DO

      IF (.NOT. PRESENT(rr_work)) DEALLOCATE (rr)
      NULLIFY (rr)

   END SUBROUTINE overlap_ab

! **************************************************************************************************
!> \brief   Calculation of the two-center overlap integrals [aa|b] over
!>          Cartesian Gaussian-type functions.
!> \param la1_max    Max L on center A (basis 1)
!> \param la1_min    Min L on center A (basis 1)
!> \param npgfa1     Number of primitives on center A (basis 1)
!> \param rpgfa1     Range of functions on A, used for screening (basis 1)
!> \param zeta1      Exponents on center A (basis 1)
!> \param la2_max    Max L on center A (basis 2)
!> \param la2_min    Min L on center A (basis 2)
!> \param npgfa2     Number of primitives on center A (basis 2)
!> \param rpgfa2     Range of functions on A, used for screening (basis 2)
!> \param zeta2      Exponents on center A (basis 2)
!> \param lb_max     Max L on center B
!> \param lb_min     Min L on center B
!> \param npgfb      Number of primitives on center B
!> \param rpgfb      Range of functions on B, used for screening
!> \param zetb       Exponents on center B
!> \param rab        Distance vector A-B
!> \param saab       Final overlap integrals
!> \param daab       First derivative overlap integrals
!> \param saba       Final overlap integrals; different order
!> \param daba       First derivative overlap integrals; different order
!> \date    01.07.2014
!> \author  JGH
! **************************************************************************************************
   SUBROUTINE overlap_aab(la1_max, la1_min, npgfa1, rpgfa1, zeta1, &
                          la2_max, la2_min, npgfa2, rpgfa2, zeta2, &
                          lb_max, lb_min, npgfb, rpgfb, zetb, &
                          rab, saab, daab, saba, daba)
      INTEGER, INTENT(IN)                                :: la1_max, la1_min, npgfa1
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rpgfa1, zeta1
      INTEGER, INTENT(IN)                                :: la2_max, la2_min, npgfa2
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rpgfa2, zeta2
      INTEGER, INTENT(IN)                                :: lb_max, lb_min, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rpgfb, zetb
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), &
         OPTIONAL                                        :: saab
      REAL(KIND=dp), DIMENSION(:, :, :, :), &
         INTENT(INOUT), OPTIONAL                         :: daab
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), &
         OPTIONAL                                        :: saba
      REAL(KIND=dp), DIMENSION(:, :, :, :), &
         INTENT(INOUT), OPTIONAL                         :: daba

      INTEGER :: ax, ax1, ax2, ay, ay1, ay2, az, az1, az2, bx, by, bz, coa1, coa2, cob, i1pgf, &
         i2pgf, ia1, ia2, ib, jpgf, la1, la2, lb, ldrr, lma, lmb, ma1, ma2, mb, na1, na2, nb, &
         ofa1, ofa2, ofb
      REAL(KIND=dp)                                      :: a, b, dumx, dumy, dumz, f0, rab2, rpgfa, &
                                                            tab, xhi, zet
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: rr
      REAL(KIND=dp), DIMENSION(3)                        :: rap, rbp

      ! Distance of the centers a and b

      rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
      tab = SQRT(rab2)

      ! Maximum l for auxiliary integrals
      CPASSERT(PRESENT(saab) .OR. PRESENT(daab) .OR. PRESENT(saba) .OR. PRESENT(daba))
      IF (PRESENT(saab) .OR. PRESENT(saba)) THEN
         lma = la1_max + la2_max
         lmb = lb_max
      END IF
      IF (PRESENT(daab) .OR. PRESENT(daba)) THEN
         lma = la1_max + la2_max + 1
         lmb = lb_max
      END IF
      ldrr = MAX(lma, lmb) + 1

      ! Allocate space for auxiliary integrals
      ALLOCATE (rr(0:ldrr - 1, 0:ldrr - 1, 3))

      ! Number of integrals, check size of arrays
      ofa1 = ncoset(la1_min - 1)
      ofa2 = ncoset(la2_min - 1)
      ofb = ncoset(lb_min - 1)
      na1 = ncoset(la1_max) - ofa1
      na2 = ncoset(la2_max) - ofa2
      nb = ncoset(lb_max) - ofb
      IF (PRESENT(saab)) THEN
         CPASSERT((SIZE(saab, 1) >= na1*npgfa1))
         CPASSERT((SIZE(saab, 2) >= na2*npgfa2))
         CPASSERT((SIZE(saab, 3) >= nb*npgfb))
      END IF
      IF (PRESENT(daab)) THEN
         CPASSERT((SIZE(daab, 1) >= na1*npgfa1))
         CPASSERT((SIZE(daab, 2) >= na2*npgfa2))
         CPASSERT((SIZE(daab, 3) >= nb*npgfb))
         CPASSERT((SIZE(daab, 4) >= 3))
      END IF
      IF (PRESENT(saba)) THEN
         CPASSERT((SIZE(saba, 1) >= na1*npgfa1))
         CPASSERT((SIZE(saba, 2) >= nb*npgfb))
         CPASSERT((SIZE(saba, 3) >= na2*npgfa2))
      END IF
      IF (PRESENT(daba)) THEN
         CPASSERT((SIZE(daba, 1) >= na1*npgfa1))
         CPASSERT((SIZE(daba, 2) >= nb*npgfb))
         CPASSERT((SIZE(daba, 3) >= na2*npgfa2))
         CPASSERT((SIZE(daba, 4) >= 3))
      END IF

      ! Loops over all primitive Gaussian-type functions
      ma1 = 0
      DO i1pgf = 1, npgfa1
         ma2 = 0
         DO i2pgf = 1, npgfa2
            rpgfa = MIN(rpgfa1(i1pgf), rpgfa2(i2pgf))
            mb = 0
            DO jpgf = 1, npgfb
               ! Distance Screening
               IF (rpgfa + rpgfb(jpgf) < tab) THEN
                  IF (PRESENT(saab)) saab(ma1 + 1:ma1 + na1, ma2 + 1:ma2 + na2, mb + 1:mb + nb) = 0.0_dp
                  IF (PRESENT(daab)) daab(ma1 + 1:ma1 + na1, ma2 + 1:ma2 + na2, mb + 1:mb + nb, 1:3) = 0.0_dp
                  IF (PRESENT(saba)) saba(ma1 + 1:ma1 + na1, mb + 1:mb + nb, ma2 + 1:ma2 + na2) = 0.0_dp
                  IF (PRESENT(daba)) daba(ma1 + 1:ma1 + na1, mb + 1:mb + nb, ma2 + 1:ma2 + na2, 1:3) = 0.0_dp
                  mb = mb + nb
                  CYCLE
               END IF

               ! Calculate some prefactors
               a = zeta1(i1pgf) + zeta2(i2pgf)
               b = zetb(jpgf)
               zet = a + b
               xhi = a*b/zet
               rap = b*rab/zet
               rbp = -a*rab/zet

               ! [ss|s] integral
               f0 = (pi/zet)**(1.5_dp)*EXP(-xhi*rab2)

               ! Calculate the recurrence relation
               CALL os_rr_ovlp(rap, lma, rbp, lmb, zet, ldrr, rr)

               DO lb = lb_min, lb_max
               DO bx = 0, lb
               DO by = 0, lb - bx
                  bz = lb - bx - by
                  cob = coset(bx, by, bz) - ofb
                  ib = mb + cob
                  DO la2 = la2_min, la2_max
                  DO ax2 = 0, la2
                  DO ay2 = 0, la2 - ax2
                     az2 = la2 - ax2 - ay2
                     coa2 = coset(ax2, ay2, az2) - ofa2
                     ia2 = ma2 + coa2
                     DO la1 = la1_min, la1_max
                     DO ax1 = 0, la1
                     DO ay1 = 0, la1 - ax1
                        az1 = la1 - ax1 - ay1
                        coa1 = coset(ax1, ay1, az1) - ofa1
                        ia1 = ma1 + coa1
                        ! integrals
                        IF (PRESENT(saab)) THEN
                           saab(ia1, ia2, ib) = f0*rr(ax1 + ax2, bx, 1)*rr(ay1 + ay2, by, 2)*rr(az1 + az2, bz, 3)
                        END IF
                        IF (PRESENT(saba)) THEN
                           saba(ia1, ib, ia2) = f0*rr(ax1 + ax2, bx, 1)*rr(ay1 + ay2, by, 2)*rr(az1 + az2, bz, 3)
                        END IF
                        ! first derivatives
                        IF (PRESENT(daab) .OR. PRESENT(daba)) THEN
                           ax = ax1 + ax2
                           ay = ay1 + ay2
                           az = az1 + az2
                           ! (da|b) = 2*a*(a+1|b) - N(a)*(a-1|b)
                           ! dx
                           dumx = 2.0_dp*a*rr(ax + 1, bx, 1)
                           IF (ax > 0) dumx = dumx - REAL(ax, dp)*rr(ax - 1, bx, 1)
                           dumx = f0*dumx*rr(ay, by, 2)*rr(az, bz, 3)
                           ! dy
                           dumy = 2.0_dp*a*rr(ay + 1, by, 2)
                           IF (ay > 0) dumy = dumy - REAL(ay, dp)*rr(ay - 1, by, 2)
                           dumy = f0*rr(ax, bx, 1)*dumy*rr(az, bz, 3)
                           ! dz
                           dumz = 2.0_dp*a*rr(az + 1, bz, 3)
                           IF (az > 0) dumz = dumz - REAL(az, dp)*rr(az - 1, bz, 3)
                           dumz = f0*rr(ax, bx, 1)*rr(ay, by, 2)*dumz
                           IF (PRESENT(daab)) THEN
                              daab(ia1, ia2, ib, 1) = dumx
                              daab(ia1, ia2, ib, 2) = dumy
                              daab(ia1, ia2, ib, 3) = dumz
                           END IF
                           IF (PRESENT(daba)) THEN
                              daba(ia1, ib, ia2, 1) = dumx
                              daba(ia1, ib, ia2, 2) = dumy
                              daba(ia1, ib, ia2, 3) = dumz
                           END IF
                        END IF
                        !
                     END DO
                     END DO
                     END DO !la1
                  END DO
                  END DO
                  END DO !la2
               END DO
               END DO
               END DO !lb

               mb = mb + nb
            END DO
            ma2 = ma2 + na2
         END DO
         ma1 = ma1 + na1
      END DO

      DEALLOCATE (rr)

   END SUBROUTINE overlap_aab

! **************************************************************************************************
!> \brief   Calculation of the two-center overlap integrals [a|bb] over
!>          Cartesian Gaussian-type functions.
!> \param la_max     Max L on center A
!> \param la_min     Min L on center A
!> \param npgfa      Number of primitives on center A
!> \param rpgfa      Range of functions on A, used for screening
!> \param zeta       Exponents on center A
!> \param lb1_max    Max L on center B (basis 1)
!> \param lb1_min    Min L on center B (basis 1)
!> \param npgfb1     Number of primitives on center B (basis 1)
!> \param rpgfb1     Range of functions on B, used for screening (basis 1)
!> \param zetb1      Exponents on center B (basis 1)
!> \param lb2_max    Max L on center B (basis 2)
!> \param lb2_min    Min L on center B (basis 2)
!> \param npgfb2     Number of primitives on center B (basis 2)
!> \param rpgfb2     Range of functions on B, used for screening (basis 2)
!> \param zetb2      Exponents on center B (basis 2)
!> \param rab        Distance vector A-B
!> \param sabb       Final overlap integrals
!> \param dabb       First derivative overlap integrals
!> \date    01.07.2014
!> \author  JGH
! **************************************************************************************************
   SUBROUTINE overlap_abb(la_max, la_min, npgfa, rpgfa, zeta, &
                          lb1_max, lb1_min, npgfb1, rpgfb1, zetb1, &
                          lb2_max, lb2_min, npgfb2, rpgfb2, zetb2, &
                          rab, sabb, dabb)
      INTEGER, INTENT(IN)                                :: la_max, la_min, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rpgfa, zeta
      INTEGER, INTENT(IN)                                :: lb1_max, lb1_min, npgfb1
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rpgfb1, zetb1
      INTEGER, INTENT(IN)                                :: lb2_max, lb2_min, npgfb2
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rpgfb2, zetb2
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), &
         OPTIONAL                                        :: sabb
      REAL(KIND=dp), DIMENSION(:, :, :, :), &
         INTENT(INOUT), OPTIONAL                         :: dabb

      INTEGER :: ax, ay, az, bx, bx1, bx2, by, by1, by2, bz, bz1, bz2, coa, cob1, cob2, ia, ib1, &
         ib2, ipgf, j1pgf, j2pgf, la, lb1, lb2, ldrr, lma, lmb, ma, mb1, mb2, na, nb1, nb2, ofa, &
         ofb1, ofb2
      REAL(KIND=dp)                                      :: a, b, dumx, dumy, dumz, f0, rab2, rpgfb, &
                                                            tab, xhi, zet
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: rr
      REAL(KIND=dp), DIMENSION(3)                        :: rap, rbp

      ! Distance of the centers a and b

      rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
      tab = SQRT(rab2)

      ! Maximum l for auxiliary integrals
      CPASSERT(PRESENT(sabb) .OR. PRESENT(dabb))
      IF (PRESENT(sabb)) THEN
         lma = la_max
         lmb = lb1_max + lb2_max
      END IF
      IF (PRESENT(dabb)) THEN
         lma = la_max + 1
         lmb = lb1_max + lb2_max
      END IF
      ldrr = MAX(lma, lmb) + 1

      ! Allocate space for auxiliary integrals
      ALLOCATE (rr(0:ldrr - 1, 0:ldrr - 1, 3))

      ! Number of integrals, check size of arrays
      ofa = ncoset(la_min - 1)
      ofb1 = ncoset(lb1_min - 1)
      ofb2 = ncoset(lb2_min - 1)
      na = ncoset(la_max) - ofa
      nb1 = ncoset(lb1_max) - ofb1
      nb2 = ncoset(lb2_max) - ofb2
      IF (PRESENT(sabb)) THEN
         CPASSERT((SIZE(sabb, 1) >= na*npgfa))
         CPASSERT((SIZE(sabb, 2) >= nb1*npgfb1))
         CPASSERT((SIZE(sabb, 3) >= nb2*npgfb2))
      END IF
      IF (PRESENT(dabb)) THEN
         CPASSERT((SIZE(dabb, 1) >= na*npgfa))
         CPASSERT((SIZE(dabb, 2) >= nb1*npgfb1))
         CPASSERT((SIZE(dabb, 3) >= nb2*npgfb2))
         CPASSERT((SIZE(dabb, 4) >= 3))
      END IF

      ! Loops over all pairs of primitive Gaussian-type functions
      ma = 0
      DO ipgf = 1, npgfa
         mb1 = 0
         DO j1pgf = 1, npgfb1
            mb2 = 0
            DO j2pgf = 1, npgfb2
               ! Distance Screening
               rpgfb = MIN(rpgfb1(j1pgf), rpgfb2(j2pgf))
               IF (rpgfa(ipgf) + rpgfb < tab) THEN
                  IF (PRESENT(sabb)) sabb(ma + 1:ma + na, mb1 + 1:mb1 + nb1, mb2 + 1:mb2 + nb2) = 0.0_dp
                  IF (PRESENT(dabb)) dabb(ma + 1:ma + na, mb1 + 1:mb1 + nb1, mb2 + 1:mb2 + nb2, 1:3) = 0.0_dp
                  mb2 = mb2 + nb2
                  CYCLE
               END IF

               ! Calculate some prefactors
               a = zeta(ipgf)
               b = zetb1(j1pgf) + zetb2(j2pgf)
               zet = a + b
               xhi = a*b/zet
               rap = b*rab/zet
               rbp = -a*rab/zet

               ! [s|s] integral
               f0 = (pi/zet)**(1.5_dp)*EXP(-xhi*rab2)

               ! Calculate the recurrence relation
               CALL os_rr_ovlp(rap, lma, rbp, lmb, zet, ldrr, rr)

               DO lb2 = lb2_min, lb2_max
               DO bx2 = 0, lb2
               DO by2 = 0, lb2 - bx2
                  bz2 = lb2 - bx2 - by2
                  cob2 = coset(bx2, by2, bz2) - ofb2
                  ib2 = mb2 + cob2
                  DO lb1 = lb1_min, lb1_max
                  DO bx1 = 0, lb1
                  DO by1 = 0, lb1 - bx1
                     bz1 = lb1 - bx1 - by1
                     cob1 = coset(bx1, by1, bz1) - ofb1
                     ib1 = mb1 + cob1
                     DO la = la_min, la_max
                     DO ax = 0, la
                     DO ay = 0, la - ax
                        az = la - ax - ay
                        coa = coset(ax, ay, az) - ofa
                        ia = ma + coa
                        ! integrals
                        IF (PRESENT(sabb)) THEN
                           sabb(ia, ib1, ib2) = f0*rr(ax, bx1 + bx2, 1)*rr(ay, by1 + by2, 2)*rr(az, bz1 + bz2, 3)
                        END IF
                        ! first derivatives
                        IF (PRESENT(dabb)) THEN
                           bx = bx1 + bx2
                           by = by1 + by2
                           bz = bz1 + bz2
                           ! (da|b) = 2*a*(a+1|b) - N(a)*(a-1|b)
                           ! dx
                           dumx = 2.0_dp*a*rr(ax + 1, bx, 1)
                           IF (ax > 0) dumx = dumx - REAL(ax, dp)*rr(ax - 1, bx, 1)
                           dabb(ia, ib1, ib2, 1) = f0*dumx*rr(ay, by, 2)*rr(az, bz, 3)
                           ! dy
                           dumy = 2.0_dp*a*rr(ay + 1, by, 2)
                           IF (ay > 0) dumy = dumy - REAL(ay, dp)*rr(ay - 1, by, 2)
                           dabb(ia, ib1, ib2, 2) = f0*rr(ax, bx, 1)*dumy*rr(az, bz, 3)
                           ! dz
                           dumz = 2.0_dp*a*rr(az + 1, bz, 3)
                           IF (az > 0) dumz = dumz - REAL(az, dp)*rr(az - 1, bz, 3)
                           dabb(ia, ib1, ib2, 3) = f0*rr(ax, bx, 1)*rr(ay, by, 2)*dumz
                        END IF
                        !
                     END DO
                     END DO
                     END DO !la
                  END DO
                  END DO
                  END DO !lb1
               END DO
               END DO
               END DO !lb2

               mb2 = mb2 + nb2
            END DO
            mb1 = mb1 + nb1
         END DO
         ma = ma + na
      END DO

      DEALLOCATE (rr)

   END SUBROUTINE overlap_abb

! **************************************************************************************************

! **************************************************************************************************
!> \brief   Calculation of the two-center overlap integrals [a|b] over
!>          Spherical Gaussian-type functions.
!> \param la         Max L on center A
!> \param zeta       Exponents on center A
!> \param lb         Max L on center B
!> \param zetb       Exponents on center B
!> \param rab        Distance vector A-B
!> \param sab        Final overlap integrals
!> \date    01.03.2016
!> \author  JGH
! **************************************************************************************************
   SUBROUTINE overlap_ab_s(la, zeta, lb, zetb, rab, sab)
      INTEGER, INTENT(IN)                                :: la
      REAL(KIND=dp), INTENT(IN)                          :: zeta
      INTEGER, INTENT(IN)                                :: lb
      REAL(KIND=dp), INTENT(IN)                          :: zetb
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: sab

      REAL(KIND=dp), PARAMETER                           :: huge4 = HUGE(1._dp)/4._dp

      INTEGER                                            :: nca, ncb, nsa, nsb
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: cab
      REAL(KIND=dp), DIMENSION(1)                        :: rpgf, za, zb
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: c2sa, c2sb

      rpgf(1) = huge4
      za(1) = zeta
      zb(1) = zetb

      nca = nco(la)
      ncb = nco(lb)
      ALLOCATE (cab(nca, ncb))
      nsa = nso(la)
      nsb = nso(lb)

      CALL overlap_ab(la, la, 1, rpgf, za, lb, lb, 1, rpgf, zb, rab, cab)

      c2sa => orbtramat(la)%c2s
      c2sb => orbtramat(lb)%c2s
      sab(1:nsa, 1:nsb) = MATMUL(c2sa(1:nsa, 1:nca), &
                                 MATMUL(cab(1:nca, 1:ncb), TRANSPOSE(c2sb(1:nsb, 1:ncb))))

      DEALLOCATE (cab)

   END SUBROUTINE overlap_ab_s

! **************************************************************************************************
!> \brief   Calculation of the overlap integrals [a|b] over
!>          cubic periodic Spherical Gaussian-type functions.
!> \param la         Max L on center A
!> \param zeta       Exponents on center A
!> \param lb         Max L on center B
!> \param zetb       Exponents on center B
!> \param alat       Lattice constant
!> \param sab        Final overlap integrals
!> \date    01.03.2016
!> \author  JGH
! **************************************************************************************************
   SUBROUTINE overlap_ab_sp(la, zeta, lb, zetb, alat, sab)
      INTEGER, INTENT(IN)                                :: la
      REAL(KIND=dp), INTENT(IN)                          :: zeta
      INTEGER, INTENT(IN)                                :: lb
      REAL(KIND=dp), INTENT(IN)                          :: zetb, alat
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: sab

      COMPLEX(KIND=dp)                                   :: zfg
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: fun, gun
      INTEGER                                            :: ax, ay, az, bx, by, bz, i, ia, ib, l, &
                                                            l1, l2, na, nb, nca, ncb, nmax, nsa, &
                                                            nsb
      REAL(KIND=dp)                                      :: oa, ob, ovol, zm
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: fexp, gexp, gval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: cab
      REAL(KIND=dp), DIMENSION(0:3, 0:3)                 :: fgsum
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: c2sa, c2sb

      nca = nco(la)
      ncb = nco(lb)
      ALLOCATE (cab(nca, ncb))
      cab = 0.0_dp
      nsa = nso(la)
      nsb = nso(lb)

      zm = MIN(zeta, zetb)
      nmax = NINT(1.81_dp*alat*SQRT(zm) + 1.0_dp)
      ALLOCATE (fun(-nmax:nmax, 0:la), gun(-nmax:nmax, 0:lb), &
                fexp(-nmax:nmax), gexp(-nmax:nmax), gval(-nmax:nmax))

      oa = 1._dp/zeta
      ob = 1._dp/zetb
      DO i = -nmax, nmax
         gval(i) = twopi/alat*REAL(i, KIND=dp)
         fexp(i) = SQRT(oa*pi)*EXP(-0.25_dp*oa*gval(i)**2)
         gexp(i) = SQRT(ob*pi)*EXP(-0.25_dp*ob*gval(i)**2)
      END DO
      DO l = 0, la
         IF (l == 0) THEN
            fun(:, l) = z_one
         ELSE IF (l == 1) THEN
            fun(:, l) = CMPLX(0.0_dp, 0.5_dp*oa*gval(:), KIND=dp)
         ELSE IF (l == 2) THEN
            fun(:, l) = CMPLX(-(0.5_dp*oa*gval(:))**2, 0.0_dp, KIND=dp)
            fun(:, l) = fun(:, l) + CMPLX(0.5_dp*oa, 0.0_dp, KIND=dp)
         ELSE IF (l == 3) THEN
            fun(:, l) = CMPLX(0.0_dp, -(0.5_dp*oa*gval(:))**3, KIND=dp)
            fun(:, l) = fun(:, l) + CMPLX(0.0_dp, 0.75_dp*oa*oa*gval(:), KIND=dp)
         ELSE
            CPABORT("l value too high")
         END IF
      END DO
      DO l = 0, lb
         IF (l == 0) THEN
            gun(:, l) = z_one
         ELSE IF (l == 1) THEN
            gun(:, l) = CMPLX(0.0_dp, 0.5_dp*ob*gval(:), KIND=dp)
         ELSE IF (l == 2) THEN
            gun(:, l) = CMPLX(-(0.5_dp*ob*gval(:))**2, 0.0_dp, KIND=dp)
            gun(:, l) = gun(:, l) + CMPLX(0.5_dp*ob, 0.0_dp, KIND=dp)
         ELSE IF (l == 3) THEN
            gun(:, l) = CMPLX(0.0_dp, -(0.5_dp*ob*gval(:))**3, KIND=dp)
            gun(:, l) = gun(:, l) + CMPLX(0.0_dp, 0.75_dp*ob*ob*gval(:), KIND=dp)
         ELSE
            CPABORT("l value too high")
         END IF
      END DO

      fgsum = 0.0_dp
      DO l1 = 0, la
         DO l2 = 0, lb
            zfg = SUM(CONJG(fun(:, l1))*fexp(:)*gun(:, l2)*gexp(:))
            fgsum(l1, l2) = REAL(zfg, KIND=dp)
         END DO
      END DO

      na = ncoset(la - 1)
      nb = ncoset(lb - 1)
      DO ax = 0, la
         DO ay = 0, la - ax
            az = la - ax - ay
            ia = coset(ax, ay, az) - na
            DO bx = 0, lb
               DO by = 0, lb - bx
                  bz = lb - bx - by
                  ib = coset(bx, by, bz) - nb
                  cab(ia, ib) = fgsum(ax, bx)*fgsum(ay, by)*fgsum(az, bz)
               END DO
            END DO
         END DO
      END DO

      c2sa => orbtramat(la)%c2s
      c2sb => orbtramat(lb)%c2s
      sab(1:nsa, 1:nsb) = MATMUL(c2sa(1:nsa, 1:nca), &
                                 MATMUL(cab(1:nca, 1:ncb), TRANSPOSE(c2sb(1:nsb, 1:ncb))))
      ovol = 1._dp/(alat**3)
      sab(1:nsa, 1:nsb) = ovol*sab(1:nsa, 1:nsb)

      DEALLOCATE (cab, fun, gun, fexp, gexp, gval)

   END SUBROUTINE overlap_ab_sp

END MODULE ai_overlap
