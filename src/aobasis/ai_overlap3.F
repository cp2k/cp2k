!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Three-center overlap integrals over Cartesian Gaussian functions.
!> \author Matthias Krack (26.06.2001)
!> \par Literature
!>      S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
! **************************************************************************************************
MODULE ai_overlap3

   USE ai_cartesian_operators,          ONLY: cartesian_ladder,&
                                              cartesian_ladder_map
   USE ai_os_rr,                        ONLY: os_rr_ovlp3
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
   USE orbital_pointers,                ONLY: indco,&
                                              ncoset
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_overlap3'

   PUBLIC :: overlap3

CONTAINS

! **************************************************************************************************
!> \brief Calculation of three-center overlap integrals [a|b|c] over primitive
!>        Cartesian Gaussian functions using separable one-dimensional tables.
!> \param la_max_set ...
!> \param npgfa ...
!> \param zeta ...
!> \param rpgfa ...
!> \param la_min_set ...
!> \param lb_max_set ...
!> \param npgfb ...
!> \param zetb ...
!> \param rpgfb ...
!> \param lb_min_set ...
!> \param lc_max_set ...
!> \param npgfc ...
!> \param zetc ...
!> \param rpgfc ...
!> \param lc_min_set ...
!> \param rab displacement B - A
!> \param dab distance between A and B
!> \param rac displacement C - A
!> \param dac distance between A and C
!> \param rbc displacement C - B
!> \param dbc distance between B and C
!> \param sabc integrals [a|b|c]
!> \param sdabc derivative [da/dAi|b|c]
!> \param sabdc derivative [a|b|dc/dCi]
!> \param int_abc_ext largest absolute integral in the requested shell ranges
!> \par History
!>      05.2014 created (Dorothea Golze)
!> \author Dorothea Golze
! **************************************************************************************************
   SUBROUTINE overlap3(la_max_set, npgfa, zeta, rpgfa, la_min_set, &
                       lb_max_set, npgfb, zetb, rpgfb, lb_min_set, &
                       lc_max_set, npgfc, zetc, rpgfc, lc_min_set, &
                       rab, dab, rac, dac, rbc, dbc, sabc, &
                       sdabc, sabdc, int_abc_ext)

      INTEGER, INTENT(IN)                                :: la_max_set, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta, rpgfa
      INTEGER, INTENT(IN)                                :: la_min_set, lb_max_set, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb, rpgfb
      INTEGER, INTENT(IN)                                :: lb_min_set, lc_max_set, npgfc
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetc, rpgfc
      INTEGER, INTENT(IN)                                :: lc_min_set
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), INTENT(IN)                          :: dab
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rac
      REAL(KIND=dp), INTENT(IN)                          :: dac
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rbc
      REAL(KIND=dp), INTENT(IN)                          :: dbc
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: sabc
      REAL(KIND=dp), DIMENSION(:, :, :, :), &
         INTENT(INOUT), OPTIONAL                         :: sdabc, sabdc
      REAL(dp), INTENT(OUT), OPTIONAL                    :: int_abc_ext

      CHARACTER(len=*), PARAMETER                        :: routineN = 'overlap3'

      INTEGER :: a_first, a_last, avec(3), axis, b_first, b_last, bvec(3), c_first, c_last, coa, &
         cob, coc, cvec(3), handle, ia, ib, ic, ipgf, jaxis, jpgf, kaxis, kpgf, la_max, lc_max, &
         na, nb, nc, map_c(4, lc_max_set + 1), map_a(4, la_max_set + 1)
      REAL(KIND=dp)                                      :: f0, f1, factor, rcp2, value, zet, zetg, &
                                                            zetp
      REAL(KIND=dp), ALLOCATABLE                         :: derivative_a(:, :, :, :), &
                                                            derivative_c(:, :, :, :), &
                                                            rr(:, :, :, :)
      REAL(KIND=dp), DIMENSION(3)                        :: rap, rbp, rcg, rcp

      CALL timeset(routineN, handle)

      la_max = la_max_set
      lc_max = lc_max_set
      IF (PRESENT(sdabc)) la_max = la_max + 1
      IF (PRESENT(sabdc)) lc_max = lc_max + 1
      ALLOCATE (rr(0:la_max, 0:lb_max_set, 0:lc_max, 3))

      a_first = ncoset(la_min_set - 1) + 1
      a_last = ncoset(la_max_set)
      b_first = ncoset(lb_min_set - 1) + 1
      b_last = ncoset(lb_max_set)
      c_first = ncoset(lc_min_set - 1) + 1
      c_last = ncoset(lc_max_set)

      ! Clear the requested derivative extents, including screened blocks and lower-shell padding.
      IF (PRESENT(sdabc)) THEN
         sdabc(1:npgfa*a_last, 1:npgfb*b_last, 1:npgfc*c_last, 1:3) = 0.0_dp
         ALLOCATE (derivative_a(0:la_max_set, 0:lb_max_set, 0:lc_max_set, 3))
         CALL cartesian_ladder_map(la_max_set, 0, 0, map_a)
      END IF
      IF (PRESENT(sabdc)) THEN
         sabdc(1:npgfa*a_last, 1:npgfb*b_last, 1:npgfc*c_last, 1:3) = 0.0_dp
         ALLOCATE (derivative_c(0:la_max_set, 0:lb_max_set, 0:lc_max_set, 3))
         CALL cartesian_ladder_map(lc_max_set, 0, 0, map_c)
      END IF
      IF (PRESENT(int_abc_ext)) int_abc_ext = 0.0_dp

      DO ipgf = 1, npgfa
         na = (ipgf - 1)*a_last
         DO jpgf = 1, npgfb
            nb = (jpgf - 1)*b_last
            IF (rpgfa(ipgf) + rpgfb(jpgf) < dab) THEN
               sabc(na + 1:na + a_last, nb + 1:nb + b_last, 1:npgfc*c_last) = 0.0_dp
               CYCLE
            END IF
            zetp = 1.0_dp/(zeta(ipgf) + zetb(jpgf))
            f1 = zetb(jpgf)*zetp
            rcp = f1*rab - rac
            rcp2 = DOT_PRODUCT(rcp, rcp)
            DO kpgf = 1, npgfc
               nc = (kpgf - 1)*c_last
               IF ((rpgfb(jpgf) + rpgfc(kpgf) < dbc) .OR. &
                   (rpgfa(ipgf) + rpgfc(kpgf) < dac)) THEN
                  sabc(na + 1:na + a_last, nb + 1:nb + b_last, nc + 1:nc + c_last) = 0.0_dp
                  CYCLE
               END IF

               zet = zeta(ipgf) + zetb(jpgf) + zetc(kpgf)
               zetg = 1.0_dp/zet
               f0 = (pi*zetg)**1.5_dp
               f0 = f0*EXP(-(zeta(ipgf)*f1*dab*dab + zetc(kpgf)*zetg*rcp2/zetp))

               ! G is the product center of all three Gaussians. The 1D tables are normalized;
               ! the three-dimensional [s|s|s] prefactor f0 is applied only on reconstruction.
               rap = zetg*(zetb(jpgf)*rab + zetc(kpgf)*rac)
               rbp = zetg*(-zeta(ipgf)*rab + zetc(kpgf)*rbc)
               rcg = -zetg*(zeta(ipgf)*rac + zetb(jpgf)*rbc)
               CALL os_rr_ovlp3(rap, la_max, rbp, lb_max_set, rcg, lc_max, zet, rr)

               ! Differentiate the original A and C basis functions after three-center table construction.
               IF (PRESENT(sdabc)) THEN
                  DO axis = 1, 3
                     DO ic = 0, lc_max_set
                        CALL cartesian_ladder(rr(:, :, ic, axis), derivative_a(:, :, ic, axis), &
                                              map_a, .TRUE., zeta(ipgf), 0.0_dp, .FALSE.)
                     END DO
                  END DO
               END IF
               IF (PRESENT(sabdc)) THEN
                  DO axis = 1, 3
                     DO ib = 0, lb_max_set
                        CALL cartesian_ladder(rr(0:la_max_set, ib, :, axis), derivative_c(:, ib, :, axis), &
                                              map_c, .TRUE., zetc(kpgf), 0.0_dp, .TRUE.)
                     END DO
                  END DO
               END IF

               DO coc = c_first, c_last
                  cvec = indco(:, coc)
                  ic = nc + coc
                  DO cob = b_first, b_last
                     bvec = indco(:, cob)
                     ib = nb + cob
                     DO coa = a_first, a_last
                        avec = indco(:, coa)
                        ia = na + coa
                        value = f0*rr(avec(1), bvec(1), cvec(1), 1)* &
                                rr(avec(2), bvec(2), cvec(2), 2)*rr(avec(3), bvec(3), cvec(3), 3)
                        sabc(ia, ib, ic) = value
                        IF (PRESENT(int_abc_ext)) int_abc_ext = MAX(int_abc_ext, ABS(value))
                        IF (PRESENT(sdabc) .OR. PRESENT(sabdc)) THEN
                           DO axis = 1, 3
                              jaxis = MOD(axis, 3) + 1
                              kaxis = MOD(axis + 1, 3) + 1
                              factor = f0*rr(avec(jaxis), bvec(jaxis), cvec(jaxis), jaxis)* &
                                       rr(avec(kaxis), bvec(kaxis), cvec(kaxis), kaxis)
                              IF (PRESENT(sdabc)) THEN
                                 sdabc(ia, ib, ic, axis) = factor*derivative_a(avec(axis), bvec(axis), cvec(axis), axis)
                              END IF
                              IF (PRESENT(sabdc)) THEN
                                 sabdc(ia, ib, ic, axis) = factor*derivative_c(avec(axis), bvec(axis), cvec(axis), axis)
                              END IF
                           END DO
                        END IF
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE overlap3

END MODULE ai_overlap3
