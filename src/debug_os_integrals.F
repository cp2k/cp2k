!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Debugs Obara-Saika integral matrices
!> \par History
!>      created [07.2014]
!> \authors Dorothea Golze
! **************************************************************************************************
MODULE debug_os_integrals

   USE ai_overlap3_debug,               ONLY: init_os_overlap3,&
                                              os_overlap3
   USE ai_overlap_debug,                ONLY: init_os_overlap2,&
                                              os_overlap2
   USE kinds,                           ONLY: dp
   USE orbital_pointers,                ONLY: coset,&
                                              indco,&
                                              ncoset
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! **************************************************************************************************

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'debug_os_integrals'

   PUBLIC :: overlap_ab_test, overlap_abc_test, overlap_aabb_test

! **************************************************************************************************

CONTAINS

! ***************************************************************************************************
!> \brief recursive test routines for integral (a,b)
!> \param la_max ...
!> \param la_min ...
!> \param npgfa ...
!> \param zeta ...
!> \param lb_max ...
!> \param lb_min ...
!> \param npgfb ...
!> \param zetb ...
!> \param ra ...
!> \param rb ...
!> \param sab ...
!> \param dmax ...
! **************************************************************************************************
   SUBROUTINE overlap_ab_test(la_max, la_min, npgfa, zeta, lb_max, lb_min, npgfb, zetb, &
                              ra, rb, sab, dmax)

      INTEGER, INTENT(IN)                                :: la_max, la_min, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta
      INTEGER, INTENT(IN)                                :: lb_max, lb_min, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: ra, rb
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: sab
      REAL(KIND=dp), INTENT(INOUT)                       :: dmax

      INTEGER                                            :: coa, cob, ia1, iax, iay, iaz, ib1, ibx, &
                                                            iby, ibz, ipgf, jpgf, ma, mb
      INTEGER, DIMENSION(3)                              :: na, nb
      REAL(KIND=dp)                                      :: res1, res2, xa, xb
      REAL(KIND=dp), DIMENSION(3)                        :: A, B

      coa = 0
      DO ipgf = 1, npgfa
         cob = 0
         DO jpgf = 1, npgfb
            xa = zeta(ipgf) !exponents
            xb = zetb(jpgf)
            A = ra !positions
            B = rb
            CALL init_os_overlap2(xa, xb, A, B)
            DO ma = la_min, la_max
               DO mb = lb_min, lb_max
                  DO iax = 0, ma
                     DO iay = 0, ma - iax
                        iaz = ma - iax - iay
                        na(1) = iax; na(2) = iay; na(3) = iaz
                        ia1 = coset(iax, iay, iaz)
                        DO ibx = 0, mb
                           DO iby = 0, mb - ibx
                              ibz = mb - ibx - iby
                              nb(1) = ibx; nb(2) = iby; nb(3) = ibz
                              ib1 = coset(ibx, iby, ibz)
                              res1 = os_overlap2(na, nb)
                              res2 = sab(coa + ia1, cob + ib1)
                              dmax = MAX(dmax, ABS(res1 - res2))
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
            cob = cob + ncoset(lb_max)
         END DO
         coa = coa + ncoset(la_max)
      END DO
      !WRITE(*,*) "dmax overlap_ab_test", dmax

   END SUBROUTINE overlap_ab_test

! ***************************************************************************************************
!> \brief recursive test routines for integral (a,b,c)
!> \param la_max ...
!> \param npgfa ...
!> \param zeta ...
!> \param la_min ...
!> \param lb_max ...
!> \param npgfb ...
!> \param zetb ...
!> \param lb_min ...
!> \param lc_max ...
!> \param npgfc ...
!> \param zetc ...
!> \param lc_min ...
!> \param ra ...
!> \param rb ...
!> \param rc ...
!> \param sabc ...
!> \param dmax ...
! **************************************************************************************************
   SUBROUTINE overlap_abc_test(la_max, npgfa, zeta, la_min, &
                               lb_max, npgfb, zetb, lb_min, &
                               lc_max, npgfc, zetc, lc_min, &
                               ra, rb, rc, sabc, dmax)

      INTEGER, INTENT(IN)                                :: la_max, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta
      INTEGER, INTENT(IN)                                :: la_min, lb_max, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb
      INTEGER, INTENT(IN)                                :: lb_min, lc_max, npgfc
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetc
      INTEGER, INTENT(IN)                                :: lc_min
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: ra, rb, rc
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: sabc
      REAL(KIND=dp), INTENT(INOUT)                       :: dmax

      INTEGER                                            :: coa, cob, coc, ia1, iax, iay, iaz, ib1, &
                                                            ibx, iby, ibz, ic1, icx, icy, icz, &
                                                            ipgf, jpgf, kpgf, ma, mb, mc
      INTEGER, DIMENSION(3)                              :: na, nb, nc
      REAL(KIND=dp)                                      :: res1, res2, xa, xb, xc
      REAL(KIND=dp), DIMENSION(3)                        :: A, B, C

      coa = 0
      DO ipgf = 1, npgfa
         cob = 0
         DO jpgf = 1, npgfb
            coc = 0
            DO kpgf = 1, npgfc

               xa = zeta(ipgf) ! exponents
               xb = zetb(jpgf)
               xc = zetc(kpgf)

               A = Ra !positions
               B = Rb
               C = Rc

               CALL init_os_overlap3(xa, xb, xc, A, B, C)

               DO ma = la_min, la_max
                  DO mc = lc_min, lc_max
                     DO mb = lb_min, lb_max
                        DO iax = 0, ma
                           DO iay = 0, ma - iax
                              iaz = ma - iax - iay
                              na(1) = iax; na(2) = iay; na(3) = iaz
                              ia1 = coset(iax, iay, iaz)
                              DO icx = 0, mc
                                 DO icy = 0, mc - icx
                                    icz = mc - icx - icy
                                    nc(1) = icx; nc(2) = icy; nc(3) = icz
                                    ic1 = coset(icx, icy, icz)
                                    DO ibx = 0, mb
                                       DO iby = 0, mb - ibx
                                          ibz = mb - ibx - iby
                                          nb(1) = ibx; nb(2) = iby; nb(3) = ibz
                                          ib1 = coset(ibx, iby, ibz)
                                          res1 = os_overlap3(na, nc, nb)
                                          res2 = sabc(coa + ia1, cob + ib1, coc + ic1)
                                          dmax = MAX(dmax, ABS(res1 - res2))
                                          !IF(dmax > 1.E-10) WRITE(*,*) "dmax in loop", dmax
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
               coc = coc + ncoset(lc_max)
            END DO
            cob = cob + ncoset(lb_max)
         END DO
         coa = coa + ncoset(la_max)
      END DO
      !WRITE(*,*) "dmax abc", dmax

   END SUBROUTINE overlap_abc_test

! ***************************************************************************************************
!> \brief recursive test routines for integral (aa,bb)
!> \param la_max1 ...
!> \param la_min1 ...
!> \param npgfa1 ...
!> \param zeta1 ...
!> \param la_max2 ...
!> \param la_min2 ...
!> \param npgfa2 ...
!> \param zeta2 ...
!> \param lb_max1 ...
!> \param lb_min1 ...
!> \param npgfb1 ...
!> \param zetb1 ...
!> \param lb_max2 ...
!> \param lb_min2 ...
!> \param npgfb2 ...
!> \param zetb2 ...
!> \param ra ...
!> \param rb ...
!> \param saabb ...
!> \param dmax ...
! **************************************************************************************************
   SUBROUTINE overlap_aabb_test(la_max1, la_min1, npgfa1, zeta1, &
                                la_max2, la_min2, npgfa2, zeta2, &
                                lb_max1, lb_min1, npgfb1, zetb1, &
                                lb_max2, lb_min2, npgfb2, zetb2, &
                                ra, rb, saabb, dmax)

      INTEGER, INTENT(IN)                                :: la_max1, la_min1, npgfa1
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta1
      INTEGER, INTENT(IN)                                :: la_max2, la_min2, npgfa2
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta2
      INTEGER, INTENT(IN)                                :: lb_max1, lb_min1, npgfb1
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb1
      INTEGER, INTENT(IN)                                :: lb_max2, lb_min2, npgfb2
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb2
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: ra, rb
      REAL(KIND=dp), DIMENSION(:, :, :, :), INTENT(IN)   :: saabb
      REAL(KIND=dp), INTENT(INOUT)                       :: dmax

      INTEGER                                            :: coa1, coa2, cob1, cob2, i, iax, iay, &
                                                            iaz, ibx, iby, ibz, ipgf, j, jpgf, k, &
                                                            kpgf, l, la_max, la_min, lb_max, &
                                                            lb_min, lpgf, ma, mb
      INTEGER, DIMENSION(3)                              :: na, naa, nb, nbb
      REAL(KIND=dp)                                      :: res1, xa, xb
      REAL(KIND=dp), DIMENSION(3)                        :: A, B

      coa1 = 0
      DO ipgf = 1, npgfa1
         coa2 = 0
         DO jpgf = 1, npgfa2
            cob1 = 0
            DO kpgf = 1, npgfb1
               cob2 = 0
               DO lpgf = 1, npgfb2

                  xa = zeta1(ipgf) + zeta2(jpgf) ! exponents
                  xb = zetb1(kpgf) + zetb2(lpgf) ! exponents
                  la_max = la_max1 + la_max2
                  lb_max = lb_max1 + lb_max2
                  la_min = la_min1 + la_min2
                  lb_min = lb_min1 + lb_min2

                  A = ra !positions
                  B = rb

                  CALL init_os_overlap2(xa, xb, A, B)

                  DO ma = la_min, la_max
                     DO mb = lb_min, lb_max
                        DO iax = 0, ma
                           DO iay = 0, ma - iax
                              iaz = ma - iax - iay
                              na(1) = iax; na(2) = iay; na(3) = iaz
                              DO ibx = 0, mb
                                 DO iby = 0, mb - ibx
                                    ibz = mb - ibx - iby
                                    nb(1) = ibx; nb(2) = iby; nb(3) = ibz
                                    res1 = os_overlap2(na, nb)
                                    DO i = ncoset(la_min1 - 1) + 1, ncoset(la_max1)
                                       DO j = ncoset(la_min2 - 1) + 1, ncoset(la_max2)
                                          naa = indco(1:3, i) + indco(1:3, j)
                                          DO k = ncoset(lb_min1 - 1) + 1, ncoset(lb_max1)
                                             DO l = ncoset(lb_min2 - 1) + 1, ncoset(lb_max2)
                                                nbb = indco(1:3, k) + indco(1:3, l)
                                                IF (ALL(na == naa) .AND. ALL(nb == nbb)) THEN
                                                   dmax = MAX(dmax, ABS(res1 - saabb(coa1 + i, coa2 + j, cob1 + k, cob2 + l)))
                                                END IF
                                             END DO
                                          END DO
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
                  cob2 = cob2 + ncoset(lb_max2)
               END DO
               cob1 = cob1 + ncoset(lb_max1)
            END DO
            coa2 = coa2 + ncoset(la_max2)
         END DO
         coa1 = coa1 + ncoset(la_max1)
      END DO

      !WRITE(*,*) "dmax aabb", dmax

   END SUBROUTINE overlap_aabb_test

END MODULE debug_os_integrals
