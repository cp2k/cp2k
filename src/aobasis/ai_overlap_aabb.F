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
!> \par Parameters
!>      - ax,ay,az  : Angular momentum index numbers of orbital a.
!>      - bx,by,bz  : Angular momentum index numbers of orbital b.
!>      - coset     : Cartesian orbital set pointer.
!>      - dab       : Distance between the atomic centers a and b.
!>      - l{a,b}    : Angular momentum quantum number of shell a or b.
!>      - l{a,b}_max: Maximum angular momentum quantum number of shell a or b.
!>      - l{a,b}_min: Minimum angular momentum quantum number of shell a or b.
!>      - rab       : Distance vector between the atomic centers a and b.
!>      - rpgf{a,b} : Radius of the primitive Gaussian-type function a or b.
!>      - sab       : Shell set of overlap integrals.
!>      - zet{a,b}  : Exponents of the Gaussian-type functions a or b.
!>      - zetp      : Reciprocal of the sum of the exponents of orbital a and b.
! **************************************************************************************************
MODULE ai_overlap_aabb

   USE ai_os_rr,                        ONLY: os_rr_ovlp
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
   USE orbital_pointers,                ONLY: coset,&
                                              indco,&
                                              ncoset
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_overlap_aabb'

! *** Public subroutines ***
   PUBLIC :: overlap_aabb

CONTAINS

! **************************************************************************************************
!> \brief   Purpose: Calculation of the two-center overlap integrals [aa|bb]
!>          over Cartesian Gaussian-type functions.
!> \param la_max_set1 ...
!> \param la_min_set1 ...
!> \param npgfa1 ...
!> \param rpgfa1 ...
!> \param zeta1 ...
!> \param la_max_set2 ...
!> \param la_min_set2 ...
!> \param npgfa2 ...
!> \param rpgfa2 ...
!> \param zeta2 ...
!> \param lb_max_set1 ...
!> \param lb_min_set1 ...
!> \param npgfb1 ...
!> \param rpgfb1 ...
!> \param zetb1 ...
!> \param lb_max_set2 ...
!> \param lb_min_set2 ...
!> \param npgfb2 ...
!> \param rpgfb2 ...
!> \param zetb2 ...
!> \param asets_equal ...
!> \param bsets_equal ...
!> \param rab ...
!> \param dab ...
!> \param saabb ...
!> \param s ...
!> \param lds ...
!> \date    06.2014
!> \author  Dorothea Golze
! **************************************************************************************************
   SUBROUTINE overlap_aabb(la_max_set1, la_min_set1, npgfa1, rpgfa1, zeta1, &
                           la_max_set2, la_min_set2, npgfa2, rpgfa2, zeta2, &
                           lb_max_set1, lb_min_set1, npgfb1, rpgfb1, zetb1, &
                           lb_max_set2, lb_min_set2, npgfb2, rpgfb2, zetb2, &
                           asets_equal, bsets_equal, rab, dab, saabb, s, lds)

      INTEGER, INTENT(IN)                                :: la_max_set1, la_min_set1, npgfa1
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rpgfa1, zeta1
      INTEGER, INTENT(IN)                                :: la_max_set2, la_min_set2, npgfa2
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rpgfa2, zeta2
      INTEGER, INTENT(IN)                                :: lb_max_set1, lb_min_set1, npgfb1
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rpgfb1, zetb1
      INTEGER, INTENT(IN)                                :: lb_max_set2, lb_min_set2, npgfb2
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rpgfb2, zetb2
      LOGICAL, INTENT(IN)                                :: asets_equal, bsets_equal
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), INTENT(IN)                          :: dab
      REAL(KIND=dp), DIMENSION(:, :, :, :), &
         INTENT(INOUT)                                   :: saabb
      INTEGER, INTENT(IN)                                :: lds
      REAL(KIND=dp), DIMENSION(lds, lds), INTENT(INOUT)  :: s

      CHARACTER(len=*), PARAMETER                        :: routineN = 'overlap_aabb'

      INTEGER :: ax, ay, az, bx, by, bz, coa, cob, handle, i, ia, ib, ipgf, j, ja, jb, jpgf, &
         jpgf_start, kpgf, la, la_max, la_min, lb, lb_max, lb_min, ldrr, lpgf, lpgf_start, ncoa1, &
         ncoa2, ncob1, ncob2
      INTEGER, DIMENSION(3)                              :: na, naa, nb, nbb, nia, nib, nja, njb
      REAL(KIND=dp)                                      :: f0, zeta, zetb, zetp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: rr
      REAL(KIND=dp), DIMENSION(3)                        :: rap, rbp

      CALL timeset(routineN, handle)
      ldrr = MAX(la_max_set1 + la_max_set2, lb_max_set1 + lb_max_set2) + 1
      ALLOCATE (rr(0:ldrr - 1, 0:ldrr - 1, 3))

!   *** Loop over all pairs of primitive Gaussian-type functions ***

      ncoa1 = 0
      ncoa2 = 0
      ncob1 = 0
      ncob2 = 0

      DO ipgf = 1, npgfa1

         ncoa2 = 0

         IF (asets_equal) THEN
            jpgf_start = ipgf
            DO i = 1, jpgf_start - 1
               ncoa2 = ncoa2 + ncoset(la_max_set2)
            END DO
         ELSE
            jpgf_start = 1
         END IF

         DO jpgf = jpgf_start, npgfa2

            ncob1 = 0
            zeta = zeta1(ipgf) + zeta2(jpgf)
            la_max = la_max_set1 + la_max_set2
            la_min = la_min_set1 + la_min_set2

            DO kpgf = 1, npgfb1

               ncob2 = 0

               IF (bsets_equal) THEN
                  lpgf_start = kpgf
                  DO i = 1, lpgf_start - 1
                     ncob2 = ncob2 + ncoset(lb_max_set2)
                  END DO
               ELSE
                  lpgf_start = 1
               END IF

               DO lpgf = lpgf_start, npgfb2

                  ! *** Screening ***
                  IF ((rpgfa1(ipgf) + rpgfb1(kpgf) < dab) .OR. &
                      (rpgfa2(jpgf) + rpgfb1(kpgf) < dab) .OR. &
                      (rpgfa1(ipgf) + rpgfb2(lpgf) < dab) .OR. &
                      (rpgfa2(jpgf) + rpgfb2(lpgf) < dab)) THEN
                     DO jb = ncoset(lb_min_set2 - 1) + 1, ncoset(lb_max_set2)
                        DO ib = ncoset(lb_min_set1 - 1) + 1, ncoset(lb_max_set1)
                           DO ja = ncoset(la_min_set2 - 1) + 1, ncoset(la_max_set2)
                              DO ia = ncoset(la_min_set1 - 1) + 1, ncoset(la_max_set1)
                                 saabb(ncoa1 + ia, ncoa2 + ja, ncob1 + ib, ncob2 + jb) = 0._dp
                                 IF (asets_equal) saabb(ncoa2 + ja, ncoa1 + ia, ncob1 + ib, ncob2 + jb) = 0._dp
                                 IF (bsets_equal) saabb(ncoa1 + ia, ncoa2 + ja, ncob2 + jb, ncob1 + ib) = 0._dp
                                 IF (asets_equal .AND. bsets_equal) THEN
                                    saabb(ncoa2 + ja, ncoa1 + ia, ncob2 + jb, ncob1 + ib) = 0._dp
                                 END IF
                              END DO
                           END DO
                        END DO
                     END DO
                     ncob2 = ncob2 + ncoset(lb_max_set2)
                     CYCLE
                  END IF

                  zetb = zetb1(kpgf) + zetb2(lpgf)
                  lb_max = lb_max_set1 + lb_max_set2
                  lb_min = lb_min_set1 + lb_min_set2

!           *** Calculate some prefactors ***

                  zetp = 1.0_dp/(zeta + zetb)

                  f0 = SQRT((pi*zetp)**3)*EXP(-zeta*zetb*zetp*dab*dab)
                  rap(:) = zetb*zetp*rab(:)
                  rbp(:) = -zeta*zetp*rab(:)

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
                                    s(coa, cob) = f0*rr(ax, bx, 1)*rr(ay, by, 2)*rr(az, bz, 3)
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO

!           *** Store the primitive overlap integrals ***
                  DO jb = ncoset(lb_min_set2 - 1) + 1, ncoset(lb_max_set2)
                     njb(1:3) = indco(1:3, jb)
                     DO ib = ncoset(lb_min_set1 - 1) + 1, ncoset(lb_max_set1)
                        nib(1:3) = indco(1:3, ib)
                        nbb(1:3) = nib + njb
                        DO ja = ncoset(la_min_set2 - 1) + 1, ncoset(la_max_set2)
                           nja(1:3) = indco(1:3, ja)
                           DO ia = ncoset(la_min_set1 - 1) + 1, ncoset(la_max_set1)
                              nia(1:3) = indco(1:3, ia)
                              naa(1:3) = nia + nja
                              ! now loop over all elements of s
                              DO j = ncoset(lb_min - 1) + 1, ncoset(lb_max)
                                 nb(1:3) = indco(1:3, j)
                                 DO i = ncoset(la_min - 1) + 1, ncoset(la_max)
                                    na(1:3) = indco(1:3, i)
                                    IF (ALL(na == naa) .AND. ALL(nb == nbb)) THEN
                                       saabb(ncoa1 + ia, ncoa2 + ja, ncob1 + ib, ncob2 + jb) = s(i, j)
                                       IF (asets_equal) saabb(ncoa2 + ja, ncoa1 + ia, ncob1 + ib, ncob2 + jb) = s(i, j)
                                       IF (bsets_equal) saabb(ncoa1 + ia, ncoa2 + ja, ncob2 + jb, ncob1 + ib) = s(i, j)
                                       IF (asets_equal .AND. bsets_equal) THEN
                                          saabb(ncoa2 + ja, ncoa1 + ia, ncob2 + jb, ncob1 + ib) = s(i, j)
                                       END IF
                                    END IF
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO

                  ncob2 = ncob2 + ncoset(lb_max_set2)

               END DO

               ncob1 = ncob1 + ncoset(lb_max_set1)

            END DO

            ncoa2 = ncoa2 + ncoset(la_max_set2)

         END DO

         ncoa1 = ncoa1 + ncoset(la_max_set1)

      END DO

      DEALLOCATE (rr)
      CALL timestop(handle)

   END SUBROUTINE overlap_aabb

END MODULE ai_overlap_aabb
