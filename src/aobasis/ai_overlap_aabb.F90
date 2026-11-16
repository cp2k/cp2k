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
!> \param saabb primitive integrals, including zero screened blocks and lower-shell padding
!> \date    06.2014
!> \author  Dorothea Golze
! **************************************************************************************************
   SUBROUTINE overlap_aabb(la_max_set1, la_min_set1, npgfa1, rpgfa1, zeta1, &
                           la_max_set2, la_min_set2, npgfa2, rpgfa2, zeta2, &
                           lb_max_set1, lb_min_set1, npgfb1, rpgfb1, zetb1, &
                           lb_max_set2, lb_min_set2, npgfb2, rpgfb2, zetb2, &
                           asets_equal, bsets_equal, rab, dab, saabb)

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
      REAL(KIND=dp), DIMENSION(:, :, :, :), INTENT(OUT)  :: saabb

      CHARACTER(len=*), PARAMETER                        :: routineN = 'overlap_aabb'

      INTEGER :: coa, cob, handle, ia, ib, ipgf, ja, jb, jpgf, jpgf_start, kpgf, la_max, la_min, &
         lb_max, lb_min, ldrr, lpgf, lpgf_start, ncoa1, ncoa2, ncob1, ncob2
      INTEGER, ALLOCATABLE                               :: map_a(:, :), map_b(:, :)
      INTEGER, DIMENSION(3)                              :: na, nb
      REAL(KIND=dp)                                      :: f0, zeta, zetb, zetp
      REAL(KIND=dp), ALLOCATABLE                         :: s(:, :)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: rr
      REAL(KIND=dp), DIMENSION(3)                        :: rap, rbp

      CALL timeset(routineN, handle)
      saabb = 0.0_dp
      la_max = la_max_set1 + la_max_set2
      la_min = la_min_set1 + la_min_set2
      lb_max = lb_max_set1 + lb_max_set2
      lb_min = lb_min_set1 + lb_min_set2
      ldrr = MAX(la_max, lb_max) + 1
      ALLOCATE (rr(0:ldrr - 1, 0:ldrr - 1, 3))
      ALLOCATE (s(ncoset(la_min - 1) + 1:ncoset(la_max), ncoset(lb_min - 1) + 1:ncoset(lb_max)))
      ALLOCATE (map_a(ncoset(la_min_set1 - 1) + 1:ncoset(la_max_set1), &
                      ncoset(la_min_set2 - 1) + 1:ncoset(la_max_set2)), &
                map_b(ncoset(lb_min_set1 - 1) + 1:ncoset(lb_max_set1), &
                      ncoset(lb_min_set2 - 1) + 1:ncoset(lb_max_set2)))

      ! Same-center products add Cartesian degrees; the maps are shared by all primitives.
      DO ja = LBOUND(map_a, 2), UBOUND(map_a, 2)
         DO ia = LBOUND(map_a, 1), UBOUND(map_a, 1)
            na = indco(:, ia) + indco(:, ja)
            map_a(ia, ja) = coset(na(1), na(2), na(3))
         END DO
      END DO
      DO jb = LBOUND(map_b, 2), UBOUND(map_b, 2)
         DO ib = LBOUND(map_b, 1), UBOUND(map_b, 1)
            nb = indco(:, ib) + indco(:, jb)
            map_b(ib, jb) = coset(nb(1), nb(2), nb(3))
         END DO
      END DO

!   *** Loop over all pairs of primitive Gaussian-type functions ***

      DO ipgf = 1, npgfa1
         ncoa1 = (ipgf - 1)*ncoset(la_max_set1)
         jpgf_start = MERGE(ipgf, 1, asets_equal)

         DO jpgf = jpgf_start, npgfa2

            ncoa2 = (jpgf - 1)*ncoset(la_max_set2)
            zeta = zeta1(ipgf) + zeta2(jpgf)

            DO kpgf = 1, npgfb1

               ncob1 = (kpgf - 1)*ncoset(lb_max_set1)
               lpgf_start = MERGE(kpgf, 1, bsets_equal)

               DO lpgf = lpgf_start, npgfb2
                  ncob2 = (lpgf - 1)*ncoset(lb_max_set2)

                  ! *** Screening ***
                  IF (MIN(rpgfa1(ipgf), rpgfa2(jpgf)) + &
                      MIN(rpgfb1(kpgf), rpgfb2(lpgf)) < dab) CYCLE

                  zetb = zetb1(kpgf) + zetb2(lpgf)

!           *** Calculate some prefactors ***

                  zetp = 1.0_dp/(zeta + zetb)

                  f0 = SQRT((pi*zetp)**3)*EXP(-zeta*zetb*zetp*dab*dab)
                  rap(:) = zetb*zetp*rab(:)
                  rbp(:) = -zeta*zetp*rab(:)

                  CALL os_rr_ovlp(rap, la_max, rbp, lb_max, 1.0_dp/zetp, ldrr, rr)

                  DO cob = ncoset(lb_min - 1) + 1, ncoset(lb_max)
                     nb = indco(:, cob)
                     DO coa = ncoset(la_min - 1) + 1, ncoset(la_max)
                        na = indco(:, coa)
                        s(coa, cob) = f0*rr(na(1), nb(1), 1)*rr(na(2), nb(2), 2)*rr(na(3), nb(3), 3)
                     END DO
                  END DO

!           *** Store the primitive overlap integrals ***
                  DO jb = ncoset(lb_min_set2 - 1) + 1, ncoset(lb_max_set2)
                     DO ib = ncoset(lb_min_set1 - 1) + 1, ncoset(lb_max_set1)
                        cob = map_b(ib, jb)
                        DO ja = ncoset(la_min_set2 - 1) + 1, ncoset(la_max_set2)
                           DO ia = ncoset(la_min_set1 - 1) + 1, ncoset(la_max_set1)
                              coa = map_a(ia, ja)
                              saabb(ncoa1 + ia, ncoa2 + ja, ncob1 + ib, ncob2 + jb) = s(coa, cob)
                              IF (asets_equal) saabb(ncoa2 + ja, ncoa1 + ia, ncob1 + ib, ncob2 + jb) = s(coa, cob)
                              IF (bsets_equal) saabb(ncoa1 + ia, ncoa2 + ja, ncob2 + jb, ncob1 + ib) = s(coa, cob)
                              IF (asets_equal .AND. bsets_equal) THEN
                                 saabb(ncoa2 + ja, ncoa1 + ia, ncob2 + jb, ncob1 + ib) = s(coa, cob)
                              END IF
                           END DO
                        END DO
                     END DO
                  END DO

               END DO
            END DO
         END DO
      END DO

      DEALLOCATE (rr, s, map_a, map_b)
      CALL timestop(handle)

   END SUBROUTINE overlap_aabb

END MODULE ai_overlap_aabb
