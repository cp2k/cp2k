!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Expansion of the product of two real spherical harmonics in real
!>        spherical harmonics
!> \note
!>      The expansion coefficients are the Gaunt coefficients of the real
!>      spherical harmonics; they are obtained from libwignernj, which evaluates
!>      them in exact integer arithmetic. Historically CP2K called these
!>      "Clebsch-Gordon coefficients", which they are not, and built them from a
!>      table of factorial ratios that lost about one digit per two units of l.
!> \par History
!>      09.2026 created out of the five copies of this loop that lived in
!>              qs_rho_atom_methods, qs_tddfpt2_soc, qs_cneo_utils,
!>              xas_tdp_atom and generic_shg_integrals_init [Susi Lehtola]
!> \author Susi Lehtola
! **************************************************************************************************
MODULE gaunt_coefficients

   USE kinds,                           ONLY: dp
   USE memory_utilities,                ONLY: reallocate
   USE orbital_pointers,                ONLY: indso,&
                                              nsoset
   USE wignernj_interface,              ONLY: real_gaunt
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: real_gaunt_coefficients

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gaunt_coefficients'

   ! Coefficients below this threshold are not entered in the non-zero list
   REAL(KIND=dp), PARAMETER, PRIVATE :: cg_eps = 1.0E-8_dp

CONTAINS

! **************************************************************************************************
!> \brief Calculate the expansion of the product of two real spherical harmonics
!>        in real spherical harmonics,
!>           S(l1,m1) S(l2,m2) = sum_lm my_cg(iso1,iso2,iso) S(l,m)
!> \param my_cg the expansion coefficients, allocated here; the pointer is
!>        INTENT(OUT), so the caller must not rely on its previous contents
!> \param maxl1 maximal angular momentum of the first harmonic
!> \param maxl2 maximal angular momentum of the second harmonic
!> \param llmax maximal angular momentum retained in the expansion
!> \param cg_none0_list list of the harmonics with a non-zero coefficient,
!>        allocated here if present
!> \param ncg_none0 number of harmonics with a non-zero coefficient,
!>        allocated here if present
! **************************************************************************************************
   SUBROUTINE real_gaunt_coefficients(my_cg, maxl1, maxl2, llmax, cg_none0_list, ncg_none0)

      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT), &
         POINTER                                         :: my_cg
      INTEGER, INTENT(IN)                                :: maxl1, maxl2, llmax
      INTEGER, DIMENSION(:, :, :), INTENT(OUT), &
         OPTIONAL, POINTER                               :: cg_none0_list
      INTEGER, DIMENSION(:, :), INTENT(OUT), OPTIONAL, &
         POINTER                                         :: ncg_none0

      CHARACTER(len=*), PARAMETER :: routineN = 'real_gaunt_coefficients'

      INTEGER                                            :: handle, ilist, im, iso, iso1, iso2, l1, &
                                                            l2, lp, m1, m2, mp(4), nlist, &
                                                            nlist_max, nsfunc, nsfunc1, nsfunc2
      LOGICAL                                            :: do_list

      CALL timeset(routineN, handle)

      CPASSERT(maxl1 >= 0)
      CPASSERT(maxl2 >= 0)
      CPASSERT(llmax >= 0)

      nsfunc1 = nsoset(maxl1)
      nsfunc2 = nsoset(maxl2)
      nsfunc = nsoset(llmax)

      ALLOCATE (my_cg(nsfunc1, nsfunc2, nsfunc))
      my_cg = 0.0_dp

      DO iso1 = 1, nsfunc1
         l1 = indso(1, iso1)
         m1 = indso(2, iso1)
         DO iso2 = 1, nsfunc2
            l2 = indso(1, iso2)
            m2 = indso(2, iso2)
            ! The product of two real harmonics only contains the azimuthal
            ! frequencies |m1+m2| and |m1-m2|; the cosine part goes to m >= 0
            ! and the sine part to m < 0, so of the four candidates below at
            ! most two have a non-vanishing coefficient. Rather than working
            ! out which one it is, we evaluate all four and let the selection
            ! rules zero out the rest.
            mp(1) = ABS(m1 + m2)
            mp(2) = -mp(1)
            mp(3) = ABS(m1 - m2)
            mp(4) = -mp(3)
            DO lp = MOD(l1 + l2, 2), MIN(l1 + l2, llmax), 2
               DO im = 1, 4
                  IF (ABS(mp(im)) > lp) CYCLE
                  iso = nsoset(lp - 1) + lp + 1 + mp(im)
                  my_cg(iso1, iso2, iso) = real_gaunt(l1, m1, l2, m2, lp, mp(im))
               END DO
            END DO
         END DO
      END DO

      do_list = PRESENT(cg_none0_list) .OR. PRESENT(ncg_none0)
      IF (do_list) THEN
         CPASSERT(PRESENT(cg_none0_list))
         CPASSERT(PRESENT(ncg_none0))
         nlist_max = 6
         ALLOCATE (ncg_none0(nsfunc1, nsfunc2))
         ALLOCATE (cg_none0_list(nsfunc1, nsfunc2, nlist_max))
         ncg_none0 = 0
         cg_none0_list = 0
         DO iso1 = 1, nsfunc1
            DO iso2 = 1, nsfunc2
               nlist = 0
               DO ilist = 1, nsfunc
                  IF (ABS(my_cg(iso1, iso2, ilist)) > cg_eps) THEN
                     nlist = nlist + 1
                     IF (nlist > nlist_max) THEN
                        CALL reallocate(cg_none0_list, 1, nsfunc1, 1, nsfunc2, 1, nlist)
                        nlist_max = nlist
                     END IF
                     cg_none0_list(iso1, iso2, nlist) = ilist
                  END IF
               END DO
               ncg_none0(iso1, iso2) = nlist
            END DO
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE real_gaunt_coefficients

END MODULE gaunt_coefficients
