!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Local and semi-local ECP integrals using the libgrpp library
! **************************************************************************************************

MODULE libgrpp_integrals
   USE kinds, ONLY: dp
   USE mathconstants, ONLY: pi
   USE orbital_pointers, ONLY: nco, &
                               ncoset
#if defined(__LIBGRPP)
   USE libgrpp, ONLY: libgrpp_type1_integrals, libgrpp_type2_integrals
#endif
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'libgrpp_integrals'

   PUBLIC :: libgrpp_semilocal_integral, libgrpp_local_integral

CONTAINS

! **************************************************************************************************
!> \brief Local ECP integrals using libgrpp
!> \param la_max_set ...
!> \param la_min_set ...
!> \param npgfa ...
!> \param rpgfa ...
!> \param zeta ...
!> \param lb_max_set ...
!> \param lb_min_set ...
!> \param npgfb ...
!> \param rpgfb ...
!> \param zetb ...
!> \param npot_ecp ...
!> \param alpha_ecp ...
!> \param coeffs_ecp ...
!> \param nrpot_ecp ...
!> \param rpgfc ...
!> \param rab ...
!> \param dab ...
!> \param rac ...
!> \param dac ...
!> \param dbc ...
!> \param vab ...
! **************************************************************************************************
   SUBROUTINE libgrpp_local_integral(la_max_set, la_min_set, npgfa, rpgfa, zeta, &
                                     lb_max_set, lb_min_set, npgfb, rpgfb, zetb, &
                                     npot_ecp, alpha_ecp, coeffs_ecp, nrpot_ecp, &
                                     rpgfc, rab, dab, rac, dac, dbc, vab)

      INTEGER, INTENT(IN)                                :: la_max_set, la_min_set, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rpgfa, zeta
      INTEGER, INTENT(IN)                                :: lb_max_set, lb_min_set, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rpgfb, zetb
      INTEGER, INTENT(IN)                                :: npot_ecp
      REAL(KIND=dp), DIMENSION(1:npot_ecp), INTENT(IN)   :: alpha_ecp, coeffs_ecp
      INTEGER, DIMENSION(1:npot_ecp), INTENT(IN)         :: nrpot_ecp
      REAL(KIND=dp), INTENT(IN)                          :: rpgfc
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), INTENT(IN)                          :: dab
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rac
      REAL(KIND=dp), INTENT(IN)                          :: dac
      REAL(KIND=dp), INTENT(IN)                          :: dbc
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: vab

#if defined(__LIBGRPP)
      INTEGER                                            :: a_offset, a_start, b_offset, b_start, i, &
                                                            ipgf, j, jpgf, li, lj, ncoa, ncob
      REAL(dp)                                           :: expi, expj, normi, normj, prefi, prefj, &
                                                            zeti, zetj
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: tmp
      REAL(dp), DIMENSION(3)                             :: ra, rb, rc

      !libgrpp requires absolute positions, not relative ones
      ra(:) = 0.0_dp
      rb(:) = rab(:)
      rc(:) = rac(:)

      DO ipgf = 1, npgfa
         IF (rpgfa(ipgf) + rpgfc < dac) CYCLE
         zeti = zeta(ipgf)
         a_start = (ipgf - 1)*ncoset(la_max_set)

         DO jpgf = 1, npgfb
            IF (rpgfb(jpgf) + rpgfc < dbc) CYCLE
            IF (rpgfa(ipgf) + rpgfb(jpgf) < dab) CYCLE
            zetj = zetb(jpgf)
            b_start = (jpgf - 1)*ncoset(lb_max_set)

            DO li = la_min_set, la_max_set
               a_offset = a_start + ncoset(li - 1)
               ncoa = nco(li)
               prefi = 2.0_dp**li*(2.0_dp/pi)**0.75_dp
               expi = 0.25_dp*REAL(2*li + 3, dp)
               normi = 1.0_dp/(prefi*zeti**expi)

               DO lj = lb_min_set, lb_max_set
                  b_offset = b_start + ncoset(lj - 1)
                  ncob = nco(lj)
                  prefj = 2.0_dp**lj*(2.0_dp/pi)**0.75_dp
                  expj = 0.25_dp*REAL(2*lj + 3, dp)
                  normj = 1.0_dp/(prefj*zetj**expj)

                  ALLOCATE (tmp(ncoa*ncob))
                  tmp = 0.0_dp
                  !libgrpp implicitely normalizes cartesian Gaussian. In CP2K, we do not, hence
                  !the 1/norm coefficients for PGFi and PGFj
                  CALL libgrpp_type1_integrals(ra, li, 1, [normi], [zeti], &
                                               rb, lj, 1, [normj], [zetj], &
                                               rc, [npot_ecp], nrpot_ecp, &
                                               coeffs_ecp, alpha_ecp, tmp)

                  !note: tmp array is in C row-major ordering
                  DO j = 1, ncob
                     DO i = 1, ncoa
                        vab(a_offset + i, b_offset + j) = vab(a_offset + i, b_offset + j) + tmp((i - 1)*ncob + j)
                     END DO
                  END DO
                  DEALLOCATE (tmp)

               END DO !lj
            END DO !li

         END DO !jpgf
      END DO !ipgf
#else

      MARK_USED(la_max_set)
      MARK_USED(la_min_set)
      MARK_USED(npgfa)
      MARK_USED(rpgfa)
      MARK_USED(zeta)
      MARK_USED(lb_max_set)
      MARK_USED(lb_min_set)
      MARK_USED(npgfb)
      MARK_USED(rpgfb)
      MARK_USED(zetb)
      MARK_USED(npot_ecp)
      MARK_USED(alpha_ecp)
      MARK_USED(coeffs_ecp)
      MARK_USED(nrpot_ecp)
      MARK_USED(rpgfc)
      MARK_USED(rab)
      MARK_USED(dab)
      MARK_USED(rac)
      MARK_USED(dac)
      MARK_USED(dbc)
      MARK_USED(vab)

      CPABORT("Please compile CP2K with libgrpp support for calculations with ECPs")
#endif

   END SUBROUTINE libgrpp_local_integral

! **************************************************************************************************
!> \brief Semi-local ECP integrals using libgrpp.
!> \param la_max_set ...
!> \param la_min_set ...
!> \param npgfa ...
!> \param rpgfa ...
!> \param zeta ...
!> \param lb_max_set ...
!> \param lb_min_set ...
!> \param npgfb ...
!> \param rpgfb ...
!> \param zetb ...
!> \param lmax_ecp ...
!> \param npot_ecp ...
!> \param alpha_ecp ...
!> \param coeffs_ecp ...
!> \param nrpot_ecp ...
!> \param rpgfc ...
!> \param rab ...
!> \param dab ...
!> \param rac ...
!> \param dac ...
!> \param dbc ...
!> \param vab ...
! **************************************************************************************************
   SUBROUTINE libgrpp_semilocal_integral(la_max_set, la_min_set, npgfa, rpgfa, zeta, &
                                         lb_max_set, lb_min_set, npgfb, rpgfb, zetb, &
                                         lmax_ecp, npot_ecp, alpha_ecp, coeffs_ecp, nrpot_ecp, &
                                         rpgfc, rab, dab, rac, dac, dbc, vab)

      INTEGER, INTENT(IN)                                :: la_max_set, la_min_set, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rpgfa, zeta
      INTEGER, INTENT(IN)                                :: lb_max_set, lb_min_set, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rpgfb, zetb
      INTEGER, INTENT(IN)                                :: lmax_ecp
      INTEGER, DIMENSION(0:10), INTENT(IN)               :: npot_ecp
      REAL(KIND=dp), DIMENSION(1:15, 0:10), INTENT(IN)   :: alpha_ecp, coeffs_ecp
      INTEGER, DIMENSION(1:15, 0:10), INTENT(IN)         :: nrpot_ecp
      REAL(KIND=dp), INTENT(IN)                          :: rpgfc
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), INTENT(IN)                          :: dab
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rac
      REAL(KIND=dp), INTENT(IN)                          :: dac
      REAL(KIND=dp), INTENT(IN)                          :: dbc
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: vab

#if defined(__LIBGRPP)
      INTEGER                                            :: a_offset, a_start, b_offset, b_start, i, &
                                                            ipgf, j, jpgf, li, lj, lk, ncoa, ncob
      REAL(dp)                                           :: expi, expj, normi, normj, prefi, prefj, &
                                                            zeti, zetj
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: tmp
      REAL(dp), DIMENSION(3)                             :: ra, rb, rc

      !libgrpp requires absolute positions, not relative ones
      ra(:) = 0.0_dp
      rb(:) = rab(:)
      rc(:) = rac(:)

      DO ipgf = 1, npgfa
         IF (rpgfa(ipgf) + rpgfc < dac) CYCLE
         zeti = zeta(ipgf)
         a_start = (ipgf - 1)*ncoset(la_max_set)

         DO jpgf = 1, npgfb
            IF (rpgfb(jpgf) + rpgfc < dbc) CYCLE
            IF (rpgfa(ipgf) + rpgfb(jpgf) < dab) CYCLE
            zetj = zetb(jpgf)
            b_start = (jpgf - 1)*ncoset(lb_max_set)

            DO li = la_min_set, la_max_set
               a_offset = a_start + ncoset(li - 1)
               ncoa = nco(li)
               prefi = 2.0_dp**li*(2.0_dp/pi)**0.75_dp
               expi = 0.25_dp*REAL(2*li + 3, dp)
               normi = 1.0_dp/(prefi*zeti**expi)

               DO lj = lb_min_set, lb_max_set
                  b_offset = b_start + ncoset(lj - 1)
                  ncob = nco(lj)
                  prefj = 2.0_dp**lj*(2.0_dp/pi)**0.75_dp
                  expj = 0.25_dp*REAL(2*lj + 3, dp)
                  normj = 1.0_dp/(prefj*zetj**expj)

                  ALLOCATE (tmp(ncoa*ncob))
                  !Loop over ECP angular momentum
                  DO lk = 0, lmax_ecp
                     tmp = 0.0_dp
                     !libgrpp implicitely normalizes cartesian Gaussian. In CP2K, we do not, hence
                     !the 1/norm coefficients for PGFi and PGFj
                     CALL libgrpp_type2_integrals(ra, li, 1, [normi], [zeti], &
                                                  rb, lj, 1, [normj], [zetj], &
                                                  rc, lk, [npot_ecp(lk)], nrpot_ecp(:, lk), &
                                                  coeffs_ecp(:, lk), alpha_ecp(:, lk), tmp)

                     !note: tmp array is in C row-major ordering
                     DO j = 1, ncob
                        DO i = 1, ncoa
                           vab(a_offset + i, b_offset + j) = vab(a_offset + i, b_offset + j) + tmp((i - 1)*ncob + j)
                        END DO
                     END DO
                  END DO !lk
                  DEALLOCATE (tmp)

               END DO !lj
            END DO !li

         END DO !jpgf
      END DO !ipgf

#else

      MARK_USED(la_max_set)
      MARK_USED(la_min_set)
      MARK_USED(npgfa)
      MARK_USED(rpgfa)
      MARK_USED(zeta)
      MARK_USED(lb_max_set)
      MARK_USED(lb_min_set)
      MARK_USED(npgfb)
      MARK_USED(rpgfb)
      MARK_USED(zetb)
      MARK_USED(lmax_ecp)
      MARK_USED(npot_ecp)
      MARK_USED(alpha_ecp)
      MARK_USED(coeffs_ecp)
      MARK_USED(nrpot_ecp)
      MARK_USED(rpgfc)
      MARK_USED(rab)
      MARK_USED(dab)
      MARK_USED(rac)
      MARK_USED(dac)
      MARK_USED(dbc)
      MARK_USED(vab)

      CPABORT("Please compile CP2K with libgrpp support for calculations with ECPs")
#endif

   END SUBROUTINE libgrpp_semilocal_integral

END MODULE libgrpp_integrals
