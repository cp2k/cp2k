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
   USE ai_derivatives, ONLY: dabdr_noscreen, adbdr, dabdr
   USE orbital_pointers, ONLY: nco, &
                               ncoset
#if defined(__LIBGRPP)
   USE libgrpp, ONLY: libgrpp_init, libgrpp_type1_integrals, libgrpp_type2_integrals, &
                      libgrpp_type1_integrals_gradient, libgrpp_type2_integrals_gradient
#endif
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'libgrpp_integrals'

   PUBLIC :: libgrpp_semilocal_integrals, libgrpp_local_integrals, &
             libgrpp_local_forces_ref, libgrpp_semilocal_forces_ref

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
!> \param pab ...
!> \param force_a ...
!> \param force_b ...
! **************************************************************************************************
   SUBROUTINE libgrpp_local_integrals(la_max_set, la_min_set, npgfa, rpgfa, zeta, &
                                      lb_max_set, lb_min_set, npgfb, rpgfb, zetb, &
                                      npot_ecp, alpha_ecp, coeffs_ecp, nrpot_ecp, &
                                      rpgfc, rab, dab, rac, dac, dbc, vab, pab, force_a, force_b)

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
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: pab
      REAL(KIND=dp), DIMENSION(3), INTENT(INOUT), &
         OPTIONAL                                        :: force_a, force_b

#if defined(__LIBGRPP)
      INTEGER                                            :: a_offset, a_start, b_offset, b_start, i, &
                                                            ipgf, j, jpgf, li, lj, ncoa, ncob
      LOGICAL                                            :: calc_forces
      REAL(dp)                                           :: expi, expj, normi, normj, prefi, prefj, &
                                                            zeti, zetj, mindist, fac_a, fac_b
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: tmp, tmpx, tmpy, tmpz
      REAL(dp), DIMENSION(3)                             :: ra, rb, rc

      CALL libgrpp_init()

      calc_forces = .FALSE.
      IF (PRESENT(pab) .AND. PRESENT(force_a) .AND. PRESENT(force_b)) calc_forces = .TRUE.

      IF (calc_forces) THEN

         !Note: warning against numerical stability of libgrpp gradients. The day the library becomes
         !      stable, this routine can be used immediatly as is, and the warning removed.
         CALL cp_warn(__LOCATION__, &
                      "ECP gradients calculated with the libgrpp library are, to this date, not numerically stable. "// &
                      "Please use the reference routine 'libgrpp_local_forces_ref' instead.")

         !there is a weird feature of libgrpp gradients, which is such that the gradient is calculated
         !for a point in space, and not with respect to an atomic center. For example, if atoms A and
         !B are the same (and C is different), then d<A | U_C | B>/dPx = d<A | U_C | B>/dAx + d<A | U_C | B>/dBx
         !Because we want the forces on centers A and B seprately, we need a case study on atomic positions
         !We always calculate the gradient wrt to atomic position of A and B, and we scale accordingly

         mindist = 1.0E-6_dp
         !If ra != rb != rc
         IF (dab >= mindist .AND. dbc >= mindist .AND. dac >= mindist) THEN
            fac_a = 1.0_dp
            fac_b = 1.0_dp

            !If ra = rb, but ra != rc
         ELSE IF (dab < mindist .AND. dac >= mindist) THEN
            fac_a = 0.5_dp
            fac_b = 0.5_dp

            !IF ra != rb but ra = rc
         ELSE IF (dab >= mindist .AND. dac < mindist) THEN
            fac_a = 0.5_dp
            fac_b = 1.0_dp

            !IF ra != rb but rb = rc
         ELSE IF (dab >= mindist .AND. dbc < mindist) THEN
            fac_a = 1.0_dp
            fac_b = 0.5_dp

            !If all atoms the same --> no force
         ELSE
            calc_forces = .FALSE.
         END IF
      END IF

      !libgrpp requires absolute positions, not relative ones
      ra(:) = 0.0_dp
      rb(:) = rab(:)
      rc(:) = rac(:)

      ALLOCATE (tmp(nco(la_max_set)*nco(lb_max_set)))
      IF (calc_forces) THEN
         ALLOCATE (tmpx(nco(la_max_set)*nco(lb_max_set)))
         ALLOCATE (tmpy(nco(la_max_set)*nco(lb_max_set)))
         ALLOCATE (tmpz(nco(la_max_set)*nco(lb_max_set)))
      END IF

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

                  tmp(1:ncoa*ncob) = 0.0_dp
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

                  IF (calc_forces) THEN
                     tmpx(1:ncoa*ncob) = 0.0_dp
                     tmpy(1:ncoa*ncob) = 0.0_dp
                     tmpz(1:ncoa*ncob) = 0.0_dp

                     !force wrt to atomic position A
                     CALL libgrpp_type1_integrals_gradient(ra, li, 1, [normi], [zeti], &
                                                           rb, lj, 1, [normj], [zetj], &
                                                           rc, [npot_ecp], nrpot_ecp, &
                                                           coeffs_ecp, alpha_ecp, ra, &
                                                           tmpx, tmpy, tmpz)

                     !note: tmp array is in C row-major ordering
                     !note: zero-gradients sometime comes out as NaN, hence tampval==tmpval check
                     DO j = 1, ncob
                        DO i = 1, ncoa
                           force_a(1) = force_a(1) + fac_a*pab(a_offset + i, b_offset + j)*tmpx((i - 1)*ncob + j)
                           force_a(2) = force_a(2) + fac_a*pab(a_offset + i, b_offset + j)*tmpy((i - 1)*ncob + j)
                           force_a(3) = force_a(3) + fac_a*pab(a_offset + i, b_offset + j)*tmpz((i - 1)*ncob + j)
                        END DO
                     END DO

                     tmpx(1:ncoa*ncob) = 0.0_dp
                     tmpy(1:ncoa*ncob) = 0.0_dp
                     tmpz(1:ncoa*ncob) = 0.0_dp

                     !force wrt to atomic position B
                     CALL libgrpp_type1_integrals_gradient(ra, li, 1, [normi], [zeti], &
                                                           rb, lj, 1, [normj], [zetj], &
                                                           rc, [npot_ecp], nrpot_ecp, &
                                                           coeffs_ecp, alpha_ecp, rb, &
                                                           tmpx, tmpy, tmpz)

                     !note: tmp array is in C row-major ordering
                     !note: zero-gradients sometime comes out as NaN, hence tampval==tmpval check
                     DO j = 1, ncob
                        DO i = 1, ncoa
                           force_b(1) = force_b(1) + fac_b*pab(a_offset + i, b_offset + j)*tmpx((i - 1)*ncob + j)
                           force_b(2) = force_b(2) + fac_b*pab(a_offset + i, b_offset + j)*tmpy((i - 1)*ncob + j)
                           force_b(3) = force_b(3) + fac_b*pab(a_offset + i, b_offset + j)*tmpz((i - 1)*ncob + j)
                        END DO
                     END DO
                  END IF

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
      MARK_USED(pab)
      MARK_USED(force_a)
      MARK_USED(force_b)

      CPABORT("Please compile CP2K with libgrpp support for calculations with ECPs")
#endif

   END SUBROUTINE libgrpp_local_integrals

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
!> \param pab ...
!> \param force_a ...
!> \param force_b ...
! **************************************************************************************************
   SUBROUTINE libgrpp_semilocal_integrals(la_max_set, la_min_set, npgfa, rpgfa, zeta, &
                                          lb_max_set, lb_min_set, npgfb, rpgfb, zetb, &
                                          lmax_ecp, npot_ecp, alpha_ecp, coeffs_ecp, nrpot_ecp, &
                                          rpgfc, rab, dab, rac, dac, dbc, vab, pab, force_a, force_b)

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
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: pab
      REAL(KIND=dp), DIMENSION(3), INTENT(INOUT), &
         OPTIONAL                                        :: force_a, force_b

#if defined(__LIBGRPP)
      INTEGER                                            :: a_offset, a_start, b_offset, b_start, i, &
                                                            ipgf, j, jpgf, li, lj, lk, ncoa, ncob
      LOGICAL                                            :: calc_forces
      REAL(dp)                                           :: expi, expj, normi, normj, prefi, prefj, &
                                                            zeti, zetj, mindist, fac_a, fac_b
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: tmp, tmpx, tmpz, tmpy
      REAL(dp), DIMENSION(3)                             :: ra, rb, rc

      CALL libgrpp_init()

      calc_forces = .FALSE.
      IF (PRESENT(pab) .AND. PRESENT(force_a) .AND. PRESENT(force_b)) calc_forces = .TRUE.

      IF (calc_forces) THEN

         !Note: warning against numerical stability of libgrpp gradients. The day the library becomes
         !      stable, this routine can be used immediatly as is, and the warning removed.
         CALL cp_warn(__LOCATION__, &
                      "ECP gradients calculated with the libgrpp library are, to this date, not numerically stable. "// &
                      "Please use the reference routine 'libgrpp_semilocal_forces_ref' instead.")

         !there is a weird feature of libgrpp gradients, which is such that the gradient is calculated
         !for a point in space, and not with respect to an atomic center. For example, if atoms A and
         !B are the same (and C is different), then d<A | U_C | B>/dPx = d<A | U_C | B>/dAx + d<A | U_C | B>/dBx
         !Because we want the forces on centers A and B seprately, we need a case study on atomic positions
         !We always calculate the gradient wrt to atomic position of A and B, and we scale accordingly

         mindist = 1.0E-6_dp
         !If ra != rb != rc
         IF (dab >= mindist .AND. dbc >= mindist .AND. dac >= mindist) THEN
            fac_a = 1.0_dp
            fac_b = 1.0_dp

            !If ra = rb, but ra != rc
         ELSE IF (dab < mindist .AND. dac >= mindist) THEN
            fac_a = 0.5_dp
            fac_b = 0.5_dp

            !IF ra != rb but ra = rc
         ELSE IF (dab >= mindist .AND. dac < mindist) THEN
            fac_a = 0.5_dp
            fac_b = 1.0_dp

            !IF ra != rb but rb = rc
         ELSE IF (dab >= mindist .AND. dbc < mindist) THEN
            fac_a = 1.0_dp
            fac_b = 0.5_dp

            !If all atoms the same --> no force
         ELSE
            calc_forces = .FALSE.
         END IF
      END IF

      !libgrpp requires absolute positions, not relative ones
      ra(:) = 0.0_dp
      rb(:) = rab(:)
      rc(:) = rac(:)

      ALLOCATE (tmp(nco(la_max_set)*nco(lb_max_set)))
      IF (calc_forces) THEN
         ALLOCATE (tmpx(nco(la_max_set)*nco(lb_max_set)))
         ALLOCATE (tmpy(nco(la_max_set)*nco(lb_max_set)))
         ALLOCATE (tmpz(nco(la_max_set)*nco(lb_max_set)))
      END IF

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

                  !Loop over ECP angular momentum
                  DO lk = 0, lmax_ecp
                     tmp(1:ncoa*ncob) = 0.0_dp
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

                     IF (calc_forces) THEN

                        tmpx(1:ncoa*ncob) = 0.0_dp
                        tmpy(1:ncoa*ncob) = 0.0_dp
                        tmpz(1:ncoa*ncob) = 0.0_dp

                        !force wrt to atomic position A
                        CALL libgrpp_type2_integrals_gradient(ra, li, 1, [normi], [zeti], &
                                                              rb, lj, 1, [normj], [zetj], &
                                                              rc, lk, [npot_ecp(lk)], nrpot_ecp(:, lk), &
                                                              coeffs_ecp(:, lk), alpha_ecp(:, lk), ra, &
                                                              tmpx, tmpy, tmpz)

                        !note: tmp array is in C row-major ordering
                        !note: zero-gradients sometime comes out as NaN, hence tampval==tmpval check
                        DO j = 1, ncob
                           DO i = 1, ncoa
                              force_a(1) = force_a(1) + fac_a*pab(a_offset + i, b_offset + j)*tmpx((i - 1)*ncob + j)
                              force_a(2) = force_a(2) + fac_a*pab(a_offset + i, b_offset + j)*tmpy((i - 1)*ncob + j)
                              force_a(3) = force_a(3) + fac_a*pab(a_offset + i, b_offset + j)*tmpz((i - 1)*ncob + j)
                           END DO
                        END DO

                        tmpx(1:ncoa*ncob) = 0.0_dp
                        tmpy(1:ncoa*ncob) = 0.0_dp
                        tmpz(1:ncoa*ncob) = 0.0_dp

                        !force wrt to atomic position B
                        CALL libgrpp_type2_integrals_gradient(ra, li, 1, [normi], [zeti], &
                                                              rb, lj, 1, [normj], [zetj], &
                                                              rc, lk, [npot_ecp(lk)], nrpot_ecp(:, lk), &
                                                              coeffs_ecp(:, lk), alpha_ecp(:, lk), rb, &
                                                              tmpx, tmpy, tmpz)
                        !note: tmp array is in C row-major ordering
                        !note: zero-gradients sometime comes out as NaN, hence tampval==tmpval check
                        DO j = 1, ncob
                           DO i = 1, ncoa
                              force_b(1) = force_b(1) + fac_b*pab(a_offset + i, b_offset + j)*tmpx((i - 1)*ncob + j)
                              force_b(2) = force_b(2) + fac_b*pab(a_offset + i, b_offset + j)*tmpy((i - 1)*ncob + j)
                              force_b(3) = force_b(3) + fac_b*pab(a_offset + i, b_offset + j)*tmpz((i - 1)*ncob + j)
                           END DO
                        END DO

                     END IF !calc_forces

                  END DO !lk

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
      MARK_USED(pab)
      MARK_USED(force_a)
      MARK_USED(force_b)

      CPABORT("Please compile CP2K with libgrpp support for calculations with ECPs")
#endif

   END SUBROUTINE libgrpp_semilocal_integrals

! **************************************************************************************************
!> \brief Reference local ECP force routine using l+-1 integrals. No call is made to the numerically
!>        unstable gradient routine of libgrpp. Calculates both the integrals and the forces.
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
!> \param pab ...
!> \param force_a ...
!> \param force_b ...
!> \note: this is a reference routine, which has no reason to be used once the libgrpp gradients
!>        become numerically stable
! **************************************************************************************************
   SUBROUTINE libgrpp_local_forces_ref(la_max_set, la_min_set, npgfa, rpgfa, zeta, &
                                       lb_max_set, lb_min_set, npgfb, rpgfb, zetb, &
                                       npot_ecp, alpha_ecp, coeffs_ecp, nrpot_ecp, &
                                       rpgfc, rab, dab, rac, dac, dbc, vab, pab, force_a, force_b)

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
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: pab
      REAL(KIND=dp), DIMENSION(3), INTENT(INOUT)         :: force_a, force_b

#if defined(__LIBGRPP)
      INTEGER                                            :: a_offset, a_start, b_offset, b_start, i, &
                                                            ipgf, j, jpgf, li, lj, ncoa, ncob, a_offset_f, &
                                                            b_offset_f, a_start_f, b_start_f
      REAL(dp)                                           :: expi, expj, normi, normj, prefi, prefj, &
                                                            zeti, zetj
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: tmp
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: vab_f, tmpx, tmpy, tmpz
      REAL(dp), DIMENSION(3)                             :: ra, rb, rc

      CALL libgrpp_init()

      !Contains the integrals necessary for the forces, with angular momenta from lmin-1 to lmax+1
      ALLOCATE (vab_f(npgfa*ncoset(la_max_set + 1), npgfb*ncoset(lb_max_set + 1)))
      vab_f(:, :) = 0.0_dp

      !libgrpp requires absolute positions, not relative ones
      ra(:) = 0.0_dp
      rb(:) = rab(:)
      rc(:) = rac(:)

      ALLOCATE (tmp(nco(la_max_set + 1)*nco(lb_max_set + 1)))

      DO ipgf = 1, npgfa
         IF (rpgfa(ipgf) + rpgfc < dac) CYCLE
         zeti = zeta(ipgf)
         a_start = (ipgf - 1)*ncoset(la_max_set)
         a_start_f = (ipgf - 1)*ncoset(la_max_set + 1)

         DO jpgf = 1, npgfb
            IF (rpgfb(jpgf) + rpgfc < dbc) CYCLE
            IF (rpgfa(ipgf) + rpgfb(jpgf) < dab) CYCLE
            zetj = zetb(jpgf)
            b_start = (jpgf - 1)*ncoset(lb_max_set)
            b_start_f = (jpgf - 1)*ncoset(lb_max_set + 1)

            DO li = MAX(0, la_min_set - 1), la_max_set + 1
               a_offset = a_start + ncoset(li - 1)
               a_offset_f = a_start_f + ncoset(li - 1)
               ncoa = nco(li)
               prefi = 2.0_dp**li*(2.0_dp/pi)**0.75_dp
               expi = 0.25_dp*REAL(2*li + 3, dp)
               normi = 1.0_dp/(prefi*zeti**expi)

               DO lj = MAX(0, lb_min_set - 1), lb_max_set + 1
                  b_offset = b_start + ncoset(lj - 1)
                  b_offset_f = b_start_f + ncoset(lj - 1)
                  ncob = nco(lj)
                  prefj = 2.0_dp**lj*(2.0_dp/pi)**0.75_dp
                  expj = 0.25_dp*REAL(2*lj + 3, dp)
                  normj = 1.0_dp/(prefj*zetj**expj)

                  tmp(1:ncoa*ncob) = 0.0_dp
                  !libgrpp implicitely normalizes cartesian Gaussian. In CP2K, we do not, hence
                  !the 1/norm coefficients for PGFi and PGFj
                  CALL libgrpp_type1_integrals(ra, li, 1, [normi], [zeti], &
                                               rb, lj, 1, [normj], [zetj], &
                                               rc, [npot_ecp], nrpot_ecp, &
                                               coeffs_ecp, alpha_ecp, tmp)

                  !the l+-1 integrals for gradient calculation
                  DO j = 1, ncob
                     DO i = 1, ncoa
                        vab_f(a_offset_f + i, b_offset_f + j) = &
                           vab_f(a_offset_f + i, b_offset_f + j) + tmp((i - 1)*ncob + j)
                     END DO
                  END DO

                  !the actual integrals
                  IF (li >= la_min_set .AND. li <= la_max_set .AND. lj >= lb_min_set .AND. lj <= lb_max_set) THEN
                     DO j = 1, ncob
                        DO i = 1, ncoa
                           vab(a_offset + i, b_offset + j) = vab(a_offset + i, b_offset + j) + tmp((i - 1)*ncob + j)
                        END DO
                     END DO
                  END IF

               END DO !lj
            END DO !li

         END DO !jpgf
      END DO !ipgf

      ALLOCATE (tmpx(npgfa*ncoset(la_max_set), npgfb*ncoset(lb_max_set)))
      ALLOCATE (tmpy(npgfa*ncoset(la_max_set), npgfb*ncoset(lb_max_set)))
      ALLOCATE (tmpz(npgfa*ncoset(la_max_set), npgfb*ncoset(lb_max_set)))

      !Derivative wrt to center A
      tmpx(:, :) = 0.0_dp
      tmpy(:, :) = 0.0_dp
      tmpz(:, :) = 0.0_dp
      CALL dabdr(la_max_set, npgfa, zeta, rpgfa, la_min_set, lb_max_set, npgfb, rpgfb, lb_min_set, &
                 dab, vab_f, tmpx, tmpy, tmpz)
      DO j = 1, npgfb*ncoset(lb_max_set)
         DO i = 1, npgfa*ncoset(la_max_set)
            force_a(1) = force_a(1) + tmpx(i, j)*pab(i, j)
            force_a(2) = force_a(2) + tmpy(i, j)*pab(i, j)
            force_a(3) = force_a(3) + tmpz(i, j)*pab(i, j)
         END DO
      END DO

      !Derivative wrt to center B
      tmpx(:, :) = 0.0_dp
      tmpy(:, :) = 0.0_dp
      tmpz(:, :) = 0.0_dp
      CALL adbdr(la_max_set, npgfa, rpgfa, la_min_set, lb_max_set, npgfb, zetb, rpgfb, lb_min_set, &
                 dab, vab_f, tmpx, tmpy, tmpz)
      DO j = 1, npgfb*ncoset(lb_max_set)
         DO i = 1, npgfa*ncoset(la_max_set)
            force_b(1) = force_b(1) + tmpx(i, j)*pab(i, j)
            force_b(2) = force_b(2) + tmpy(i, j)*pab(i, j)
            force_b(3) = force_b(3) + tmpz(i, j)*pab(i, j)
         END DO
      END DO
      DEALLOCATE (tmpx, tmpy, tmpz)

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
      MARK_USED(pab)
      MARK_USED(vab)
      MARK_USED(force_a)
      MARK_USED(force_b)

      CPABORT("Please compile CP2K with libgrpp support for calculations with ECPs")
#endif

   END SUBROUTINE libgrpp_local_forces_ref

! **************************************************************************************************
!> \brief Reference semi-local ECP forces using l+-1 integrals. No call is made to the numerically
!>        unstable gradient routine of libgrpp. Calculates both the integrals and the forces.
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
!> \param pab ...
!> \param force_a ...
!> \param force_b ...
!> \note: this is a reference routine, which has no reason to be used once the libgrpp gradients
!>        become numerically stable
! **************************************************************************************************
   SUBROUTINE libgrpp_semilocal_forces_ref(la_max_set, la_min_set, npgfa, rpgfa, zeta, &
                                           lb_max_set, lb_min_set, npgfb, rpgfb, zetb, &
                                           lmax_ecp, npot_ecp, alpha_ecp, coeffs_ecp, nrpot_ecp, &
                                           rpgfc, rab, dab, rac, dac, dbc, vab, pab, force_a, force_b)

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
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: pab
      REAL(KIND=dp), DIMENSION(3), INTENT(INOUT)         :: force_a, force_b

#if defined(__LIBGRPP)
      INTEGER                                            :: a_offset, a_start, b_offset, b_start, i, &
                                                            ipgf, j, jpgf, li, lj, lk, ncoa, ncob, &
                                                            a_start_f, b_start_f, a_offset_f, b_offset_f
      REAL(dp)                                           :: expi, expj, normi, normj, prefi, prefj, &
                                                            zeti, zetj
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: tmp
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: vab_f, tmpx, tmpy, tmpz
      REAL(dp), DIMENSION(3)                             :: ra, rb, rc

      CALL libgrpp_init()

      !Contains the integrals necessary for the forces, with angular momenta from lmin-1 to lmax+1
      ALLOCATE (vab_f(npgfa*ncoset(la_max_set + 1), npgfb*ncoset(lb_max_set + 1)))
      vab_f(:, :) = 0.0_dp

      !libgrpp requires absolute positions, not relative ones
      ra(:) = 0.0_dp
      rb(:) = rab(:)
      rc(:) = rac(:)

      ALLOCATE (tmp(nco(la_max_set + 1)*nco(lb_max_set + 1)))

      DO ipgf = 1, npgfa
         IF (rpgfa(ipgf) + rpgfc < dac) CYCLE
         zeti = zeta(ipgf)
         a_start = (ipgf - 1)*ncoset(la_max_set)
         a_start_f = (ipgf - 1)*ncoset(la_max_set + 1)

         DO jpgf = 1, npgfb
            IF (rpgfb(jpgf) + rpgfc < dbc) CYCLE
            IF (rpgfa(ipgf) + rpgfb(jpgf) < dab) CYCLE
            zetj = zetb(jpgf)
            b_start = (jpgf - 1)*ncoset(lb_max_set)
            b_start_f = (jpgf - 1)*ncoset(lb_max_set + 1)

            DO li = MAX(0, la_min_set - 1), la_max_set + 1
               a_offset = a_start + ncoset(li - 1)
               a_offset_f = a_start_f + ncoset(li - 1)
               ncoa = nco(li)
               prefi = 2.0_dp**li*(2.0_dp/pi)**0.75_dp
               expi = 0.25_dp*REAL(2*li + 3, dp)
               normi = 1.0_dp/(prefi*zeti**expi)

               DO lj = MAX(0, lb_min_set - 1), lb_max_set + 1
                  b_offset = b_start + ncoset(lj - 1)
                  b_offset_f = b_start_f + ncoset(lj - 1)
                  ncob = nco(lj)
                  prefj = 2.0_dp**lj*(2.0_dp/pi)**0.75_dp
                  expj = 0.25_dp*REAL(2*lj + 3, dp)
                  normj = 1.0_dp/(prefj*zetj**expj)

                  !Loop over ECP angular momentum
                  DO lk = 0, lmax_ecp
                     tmp(1:ncoa*ncob) = 0.0_dp
                     !libgrpp implicitely normalizes cartesian Gaussian. In CP2K, we do not, hence
                     !the 1/norm coefficients for PGFi and PGFj
                     CALL libgrpp_type2_integrals(ra, li, 1, [normi], [zeti], &
                                                  rb, lj, 1, [normj], [zetj], &
                                                  rc, lk, [npot_ecp(lk)], nrpot_ecp(:, lk), &
                                                  coeffs_ecp(:, lk), alpha_ecp(:, lk), tmp)

                     !the l+-1 integrals for gradient calculation
                     DO j = 1, ncob
                        DO i = 1, ncoa
                           vab_f(a_offset_f + i, b_offset_f + j) = &
                              vab_f(a_offset_f + i, b_offset_f + j) + tmp((i - 1)*ncob + j)
                        END DO
                     END DO

                     !the actual integrals
                     IF (li >= la_min_set .AND. li <= la_max_set .AND. lj >= lb_min_set .AND. lj <= lb_max_set) THEN
                        DO j = 1, ncob
                           DO i = 1, ncoa
                              vab(a_offset + i, b_offset + j) = vab(a_offset + i, b_offset + j) + tmp((i - 1)*ncob + j)
                           END DO
                        END DO
                     END IF

                  END DO !lk

               END DO !lj
            END DO !li

         END DO !jpgf
      END DO !ipgf

      ALLOCATE (tmpx(npgfa*ncoset(la_max_set), npgfb*ncoset(lb_max_set)))
      ALLOCATE (tmpy(npgfa*ncoset(la_max_set), npgfb*ncoset(lb_max_set)))
      ALLOCATE (tmpz(npgfa*ncoset(la_max_set), npgfb*ncoset(lb_max_set)))

      !Derivative wrt to center A
      tmpx(:, :) = 0.0_dp
      tmpy(:, :) = 0.0_dp
      tmpz(:, :) = 0.0_dp
      CALL dabdr(la_max_set, npgfa, zeta, rpgfa, la_min_set, lb_max_set, npgfb, rpgfb, lb_min_set, &
                 0.0_dp, vab_f, tmpx, tmpy, tmpz)
      DO j = 1, npgfb*ncoset(lb_max_set)
         DO i = 1, npgfa*ncoset(la_max_set)
            force_a(1) = force_a(1) + tmpx(i, j)*pab(i, j)
            force_a(2) = force_a(2) + tmpy(i, j)*pab(i, j)
            force_a(3) = force_a(3) + tmpz(i, j)*pab(i, j)
         END DO
      END DO

      !Derivative wrt to center B
      tmpx(:, :) = 0.0_dp
      tmpy(:, :) = 0.0_dp
      tmpz(:, :) = 0.0_dp
      CALL adbdr(la_max_set, npgfa, rpgfa, la_min_set, lb_max_set, npgfb, zetb, rpgfb, lb_min_set, &
                 0.0_dp, vab_f, tmpx, tmpy, tmpz)
      DO j = 1, npgfb*ncoset(lb_max_set)
         DO i = 1, npgfa*ncoset(la_max_set)
            force_b(1) = force_b(1) + tmpx(i, j)*pab(i, j)
            force_b(2) = force_b(2) + tmpy(i, j)*pab(i, j)
            force_b(3) = force_b(3) + tmpz(i, j)*pab(i, j)
         END DO
      END DO
      DEALLOCATE (tmpx, tmpy, tmpz)

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
      MARK_USED(pab)
      MARK_USED(vab)
      MARK_USED(force_a)
      MARK_USED(force_b)

      CPABORT("Please compile CP2K with libgrpp support for calculations with ECPs")
#endif

   END SUBROUTINE libgrpp_semilocal_forces_ref

END MODULE libgrpp_integrals
