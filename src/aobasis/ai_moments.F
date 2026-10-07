!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of the moment integrals over Cartesian Gaussian-type
!>      functions.
!> \par Literature
!>      S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!> \par History
!>      none
!> \author J. Hutter (16.02.2005)
! **************************************************************************************************
MODULE ai_moments

! ax,ay,az  : Angular momentum index numbers of orbital a.
! bx,by,bz  : Angular momentum index numbers of orbital b.
! coset     : Cartesian orbital set pointer.
! dab       : Distance between the atomic centers a and b.
! l{a,b}    : Angular momentum quantum number of shell a or b.
! l{a,b}_max: Maximum angular momentum quantum number of shell a or b.
! l{a,b}_min: Minimum angular momentum quantum number of shell a or b.
! rac       : Distance vector between the atomic center a and reference point c.
! rbc       : Distance vector between the atomic center b and reference point c.
! rpgf{a,b} : Radius of the primitive Gaussian-type function a or b.
! zet{a,b}  : Exponents of the Gaussian-type functions a or b.
! zetp      : Reciprocal of the sum of the exponents of orbital a and b.

   USE ai_cartesian_operators,          ONLY: apply_cartesian_operator,&
                                              cartesian_ladder,&
                                              cartesian_ladder_map,&
                                              cartesian_polynomial_1d
   USE ai_os_rr,                        ONLY: os_rr_ovlp
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
   USE orbital_pointers,                ONLY: coset,&
                                              indco,&
                                              ncoset
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_moments'

   PUBLIC :: cossin, moment, diff_momop, contract_cossin, dipole_force

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param cos_block ...
!> \param sin_block ...
!> \param iatom ...
!> \param ncoa ...
!> \param nsgfa ...
!> \param sgfa ...
!> \param sphi_a ...
!> \param ldsa ...
!> \param jatom ...
!> \param ncob ...
!> \param nsgfb ...
!> \param sgfb ...
!> \param sphi_b ...
!> \param ldsb ...
!> \param cosab ...
!> \param sinab ...
!> \param ldab ...
!> \param work ...
!> \param ldwork ...
!> \param ordered retain the ordered atom pair instead of mapping to an upper block
! **************************************************************************************************
   SUBROUTINE contract_cossin(cos_block, sin_block, &
                              iatom, ncoa, nsgfa, sgfa, sphi_a, ldsa, &
                              jatom, ncob, nsgfb, sgfb, sphi_b, ldsb, &
                              cosab, sinab, ldab, work, ldwork, ordered)

      REAL(dp), DIMENSION(:, :), POINTER                 :: cos_block, sin_block
      INTEGER, INTENT(IN)                                :: iatom, ncoa, nsgfa, sgfa
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: sphi_a
      INTEGER, INTENT(IN)                                :: ldsa, jatom, ncob, nsgfb, sgfb
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: sphi_b
      INTEGER, INTENT(IN)                                :: ldsb
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: cosab, sinab
      INTEGER, INTENT(IN)                                :: ldab
      REAL(dp), DIMENSION(:, :)                          :: work
      INTEGER, INTENT(IN)                                :: ldwork
      LOGICAL, INTENT(IN), OPTIONAL                      :: ordered

      LOGICAL                                            :: direct_block

      direct_block = iatom <= jatom
      IF (PRESENT(ordered)) direct_block = direct_block .OR. ordered

! Calculate cosine

      CALL dgemm("N", "N", ncoa, nsgfb, ncob, &
                 1.0_dp, cosab(1, 1), ldab, &
                 sphi_b(1, sgfb), ldsb, &
                 0.0_dp, work(1, 1), ldwork)

      IF (direct_block) THEN
         CALL dgemm("T", "N", nsgfa, nsgfb, ncoa, &
                    1.0_dp, sphi_a(1, sgfa), ldsa, &
                    work(1, 1), ldwork, &
                    1.0_dp, cos_block(sgfa, sgfb), &
                    SIZE(cos_block, 1))
      ELSE
         CALL dgemm("T", "N", nsgfb, nsgfa, ncoa, &
                    1.0_dp, work(1, 1), ldwork, &
                    sphi_a(1, sgfa), ldsa, &
                    1.0_dp, cos_block(sgfb, sgfa), &
                    SIZE(cos_block, 1))
      END IF

      ! Calculate sine
      CALL dgemm("N", "N", ncoa, nsgfb, ncob, &
                 1.0_dp, sinab(1, 1), ldab, &
                 sphi_b(1, sgfb), ldsb, &
                 0.0_dp, work(1, 1), ldwork)

      IF (direct_block) THEN
         CALL dgemm("T", "N", nsgfa, nsgfb, ncoa, &
                    1.0_dp, sphi_a(1, sgfa), ldsa, &
                    work(1, 1), ldwork, &
                    1.0_dp, sin_block(sgfa, sgfb), &
                    SIZE(sin_block, 1))
      ELSE
         CALL dgemm("T", "N", nsgfb, nsgfa, ncoa, &
                    1.0_dp, work(1, 1), ldwork, &
                    sphi_a(1, sgfa), ldsa, &
                    1.0_dp, sin_block(sgfb, sgfa), &
                    SIZE(sin_block, 1))
      END IF

   END SUBROUTINE contract_cossin

! **************************************************************************************************
!> \brief ...
!> \param la_max_set ...
!> \param npgfa ...
!> \param zeta ...
!> \param rpgfa ...
!> \param la_min_set ...
!> \param lb_max ...
!> \param npgfb ...
!> \param zetb ...
!> \param rpgfb ...
!> \param lb_min ...
!> \param rac ...
!> \param rbc ...
!> \param kvec ...
!> \param cosab ...
!> \param sinab ...
!> \param dcosab ...
!> \param dsinab ...
! **************************************************************************************************
   SUBROUTINE cossin(la_max_set, npgfa, zeta, rpgfa, la_min_set, &
                     lb_max, npgfb, zetb, rpgfb, lb_min, &
                     rac, rbc, kvec, cosab, sinab, dcosab, dsinab)

      INTEGER, INTENT(IN)                                :: la_max_set, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta, rpgfa
      INTEGER, INTENT(IN)                                :: la_min_set, lb_max, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb, rpgfb
      INTEGER, INTENT(IN)                                :: lb_min
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rac, rbc, kvec
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: cosab, sinab
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), &
         OPTIONAL                                        :: dcosab, dsinab

      INTEGER                                            :: ax, ay, az, bx, by, bz, i, ipgf, j, &
                                                            jpgf, k, la, la_max, la_min, la_start, &
                                                            lb, na, nb
      INTEGER, ALLOCATABLE                               :: derivative_map(:, :, :)
      REAL(KIND=dp)                                      :: dab, f0, f1, f2, f3, fx, fy, fz, k2, &
                                                            kdp, rab2, s, zetp
      REAL(KIND=dp), DIMENSION(3)                        :: rab, rap, rbp
      REAL(KIND=dp), &
         DIMENSION(ncoset(la_max_set+1), ncoset(lb_max)) :: sc, ss

      rab = rbc - rac
      rab2 = SUM(rab**2)
      dab = SQRT(rab2)
      k2 = kvec(1)*kvec(1) + kvec(2)*kvec(2) + kvec(3)*kvec(3)

      IF (PRESENT(dcosab)) THEN
         la_max = la_max_set + 1
         la_min = MAX(0, la_min_set - 1)
         ALLOCATE (derivative_map(4, ncoset(la_max_set), 3))
         DO k = 1, 3
            CALL cartesian_ladder_map(la_max_set, 0, k, derivative_map(:, :, k))
         END DO
      ELSE
         la_max = la_max_set
         la_min = la_min_set
      END IF

      ! initialize all matrix elements to zero
      na = ncoset(la_max_set)*npgfa
      nb = ncoset(lb_max)*npgfb
      cosab(1:na, 1:nb) = 0.0_dp
      sinab(1:na, 1:nb) = 0.0_dp
      IF (PRESENT(dcosab)) THEN
         dcosab(1:na, 1:nb, :) = 0.0_dp
         dsinab(1:na, 1:nb, :) = 0.0_dp
      END IF
!   *** Loop over all pairs of primitive Gaussian-type functions ***

      na = 0
      DO ipgf = 1, npgfa

         nb = 0

         DO jpgf = 1, npgfb

            ss = 0.0_dp
            sc = 0.0_dp

!       *** Screening ***
            IF (rpgfa(ipgf) + rpgfb(jpgf) < dab) THEN
               nb = nb + ncoset(lb_max)
               CYCLE
            END IF

!       *** Calculate some prefactors ***

            zetp = 1.0_dp/(zeta(ipgf) + zetb(jpgf))

            f0 = (pi*zetp)**1.5_dp
            f1 = zetb(jpgf)*zetp
            f2 = 0.5_dp*zetp

            kdp = zetp*DOT_PRODUCT(kvec, zeta(ipgf)*rac + zetb(jpgf)*rbc)

!       *** Calculate the basic two-center cos/sin integral [s|cos/sin|s] ***

            s = f0*EXP(-zeta(ipgf)*f1*rab2)*EXP(-0.25_dp*k2*zetp)
            sc(1, 1) = s*COS(kdp)
            ss(1, 1) = s*SIN(kdp)

!       *** Recurrence steps: [s|O|s] -> [a|O|b] ***

            IF (la_max > 0) THEN

!         *** Vertical recurrence steps: [s|O|s] -> [a|O|s] ***

               rap(:) = f1*rab(:)

!         *** [p|O|s] = (Pi - Ai)*[s|O|s] +[s|dO|s]  (i = x,y,z) ***

               sc(2, 1) = rap(1)*sc(1, 1) - f2*kvec(1)*ss(1, 1)
               sc(3, 1) = rap(2)*sc(1, 1) - f2*kvec(2)*ss(1, 1)
               sc(4, 1) = rap(3)*sc(1, 1) - f2*kvec(3)*ss(1, 1)
               ss(2, 1) = rap(1)*ss(1, 1) + f2*kvec(1)*sc(1, 1)
               ss(3, 1) = rap(2)*ss(1, 1) + f2*kvec(2)*sc(1, 1)
               ss(4, 1) = rap(3)*ss(1, 1) + f2*kvec(3)*sc(1, 1)

!         *** [a|O|s] = (Pi - Ai)*[a-1i|O|s] + f2*Ni(a-1i)*[a-2i|s] ***
!         ***           + [a-1i|dO|s]                               ***

               DO la = 2, la_max

!           *** Increase the angular momentum component z of function a ***

                  sc(coset(0, 0, la), 1) = rap(3)*sc(coset(0, 0, la - 1), 1) + &
                                           f2*REAL(la - 1, dp)*sc(coset(0, 0, la - 2), 1) - &
                                           f2*kvec(3)*ss(coset(0, 0, la - 1), 1)
                  ss(coset(0, 0, la), 1) = rap(3)*ss(coset(0, 0, la - 1), 1) + &
                                           f2*REAL(la - 1, dp)*ss(coset(0, 0, la - 2), 1) + &
                                           f2*kvec(3)*sc(coset(0, 0, la - 1), 1)

!           *** Increase the angular momentum component y of function a ***

                  az = la - 1
                  sc(coset(0, 1, az), 1) = rap(2)*sc(coset(0, 0, az), 1) - &
                                           f2*kvec(2)*ss(coset(0, 0, az), 1)
                  ss(coset(0, 1, az), 1) = rap(2)*ss(coset(0, 0, az), 1) + &
                                           f2*kvec(2)*sc(coset(0, 0, az), 1)

                  DO ay = 2, la
                     az = la - ay
                     sc(coset(0, ay, az), 1) = rap(2)*sc(coset(0, ay - 1, az), 1) + &
                                               f2*REAL(ay - 1, dp)*sc(coset(0, ay - 2, az), 1) - &
                                               f2*kvec(2)*ss(coset(0, ay - 1, az), 1)
                     ss(coset(0, ay, az), 1) = rap(2)*ss(coset(0, ay - 1, az), 1) + &
                                               f2*REAL(ay - 1, dp)*ss(coset(0, ay - 2, az), 1) + &
                                               f2*kvec(2)*sc(coset(0, ay - 1, az), 1)
                  END DO

!           *** Increase the angular momentum component x of function a ***

                  DO ay = 0, la - 1
                     az = la - 1 - ay
                     sc(coset(1, ay, az), 1) = rap(1)*sc(coset(0, ay, az), 1) - &
                                               f2*kvec(1)*ss(coset(0, ay, az), 1)
                     ss(coset(1, ay, az), 1) = rap(1)*ss(coset(0, ay, az), 1) + &
                                               f2*kvec(1)*sc(coset(0, ay, az), 1)
                  END DO

                  DO ax = 2, la
                     f3 = f2*REAL(ax - 1, dp)
                     DO ay = 0, la - ax
                        az = la - ax - ay
                        sc(coset(ax, ay, az), 1) = rap(1)*sc(coset(ax - 1, ay, az), 1) + &
                                                   f3*sc(coset(ax - 2, ay, az), 1) - &
                                                   f2*kvec(1)*ss(coset(ax - 1, ay, az), 1)
                        ss(coset(ax, ay, az), 1) = rap(1)*ss(coset(ax - 1, ay, az), 1) + &
                                                   f3*ss(coset(ax - 2, ay, az), 1) + &
                                                   f2*kvec(1)*sc(coset(ax - 1, ay, az), 1)
                     END DO
                  END DO

               END DO

!         *** Recurrence steps: [a|O|s] -> [a|O|b] ***

               IF (lb_max > 0) THEN

                  DO j = 2, ncoset(lb_max)
                     DO i = 1, ncoset(la_max)
                        sc(i, j) = 0.0_dp
                        ss(i, j) = 0.0_dp
                     END DO
                  END DO

!           *** Horizontal recurrence steps ***

                  rbp(:) = rap(:) - rab(:)

!           *** [a|O|p] = [a+1i|O|s] - (Bi - Ai)*[a|O|s] ***

                  IF (lb_max == 1) THEN
                     la_start = la_min
                  ELSE
                     la_start = MAX(0, la_min - 1)
                  END IF

                  DO la = la_start, la_max - 1
                     DO ax = 0, la
                        DO ay = 0, la - ax
                           az = la - ax - ay
                           sc(coset(ax, ay, az), 2) = sc(coset(ax + 1, ay, az), 1) - &
                                                      rab(1)*sc(coset(ax, ay, az), 1)
                           sc(coset(ax, ay, az), 3) = sc(coset(ax, ay + 1, az), 1) - &
                                                      rab(2)*sc(coset(ax, ay, az), 1)
                           sc(coset(ax, ay, az), 4) = sc(coset(ax, ay, az + 1), 1) - &
                                                      rab(3)*sc(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 2) = ss(coset(ax + 1, ay, az), 1) - &
                                                      rab(1)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 3) = ss(coset(ax, ay + 1, az), 1) - &
                                                      rab(2)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 4) = ss(coset(ax, ay, az + 1), 1) - &
                                                      rab(3)*ss(coset(ax, ay, az), 1)
                        END DO
                     END DO
                  END DO

!           *** Vertical recurrence step ***

!           *** [a|O|p] = (Pi - Bi)*[a|O|s] + f2*Ni(a)*[a-1i|O|s] ***
!           ***           + [a|dO|s]                              ***

                  DO ax = 0, la_max
                     fx = f2*REAL(ax, dp)
                     DO ay = 0, la_max - ax
                        fy = f2*REAL(ay, dp)
                        az = la_max - ax - ay
                        fz = f2*REAL(az, dp)
                        IF (ax == 0) THEN
                           sc(coset(ax, ay, az), 2) = rbp(1)*sc(coset(ax, ay, az), 1) - &
                                                      f2*kvec(1)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 2) = rbp(1)*ss(coset(ax, ay, az), 1) + &
                                                      f2*kvec(1)*sc(coset(ax, ay, az), 1)
                        ELSE
                           sc(coset(ax, ay, az), 2) = rbp(1)*sc(coset(ax, ay, az), 1) + &
                                                      fx*sc(coset(ax - 1, ay, az), 1) - &
                                                      f2*kvec(1)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 2) = rbp(1)*ss(coset(ax, ay, az), 1) + &
                                                      fx*ss(coset(ax - 1, ay, az), 1) + &
                                                      f2*kvec(1)*sc(coset(ax, ay, az), 1)
                        END IF
                        IF (ay == 0) THEN
                           sc(coset(ax, ay, az), 3) = rbp(2)*sc(coset(ax, ay, az), 1) - &
                                                      f2*kvec(2)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 3) = rbp(2)*ss(coset(ax, ay, az), 1) + &
                                                      f2*kvec(2)*sc(coset(ax, ay, az), 1)
                        ELSE
                           sc(coset(ax, ay, az), 3) = rbp(2)*sc(coset(ax, ay, az), 1) + &
                                                      fy*sc(coset(ax, ay - 1, az), 1) - &
                                                      f2*kvec(2)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 3) = rbp(2)*ss(coset(ax, ay, az), 1) + &
                                                      fy*ss(coset(ax, ay - 1, az), 1) + &
                                                      f2*kvec(2)*sc(coset(ax, ay, az), 1)
                        END IF
                        IF (az == 0) THEN
                           sc(coset(ax, ay, az), 4) = rbp(3)*sc(coset(ax, ay, az), 1) - &
                                                      f2*kvec(3)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 4) = rbp(3)*ss(coset(ax, ay, az), 1) + &
                                                      f2*kvec(3)*sc(coset(ax, ay, az), 1)
                        ELSE
                           sc(coset(ax, ay, az), 4) = rbp(3)*sc(coset(ax, ay, az), 1) + &
                                                      fz*sc(coset(ax, ay, az - 1), 1) - &
                                                      f2*kvec(3)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 4) = rbp(3)*ss(coset(ax, ay, az), 1) + &
                                                      fz*ss(coset(ax, ay, az - 1), 1) + &
                                                      f2*kvec(3)*sc(coset(ax, ay, az), 1)
                        END IF
                     END DO
                  END DO

!           *** Recurrence steps: [a|O|p] -> [a|O|b] ***

                  DO lb = 2, lb_max

!             *** Horizontal recurrence steps ***

!             *** [a|O|b] = [a+1i|O|b-1i] - (Bi - Ai)*[a|O|b-1i] ***

                     IF (lb == lb_max) THEN
                        la_start = la_min
                     ELSE
                        la_start = MAX(0, la_min - 1)
                     END IF

                     DO la = la_start, la_max - 1
                        DO ax = 0, la
                           DO ay = 0, la - ax
                              az = la - ax - ay

!                   *** Shift of angular momentum component z from a to b ***

                              sc(coset(ax, ay, az), coset(0, 0, lb)) = &
                                 sc(coset(ax, ay, az + 1), coset(0, 0, lb - 1)) - &
                                 rab(3)*sc(coset(ax, ay, az), coset(0, 0, lb - 1))
                              ss(coset(ax, ay, az), coset(0, 0, lb)) = &
                                 ss(coset(ax, ay, az + 1), coset(0, 0, lb - 1)) - &
                                 rab(3)*ss(coset(ax, ay, az), coset(0, 0, lb - 1))

!                   *** Shift of angular momentum component y from a to b ***

                              DO by = 1, lb
                                 bz = lb - by
                                 sc(coset(ax, ay, az), coset(0, by, bz)) = &
                                    sc(coset(ax, ay + 1, az), coset(0, by - 1, bz)) - &
                                    rab(2)*sc(coset(ax, ay, az), coset(0, by - 1, bz))
                                 ss(coset(ax, ay, az), coset(0, by, bz)) = &
                                    ss(coset(ax, ay + 1, az), coset(0, by - 1, bz)) - &
                                    rab(2)*ss(coset(ax, ay, az), coset(0, by - 1, bz))
                              END DO

!                   *** Shift of angular momentum component x from a to b ***

                              DO bx = 1, lb
                                 DO by = 0, lb - bx
                                    bz = lb - bx - by
                                    sc(coset(ax, ay, az), coset(bx, by, bz)) = &
                                       sc(coset(ax + 1, ay, az), coset(bx - 1, by, bz)) - &
                                       rab(1)*sc(coset(ax, ay, az), coset(bx - 1, by, bz))
                                    ss(coset(ax, ay, az), coset(bx, by, bz)) = &
                                       ss(coset(ax + 1, ay, az), coset(bx - 1, by, bz)) - &
                                       rab(1)*ss(coset(ax, ay, az), coset(bx - 1, by, bz))
                                 END DO
                              END DO

                           END DO
                        END DO
                     END DO

!             *** Vertical recurrence step ***

!             *** [a|O|b] = (Pi - Bi)*[a|O|b-1i] + f2*Ni(a)*[a-1i|O|b-1i] + ***
!             ***           f2*Ni(b-1i)*[a|O|b-2i] + [a|dO|b-1i]            ***

                     DO ax = 0, la_max
                        fx = f2*REAL(ax, dp)
                        DO ay = 0, la_max - ax
                           fy = f2*REAL(ay, dp)
                           az = la_max - ax - ay
                           fz = f2*REAL(az, dp)

!                 *** Increase the angular momentum component z of function b ***

                           f3 = f2*REAL(lb - 1, dp)

                           IF (az == 0) THEN
                              sc(coset(ax, ay, az), coset(0, 0, lb)) = &
                                 rbp(3)*sc(coset(ax, ay, az), coset(0, 0, lb - 1)) + &
                                 f3*sc(coset(ax, ay, az), coset(0, 0, lb - 2)) - &
                                 f2*kvec(3)*ss(coset(ax, ay, az), coset(0, 0, lb - 1))
                              ss(coset(ax, ay, az), coset(0, 0, lb)) = &
                                 rbp(3)*ss(coset(ax, ay, az), coset(0, 0, lb - 1)) + &
                                 f3*ss(coset(ax, ay, az), coset(0, 0, lb - 2)) + &
                                 f2*kvec(3)*sc(coset(ax, ay, az), coset(0, 0, lb - 1))
                           ELSE
                              sc(coset(ax, ay, az), coset(0, 0, lb)) = &
                                 rbp(3)*sc(coset(ax, ay, az), coset(0, 0, lb - 1)) + &
                                 fz*sc(coset(ax, ay, az - 1), coset(0, 0, lb - 1)) + &
                                 f3*sc(coset(ax, ay, az), coset(0, 0, lb - 2)) - &
                                 f2*kvec(3)*ss(coset(ax, ay, az), coset(0, 0, lb - 1))
                              ss(coset(ax, ay, az), coset(0, 0, lb)) = &
                                 rbp(3)*ss(coset(ax, ay, az), coset(0, 0, lb - 1)) + &
                                 fz*ss(coset(ax, ay, az - 1), coset(0, 0, lb - 1)) + &
                                 f3*ss(coset(ax, ay, az), coset(0, 0, lb - 2)) + &
                                 f2*kvec(3)*sc(coset(ax, ay, az), coset(0, 0, lb - 1))
                           END IF

!                 *** Increase the angular momentum component y of function b ***

                           IF (ay == 0) THEN
                              bz = lb - 1
                              sc(coset(ax, ay, az), coset(0, 1, bz)) = &
                                 rbp(2)*sc(coset(ax, ay, az), coset(0, 0, bz)) - &
                                 f2*kvec(2)*ss(coset(ax, ay, az), coset(0, 0, bz))
                              ss(coset(ax, ay, az), coset(0, 1, bz)) = &
                                 rbp(2)*ss(coset(ax, ay, az), coset(0, 0, bz)) + &
                                 f2*kvec(2)*sc(coset(ax, ay, az), coset(0, 0, bz))
                              DO by = 2, lb
                                 bz = lb - by
                                 f3 = f2*REAL(by - 1, dp)
                                 sc(coset(ax, ay, az), coset(0, by, bz)) = &
                                    rbp(2)*sc(coset(ax, ay, az), coset(0, by - 1, bz)) + &
                                    f3*sc(coset(ax, ay, az), coset(0, by - 2, bz)) - &
                                    f2*kvec(2)*ss(coset(ax, ay, az), coset(0, by - 1, bz))
                                 ss(coset(ax, ay, az), coset(0, by, bz)) = &
                                    rbp(2)*ss(coset(ax, ay, az), coset(0, by - 1, bz)) + &
                                    f3*ss(coset(ax, ay, az), coset(0, by - 2, bz)) + &
                                    f2*kvec(2)*sc(coset(ax, ay, az), coset(0, by - 1, bz))
                              END DO
                           ELSE
                              bz = lb - 1
                              sc(coset(ax, ay, az), coset(0, 1, bz)) = &
                                 rbp(2)*sc(coset(ax, ay, az), coset(0, 0, bz)) + &
                                 fy*sc(coset(ax, ay - 1, az), coset(0, 0, bz)) - &
                                 f2*kvec(2)*ss(coset(ax, ay, az), coset(0, 0, bz))
                              ss(coset(ax, ay, az), coset(0, 1, bz)) = &
                                 rbp(2)*ss(coset(ax, ay, az), coset(0, 0, bz)) + &
                                 fy*ss(coset(ax, ay - 1, az), coset(0, 0, bz)) + &
                                 f2*kvec(2)*sc(coset(ax, ay, az), coset(0, 0, bz))
                              DO by = 2, lb
                                 bz = lb - by
                                 f3 = f2*REAL(by - 1, dp)
                                 sc(coset(ax, ay, az), coset(0, by, bz)) = &
                                    rbp(2)*sc(coset(ax, ay, az), coset(0, by - 1, bz)) + &
                                    fy*sc(coset(ax, ay - 1, az), coset(0, by - 1, bz)) + &
                                    f3*sc(coset(ax, ay, az), coset(0, by - 2, bz)) - &
                                    f2*kvec(2)*ss(coset(ax, ay, az), coset(0, by - 1, bz))
                                 ss(coset(ax, ay, az), coset(0, by, bz)) = &
                                    rbp(2)*ss(coset(ax, ay, az), coset(0, by - 1, bz)) + &
                                    fy*ss(coset(ax, ay - 1, az), coset(0, by - 1, bz)) + &
                                    f3*ss(coset(ax, ay, az), coset(0, by - 2, bz)) + &
                                    f2*kvec(2)*sc(coset(ax, ay, az), coset(0, by - 1, bz))
                              END DO
                           END IF

!                 *** Increase the angular momentum component x of function b ***

                           IF (ax == 0) THEN
                              DO by = 0, lb - 1
                                 bz = lb - 1 - by
                                 sc(coset(ax, ay, az), coset(1, by, bz)) = &
                                    rbp(1)*sc(coset(ax, ay, az), coset(0, by, bz)) - &
                                    f2*kvec(1)*ss(coset(ax, ay, az), coset(0, by, bz))
                                 ss(coset(ax, ay, az), coset(1, by, bz)) = &
                                    rbp(1)*ss(coset(ax, ay, az), coset(0, by, bz)) + &
                                    f2*kvec(1)*sc(coset(ax, ay, az), coset(0, by, bz))
                              END DO
                              DO bx = 2, lb
                                 f3 = f2*REAL(bx - 1, dp)
                                 DO by = 0, lb - bx
                                    bz = lb - bx - by
                                    sc(coset(ax, ay, az), coset(bx, by, bz)) = &
                                       rbp(1)*sc(coset(ax, ay, az), &
                                                 coset(bx - 1, by, bz)) + &
                                       f3*sc(coset(ax, ay, az), coset(bx - 2, by, bz)) - &
                                       f2*kvec(1)*ss(coset(ax, ay, az), coset(bx - 1, by, bz))
                                    ss(coset(ax, ay, az), coset(bx, by, bz)) = &
                                       rbp(1)*ss(coset(ax, ay, az), &
                                                 coset(bx - 1, by, bz)) + &
                                       f3*ss(coset(ax, ay, az), coset(bx - 2, by, bz)) + &
                                       f2*kvec(1)*sc(coset(ax, ay, az), coset(bx - 1, by, bz))
                                 END DO
                              END DO
                           ELSE
                              DO by = 0, lb - 1
                                 bz = lb - 1 - by
                                 sc(coset(ax, ay, az), coset(1, by, bz)) = &
                                    rbp(1)*sc(coset(ax, ay, az), coset(0, by, bz)) + &
                                    fx*sc(coset(ax - 1, ay, az), coset(0, by, bz)) - &
                                    f2*kvec(1)*ss(coset(ax, ay, az), coset(0, by, bz))
                                 ss(coset(ax, ay, az), coset(1, by, bz)) = &
                                    rbp(1)*ss(coset(ax, ay, az), coset(0, by, bz)) + &
                                    fx*ss(coset(ax - 1, ay, az), coset(0, by, bz)) + &
                                    f2*kvec(1)*sc(coset(ax, ay, az), coset(0, by, bz))
                              END DO
                              DO bx = 2, lb
                                 f3 = f2*REAL(bx - 1, dp)
                                 DO by = 0, lb - bx
                                    bz = lb - bx - by
                                    sc(coset(ax, ay, az), coset(bx, by, bz)) = &
                                       rbp(1)*sc(coset(ax, ay, az), &
                                                 coset(bx - 1, by, bz)) + &
                                       fx*sc(coset(ax - 1, ay, az), coset(bx - 1, by, bz)) + &
                                       f3*sc(coset(ax, ay, az), coset(bx - 2, by, bz)) - &
                                       f2*kvec(1)*ss(coset(ax, ay, az), coset(bx - 1, by, bz))
                                    ss(coset(ax, ay, az), coset(bx, by, bz)) = &
                                       rbp(1)*ss(coset(ax, ay, az), &
                                                 coset(bx - 1, by, bz)) + &
                                       fx*ss(coset(ax - 1, ay, az), coset(bx - 1, by, bz)) + &
                                       f3*ss(coset(ax, ay, az), coset(bx - 2, by, bz)) + &
                                       f2*kvec(1)*sc(coset(ax, ay, az), coset(bx - 1, by, bz))
                                 END DO
                              END DO
                           END IF

                        END DO
                     END DO

                  END DO

               END IF

            ELSE

               IF (lb_max > 0) THEN

!           *** Vertical recurrence steps: [s|O|s] -> [s|O|b] ***

                  rbp(:) = (f1 - 1.0_dp)*rab(:)

!           *** [s|O|p] = (Pi - Bi)*[s|O|s] + [s|dO|s] ***

                  sc(1, 2) = rbp(1)*sc(1, 1) - f2*kvec(1)*ss(1, 1)
                  sc(1, 3) = rbp(2)*sc(1, 1) - f2*kvec(2)*ss(1, 1)
                  sc(1, 4) = rbp(3)*sc(1, 1) - f2*kvec(3)*ss(1, 1)
                  ss(1, 2) = rbp(1)*ss(1, 1) + f2*kvec(1)*sc(1, 1)
                  ss(1, 3) = rbp(2)*ss(1, 1) + f2*kvec(2)*sc(1, 1)
                  ss(1, 4) = rbp(3)*ss(1, 1) + f2*kvec(3)*sc(1, 1)

!           *** [s|O|b] = (Pi - Bi)*[s|O|b-1i] + f2*Ni(b-1i)*[s|O|b-2i] ***
!           ***           + [s|dO|b-1i]                                 ***

                  DO lb = 2, lb_max

!             *** Increase the angular momentum component z of function b ***

                     sc(1, coset(0, 0, lb)) = rbp(3)*sc(1, coset(0, 0, lb - 1)) + &
                                              f2*REAL(lb - 1, dp)*sc(1, coset(0, 0, lb - 2)) - &
                                              f2*kvec(3)*ss(1, coset(0, 0, lb - 1))
                     ss(1, coset(0, 0, lb)) = rbp(3)*ss(1, coset(0, 0, lb - 1)) + &
                                              f2*REAL(lb - 1, dp)*ss(1, coset(0, 0, lb - 2)) + &
                                              f2*kvec(3)*sc(1, coset(0, 0, lb - 1))

!             *** Increase the angular momentum component y of function b ***

                     bz = lb - 1
                     sc(1, coset(0, 1, bz)) = rbp(2)*sc(1, coset(0, 0, bz)) - &
                                              f2*kvec(2)*ss(1, coset(0, 0, bz))
                     ss(1, coset(0, 1, bz)) = rbp(2)*ss(1, coset(0, 0, bz)) + &
                                              f2*kvec(2)*sc(1, coset(0, 0, bz))

                     DO by = 2, lb
                        bz = lb - by
                        sc(1, coset(0, by, bz)) = rbp(2)*sc(1, coset(0, by - 1, bz)) + &
                                                  f2*REAL(by - 1, dp)*sc(1, coset(0, by - 2, bz)) - &
                                                  f2*kvec(2)*ss(1, coset(0, by - 1, bz))
                        ss(1, coset(0, by, bz)) = rbp(2)*ss(1, coset(0, by - 1, bz)) + &
                                                  f2*REAL(by - 1, dp)*ss(1, coset(0, by - 2, bz)) + &
                                                  f2*kvec(2)*sc(1, coset(0, by - 1, bz))
                     END DO

!             *** Increase the angular momentum component x of function b ***

                     DO by = 0, lb - 1
                        bz = lb - 1 - by
                        sc(1, coset(1, by, bz)) = rbp(1)*sc(1, coset(0, by, bz)) - &
                                                  f2*kvec(1)*ss(1, coset(0, by, bz))
                        ss(1, coset(1, by, bz)) = rbp(1)*ss(1, coset(0, by, bz)) + &
                                                  f2*kvec(1)*sc(1, coset(0, by, bz))
                     END DO

                     DO bx = 2, lb
                        f3 = f2*REAL(bx - 1, dp)
                        DO by = 0, lb - bx
                           bz = lb - bx - by
                           sc(1, coset(bx, by, bz)) = rbp(1)*sc(1, coset(bx - 1, by, bz)) + &
                                                      f3*sc(1, coset(bx - 2, by, bz)) - &
                                                      f2*kvec(1)*ss(1, coset(bx - 1, by, bz))
                           ss(1, coset(bx, by, bz)) = rbp(1)*ss(1, coset(bx - 1, by, bz)) + &
                                                      f3*ss(1, coset(bx - 2, by, bz)) + &
                                                      f2*kvec(1)*sc(1, coset(bx - 1, by, bz))
                        END DO
                     END DO

                  END DO

               END IF

            END IF

            DO j = ncoset(lb_min - 1) + 1, ncoset(lb_max)
               DO i = ncoset(la_min_set - 1) + 1, ncoset(la_max_set)
                  cosab(na + i, nb + j) = sc(i, j)
                  sinab(na + i, nb + j) = ss(i, j)
               END DO
            END DO

            IF (PRESENT(dcosab)) THEN
               DO k = 1, 3
                  CALL cartesian_ladder(sc, dcosab(na + 1:na + ncoset(la_max_set), nb + 1:nb + ncoset(lb_max), k), &
                                        derivative_map(:, :, k), .TRUE., zeta(ipgf), 0.0_dp, .FALSE.)
                  CALL cartesian_ladder(ss, dsinab(na + 1:na + ncoset(la_max_set), nb + 1:nb + ncoset(lb_max), k), &
                                        derivative_map(:, :, k), .TRUE., zeta(ipgf), 0.0_dp, .FALSE.)
               END DO
            END IF

            nb = nb + ncoset(lb_max)

         END DO

         na = na + ncoset(la_max_set)

      END DO

   END SUBROUTINE cossin

! **************************************************************************************************
!> \brief ...
!> \param la_max ...
!> \param npgfa ...
!> \param zeta ...
!> \param rpgfa ...
!> \param la_min ...
!> \param lb_max ...
!> \param npgfb ...
!> \param zetb ...
!> \param rpgfb ...
!> \param lc_max ...
!> \param rac ...
!> \param rbc ...
!> \param mab ...
! **************************************************************************************************
   SUBROUTINE moment(la_max, npgfa, zeta, rpgfa, la_min, &
                     lb_max, npgfb, zetb, rpgfb, lc_max, rac, rbc, mab)
      INTEGER, INTENT(IN)                                :: la_max, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta, rpgfa
      INTEGER, INTENT(IN)                                :: la_min, lb_max, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb, rpgfb
      INTEGER, INTENT(IN)                                :: lc_max
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rac, rbc
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: mab

      INTEGER                                            :: ax, ay, az, bx, by, bz, coa, cob, i, ia, &
                                                            ib, imom, ipgf, jpgf, ldrr, lx, ly, &
                                                            lz, ma, mb, powers(0:lc_max)
      REAL(KIND=dp)                                      :: a, b, dab, f0, rab(3), rab2, rap(3), &
                                                            rbp(3), zet
      REAL(KIND=dp), ALLOCATABLE                         :: polynomial(:, :, :, :), rr(:, :, :)

      MARK_USED(la_min)
      rab = rbc - rac
      dab = NORM2(rab)
      rab2 = dab**2
      ldrr = MAX(la_max + lc_max, lb_max) + 1
      ALLOCATE (rr(0:ldrr - 1, 0:ldrr - 1, 3), polynomial(0:la_max, 0:lb_max, 0:lc_max, 3))
      powers = [(i, i=0, lc_max)]
      mab(:npgfa*ncoset(la_max), :npgfb*ncoset(lb_max), :ncoset(lc_max) - 1) = 0.0_dp
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
            CALL os_rr_ovlp(rap, la_max + lc_max, rbp, lb_max, zet, ldrr, rr)
            CALL cartesian_polynomial_1d(rr, polynomial, la_max, lb_max, powers, rac)

            DO cob = 1, ncoset(lb_max)
               bx = indco(1, cob)
               by = indco(2, cob)
               bz = indco(3, cob)
               ib = mb + cob
               DO coa = 1, ncoset(la_max)
                  ax = indco(1, coa)
                  ay = indco(2, coa)
                  az = indco(3, coa)
                  ia = ma + coa
                  ! Preserve the coset ordering and omission of the zeroth moment.
                  DO imom = 1, ncoset(lc_max) - 1
                     lx = indco(1, imom + 1)
                     ly = indco(2, imom + 1)
                     lz = indco(3, imom + 1)
                     mab(ia, ib, imom) = f0*polynomial(ax, bx, lx, 1)* &
                                         polynomial(ay, by, ly, 2)*polynomial(az, bz, lz, 3)
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE moment

! **************************************************************************************************
!> \brief This returns the derivative of the moment integrals [a|\mu|b].
!>       By default, it differentiates the primitive on the right:
!>       [a|\mu|d/dR_bi] =  2*zetb*[a|\mu|b+1i] - Ni(b)[a|\mu|b-1i]
!>       A weighted derivative combines the left and right primitive derivatives
!>       using deltaR for the corresponding atom centers.
!>       order indicates the max order of the moment operator to be calculated
!>       1: dipole
!>       2: quadrupole
!>       ...
!> \param la_max ...
!> \param npgfa ...
!> \param zeta ...
!> \param rpgfa ...
!> \param la_min ...
!> \param lb_max ...
!> \param npgfb ...
!> \param zetb ...
!> \param rpgfb ...
!> \param lb_min ...
!> \param order ...
!> \param rac ...
!> \param rbc ...
!> \param difmab ...
!> \param mab_ext optional precomputed moments with zero screened primitive blocks, including raised shells
!> \param deltaR optional weights for the left and right primitive derivatives
!> \param lambda optional atom selector for the factor (iatom == lambda) - (jatom == lambda)
!> \param iatom atom associated with the left basis function
!> \param jatom atom associated with the right basis function
!> \note
! **************************************************************************************************
   SUBROUTINE diff_momop(la_max, npgfa, zeta, rpgfa, la_min, &
                         lb_max, npgfb, zetb, rpgfb, lb_min, &
                         order, rac, rbc, difmab, mab_ext, deltaR, lambda, iatom, jatom)

      INTEGER, INTENT(IN)                                :: la_max, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta, rpgfa
      INTEGER, INTENT(IN)                                :: la_min, lb_max, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb, rpgfb
      INTEGER, INTENT(IN)                                :: lb_min, order
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rac, rbc
      REAL(KIND=dp), DIMENSION(:, :, :, :), INTENT(OUT)  :: difmab
      REAL(KIND=dp), DIMENSION(:, :, :), OPTIONAL, &
         POINTER                                         :: mab_ext
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: deltaR
      INTEGER, INTENT(IN), OPTIONAL                      :: lambda, iatom, jatom

      INTEGER                                            :: ider, imom, lda, lda_min, ldb, ldb_min
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: difmab_tmp
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: mab

      lda_min = MAX(0, la_min - 1)
      ldb_min = MAX(0, lb_min - 1)
      lda = ncoset(la_max)*npgfa
      ldb = ncoset(lb_max)*npgfb
      ALLOCATE (difmab_tmp(lda, ldb, 3))

      IF (PRESENT(mab_ext)) THEN
         mab => mab_ext
      ELSE
         ALLOCATE (mab(npgfa*ncoset(la_max + 1), npgfb*ncoset(lb_max + 1), &
                       ncoset(order) - 1))
         mab = 0.0_dp
!     *** Calculate the primitive overlap integrals ***
         CALL moment(la_max + 1, npgfa, zeta, rpgfa, lda_min, &
                     lb_max + 1, npgfb, zetb, rpgfb, &
                     order, rac, rbc, mab)

      END IF
      DO imom = 1, ncoset(order) - 1
         CALL apply_cartesian_operator(la_max, 0, npgfa, lb_max, 0, npgfb, mab(:, :, imom), &
                                       difmab_tmp, right=.TRUE., derivative=.TRUE., exponents=zetb)
         IF (PRESENT(deltaR)) THEN
            CPASSERT(ASSOCIATED(deltaR))
            CPASSERT(PRESENT(iatom) .AND. PRESENT(jatom))
            DO ider = 1, 3
               difmab(1:lda, 1:ldb, imom, ider) = &
                  difmab_tmp(1:lda, 1:ldb, ider)*deltaR(ider, jatom)
            END DO

            CALL apply_cartesian_operator(la_max, 0, npgfa, lb_max, 0, npgfb, mab(:, :, imom), &
                                          difmab_tmp, right=.FALSE., derivative=.TRUE., exponents=zeta)
            DO ider = 1, 3
               difmab(1:lda, 1:ldb, imom, ider) = difmab(1:lda, 1:ldb, imom, ider) &
                                                  + difmab_tmp(1:lda, 1:ldb, ider)*deltaR(ider, iatom)
            END DO
         ELSE
            difmab(1:lda, 1:ldb, imom, :) = difmab_tmp(1:lda, 1:ldb, :)
         END IF
      END DO

      IF (PRESENT(lambda)) THEN
         CPASSERT(.NOT. PRESENT(deltaR))
         CPASSERT(PRESENT(iatom) .AND. PRESENT(jatom))
         IF (iatom == lambda .AND. jatom == lambda) THEN
            difmab = 0.0_dp
         ELSE IF (iatom == lambda) THEN
            ! The right-hand derivative is selected with a positive sign.
         ELSE IF (jatom == lambda) THEN
            difmab = -difmab
         ELSE
            difmab = 0.0_dp
         END IF
      END IF

      IF (PRESENT(mab_ext)) THEN
         NULLIFY (mab)
      ELSE
         DEALLOCATE (mab)
      END IF
      DEALLOCATE (difmab_tmp)

   END SUBROUTINE diff_momop

! **************************************************************************************************
!> \brief This returns the derivative of the dipole integrals [a|x|b], with respect
!>       to the position of the primitive on the left and right, i.e.
!>       [da/dR_ai|\mu|b] =  2*zeta*[a+1i|\mu|b] - Ni(a)[a-1i|\mu|b]
!> \param la_max ...
!> \param npgfa ...
!> \param zeta ...
!> \param rpgfa ...
!> \param la_min ...
!> \param lb_max ...
!> \param npgfb ...
!> \param zetb ...
!> \param rpgfb ...
!> \param lb_min ...
!> \param order ...
!> \param rac ...
!> \param rbc ...
!> \param pab ...
!> \param forcea ...
!> \param forceb ...
!> \note
! **************************************************************************************************
   SUBROUTINE dipole_force(la_max, npgfa, zeta, rpgfa, la_min, &
                           lb_max, npgfb, zetb, rpgfb, lb_min, &
                           order, rac, rbc, pab, forcea, forceb)

      INTEGER, INTENT(IN)                                :: la_max, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta, rpgfa
      INTEGER, INTENT(IN)                                :: la_min, lb_max, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb, rpgfb
      INTEGER, INTENT(IN)                                :: lb_min, order
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rac, rbc
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: pab
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: forcea, forceb

      INTEGER                                            :: i, imom, ipgf, j, jpgf, lda, lda_min, &
                                                            ldb, ldb_min, na, nb
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: difmab, mab

      CPASSERT(order == 1)
      MARK_USED(order)

      lda_min = MAX(0, la_min - 1)
      ldb_min = MAX(0, lb_min - 1)
      lda = ncoset(la_max)*npgfa
      ldb = ncoset(lb_max)*npgfb
      ALLOCATE (difmab(lda, ldb, 3))
      ALLOCATE (mab(npgfa*ncoset(la_max + 1), npgfb*ncoset(lb_max + 1), 3))
      mab = 0.0_dp
      CALL moment(la_max + 1, npgfa, zeta, rpgfa, lda_min, &
                  lb_max + 1, npgfb, zetb, rpgfb, 1, rac, rbc, mab)

      DO imom = 1, 3
         difmab = 0.0_dp
         CALL apply_cartesian_operator(la_max, 0, npgfa, lb_max, 0, npgfb, mab(:, :, imom), &
                                       difmab, right=.TRUE., derivative=.TRUE., exponents=zetb)
         na = 0
         DO ipgf = 1, npgfa
            nb = 0
            DO jpgf = 1, npgfb
               DO j = nb + ncoset(lb_min - 1) + 1, nb + ncoset(lb_max)
                  DO i = na + ncoset(la_min - 1) + 1, na + ncoset(la_max)
                     forceb(imom, 1) = forceb(imom, 1) + pab(i, j)*difmab(i, j, 1)
                     forceb(imom, 2) = forceb(imom, 2) + pab(i, j)*difmab(i, j, 2)
                     forceb(imom, 3) = forceb(imom, 3) + pab(i, j)*difmab(i, j, 3)
                  END DO
               END DO
               nb = nb + ncoset(lb_max)
            END DO
            na = na + ncoset(la_max)
         END DO

         difmab = 0.0_dp
         CALL apply_cartesian_operator(la_max, 0, npgfa, lb_max, 0, npgfb, mab(:, :, imom), &
                                       difmab, right=.FALSE., derivative=.TRUE., exponents=zeta)
         na = 0
         DO ipgf = 1, npgfa
            nb = 0
            DO jpgf = 1, npgfb
               DO j = nb + ncoset(lb_min - 1) + 1, nb + ncoset(lb_max)
                  DO i = na + ncoset(la_min - 1) + 1, na + ncoset(la_max)
                     forcea(imom, 1) = forcea(imom, 1) + pab(i, j)*difmab(i, j, 1)
                     forcea(imom, 2) = forcea(imom, 2) + pab(i, j)*difmab(i, j, 2)
                     forcea(imom, 3) = forcea(imom, 3) + pab(i, j)*difmab(i, j, 3)
                  END DO
               END DO
               nb = nb + ncoset(lb_max)
            END DO
            na = na + ncoset(la_max)
         END DO
      END DO

      DEALLOCATE (mab, difmab)

   END SUBROUTINE dipole_force

END MODULE ai_moments
