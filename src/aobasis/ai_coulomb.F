!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of Coulomb integrals over Cartesian Gaussian-type functions
!>      (electron repulsion integrals, ERIs).
!> \par Literature
!>      S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!> \par History
!>      none
!> \par Parameters
!>       - ax,ay,az    : Angular momentum index numbers of orbital a.
!>       - bx,by,bz    : Angular momentum index numbers of orbital b.
!>       - cx,cy,cz    : Angular momentum index numbers of orbital c.
!>       - coset       : Cartesian orbital set pointer.
!>       - dab         : Distance between the atomic centers a and b.
!>       - dac         : Distance between the atomic centers a and c.
!>       - dbc         : Distance between the atomic centers b and c.
!>       - gccc        : Prefactor of the primitive Gaussian function c.
!>       - l{a,b,c}    : Angular momentum quantum number of shell a, b or c.
!>       - l{a,b,c}_max: Maximum angular momentum quantum number of shell a, b or c.
!>       - l{a,b,c}_min: Minimum angular momentum quantum number of shell a, b or c.
!>       - ncoset      : Number of orbitals in a Cartesian orbital set.
!>       - npgf{a,b}   : Degree of contraction of shell a or b.
!>       - rab         : Distance vector between the atomic centers a and b.
!>       - rab2        : Square of the distance between the atomic centers a and b.
!>       - rac         : Distance vector between the atomic centers a and c.
!>       - rac2        : Square of the distance between the atomic centers a and c.
!>       - rbc         : Distance vector between the atomic centers b and c.
!>       - rbc2        : Square of the distance between the atomic centers b and c.
!>       - rpgf{a,b,c} : Radius of the primitive Gaussian-type function a, b or c.
!>       - zet{a,b,c}  : Exponents of the Gaussian-type functions a, b or c.
!>       - zetp        : Reciprocal of the sum of the exponents of orbital a and b.
!>       - zetw        : Reciprocal of the sum of the exponents of orbital a, b and c.
!> \author Matthias Krack (22.08.2000)
! **************************************************************************************************
MODULE ai_coulomb

   USE ai_operators_r12,                ONLY: operator2_recurrence
   USE ai_os_rr,                        ONLY: os_rr_coulomb3
   USE gamma,                           ONLY: fgamma => fgamma_0
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
   USE orbital_pointers,                ONLY: ncoset
#include "../base/base_uses.f90"

   IMPLICIT NONE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_coulomb'
   PRIVATE

   ! *** Public subroutines ***

   PUBLIC :: coulomb2, coulomb3

CONTAINS

! **************************************************************************************************
!> \brief   Calculation of the primitive two-center Coulomb integrals over
!>          Cartesian Gaussian-type functions.
!> \param la_max ...
!> \param npgfa ...
!> \param zeta ...
!> \param rpgfa ...
!> \param la_min ...
!> \param lc_max ...
!> \param npgfc ...
!> \param zetc ...
!> \param rpgfc ...
!> \param lc_min ...
!> \param rac ...
!> \param rac2 ...
!> \param vac ...
!> \param v ...
!> \param f ...
!> \param screening optional primitive-pair screening switch
!> \date    05.12.2000
!> \author  Matthias Krack
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE coulomb2(la_max, npgfa, zeta, rpgfa, la_min, lc_max, npgfc, zetc, rpgfc, lc_min, &
                       rac, rac2, vac, v, f, screening)
      INTEGER, INTENT(IN)                                :: la_max, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta, rpgfa
      INTEGER, INTENT(IN)                                :: la_min, lc_max, npgfc
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetc, rpgfc
      INTEGER, INTENT(IN)                                :: lc_min
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rac
      REAL(KIND=dp), INTENT(IN)                          :: rac2
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: vac
      REAL(KIND=dp), DIMENSION(:, :, :)                  :: v
      REAL(KIND=dp), DIMENSION(0:)                       :: f
      LOGICAL, INTENT(IN), OPTIONAL                      :: screening

      INTEGER                                            :: i, ipgf, j, jpgf, n, na, nap, nc, ncp, &
                                                            nmax
      LOGICAL                                            :: do_screening
      REAL(KIND=dp)                                      :: dac, f0, rho, t, zetp, zetq, zetw

      do_screening = .TRUE.
      IF (PRESENT(screening)) do_screening = screening

      v = 0.0_dp

      nmax = la_max + lc_max + 1

      dac = SQRT(rac2)

      na = 0
      nap = 0
      DO ipgf = 1, npgfa

         nc = 0
         ncp = 0

         DO jpgf = 1, npgfc

            IF (do_screening .AND. rpgfa(ipgf) + rpgfc(jpgf) < dac) THEN
               DO j = nc + ncoset(lc_min - 1) + 1, nc + ncoset(lc_max)
                  DO i = na + ncoset(la_min - 1) + 1, na + ncoset(la_max)
                     vac(i, j) = 0.0_dp
                  END DO
               END DO
               nc = nc + ncoset(lc_max)
               CYCLE
            END IF

            zetp = 1.0_dp/zeta(ipgf)
            zetq = 1.0_dp/zetc(jpgf)
            zetw = 1.0_dp/(zeta(ipgf) + zetc(jpgf))

            rho = zeta(ipgf)*zetc(jpgf)*zetw

            f0 = 2.0_dp*SQRT(pi**5*zetw)*zetp*zetq

            t = rho*rac2

            CALL fgamma(nmax - 1, t, f)

            DO n = 1, nmax
               v(1, 1, n) = f0*f(n - 1)
            END DO

            CALL operator2_recurrence(la_max, la_min, lc_max, lc_min, zeta(ipgf), zetc(jpgf), &
                                      zetp, zetq, zetw, rho, rac, vac, v, na, nc, nap, ncp)

         END DO

         na = na + ncoset(la_max)
         nap = nap + ncoset(la_max)

      END DO

   END SUBROUTINE coulomb2
! **************************************************************************************************
!> \brief   Calculation of the primitive three-center Coulomb integrals over
!>          Cartesian Gaussian-type functions (electron repulsion integrals,
!>          ERIs).
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
!> \param lc_max ...
!> \param zetc ...
!> \param rpgfc ...
!> \param lc_min ...
!> \param gccc ...
!> \param rab ...
!> \param rab2 ...
!> \param rac ...
!> \param rac2 ...
!> \param rbc2 ...
!> \param vabc ...
!> \param int_abc ...
!> \param v ...
!> \param f ...
!> \param maxder ...
!> \param vabc_plus ...
!> \date    06.11.2000
!> \author  Matthias Krack
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE coulomb3(la_max, npgfa, zeta, rpgfa, la_min, lb_max, npgfb, zetb, rpgfb, lb_min, &
                       lc_max, zetc, rpgfc, lc_min, gccc, rab, rab2, rac, rac2, rbc2, vabc, int_abc, &
                       v, f, maxder, vabc_plus)
      INTEGER, INTENT(IN)                                :: la_max, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta, rpgfa
      INTEGER, INTENT(IN)                                :: la_min, lb_max, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb, rpgfb
      INTEGER, INTENT(IN)                                :: lb_min, lc_max
      REAL(KIND=dp), INTENT(IN)                          :: zetc, rpgfc
      INTEGER, INTENT(IN)                                :: lc_min
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: gccc
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), INTENT(IN)                          :: rab2
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rac
      REAL(KIND=dp), INTENT(IN)                          :: rac2, rbc2
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: vabc
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)     :: int_abc
      REAL(KIND=dp), DIMENSION(:, :, :, :)               :: v
      REAL(KIND=dp), DIMENSION(0:)                       :: f
      INTEGER, INTENT(IN), OPTIONAL                      :: maxder
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL           :: vabc_plus

      INTEGER                                            :: i, ipgf, j, jpgf, k, kk, maxder_local, &
                                                            n, na, nap, nb, nmax
      REAL(dp)                                           :: dab, dac, dbc, f0, f1, f4, rap(3), &
                                                            rbp(3), rcp(3), t, zetp, zetq, zetw

      v = 0.0_dp

      maxder_local = 0
      IF (PRESENT(maxder)) THEN
         maxder_local = maxder
      END IF

      nmax = la_max + lb_max + lc_max + 1

      ! *** Calculate the distances of the centers a, b and c ***

      dab = SQRT(rab2)
      dac = SQRT(rac2)
      dbc = SQRT(rbc2)

      ! *** Initialize integrals array
      int_abc = 0.0_dp

      ! *** Loop over all pairs of primitive Gaussian-type functions ***

      na = 0
      nap = 0

      DO ipgf = 1, npgfa

         ! *** Screening ***
         IF (rpgfa(ipgf) + rpgfc < dac) THEN
            na = na + ncoset(la_max - maxder_local)
            nap = nap + ncoset(la_max)
            CYCLE
         END IF

         nb = 0

         DO jpgf = 1, npgfb

            ! *** Screening ***
            IF ( &
               (rpgfb(jpgf) + rpgfc < dbc) .OR. &
               (rpgfa(ipgf) + rpgfb(jpgf) < dab)) THEN
               nb = nb + ncoset(lb_max)
               CYCLE
            END IF

            ! *** Calculate some prefactors ***

            zetp = 1.0_dp/(zeta(ipgf) + zetb(jpgf))
            zetq = 1.0_dp/zetc
            zetw = 1.0_dp/(zeta(ipgf) + zetb(jpgf) + zetc)

            f0 = 2.0_dp*SQRT(pi**5*zetw)*zetp*zetq
            f1 = zetb(jpgf)*zetp
            f4 = -zetc*zetw

            f0 = f0*EXP(-zeta(ipgf)*f1*rab2)

            rap(:) = f1*rab(:)
            rcp(:) = rap(:) - rac(:)

            ! *** Calculate the incomplete Gamma function ***

            t = -f4*(rcp(1)*rcp(1) + rcp(2)*rcp(2) + rcp(3)*rcp(3))/zetp

            CALL fgamma(nmax - 1, t, f)

            ! *** Calculate the basic three-center Coulomb integrals [ss||s]{n} ***

            DO n = 1, nmax
               v(1, 1, 1, n) = f0*f(n - 1)
            END DO

            rbp = rap - rab
            CALL os_rr_coulomb3(rap, la_max, rbp, lb_max, rcp, lc_max, zetp, zetq, zetw, v)

            ! *** Add the contribution of the current pair ***
            ! *** of primitive Gaussian-type functions     ***

            DO k = ncoset(lc_min - 1) + 1, ncoset(lc_max)
               kk = k - ncoset(lc_min - 1)
               DO j = ncoset(lb_min - 1) + 1, ncoset(lb_max)
                  DO i = ncoset(la_min - 1) + 1, ncoset(la_max - maxder_local)
                     vabc(na + i, nb + j) = vabc(na + i, nb + j) + gccc(kk)*v(i, j, k, 1)
                     int_abc(na + i, nb + j, kk) = v(i, j, k, 1)
                  END DO
               END DO
            END DO

            IF (PRESENT(maxder)) THEN
               DO k = ncoset(lc_min - 1) + 1, ncoset(lc_max)
                  kk = k - ncoset(lc_min - 1)
                  DO j = 1, ncoset(lb_max)
                     DO i = 1, ncoset(la_max)
                        vabc_plus(nap + i, nb + j) = vabc_plus(nap + i, nb + j) + gccc(kk)*v(i, j, k, 1)
                     END DO
                  END DO
               END DO
            END IF

            nb = nb + ncoset(lb_max)

         END DO

         na = na + ncoset(la_max - maxder_local)
         nap = nap + ncoset(la_max)

      END DO

   END SUBROUTINE coulomb3

END MODULE ai_coulomb
