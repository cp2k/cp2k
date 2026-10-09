!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of integrals over Cartesian Gaussian-type functions for [a|(r-Ra)^(2m)|b]
!>        Ra is the position of center a
!> \par Literature
!>      S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!> \par History
!>      none
!> \par Parameters
!>       - ax,ay,az    : Angular momentum index numbers of orbital a.
!>       - cx,cy,cz    : Angular momentum index numbers of orbital c.
!>       - coset       : Cartesian orbital set pointer.
!>       - dac         : Distance between the atomic centers a and c.
!>       - l{a,c}      : Angular momentum quantum number of shell a or c.
!>       - l{a,c}_max  : Maximum angular momentum quantum number of shell a or c.
!>       - l{a,c}_min  : Minimum angular momentum quantum number of shell a or c.
!>       - ncoset      : Number of orbitals in a Cartesian orbital set.
!>       - npgf{a,c}   : Degree of contraction of shell a or c.
!>       - rac         : Distance vector between the atomic centers a and c.
!>       - rac2        : Square of the distance between the atomic centers a and c.
!>       - zet{a,c}    : Exponents of the Gaussian-type functions a or c.
!>       - zetp        : Reciprocal of the sum of the exponents of orbital a and b.
!>       - zetw        : Reciprocal of the sum of the exponents of orbital a and c.
!> \author Dorothea Golze (08.2016)
! **************************************************************************************************
MODULE ai_operator_ra2m

   USE ai_cartesian_operators,          ONLY: apply_cartesian_1d,&
                                              cartesian_polynomial_1d
   USE ai_os_rr,                        ONLY: os_rr_ovlp
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fac,&
                                              pi
   USE orbital_pointers,                ONLY: indco,&
                                              ncoset
#include "../base/base_uses.f90"

   IMPLICIT NONE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_operator_ra2m'

   PRIVATE

   ! *** Public subroutines ***

   PUBLIC :: operator_ra2m

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Calculation of the primitive two-center [a|(r-Ra)^(2m)|b] integrals over Cartesian
!>        Gaussian-type functions; operator is here (r-Ra)^(2m)
!> \param la_max ...
!> \param la_min ...
!> \param npgfa ...
!> \param zeta ...
!> \param lb_max ...
!> \param lb_min ...
!> \param npgfb ...
!> \param zetb ...
!> \param m exponent in (r-Ra)^(2m) operator
!> \param rab ...
!> \param sab ...
!> \param dsab ...
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE operator_ra2m(la_max, la_min, npgfa, zeta, &
                            lb_max, lb_min, npgfb, zetb, m, rab, sab, dsab, calculate_forces)
      INTEGER, INTENT(IN)                                :: la_max, la_min, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta
      INTEGER, INTENT(IN)                                :: lb_max, lb_min, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb
      INTEGER, INTENT(IN)                                :: m
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: sab
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: dsab
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'operator_ra2m'

      INTEGER                                            :: ax, ay, az, bx, by, bz, coa, cob, &
                                                            handle, i, ia, ib, ipgf, j, jpgf, k, &
                                                            ldrr, lma, ma, mb, na, nb, ofa, ofb, &
                                                            powers(0:m)
      REAL(KIND=dp)                                      :: a, b, f0, prefac, rab2, rap(3), rbp(3), &
                                                            weights(0:m, 0:m), zet
      REAL(KIND=dp), ALLOCATABLE                         :: dpolynomial(:, :, :, :), drr(:, :, :), &
                                                            polynomial(:, :, :, :), rr(:, :, :)

      CALL timeset(routineN, handle)
      sab = 0.0_dp
      IF (calculate_forces) dsab = 0.0_dp
      rab2 = NORM2(rab)**2
      lma = la_max + 2*m
      IF (calculate_forces) lma = lma + 1
      ldrr = MAX(lma, lb_max) + 1
      ALLOCATE (rr(0:ldrr - 1, 0:ldrr - 1, 3), polynomial(0:la_max, 0:lb_max, 0:m, 3))
      IF (calculate_forces) THEN
         ALLOCATE (drr(0:ldrr - 1, 0:ldrr - 1, 3), dpolynomial(0:la_max, 0:lb_max, 0:m, 3))
      END IF
      powers = [(2*i, i=0, m)]
      DO i = 0, m
         DO j = 0, m - i
            k = m - i - j
            weights(i, j) = fac(m)/fac(i)/fac(j)/fac(k)
         END DO
      END DO
      ofa = ncoset(la_min - 1)
      ofb = ncoset(lb_min - 1)
      na = ncoset(la_max) - ofa
      nb = ncoset(lb_max) - ofb
      CPASSERT(SIZE(sab, 1) >= na*npgfa .AND. SIZE(sab, 2) >= nb*npgfb)

      DO ipgf = 1, npgfa
         ma = (ipgf - 1)*na
         DO jpgf = 1, npgfb
            mb = (jpgf - 1)*nb
            a = zeta(ipgf)
            b = zetb(jpgf)
            zet = a + b
            rap = b*rab/zet
            rbp = -a*rab/zet
            f0 = (pi/zet)**1.5_dp*EXP(-a*b*rab2/zet)
            CALL os_rr_ovlp(rap, lma, rbp, lb_max, zet, ldrr, rr)
            CALL cartesian_polynomial_1d(rr, polynomial, la_max, lb_max, powers, [0.0_dp, 0.0_dp, 0.0_dp])
            IF (calculate_forces) THEN
               ! The origin moves with A: differentiate the combined Gaussian degree first.
               ! Raising that derivative includes both basis and operator-origin contributions.
               CALL apply_cartesian_1d(rr, drr, la_max + 2*m, lb_max, .FALSE., .TRUE., a, &
                                       [0.0_dp, 0.0_dp, 0.0_dp])
               CALL cartesian_polynomial_1d(drr, dpolynomial, la_max, lb_max, powers, [0.0_dp, 0.0_dp, 0.0_dp])
            END IF

            DO cob = ofb + 1, ncoset(lb_max)
               bx = indco(1, cob)
               by = indco(2, cob)
               bz = indco(3, cob)
               ib = mb + cob - ofb
               DO coa = ofa + 1, ncoset(la_max)
                  ax = indco(1, coa)
                  ay = indco(2, coa)
                  az = indco(3, coa)
                  ia = ma + coa - ofa
                  DO i = 0, m
                     DO j = 0, m - i
                        k = m - i - j
                        prefac = weights(i, j)*f0
                        sab(ia, ib) = sab(ia, ib) + prefac*polynomial(ax, bx, i, 1)* &
                                      polynomial(ay, by, j, 2)*polynomial(az, bz, k, 3)
                        IF (calculate_forces) THEN
                           dsab(ia, ib, 1) = dsab(ia, ib, 1) + prefac*dpolynomial(ax, bx, i, 1)* &
                                             polynomial(ay, by, j, 2)*polynomial(az, bz, k, 3)
                           dsab(ia, ib, 2) = dsab(ia, ib, 2) + prefac*polynomial(ax, bx, i, 1)* &
                                             dpolynomial(ay, by, j, 2)*polynomial(az, bz, k, 3)
                           dsab(ia, ib, 3) = dsab(ia, ib, 3) + prefac*polynomial(ax, bx, i, 1)* &
                                             polynomial(ay, by, j, 2)*dpolynomial(az, bz, k, 3)
                        END IF
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      CALL timestop(handle)
   END SUBROUTINE operator_ra2m

END MODULE ai_operator_ra2m
