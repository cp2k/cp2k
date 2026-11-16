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
!>       - coset       : Cartesian orbital set pointer.
!>       - dab         : Distance between the atomic centers a and b.
!>       - dac         : Distance between the atomic centers a and c.
!>       - dbc         : Distance between the atomic centers b and c.
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
!> \author VW
! **************************************************************************************************
MODULE ai_elec_field
   USE ai_cartesian_operators,          ONLY: apply_cartesian_operator_pair,&
                                              cartesian_ladder_map
   USE ai_os_rr,                        ONLY: os_rr_coul
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
   USE orbital_pointers,                ONLY: ncoset
#include "../base/base_uses.f90"

   IMPLICIT NONE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_elec_field'
   PRIVATE

   ! *** Public subroutines ***

   PUBLIC :: efg

CONTAINS

! **************************************************************************************************
!> \brief   Calculation of the primitive electric field integrals over
!>          Cartesian Gaussian-type functions.
!> \param la_max ...
!> \param la_min ...
!> \param npgfa ...
!> \param rpgfa ...
!> \param zeta ...
!> \param lb_max ...
!> \param lb_min ...
!> \param npgfb ...
!> \param rpgfb ...
!> \param zetb ...
!> \param rac ...
!> \param rbc ...
!> \param rab ...
!> \param vab ...
!> \param ldrr1 ...
!> \param ldrr2 ...
!> \param rr ...
!> \date    02.03.2009
!> \author  VW
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE efg(la_max, la_min, npgfa, rpgfa, zeta, &
                  lb_max, lb_min, npgfb, rpgfb, zetb, &
                  rac, rbc, rab, vab, ldrr1, ldrr2, rr)
      INTEGER, INTENT(IN)                                :: la_max, la_min, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rpgfa, zeta
      INTEGER, INTENT(IN)                                :: lb_max, lb_min, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rpgfb, zetb
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rac, rbc, rab
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: vab
      INTEGER, INTENT(IN)                                :: ldrr1, ldrr2
      REAL(KIND=dp), DIMENSION(0:ldrr1-1, ldrr2, *), &
         INTENT(INOUT)                                   :: rr

      INTEGER, DIMENSION(2, 6), PARAMETER :: &
         symmetric_components = RESHAPE([1, 1, 1, 2, 1, 3, 2, 2, 2, 3, 3, 3], [2, 6])
      INTEGER, DIMENSION(2, 9), PARAMETER :: &
         mixed_components = RESHAPE([1, 1, 2, 1, 3, 1, 1, 2, 2, 2, 3, 2, 1, 3, 2, 3, 3, 3], [2, 9])

      INTEGER                                            :: axis, first_a, first_b, i, icomp, idir, &
                                                            ipgf, j, jdir, jpgf, na, nb, ncap, &
                                                            ncbp, ncol, nrow
      INTEGER, ALLOCATABLE                               :: map_a(:, :, :), map_ap(:, :, :), &
                                                            map_b(:, :, :), map_bp(:, :, :)
      REAL(KIND=dp)                                      :: dab, f0, rab2, trace_third, xhi, za, zb, &
                                                            zet
      REAL(KIND=dp), ALLOCATABLE                         :: hessian(:, :, :), tensor(:, :, :), &
                                                            work(:, :)
      REAL(KIND=dp), DIMENSION(3)                        :: rap, rbp, rcp

      first_a = ncoset(la_min - 1) + 1
      first_b = ncoset(lb_min - 1) + 1
      nrow = ncoset(la_max) - first_a + 1
      ncol = ncoset(lb_max) - first_b + 1
      ncap = ncoset(la_max + 1)
      ncbp = ncoset(lb_max + 1)
      ALLOCATE (map_a(4, nrow, 3), map_ap(4, ncap, 3), map_b(4, ncol, 3), map_bp(4, ncbp, 3))
      ALLOCATE (hessian(nrow, ncol, 6), tensor(nrow, ncol, 9), work(ncap, ncbp))
      DO axis = 1, 3
         CALL cartesian_ladder_map(la_max, la_min, axis, map_a(:, :, axis), source_lmin=0)
         CALL cartesian_ladder_map(la_max + 1, 0, axis, map_ap(:, :, axis))
         CALL cartesian_ladder_map(lb_max, lb_min, axis, map_b(:, :, axis), source_lmin=0)
         CALL cartesian_ladder_map(lb_max + 1, 0, axis, map_bp(:, :, axis))
      END DO

! *** Calculate the distance of the centers a and c ***

      rab2 = rab(1)**2 + rab(2)**2 + rab(3)**2
      dab = SQRT(rab2)

      ! *** Loop over all pairs of primitive Gaussian-type functions ***

      na = 0

      DO ipgf = 1, npgfa

         nb = 0

         DO jpgf = 1, npgfb

            ! *** Screening ***

            IF (rpgfa(ipgf) + rpgfb(jpgf) < dab) THEN
               DO j = nb + 1, nb + ncoset(lb_max)
                  DO i = na + 1, na + ncoset(la_max)
                     vab(i, j, 1) = 0.0_dp
                     vab(i, j, 2) = 0.0_dp
                     vab(i, j, 3) = 0.0_dp
                     vab(i, j, 4) = 0.0_dp
                     vab(i, j, 5) = 0.0_dp
                     vab(i, j, 6) = 0.0_dp
                  END DO
               END DO
               nb = nb + ncoset(lb_max)
               CYCLE
            END IF

            ! *** Calculate some prefactors ***

            za = zeta(ipgf)
            zb = zetb(jpgf)
            zet = za + zb
            xhi = za*zb/zet
            rap = zb*rab/zet
            rbp = -za*rab/zet
            rcp = -(za*rac + zb*rbc)/zet
            f0 = 2.0_dp*SQRT(zet/pi)*(pi/zet)**(1.5_dp)*EXP(-xhi*rab2)

            ! *** Calculate the recurrence relation

            CALL os_rr_coul(rap, la_max + 2, rbp, lb_max + 2, rcp, zet, ldrr1, ldrr2, rr)

            ! Translational invariance gives D_C = -(D_A + D_B).
            ! The A/B ladders differentiate the basis functions with the potential held fixed.
            ! Form its Hessian from AA, BB and both ordered AB derivatives.
            CALL apply_cartesian_operator_pair(rr(0, 1:ncoset(la_max + 2), first_b:ncoset(lb_max)), &
                                               hessian, map_ap, map_a, symmetric_components, [.TRUE., .TRUE.], [za, za], &
                                               [.FALSE., .FALSE.], work(:, 1:ncol))
            CALL apply_cartesian_operator_pair(rr(0, first_a:ncoset(la_max), 1:ncoset(lb_max + 2)), &
                                               tensor(:, :, 1:6), map_bp, map_b, symmetric_components, [.TRUE., .TRUE.], [zb, zb], &
                                               [.TRUE., .TRUE.], work(1:nrow, :))
            hessian(:, :, :) = hessian + tensor(:, :, 1:6)
            CALL apply_cartesian_operator_pair(rr(0, 1:ncap, 1:ncbp), tensor, map_a, map_b, &
                                               mixed_components, [.TRUE., .TRUE.], [za, zb], [.FALSE., .TRUE.], work(1:nrow, :))
            DO icomp = 1, 6
               idir = symmetric_components(1, icomp)
               jdir = symmetric_components(2, icomp)
               hessian(:, :, icomp) = hessian(:, :, icomp) + &
                                      tensor(:, :, idir + 3*(jdir - 1)) + tensor(:, :, jdir + 3*(idir - 1))
            END DO

            ! Preserve the symmetric, traceless EFG ordering: xx, xy, xz, yy, yz, zz.
            DO j = 1, ncol
               DO i = 1, nrow
                  trace_third = (hessian(i, j, 1) + hessian(i, j, 4) + hessian(i, j, 6))/3.0_dp
                  hessian(i, j, 1) = hessian(i, j, 1) - trace_third
                  hessian(i, j, 4) = hessian(i, j, 4) - trace_third
                  hessian(i, j, 6) = hessian(i, j, 6) - trace_third
               END DO
            END DO
            vab(na + first_a:na + ncoset(la_max), nb + first_b:nb + ncoset(lb_max), 1:6) = f0*hessian

            nb = nb + ncoset(lb_max)

         END DO

         na = na + ncoset(la_max)

      END DO

   END SUBROUTINE efg

END MODULE ai_elec_field
