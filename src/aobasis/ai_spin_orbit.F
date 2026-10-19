!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of spin orbit integrals over Cartesian Gaussian-type functions
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
MODULE ai_spin_orbit
   USE ai_cartesian_operators,          ONLY: apply_cartesian_operator_pair,&
                                              cartesian_ladder_map
   USE ai_os_rr,                        ONLY: os_rr_coul
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
   USE orbital_pointers,                ONLY: ncoset
#include "../base/base_uses.f90"

   IMPLICIT NONE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_spin_orbit'
   PRIVATE

   ! *** Public subroutines ***

   PUBLIC :: pso

CONTAINS

! **************************************************************************************************
!> \brief   Calculation of the primitive paramagnetic spin orbit integrals over
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
   SUBROUTINE pso(la_max, la_min, npgfa, rpgfa, zeta, lb_max, lb_min, npgfb, rpgfb, zetb, &
                  rac, rbc, rab, vab, ldrr1, ldrr2, rr)
      INTEGER, INTENT(IN)                                :: la_max, la_min, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rpgfa, zeta
      INTEGER, INTENT(IN)                                :: lb_max, lb_min, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rpgfb, zetb
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rac, rbc, rab
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: vab
      INTEGER, INTENT(IN)                                :: ldrr1, ldrr2
      REAL(dp), DIMENSION(0:ldrr1-1, ldrr2, *), &
         INTENT(INOUT)                                   :: rr

      INTEGER, DIMENSION(2, 6), PARAMETER :: &
         components = RESHAPE([2, 3, 3, 2, 3, 1, 1, 3, 1, 2, 2, 1], [2, 6])

      INTEGER                                            :: axis, first_a, first_b, i, ipgf, j, &
                                                            jpgf, na, nb, ncol, nrow
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: map_a, map_b
      REAL(dp)                                           :: dab, f0, rab2, xhi, zet, zetab
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: work
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: derivatives
      REAL(dp), DIMENSION(3)                             :: rap, rbp, rcp

      first_a = ncoset(la_min - 1) + 1
      first_b = ncoset(lb_min - 1) + 1
      nrow = ncoset(la_max) - first_a + 1
      ncol = ncoset(lb_max) - first_b + 1
      ALLOCATE (map_a(4, nrow, 3), map_b(4, ncol, 3))
      ALLOCATE (work(nrow, ncoset(lb_max + 1)), derivatives(nrow, ncol, 6))
      DO axis = 1, 3
         CALL cartesian_ladder_map(la_max, la_min, axis, map_a(:, :, axis), source_lmin=0)
         CALL cartesian_ladder_map(lb_max, lb_min, axis, map_b(:, :, axis), source_lmin=0)
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
                  END DO
               END DO
               nb = nb + ncoset(lb_max)
               CYCLE
            END IF

            ! *** Calculate some prefactors ***

            zetab = zeta(ipgf)*zetb(jpgf)
            zet = zeta(ipgf) + zetb(jpgf)
            xhi = zetab/zet
            rap = zetb(jpgf)*rab/zet
            rbp = -zeta(ipgf)*rab/zet
            rcp = -(zeta(ipgf)*rac + zetb(jpgf)*rbc)/zet

            f0 = 2.0_dp*SQRT(zet/pi)*(pi/zet)**(1.5_dp)*EXP(-xhi*rab2)

            ! *** Calculate the recurrence relation ***

            CALL os_rr_coul(rap, la_max + 1, rbp, lb_max + 1, rcp, zet, ldrr1, ldrr2, rr)

            ! pso = (T_yz-T_zy, T_zx-T_xz, T_xy-T_yx), T_ij = D_Ai D_Bj V.
            CALL apply_cartesian_operator_pair(rr(0, 1:ncoset(la_max + 1), 1:ncoset(lb_max + 1)), &
                                               derivatives, map_a, map_b, components, &
                                               [.TRUE., .TRUE.], [zeta(ipgf), zetb(jpgf)], [.FALSE., .TRUE.], work)
            DO axis = 1, 3
               vab(na + first_a:na + ncoset(la_max), nb + first_b:nb + ncoset(lb_max), axis) = &
                  f0*(derivatives(:, :, 2*axis - 1) - derivatives(:, :, 2*axis))
            END DO

            nb = nb + ncoset(lb_max)

         END DO

         na = na + ncoset(la_max)

      END DO

   END SUBROUTINE pso

END MODULE ai_spin_orbit
