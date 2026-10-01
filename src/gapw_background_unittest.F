!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM gapw_background_unittest
   USE hartree_local_methods,           ONLY: gapw_background_potential
   USE kinds,                           ONLY: dp
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init
   USE mathconstants,                   ONLY: fourpi,&
                                              pi

   IMPLICIT NONE

   INTEGER, PARAMETER :: nr = 1601, ng = 32
   REAL(dp), PARAMETER :: alpha = 1.4_dp, beta = 0.55_dp, box = 12.0_dp, dr = 0.00625_dp
   REAL(dp), PARAMETER :: charges(3) = [-1.0_dp, 0.3_dp, 1.0_dp]
   REAL(dp) :: rad2(nr), weights(nr), delta_rho(nr), potential(nr), shifted(nr)
   REAL(dp) :: q, g2, periodic_difference, local_difference, correction, error, h, eplus, eminus
   INTEGER :: i, ix, iy, iz

   CALL cp2k_init()
   DO i = 1, nr
      rad2(i) = ((i - 1)*dr)**2
   END DO
   weights = 2.0_dp*dr/3.0_dp
   weights(2:nr - 1:2) = 4.0_dp*dr/3.0_dp
   weights(1) = dr/3.0_dp
   weights(nr) = dr/3.0_dp
   weights = weights*rad2

   DO i = 1, SIZE(charges)
      q = charges(i)
      ! Electron-minus-nuclear l=0 density. Both normalized Gaussians carry q,
      ! so their difference has zero monopole and a nonzero second moment.
      delta_rho = q*SQRT(fourpi)*((alpha/pi)**1.5_dp*EXP(-alpha*rad2) - &
                                  (beta/pi)**1.5_dp*EXP(-beta*rad2))
      potential = gapw_background_potential(-q, box**3, rad2)
      correction = SUM(potential*delta_rho*weights)

      ! Independent reference: exact reciprocal-space Coulomb energies of the
      ! two periodic Gaussians with G=0 omitted, minus their free-space energies.
      periodic_difference = 0.0_dp
      DO iz = -ng, ng
         DO iy = -ng, ng
            DO ix = -ng, ng
               IF (ix == 0 .AND. iy == 0 .AND. iz == 0) CYCLE
               g2 = (2.0_dp*pi/box)**2*REAL(ix*ix + iy*iy + iz*iz, dp)
               periodic_difference = periodic_difference + &
                                     2.0_dp*pi*q*q/box**3/g2*(EXP(-g2/(2.0_dp*alpha)) - EXP(-g2/(2.0_dp*beta)))
            END DO
         END DO
      END DO
      local_difference = q*q*(SQRT(alpha) - SQRT(beta))/SQRT(2.0_dp*pi)
      error = ABS(periodic_difference - local_difference - correction)
      WRITE (*, '(A,F8.3,ES16.7)') "GAPW periodic Gaussian background:", q, error
      IF (error > 1.e-12_dp) ERROR STOP "GAPW background disagrees with periodic Gaussian electrostatics"

      ! The explicit volume derivative contributes E_background to each
      ! diagonal component of -dE/dstrain (local radial functions held fixed).
      h = 1.e-5_dp
      shifted = gapw_background_potential(-q, box**3*(1.0_dp + h), rad2)
      eplus = SUM(shifted*delta_rho*weights)
      shifted = gapw_background_potential(-q, box**3*(1.0_dp - h), rad2)
      eminus = SUM(shifted*delta_rho*weights)
      IF (ABS(-(eplus - eminus)/(2.0_dp*h) - correction) > 1.e-11_dp) &
         ERROR STOP "GAPW background volume derivative failed"
   END DO
   potential = gapw_background_potential(0.0_dp, box**3, rad2)
   IF (ANY(potential /= 0.0_dp)) ERROR STOP "Neutral cell has a background potential"
   CALL cp2k_finalize()
END PROGRAM gapw_background_unittest
