!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM qs_fod_unittest
   USE kinds,                           ONLY: dp
   USE qs_fod,                          ONLY: fod_weights

   IMPLICIT NONE
   REAL(dp), DIMENSION(4) :: occ, eig, weight, expected

   eig = [-1.0_dp, -0.1_dp, 0.0_dp, 0.2_dp]
   occ = [2.0_dp, 1.6_dp, 1.0_dp, 0.2_dp]
   expected = [0.0_dp, 0.4_dp, 1.0_dp, 0.2_dp]
   weight = fod_weights(occ, eig, 0.0_dp, 2.0_dp)
   IF (MAXVAL(ABS(weight - expected)) > 1.0E-14_dp) ERROR STOP "Restricted FOD weights"
   weight = fod_weights(0.5_dp*occ, eig, 0.0_dp, 1.0_dp)
   IF (MAXVAL(ABS(2.0_dp*weight - expected)) > 1.0E-14_dp) ERROR STOP "Unrestricted FOD weights"

   ! Fractional spin populations: the electron-count/index split is not valid.
   eig = [0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp]
   occ = [0.4_dp, 0.3_dp, 0.2_dp, 0.1_dp]
   weight = fod_weights(occ, eig, 0.0_dp, 1.0_dp)
   IF (MAXVAL(ABS(weight - occ)) > 1.0E-14_dp) ERROR STOP "FOD must use the chemical potential"
   weight = fod_weights(occ(4:1:-1), eig(4:1:-1), 0.0_dp, 1.0_dp)
   IF (MAXVAL(ABS(weight - occ(4:1:-1))) > 1.0E-14_dp) ERROR STOP "FOD orbital permutation"

   ! Integer occupations vanish; a shift of the energy zero changes nothing.
   eig = [-2.0_dp, -1.0_dp, 1.0_dp, 2.0_dp]
   occ = [2.0_dp, 2.0_dp, 0.0_dp, 0.0_dp]
   weight = fod_weights(occ, eig + 3.0_dp, 3.0_dp, 2.0_dp)
   IF (ANY(weight /= 0.0_dp)) ERROR STOP "Integer occupations should have no FOD"
END PROGRAM qs_fod_unittest
