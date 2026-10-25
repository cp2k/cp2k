!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! Local constants for standalone trajana programs.
! CODATA 2006: https://doi.org/10.1103/RevModPhys.80.633
MODULE trajana_constants
   USE trajana_kinds,                   ONLY: dp

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: angstrom, bohr, boltzmann, c_light, e_charge, h_planck, n_avogadro, picoseconds

   ! Exact SI value.
   REAL(dp), PARAMETER :: c_light = 299792458.0_dp

   ! Recommended CODATA table entries.
   REAL(dp), PARAMETER :: h_planck = 6.62606896E-34_dp
   REAL(dp), PARAMETER :: e_charge = 1.602176487E-19_dp
   REAL(dp), PARAMETER :: boltzmann = 1.3806504E-23_dp
   REAL(dp), PARAMETER :: n_avogadro = 6.02214179E23_dp
   REAL(dp), PARAMETER :: a_bohr = 0.52917720859E-10_dp
   REAL(dp), PARAMETER :: rydberg = 10973731.568527_dp

   ! Derived atomic-unit conversions.
   REAL(dp), PARAMETER :: angstrom = 1.0E10_dp*a_bohr
   REAL(dp), PARAMETER :: bohr = 1.0_dp/angstrom
   REAL(dp), PARAMETER :: seconds = 1.0_dp/(4.0_dp*ACOS(-1.0_dp)*rydberg*c_light)
   REAL(dp), PARAMETER :: picoseconds = 1.0E12_dp*seconds

END MODULE trajana_constants
