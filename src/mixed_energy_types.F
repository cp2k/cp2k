!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \author fschiff
!> \date   11.06
! **************************************************************************************************
MODULE mixed_energy_types

   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mixed_energy_types'
   PRIVATE

! **************************************************************************************************
   TYPE mixed_energy_type
      REAL(kind=dp) :: pot = 0.0_dp
      REAL(kind=dp) :: kin = 0.0_dp
   END TYPE mixed_energy_type

! **************************************************************************************************
   TYPE mixed_force_type
      REAL(KIND=dp), DIMENSION(:, :), POINTER         :: forces => NULL()
   END TYPE mixed_force_type

! Public data types
   PUBLIC :: mixed_energy_type, &
             mixed_force_type

END MODULE mixed_energy_types
