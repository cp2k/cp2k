!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>      none
! **************************************************************************************************
MODULE structure_factor_types

   USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: structure_factor_type

! **************************************************************************************************
   TYPE structure_factor_type
      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER :: ex => NULL(), ey => NULL(), ez => NULL()
      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER :: shell_ex => NULL(), shell_ey => NULL(), shell_ez => NULL()
      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER :: core_ex => NULL(), core_ey => NULL(), core_ez => NULL()
      INTEGER, DIMENSION(:, :), POINTER :: centre => NULL(), core_centre => NULL(), shell_centre => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: delta => NULL(), core_delta => NULL(), shell_delta => NULL()
      INTEGER :: lb(3) = -1
   END TYPE structure_factor_type

END MODULE structure_factor_types

