!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
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
      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER :: ex, ey, ez
      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER :: shell_ex, shell_ey, shell_ez
      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER :: core_ex, core_ey, core_ez
      INTEGER, DIMENSION(:, :), POINTER :: centre, core_centre, shell_centre
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: delta, core_delta, shell_delta
      INTEGER :: lb(3)
   END TYPE structure_factor_type

END MODULE structure_factor_types

