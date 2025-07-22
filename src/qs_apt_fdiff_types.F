!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Atomic Polarization Tensor calculation by dF/d(E-field) finite differences
!> \author Leo Decking, Hossam Elgabarty
! **************************************************************************************************

MODULE qs_apt_fdiff_types

   USE kinds,                           ONLY: dp
   !USE physcon, ONLY: debye
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_apt_fdiff_types'
   LOGICAL, PARAMETER, PRIVATE :: debug_this_module = .FALSE.
   PUBLIC :: apt_fdiff_point_type, apt_fdiff_points_type

   TYPE apt_fdiff_point_type
      REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE :: forces
   END TYPE apt_fdiff_point_type

   TYPE apt_fdiff_points_type
      !TYPE(apt_fdiff_point_type) :: point_no_field ! Only for forward differences
      TYPE(apt_fdiff_point_type), DIMENSION(1:3, 1:2) :: point_field = apt_fdiff_point_type()
      REAL(kind=dp)                                   :: field_strength = 0.0_dp
   END TYPE apt_fdiff_points_type

END MODULE qs_apt_fdiff_types
