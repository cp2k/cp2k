! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! tblite is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! tblite is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with tblite.  If not, see <https://www.gnu.org/licenses/>.

!> @file tblite/disp/cache.f90
!> Provides a cache for dispersion interactions

!> Reusable data container for dispersion related calculations
module tblite_disp_cache
   use mctc_env, only : wp
   use dftd4_cache, only : dftd4_cache_type => dispersion_cache
   implicit none
   private
   public :: dispersion_cache

   type :: dispersion_cache
      !> Damped R^-6 dispersion matrix: [nat, nat]
      real(wp), allocatable :: dispmat6(:, :)
      !> Damped R^-8 dispersion matrix: [nat, nat]
      real(wp), allocatable :: dispmat8(:, :)
      !> Coordination numbers: [nat]
      real(wp), allocatable :: cn(:)
      !> Derivative of CN w.r.t. coordinates: [3, nat, nat]
      real(wp), allocatable :: dcndr(:, :, :)
      !> Derivative of CN w.r.t. lattice vectors: [3, 3, nat]
      real(wp), allocatable :: dcndL(:, :, :)
      !> DFT-D4 internal cache for C6 coefficients and derivatives
      type(dftd4_cache_type) :: disp_cache
   end type dispersion_cache


end module tblite_disp_cache
