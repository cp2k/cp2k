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

!> @file tblite/basis/cache.f90
!> Provides a cache for the basis set construction

!> Data container for mutable data in basis set construction
module tblite_basis_cache
   use mctc_env, only : wp
   implicit none
   private

   public :: basis_cache, cgto_cache

   type :: cgto_cache
      !> Effective environment-dependent charge
      real(wp) :: qeff
      !> Derivative of the effective charge w.r.t. the coordination number
      real(wp) :: dqeffdcn
      !> Derivative of the effective charge w.r.t. the input charges
      real(wp) :: dqeffdq
      !> Normalization factor
      real(wp) :: norm
      !> Derivative of the normalization factor w.r.t. effective charge
      real(wp) :: dnorm
   end type cgto_cache

   type :: basis_cache
      !> Cache for CGTOs
      type(cgto_cache), allocatable :: cgto(:, :)
      !> Cache for H0 scaled CGTOs
      type(cgto_cache), allocatable :: cgto_h0(:, :)
   end type basis_cache

end module tblite_basis_cache
