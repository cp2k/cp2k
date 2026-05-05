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

!> @file tblite/acp/cache.f90
!> Provides a cache for the atomic correction potential

!> Data container for mutable data in atomic correction potential
module tblite_acp_cache
   use tblite_basis_cache, only : basis_cache
   use mctc_env, only : wp
   implicit none
   private

   public :: acp_cache

   type :: acp_cache
      !> Cache for the auxiliary basis set
      type(basis_cache) :: auxbas

      !> Scaled projector-valence overlap integral intermediate for the gradient
      real(wp), allocatable :: scaled_pv_overlap(:, :)
   end type acp_cache

end module tblite_acp_cache
