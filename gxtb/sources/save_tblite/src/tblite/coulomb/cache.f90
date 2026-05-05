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

!> @file tblite/coulomb/cache.f90
!> Provides a cache specific for all Coulomb interactions

!> Data container for mutable data in electrostatic calculations
module tblite_coulomb_cache
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_coulomb_ewald, only : get_alpha
   use tblite_wignerseitz, only : wignerseitz_cell, new_wignerseitz_cell
   implicit none
   private

   public :: coulomb_cache


   type :: coulomb_cache
      real(wp) :: alpha
      type(wignerseitz_cell) :: wsc

      ! Coordination number for Coulomb and tight-binding
      real(wp), allocatable :: cn(:)
      real(wp), allocatable :: dcndr(:, :, :)
      real(wp), allocatable :: dcndL(:, :, :)

      ! Second-order electrostatics intermediates
      real(wp), allocatable :: amat(:, :)
      real(wp), allocatable :: vvec(:)

      ! Multipole intermediates
      real(wp), allocatable :: mrad(:)
      real(wp), allocatable :: dmrdcn(:)

      real(wp), allocatable :: amat_sd(:, :, :)
      real(wp), allocatable :: amat_dd(:, :, :, :)
      real(wp), allocatable :: amat_sq(:, :, :)
      real(wp), allocatable :: amat_dq(:, :, :, :)
      real(wp), allocatable :: amat_qq(:, :, :, :)

      ! Third-order tight-binding intermediates
      real(wp), allocatable :: taumat(:, :)
   contains
      procedure :: update
   end type coulomb_cache


contains


subroutine update(self, mol)
   !> Instance of the electrostatic container
   class(coulomb_cache), intent(inout) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   if (any(mol%periodic)) then
      call new_wignerseitz_cell(self%wsc, mol)
      call get_alpha(mol%lattice, self%alpha, .false.)
   end if

end subroutine update

end module tblite_coulomb_cache
