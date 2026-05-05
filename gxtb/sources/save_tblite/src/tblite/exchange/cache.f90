
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

!> @file tblite/exchange/cache.f90
!> Provides a cache specific for all exchange interactions

!> Data container for mutable data in exchange calculations
module tblite_exchange_cache
   use mctc_env, only : wp, dp
   use mctc_io, only : structure_type
   use tblite_coulomb_ewald, only : get_alpha
   use tblite_wignerseitz, only : wignerseitz_cell, new_wignerseitz_cell
   implicit none
   private

   public :: exchange_cache

   type :: exchange_cache
      !> Ewald parameter
      real(wp) :: alpha
      !> Wigner-Seitz cell for periodic calculations
      type(wignerseitz_cell) :: wsc

      !> Mulliken approximate exchange matrix: [nsh, nsh]
      real(wp), allocatable :: g_mulliken(:, :)
      !> Integrals for the onsite exchange correction: [maxsh, maxsh, nat]
      real(wp), allocatable :: g_onsfx(:, :, :)
      !> Charge-derivative of the integrals for the onsite exchange correction: [maxsh, maxsh, nsh]
      real(wp), allocatable :: dgdq_onsfx(:, :, :)
      !> Integrals for the onsite rotational invariance correction: [maxsh, nat]
      real(wp), allocatable :: g_onsri(:, :)
      !> Charge-derivative of the onsite rotational invariance correction matrix: [maxsh, nsh]
      real(wp), allocatable :: dgdq_onsri(:, :)
      !> Integrals for the bond-order correlation correction: [nat, nat]
      real(wp), allocatable :: g_bocorr(:, :)

      !> Previously calculated Fock matrix contribution 
      real(wp), allocatable :: prev_F(:, :, :)
      !> Previously calculated shell potential contribution
      real(wp), allocatable :: prev_vsh(:, :)
   contains
      procedure :: update
   end type exchange_cache


contains


subroutine update(self, mol)
   !> Instance of the electrostatic container
   class(exchange_cache), intent(inout) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   if (any(mol%periodic)) then
      call new_wignerseitz_cell(self%wsc, mol)
      call get_alpha(mol%lattice, self%alpha, .false.)
   end if

end subroutine update

end module tblite_exchange_cache
