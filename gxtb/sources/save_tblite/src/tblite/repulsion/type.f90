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

!> @file tblite/repulsion/type.f90
!> Provides a general implementation of the screened Coulomb repulsion interaction

!> Classical repulsion interaction as used with the xTB Hamiltonian
module tblite_repulsion_type
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_classical_type, only : classical_type
   use tblite_container, only : container_cache
   use tblite_repulsion_cache, only : repulsion_cache
   implicit none
   private

   public :: repulsion_kernel

   !> Container to evaluate classical repulsion interactions for the xTB Hamiltonian
   type, public, extends(classical_type), abstract :: repulsion_type
      !> Real-space cutoff
      real(wp) :: cutoff = 25.0_wp
   end type repulsion_type

   !> Possible repulsion kernels
   type :: enum_repulsion_kernel
      !> Effective Coulomb repulsion kernel (GFNn-xTB)
      integer :: gfn = 1
      !> Effective Coulomb + Pauli repulsion kernel (g-xTB)
      integer :: gxtb = 2
   end type enum_repulsion_kernel

   !> Actual enumerator for possible repulsion kernels
   type(enum_repulsion_kernel), parameter :: repulsion_kernel = enum_repulsion_kernel()

end module tblite_repulsion_type
