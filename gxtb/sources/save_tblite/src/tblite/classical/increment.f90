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

!> @dir tblite/classical
!> Contains classical correction potentials

!> @file tblite/classical/increment.f90
!> Provides atomic energy increment for missing core electrons

!> Atomic energy increment
module tblite_classical_increment
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_container, only : container_cache
   use tblite_classical_type, only : classical_type
   implicit none
   private

   public :: new_core_increment


   !> Container for evaluating the atomic energy increment
   type, public, extends(classical_type) :: core_increment
      !> Atomic energy increment
      real(wp), allocatable :: increment(:)
   contains
      !> Entry point for evaluation of energy and gradient
      procedure :: get_engrad
   end type core_increment

   character(len=*), parameter :: label = "core increment"

contains


!> Construct new core energy increment
subroutine new_core_increment(self, mol, increment)
   !> Instance of the core increment
   type(core_increment), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Energy increment for the missing core electrons
   real(wp), intent(in) :: increment(:)

   self%label = label

   allocate(self%increment(mol%nid))
   self%increment(:) = increment

end subroutine new_core_increment

!> Evaluate classical interaction for energy and derivatives
subroutine get_engrad(self, mol, cache, energies, gradient, sigma)
   !> Instance of the core increment
   class(core_increment), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
   !> Core energy increment 
   real(wp), intent(inout) :: energies(:)
   !> Molecular gradient of the core increment 
   real(wp), contiguous, intent(inout), optional :: gradient(:, :)
   !> Strain derivatives of the core increment
   real(wp), contiguous, intent(inout), optional :: sigma(:, :)

   integer :: iat, izp

   do iat = 1, mol%nat
      izp = mol%id(iat)
      ! Add atomic energy increment
      energies(iat) = energies(iat) + self%increment(izp)
   end do 

end subroutine get_engrad

end module tblite_classical_increment
