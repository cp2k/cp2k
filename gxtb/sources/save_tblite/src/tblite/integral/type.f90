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

!> @dir tblite/integral
!> Contains the integral evaluation implementations

!> @file tblite/integral/type.f90
!> Provides an integral storage container

!> Declaration of an integral storage container to collect all overlap related integrals
module tblite_integral_type
   use mctc_env, only : wp
   implicit none
   private

   public :: new_integral

   !> Integral container to store all overlap related integrals
   type, public :: integral_type
      !> Effective one-electron Hamiltonian
      real(wp), allocatable :: hamiltonian(:, :)
      !> Overlap integrals
      real(wp), allocatable :: overlap(:, :)
      !> Dipole moment integrals, moment operator is centered on last index
      real(wp), allocatable :: dipole(:, :, :)
      !> Quadrupole moment integrals, moment operator is centered on last index
      real(wp), allocatable :: quadrupole(:, :, :)
      !> Diatomic frame scaled overlap integrals
      real(wp), allocatable :: overlap_diat(:, :)
   end type integral_type

contains

!> Create and allocate a new integral container storage
subroutine new_integral(self, nao)
   !> Instance of the integral container
   type(integral_type), intent(out) :: self
   !> Dimension of the integrals
   integer, intent(in) :: nao

   allocate(self%hamiltonian(nao, nao), source = 0.0_wp)
   allocate(self%overlap(nao, nao), source = 0.0_wp)
   allocate(self%dipole(3, nao, nao), source = 0.0_wp)
   allocate(self%quadrupole(6, nao, nao), source = 0.0_wp)
   allocate(self%overlap_diat(nao, nao), source = 0.0_wp)
end subroutine new_integral

end module tblite_integral_type
