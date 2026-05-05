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

!> @file tblite/scf/info.f90
!> Provides a density dependence information

!> Definition of information container on density dependence
module tblite_scf_info
   use mctc_env, only : sp, dp, error_type
   implicit none
   private

   public :: not_used, atom_resolved, shell_resolved, orbital_resolved, max

   integer, parameter :: not_used = 0
   integer, parameter :: atom_resolved = 1
   integer, parameter :: shell_resolved = 2
   integer, parameter :: orbital_resolved = 3

   !> Small info container about density dependence
   type, public :: scf_info
      integer :: charge = not_used
      integer :: dipole = not_used
      integer :: quadrupole = not_used
      integer :: density = not_used
   end type scf_info

   interface max
      module procedure :: max_info
   end interface max

contains

pure function max_info(lhs, rhs) result(new)
   type(scf_info), intent(in) :: lhs
   type(scf_info), intent(in) :: rhs
   type(scf_info) :: new

   new = scf_info( &
      charge=max(lhs%charge, rhs%charge), &
      dipole=max(lhs%dipole, rhs%dipole), &
      quadrupole=max(lhs%quadrupole, rhs%quadrupole), &
      density=max(lhs%density, rhs%density))

end function max_info

end module tblite_scf_info
