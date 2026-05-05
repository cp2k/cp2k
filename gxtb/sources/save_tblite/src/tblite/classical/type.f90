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

!> @file tblite/classical/type.f90
!> Provides a base class for classical interactions

!> Definition of the abstract base class for classical interactions.
module tblite_classical_type
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_container_type, only : container_type
   implicit none
   private

   !> Abstract base class for classical interactions, like repulsion interactions.
   !>
   !> This class provides a method to retrieve the contributions to the energy,
   !> gradient and virial within a given cutoff.
   type, public, extends(container_type), abstract :: classical_type
   end type classical_type

end module tblite_classical_type
