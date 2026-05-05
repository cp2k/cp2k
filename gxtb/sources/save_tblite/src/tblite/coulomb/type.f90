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

!> @file tblite/coulomb/type.f90
!> Provides a base class for all Coulombic interactions

!> Definition of the abstract base class for electrostatic and coulombic interactions
module tblite_coulomb_type
   use tblite_container_type, only : container_type
   implicit none
   private

   !> General base class for Coulombic interactions
   type, public, extends(container_type), abstract :: coulomb_type
   end type coulomb_type

end module tblite_coulomb_type
