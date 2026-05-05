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

!> @file tblite/coulomb/multipole.f90
!> Provides an implemenation of a multipole based second-order electrostatic

!> Anisotropic second-order electrostatics using a damped multipole expansion

!> @dir tblite/coulomb/multipole
!> Contains the implementation of a multipole based second-order electrostatics

!> @file tblite/coulomb/multipole.f90
!> Provides a reexport of the anisotropic second-order electrostatic implementations

!> Proxy module for defining anisotropic electrostatic interactions
module tblite_coulomb_multipole
   use tblite_coulomb_multipole_gfn2, only : gfn2_multipole, new_gfn2_multipole
   use tblite_coulomb_multipole_gxtb, only : gxtb_multipole, new_gxtb_multipole
   use tblite_coulomb_multipole_type, only : damped_multipole
   implicit none
   private

   public :: damped_multipole, multipole_damping
   public :: gfn2_multipole, new_gfn2_multipole
   public :: gxtb_multipole, new_gxtb_multipole

   !> Possible multipole damping types
   type :: enum_multipole_damping
      !> Zero-damping for GFN2-xTB
      integer :: gfn2 = 1
      !> Error-function damping for g-xTB
      integer :: gxtb = 2
   end type enum_multipole_damping

   !> Actual enumerator for possible multipole damping types
   type(enum_multipole_damping), parameter :: multipole_damping = enum_multipole_damping()

end module tblite_coulomb_multipole
