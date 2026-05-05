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

!> @dir tblite/repulsion
!> Contains the repulsive interaction implementations.

!> @file tblite/repulsion.f90
!> Provides a repulsion container base class

!> Proxy module for repulsion interactions.
module tblite_repulsion
   use tblite_repulsion_cache, only : repulsion_cache
   use tblite_repulsion_gfn, only : gfn_repulsion, new_repulsion_gfn
   use tblite_repulsion_gxtb, only : gxtb_repulsion, new_repulsion_gxtb
   use tblite_repulsion_type, only : repulsion_type, repulsion_kernel
   implicit none
   public

end module tblite_repulsion
