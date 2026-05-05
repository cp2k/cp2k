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

!> @dir tblite/acp
!> Contains the atomic correction potential.

!> @file tblite/acp.f90
!> Reexports of the atomic correction potential

!> Proxy module for reexporting atomic correction potential
module tblite_acp
   use tblite_acp_cache, only : acp_cache
   use tblite_acp_type, only : acp_type, new_acp
   implicit none
   private

   public :: acp_type, new_acp, acp_cache

end module tblite_acp
