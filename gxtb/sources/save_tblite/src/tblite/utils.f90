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

!> @dir tblite/utils
!> Contains utility functions for tblite

!> @file tblite/utils.f90
!> Reexports utility functions for tblite

!> Proxy module for rexports from utility modules
module tblite_utils
   use tblite_utils_average, only : average_type, new_average, average_id
   use tblite_utils_string, only : compact, lowercase
   implicit none
   private

   public :: compact, lowercase
   public :: average_type, new_average, average_id

end module tblite_utils
