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

!> @dir tblite/api
!> Provides C API bindings for tblite library

!> @file tblite/api/version.f90
!> Provides a stable version query with #tblite_get_version.

!> API export for the version information of the library
module tblite_api_version
   use, intrinsic :: iso_c_binding
   use tblite_version, only : get_tblite_version
   implicit none
   private

   public :: get_version_api, namespace

   character(len=*), parameter :: namespace = "tblite_"

contains


!> Obtain library version as major * 10000 + minor + 100 + patch
function get_version_api() result(version) &
      & bind(C, name=namespace//"get_version")
   integer(c_int) :: version
   integer :: major, minor, patch

   call get_tblite_version(major, minor, patch)
   version = 10000_c_int * major + 100_c_int * minor + patch

end function get_version_api


end module tblite_api_version
