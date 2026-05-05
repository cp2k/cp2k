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

!> @file tblite/version.f90
!> Provides version and feature information

!> Interfaces to query the version information of the library.
module tblite_version
   implicit none
   private

   public :: tblite_version_string, tblite_version_compact
   public :: get_tblite_version


   !> String representation of the tblite version
   character(len=*), parameter :: tblite_version_string = "0.5.0"

   !> Numeric representation of the tblite version
   integer, parameter :: tblite_version_compact(3) = [0, 5, 0]


contains


!> Getter function to retrieve tblite version
subroutine get_tblite_version(major, minor, patch, string)
   !> Major version number of the tblite version
   integer, intent(out), optional :: major
   !> Minor version number of the tblite version
   integer, intent(out), optional :: minor
   !> Patch version number of the tblite version
   integer, intent(out), optional :: patch
   !> String representation of the tblite version
   character(len=:), allocatable, intent(out), optional :: string

   if (present(major)) then
      major = tblite_version_compact(1)
   end if
   if (present(minor)) then
      minor = tblite_version_compact(2)
   end if
   if (present(patch)) then
      patch = tblite_version_compact(3)
   end if
   if (present(string)) then
      string = tblite_version_string
   end if

end subroutine get_tblite_version


end module tblite_version
