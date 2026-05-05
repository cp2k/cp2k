! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

!> Versioning information on this library.
module dftd4_version
   implicit none
   private

   public :: dftd4_version_string, dftd4_version_compact
   public :: get_dftd4_version


   !> String representation of the dftd4 version
   character(len=*), parameter :: dftd4_version_string = "4.0.2"

   !> Numeric representation of the dftd4 version
   integer, parameter :: dftd4_version_compact(3) = [4, 0, 2]


contains


!> Getter function to retrieve dftd4 version
subroutine get_dftd4_version(major, minor, patch, string)
   !DEC$ ATTRIBUTES DLLEXPORT :: get_dftd4_version

   !> Major version number of the dftd4 version
   integer, intent(out), optional :: major

   !> Minor version number of the dftd4 version
   integer, intent(out), optional :: minor

   !> Patch version number of the dftd4 version
   integer, intent(out), optional :: patch

   !> String representation of the dftd4 version
   character(len=:), allocatable, intent(out), optional :: string

   if (present(major)) then
      major = dftd4_version_compact(1)
   end if
   if (present(minor)) then
      minor = dftd4_version_compact(2)
   end if
   if (present(patch)) then
      patch = dftd4_version_compact(3)
   end if
   if (present(string)) then
      string = dftd4_version_string
   end if

end subroutine get_dftd4_version


end module dftd4_version
