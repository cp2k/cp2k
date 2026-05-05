! This file is part of s-dftd3.
! SPDX-Identifier: LGPL-3.0-or-later
!
! s-dftd3 is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! s-dftd3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.

!> Versioning information on this library.
module dftd3_version
   implicit none
   private

   public :: dftd3_version_string, dftd3_version_compact
   public :: get_dftd3_version


   !> String representation of the s-dftd3 version
   character(len=*), parameter :: dftd3_version_string = "1.3.2"

   !> Numeric representation of the s-dftd3 version
   integer, parameter :: dftd3_version_compact(3) = [1, 3, 2]


contains


!> Getter function to retrieve s-dftd3 version
subroutine get_dftd3_version(major, minor, patch, string)

   !> Major version number of the s-dftd3 version
   integer, intent(out), optional :: major

   !> Minor version number of the s-dftd3 version
   integer, intent(out), optional :: minor

   !> Patch version number of the s-dftd3 version
   integer, intent(out), optional :: patch

   !> String representation of the s-dftd3 version
   character(len=:), allocatable, intent(out), optional :: string

   if (present(major)) then
      major = dftd3_version_compact(1)
   end if
   if (present(minor)) then
      minor = dftd3_version_compact(2)
   end if
   if (present(patch)) then
      patch = dftd3_version_compact(3)
   end if
   if (present(string)) then
      string = dftd3_version_string
   end if

end subroutine get_dftd3_version


end module dftd3_version
