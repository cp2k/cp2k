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

!> @file tblite/output/format.f90
!> Provides conversion routines from intrinsic data types to string

!> Functional procedures for creating strings from intrinsic data types
module tblite_output_format
   use mctc_env, only : wp
   implicit none
   private

   public :: format_string

   interface format_string
      module procedure :: format_string_int
      module procedure :: format_string_real_dp
   end interface format_string

contains

pure function format_string_real_dp(val, format) result(str)
   real(wp), intent(in) :: val
   character(len=*), intent(in) :: format
   character(len=:), allocatable :: str

   character(len=128) :: buffer
   integer :: stat

   write(buffer, format, iostat=stat) val
   if (stat == 0) then
      str = trim(buffer)
   else
      str = "*"
   end if
end function format_string_real_dp

pure function format_string_int(val, format) result(str)
   integer, intent(in) :: val
   character(len=*), intent(in) :: format
   character(len=:), allocatable :: str

   character(len=128) :: buffer
   integer :: stat

   write(buffer, format, iostat=stat) val
   if (stat == 0) then
      str = trim(buffer)
   else
      str = "*"
   end if
end function format_string_int

end module tblite_output_format
