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

!> @file tblite/os.F90
!> Provides operating system dependent functionality.

!> Interfaces to the system libraries.
module tblite_os
   use, intrinsic :: iso_c_binding, only : c_char, c_null_char, c_int
   implicit none
   private

   public :: setenv
   public :: file_exists, delete_file


   interface
      function sys_putenv(name) result(stat) &
            & bind(c, name="putenv")
         import :: c_char, c_int
         character(len=c_char), intent(in) :: name(*)
         integer(c_int) :: stat
      end function
   end interface


contains


subroutine setenv(name, value, stat)
   character(len=*), intent(in) :: name
   character(len=*), intent(in) :: value
   integer, intent(out) :: stat

   stat = sys_putenv(as_c_char(name//"="//value))

end subroutine setenv


pure function as_c_char(str) result(res)
   character(len=*), intent(in) :: str
   character(kind=c_char) :: res(len(str)+1)
   res = transfer(str // c_null_char, res)
end function as_c_char


function file_exists(file) result(exist)
   character(len=*), intent(in) :: file
   logical :: exist
   inquire(file=file, exist=exist)
end function file_exists


subroutine delete_file(file)
   character(len=*), intent(in) :: file
   integer :: unit
   if (file_exists(file)) then
      open(file=file, newunit=unit)
      close(unit, status="delete")
   end if
end subroutine delete_file


end module tblite_os
