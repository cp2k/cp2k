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

!> @file tblite/api/utils.f90
!> Common utilities for all API exports

!> Provides utilities used for all API export modules
module tblite_api_utils
   use, intrinsic :: iso_c_binding
   implicit none
   private

   public :: f_c_character, c_f_character, strlen


   interface
      function strlen(str) result(len) bind(c)
         import :: c_ptr, c_int
         type(c_ptr), value :: str
         integer(c_int) :: len
      end function strlen
   end interface

contains


subroutine f_c_character(rhs, lhs, len)
   character(kind=c_char), intent(out) :: lhs(*)
   character(len=*), intent(in) :: rhs
   integer, intent(in) :: len
   integer :: length
   length = min(len-1, len_trim(rhs))

   lhs(1:length) = transfer(rhs(1:length), lhs(1:length))
   lhs(length+1:length+1) = c_null_char

end subroutine f_c_character

subroutine c_f_character(rhs, lhs)
   character(kind=c_char), intent(in) :: rhs(*)
   character(len=:), allocatable, intent(out) :: lhs

   integer :: ii

   do ii = 1, huge(ii) - 1
      if (rhs(ii) == c_null_char) then
         exit
      end if
   end do
   allocate(character(len=ii-1) :: lhs)
   lhs = transfer(rhs(1:ii-1), lhs)

end subroutine c_f_character

end module tblite_api_utils
