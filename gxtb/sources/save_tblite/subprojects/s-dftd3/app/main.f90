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

program dftd3_main
   use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
   use mctc_env, only : error_type
   use dftd3_app_cli, only : app_config, get_arguments
   use dftd3_app_driver, only : app_driver
   implicit none
   class(app_config), allocatable :: config
   type(error_type), allocatable :: error

   call get_arguments(config, error)
   call handle_error(error)

   call app_driver(config, error)
   call handle_error(error)

contains

subroutine handle_error(error)
   type(error_type), allocatable, intent(inout) :: error

   interface
      subroutine sys_exit(stat) bind(c, name="exit")
         use, intrinsic :: iso_c_binding, only : c_int
         integer(c_int), value :: stat
      end subroutine sys_exit
   end interface

   if (allocated(error)) then
      if (error%stat == 0) then
         write(output_unit, '(a)') error%message
         call sys_exit(0)
      else
         write(error_unit, '("[Error]", 1x, a)') error%message
         call sys_exit(1)
      end if
   end if
end subroutine handle_error

end program dftd3_main
