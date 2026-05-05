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

!> Entry point for the command line interface of dftd4
program driver
   use, intrinsic :: iso_fortran_env, only : error_unit
   use mctc_env, only : error_type
   use dftd4_cli, only : cli_config, get_arguments
   use dftd4_driver, only : main
   implicit none
   !> Configuration data deteriming the driver behaviour
   class(cli_config), allocatable :: config
   !> Error handling
   type(error_type), allocatable :: error

   call get_arguments(config, error)
   if (allocated(error)) then
      write(error_unit, '("[Error]:", 1x, a)') error%message
      error stop
   end if

   call main(config, error)
   if (allocated(error)) then
      write(error_unit, '("[Error]:", 1x, a)') error%message
      error stop
   end if

end program driver
