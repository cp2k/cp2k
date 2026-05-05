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

program main_driver
   use, intrinsic :: iso_fortran_env, only : error_unit
   use mctc_env, only : error_type
   use tblite_cli, only : get_arguments, driver_config, run_config, param_config, &
      & fit_config, tagdiff_config, guess_config
   use tblite_context_type, only : context_type
   use tblite_driver, only : main
   implicit none
   class(driver_config), allocatable :: config
   type(error_type), allocatable :: error

   call get_arguments(config, error)
   if (allocated(error)) then
      write(error_unit, '(a)') error%message
      error stop
   end if

   select type(config)
   type is (fit_config)
      call main(config, error)
   type is (param_config)
      call main(config, error)
   type is (run_config)
      call main(config, error)
   type is (tagdiff_config)
      call main(config, error)
   type is (guess_config)
      call main(config, error)
   end select
   if (allocated(error)) then
      write(error_unit, '(a)') error%message
      error stop
   end if
end program main_driver
